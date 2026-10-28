/*
Tandem Repeats Finder - Rust port of the reference C implementation
Copyright (C) 1999-2020 Gary Benson (original C program)

This file is part of the Tandem Repeats Finder (TRF) program.

TRF is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

TRF is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with TRF.  If not, see <https://www.gnu.org/licenses/>.
*/

//! Command line interface (`trf.c`).
//!
//! The reference implementation parses its options with
//! `getopt_long_only()`, which is unusually permissive:
//!
//! * long options may be introduced by a single dash (`-ngs`),
//! * long options may be abbreviated to a unique prefix (`--maxl`),
//! * short options may be clustered (`-df`) and take attached values (`-l3`),
//! * short option letters are matched case insensitively (`-D`).
//!
//! `clap` (Derive API) provides the clustering, the abbreviations and the
//! attached values; a small pre-pass spells single dash long options with
//! `--` before the vector is handed to clap.  The file name is the only
//! positional parameter; the seven search parameters are key-value options
//! with defaults.  They are kept as strings: validating them is the job of
//! the caller, which reproduces the error messages of the reference program.

use clap::Parser;

use crate::driver::print_banner;

/// Defaults of the search parameters; the values are the customary human DNA
/// set `2 5 7 80 10 50 2000`.
pub const DEFAULT_MATCH: &str = "2";
pub const DEFAULT_MISMATCH: &str = "5";
pub const DEFAULT_INDEL: &str = "7";
pub const DEFAULT_PM: &str = "80";
pub const DEFAULT_PI: &str = "10";
pub const DEFAULT_MINSCORE: &str = "50";
pub const DEFAULT_MAXPERIOD: &str = "2000";

/// The `usage()` text of the reference program, `%s` being the program name.
const USAGE_BODY: &str = "\n\nPlease use: %s File [options]\n\
\nWhere: (all weights, penalties, and scores are positive)\
\n  File = sequences input file\
\n  [options] = one or more of the following:\
\n        --match <n>       matching weight (default 2)\
\n        --mismatch <n>    mismatching penalty (default 5)\
\n        --indel <n>       indel penalty (default 7)\
\n        --pm <n>          match probability (whole number, default 80)\
\n        --pi <n>          indel probability (whole number, default 10)\
\n        --minscore <n>    minimum alignment score to report (default 50)\
\n        --maxperiod <n>   maximum period size to report. Must be between 1 and 2000,\
\n                          inclusive (default 2000)\
\n        -m        masked sequence file\
\n        -f        flanking sequence\
\n        -d        data file\
\n        -h        suppress html output\
\n        -r        no redundancy elimination\
\n        -l <n>    maximum TR length expected (in millions) (eg, -l 3 or -l=3 for 3 million)\
\n                  Human genome HG38 would need -l 6\
\n        -ngs      more compact .dat output on multisequence files, returns 0 on success.\
\n                  Output is printed to the screen, not a file. You may pipe input in with\
\n                  this option using - for file name. Short 50 flanks are appended to .dat\
\n                  output.\
\n\
\nSee more information on the TRF Unix Help web page: https://tandem.bu.edu/trf/trf.unix.help.html\
\n\
\nNote the sequence file should be in FASTA format:\
\n\
\n>Name of sequence\
\naggaaacctgccatggcctcctggtgagctgtcctcatccactgctcgctgcctctccag\
\natactctgacccatggatcccctgggtgcagccaagccacaatggccatggcgccgctgt\
\nactcccacccgccccaccctcctgatcctgctatggacatggcctttccacatccctgtg\
\n\
\n\
\nTandem Repeats Finder\
\nCopyright (C) 1999-2020 Gary Benson\
\n\
\nThis program is free software: you can redistribute it and/or modify\
\nit under the terms of the GNU Affero General Public License as\
\npublished by the Free Software Foundation, either version 3 of\
\nthe License, or (at your option) any later version.\
\n";

/// `usage()`: the help text with the program name substituted.
fn usage(program: &str) -> String {
    USAGE_BODY.replacen("%s", program, 1)
}

/// Long option names of the `getopt_long_only()` table, mapped to the
/// canonical `clap` long name.  The three spellings of `-ngs` select one
/// option.
const LONG_OPTIONS: &[(&str, &str)] = &[
    ("help", "help"),
    ("version", "version"),
    ("dat", "dat"),
    ("mask", "mask"),
    ("flank", "flank"),
    ("html-off", "html-off"),
    ("redund-off", "redund-off"),
    ("ngs", "ngs"),
    ("Ngs", "ngs"),
    ("NGS", "ngs"),
    ("maxlength", "maxlength"),
    ("match", "match"),
    ("mismatch", "mismatch"),
    ("indel", "indel"),
    ("pm", "pm"),
    ("pi", "pi"),
    ("minscore", "minscore"),
    ("maxperiod", "maxperiod"),
];

/// `getopt_long_only()` lookup: an exact match wins, otherwise a unique
/// prefix is accepted (long option names are case sensitive).
fn long_option(name: &str) -> Option<&'static str> {
    if name.is_empty() {
        return None;
    }
    if let Some((_, canonical)) = LONG_OPTIONS.iter().find(|(n, _)| *n == name) {
        return Some(canonical);
    }
    let mut matches = LONG_OPTIONS.iter().filter(|(n, _)| n.starts_with(name));
    let first = matches.next()?;
    if matches.next().is_some() {
        return None; /* ambiguous abbreviation */
    }
    Some(first.1)
}

/// Rewrites a single dash long option (`-ngs`, `-maxl=3`) into the `--`
/// spelling understood by `clap`.  Short options and short option clusters
/// are returned unchanged.
fn rewrite_long_option(arg: &str) -> Option<String> {
    let body = arg.strip_prefix('-')?;
    if body.is_empty() || body.starts_with('-') {
        return None;
    }
    let (name, value) = match body.split_once('=') {
        Some((name, value)) => (name, Some(value)),
        None => (body, None),
    };
    let canonical = long_option(name)?;
    Some(match value {
        Some(value) => format!("--{canonical}={value}"),
        None => format!("--{canonical}"),
    })
}

/// The `clap` definition of the command line (Derive API).
///
/// `File` is the only positional parameter; the seven search parameters are
/// key-value options (`--match 2`, `--match=2`, or the single dash spelling
/// `-match 2`) with defaults, so an invocation such as `trf file.fasta`
/// searches with the customary human DNA parameters.  They are kept as
/// strings, validating them is the job of the caller, which reproduces the
/// error messages of the reference program.
/// Every other behaviour mirrors `getopt_long_only()`: `-h` means
/// "suppress html output" while help lives on `-u`, `-v` prints the banner,
/// short option letters are matched case insensitively and the three
/// spellings of `-ngs` select one option.
#[derive(Parser)]
#[command(
    name = "trf",
    disable_help_flag = true, /* -h is "suppress html output" */
    disable_version_flag = true, /* -v and -V both print the banner */
    infer_long_args = true,   /* unique prefixes, as getopt_long_only() */
    args_override_self = true /* repeated options are accepted */
)]
struct Args {
    /* Optional in `clap` so that option-only command lines such as `trfrs -h`
    still reach `parse()`, which prints the usage text instead of a clap
    "missing required <File>" diagnostic. */
    #[arg(value_name = "File")]
    file: Option<String>,
    #[arg(
        long = "match",
        value_name = "Match",
        default_value = DEFAULT_MATCH
    )]
    match_weight: String,
    #[arg(
        long = "mismatch",
        value_name = "Mismatch",
        default_value = DEFAULT_MISMATCH
    )]
    mismatch: String,
    #[arg(long = "indel", value_name = "Delta", default_value = DEFAULT_INDEL)]
    indel: String,
    #[arg(long = "pm", value_name = "PM", default_value = DEFAULT_PM)]
    pm: String,
    #[arg(long = "pi", value_name = "PI", default_value = DEFAULT_PI)]
    pi: String,
    #[arg(
        long = "minscore",
        value_name = "Minscore",
        default_value = DEFAULT_MINSCORE
    )]
    minscore: String,
    #[arg(
        long = "maxperiod",
        value_name = "MaxPeriod",
        default_value = DEFAULT_MAXPERIOD
    )]
    maxperiod: String,

    #[arg(short = 'd', visible_short_alias = 'D', long = "dat")]
    dat: bool,
    #[arg(short = 'm', visible_short_alias = 'M', long = "mask")]
    mask: bool,
    #[arg(short = 'f', visible_short_alias = 'F', long = "flank")]
    flank: bool,
    #[arg(short = 'h', visible_short_alias = 'H', long = "html-off")]
    html_off: bool,
    #[arg(short = 'r', visible_short_alias = 'R', long = "redund-off")]
    redund_off: bool,
    #[arg(long = "ngs", alias = "Ngs", alias = "NGS")]
    ngs: bool,
    #[arg(
        short = 'l',
        visible_short_alias = 'L',
        long = "maxlength",
        value_name = "n",
        /* getopt reports a missing value but carries on */
        num_args = 0..=1,
        default_missing_value = ""
    )]
    maxlength: Option<String>,
    #[arg(short = 'u', visible_short_alias = 'U', long = "help")]
    help: bool,
    #[arg(short = 'v', visible_short_alias = 'V', long = "version")]
    version: bool,
}

/// The value of the `-l` / `--maxlength` option.
pub enum MaxLength {
    /// The option was not given.
    Absent,
    /// The option was given without a value; `getopt` reports it and the
    /// program continues.
    MissingValue,
    /// The option was given with a value.
    Value(String),
}

/// Command line arguments.
pub struct Cli {
    pub file: String,
    pub match_weight: String,
    pub mismatch: String,
    pub indel: String,
    pub pm: String,
    pub pi: String,
    pub minscore: String,
    pub maxperiod: String,

    pub data_file: bool,
    pub masked_file: bool,
    pub flanking_sequence: bool,
    pub html_output_off: bool,
    pub redundancy_removal_off: bool,
    pub ngs: bool,
    pub maxlength: MaxLength,
}

/// Parses `std::env::args()`.
///
/// Returns the parsed arguments, or the exit code to terminate with.  Help,
/// version and usage messages are printed here, as is the `clap` diagnostic
/// for a malformed command line.
pub fn parse() -> Result<Cli, i32> {
    let argv: Vec<String> = std::env::args().collect();

    /* The reference program inspects the raw argument vector before it looks
    at any option: a lone -v/-V prints the banner. */
    if argv.len() == 2 && (argv[1] == "-v" || argv[1] == "-V") {
        print_banner();
        return Err(0);
    }

    let args: Vec<String> = argv
        .iter()
        .enumerate()
        .map(|(index, arg)| {
            if index == 0 {
                arg.clone()
            } else {
                rewrite_long_option(arg).unwrap_or_else(|| arg.clone())
            }
        })
        .collect();

    let args = match Args::try_parse_from(args) {
        Ok(args) => args,
        Err(error) => {
            let _ = error.print();
            return Err(error.exit_code());
        }
    };

    if args.help {
        print!("{}", usage(&argv[0]));
        return Err(0);
    }
    if args.version {
        print_banner();
        return Err(0);
    }

    /* No file name: print the usage text, as the reference program does for a
    command line without arguments (this covers e.g. a lone `-h`). */
    let file = match args.file {
        Some(file) => file,
        None => {
            eprint!("{}", usage(&argv[0]));
            return Err(1);
        }
    };

    let maxlength = match args.maxlength {
        Some(value) if !value.is_empty() => MaxLength::Value(value),
        Some(_) => MaxLength::MissingValue,
        None => MaxLength::Absent,
    };

    Ok(Cli {
        file,
        match_weight: args.match_weight,
        mismatch: args.mismatch,
        indel: args.indel,
        pm: args.pm,
        pi: args.pi,
        minscore: args.minscore,
        maxperiod: args.maxperiod,
        data_file: args.dat,
        masked_file: args.mask,
        flanking_sequence: args.flank,
        html_output_off: args.html_off,
        redundancy_removal_off: args.redund_off,
        ngs: args.ngs,
        maxlength,
    })
}
