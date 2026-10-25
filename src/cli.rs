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
//! `clap` provides the clustering, the abbreviations and the attached values;
//! a small pre-pass spells single dash long options with `--` before the
//! vector is handed to `clap`.  The eight positional parameters are kept as
//! strings: validating them is the job of the caller, which reproduces the
//! error messages of the reference program.

use clap::{Arg, ArgAction, Command};

use crate::print_banner;

/// The `usage()` text of the reference program, `%s` being the program name.
const USAGE_BODY: &str =
    "\n\nPlease use: %s File Match Mismatch Delta PM PI Minscore MaxPeriod [options]\n\
\nWhere: (all weights, penalties, and scores are positive)\
\n  File = sequences input file\
\n  Match  = matching weight\
\n  Mismatch  = mismatching penalty\
\n  Delta = indel penalty\
\n  PM = match probability (whole number)\
\n  PI = indel probability (whole number)\
\n  Minscore = minimum alignment score to report\
\n  MaxPeriod = maximum period size to report. Must be between 1 and 2000, inclusive\
\n  [options] = one or more of the following:\
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
];

/// Positional parameter names, in the order `av[1]` .. `av[8]`.
const POSITIONALS: [&str; 8] = [
    "file",
    "match",
    "mismatch",
    "delta",
    "pm",
    "pi",
    "minscore",
    "maxperiod",
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

/// The `clap` definition of the command line.
fn command() -> Command {
    Command::new("trf")
        .disable_help_flag(true) /* -h is "suppress html output" */
        .disable_version_flag(true) /* -v and -V both print the banner */
        .infer_long_args(true) /* unique prefixes, as getopt_long_only() */
        .args_override_self(true) /* repeated options are accepted */
        .arg(Arg::new("file").value_name("File").index(1))
        .arg(Arg::new("match").value_name("Match").index(2))
        .arg(Arg::new("mismatch").value_name("Mismatch").index(3))
        .arg(Arg::new("delta").value_name("Delta").index(4))
        .arg(Arg::new("pm").value_name("PM").index(5))
        .arg(Arg::new("pi").value_name("PI").index(6))
        .arg(Arg::new("minscore").value_name("Minscore").index(7))
        .arg(Arg::new("maxperiod").value_name("MaxPeriod").index(8))
        .arg(
            Arg::new("dat")
                .short('d')
                .visible_short_alias('D')
                .long("dat")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("mask")
                .short('m')
                .visible_short_alias('M')
                .long("mask")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("flank")
                .short('f')
                .visible_short_alias('F')
                .long("flank")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("html-off")
                .short('h')
                .visible_short_alias('H')
                .long("html-off")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("redund-off")
                .short('r')
                .visible_short_alias('R')
                .long("redund-off")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("ngs")
                .long("ngs")
                .alias("Ngs")
                .alias("NGS")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("maxlength")
                .short('l')
                .visible_short_alias('L')
                .long("maxlength")
                .value_name("n")
                /* getopt reports a missing value but carries on */
                .num_args(0..=1)
                .default_missing_value("")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("help")
                .short('u')
                .visible_short_alias('U')
                .long("help")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("version")
                .short('v')
                .visible_short_alias('V')
                .long("version")
                .action(ArgAction::SetTrue),
        )
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

    pub datafile: bool,
    pub maskedfile: bool,
    pub flankingsequence: bool,
    pub htmloff: bool,
    pub redundoff: bool,
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
    at any option: a lone -v/-V prints the banner, a short vector the usage
    message. */
    if argv.len() == 2 && (argv[1] == "-v" || argv[1] == "-V") {
        print_banner();
        return Err(0);
    }
    if argv.len() < 9 {
        eprint!("{}", usage(&argv[0]));
        return Err(1);
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

    let matches = match command().try_get_matches_from(args) {
        Ok(matches) => matches,
        Err(error) => {
            let _ = error.print();
            return Err(error.exit_code());
        }
    };

    if matches.get_flag("help") {
        print!("{}", usage(&argv[0]));
        return Err(0);
    }
    if matches.get_flag("version") {
        print_banner();
        return Err(0);
    }

    let mut positional = Vec::with_capacity(POSITIONALS.len());
    for name in POSITIONALS {
        match matches.get_one::<String>(name) {
            Some(value) => positional.push(value.clone()),
            None => {
                eprint!("{}", usage(&argv[0]));
                return Err(1);
            }
        }
    }
    let [file, match_weight, mismatch, indel, pm, pi, minscore, maxperiod] =
        <[String; POSITIONALS.len()]>::try_from(positional).expect("eight positional values");

    let maxlength = match matches.value_source("maxlength") {
        Some(_) => match matches.get_one::<String>("maxlength") {
            Some(value) if !value.is_empty() => MaxLength::Value(value.clone()),
            _ => MaxLength::MissingValue,
        },
        None => MaxLength::Absent,
    };

    Ok(Cli {
        file,
        match_weight,
        mismatch,
        indel,
        pm,
        pi,
        minscore,
        maxperiod,
        datafile: matches.get_flag("dat"),
        maskedfile: matches.get_flag("mask"),
        flankingsequence: matches.get_flag("flank"),
        htmloff: matches.get_flag("html-off"),
        redundoff: matches.get_flag("redund-off"),
        ngs: matches.get_flag("ngs"),
        maxlength,
    })
}
