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

//! Command line interface.
//!
//! The interface is designed for this port instead of being transcribed from
//! the reference implementation's `getopt_long_only()` parsing.  It follows
//! conventional Rust/Unix command line practice:
//!
//! * long options are introduced by `--` (`--ngs`, `--no-html`),
//! * values are typed and range checked while parsing (a bad value is a
//!   command line error, not a run time failure),
//! * `--help` and `--version` behave as users expect,
//! * exit codes follow the common scheme: `0` success, `1` run time error,
//!   `2` command line error.
//!
//! The search parameters keep the meaning (and the defaults `2 5 7 80 10 50
//! 2000`) of the reference program.  `--str` is a preset that only narrows the
//! reported period range to 1–15 bp; explicit `--minperiod`/`--maxperiod`
//! take precedence.

use std::path::PathBuf;

use clap::Parser;

use crate::params::STR_MAX_PERIOD;

/// Default of `--match`: the customary human DNA set.
const DEFAULT_MATCH: u32 = 2;

/// Default of `--mismatch`.
const DEFAULT_MISMATCH: u32 = 5;

/// Default of `--indel`.
const DEFAULT_INDEL: u32 = 7;

/// Default of `--pm`.
const DEFAULT_PM: u32 = 80;

/// Default of `--pi`.
const DEFAULT_PI: u32 = 10;

/// Default of `--minscore`.
const DEFAULT_MINSCORE: u32 = 50;

/// Default of `--maxperiod` outside `--str` mode.
const DEFAULT_MAXPERIOD: u32 = 2000;

/// Default of `--max-length`, in millions of bases.
const DEFAULT_MAX_LENGTH: u32 = 2;

/// The `clap` definition of the command line (Derive API).
///
/// The values are parsed into their final types here, so a bad value is
/// reported by `clap` with exit code 2 before the program does any work.
/// Cross-field checks (`--pm` table support, the period range) live in
/// [`parse`], which is why the period options are kept optional.
#[derive(Parser, Debug)]
#[command(
    name = "trfrs",
    version,
    about = "Find tandem repeats in DNA sequences.",
    long_about = "A Rust port of Tandem Repeats Finder: searches a FASTA or FASTQ file for \
tandem repeats and writes HTML tables, alignment reports and machine readable \
data files.\n\nThe search parameter defaults are the customary human DNA set \
2 5 7 80 10 50 2000, so naming the input file is enough to run a search.",
    args_override_self = true
)]
struct Args {
    /// Input sequence file (FASTA or FASTQ); `-` reads from stdin
    #[arg(value_name = "FILE")]
    file: PathBuf,

    /// Matching weight
    #[arg(
        long = "match",
        value_name = "N",
        default_value_t = DEFAULT_MATCH,
        value_parser = clap::value_parser!(u32).range(1..=i32::MAX as i64),
        help_heading = "Search parameters"
    )]
    match_weight: u32,

    /// Mismatching penalty
    #[arg(
        long,
        value_name = "N",
        default_value_t = DEFAULT_MISMATCH,
        value_parser = clap::value_parser!(u32).range(1..=i32::MAX as i64),
        help_heading = "Search parameters"
    )]
    mismatch: u32,

    /// Insertion/deletion penalty
    #[arg(
        long,
        value_name = "N",
        default_value_t = DEFAULT_INDEL,
        value_parser = clap::value_parser!(u32).range(1..=i32::MAX as i64),
        help_heading = "Search parameters"
    )]
    indel: u32,

    /// Match probability in percent (only 75 and 80 have statistical tables)
    #[arg(
        long,
        value_name = "N",
        default_value_t = DEFAULT_PM,
        value_parser = clap::value_parser!(u32).range(0..=100),
        help_heading = "Search parameters"
    )]
    pm: u32,

    /// Insertion/deletion probability in percent
    #[arg(
        long,
        value_name = "N",
        default_value_t = DEFAULT_PI,
        value_parser = clap::value_parser!(u32).range(0..=100),
        help_heading = "Search parameters"
    )]
    pi: u32,

    /// Minimum alignment score to report
    #[arg(
        long,
        value_name = "N",
        default_value_t = DEFAULT_MINSCORE,
        value_parser = clap::value_parser!(u32).range(1..=i32::MAX as i64),
        help_heading = "Search parameters"
    )]
    minscore: u32,

    /// Minimum repeat unit length to report
    #[arg(
        long,
        value_name = "N",
        value_parser = clap::value_parser!(u32).range(1..=2000),
        help_heading = "Search parameters"
    )]
    minperiod: Option<u32>,

    /// Maximum repeat unit length to report
    #[arg(
        long,
        value_name = "N",
        value_parser = clap::value_parser!(u32).range(1..=2000),
        help_heading = "Search parameters"
    )]
    maxperiod: Option<u32>,

    /// Short tandem repeat mode: period 1-15 (explicit --minperiod/--maxperiod win)
    #[arg(long, help_heading = "Search parameters")]
    str: bool,

    /// Output directory for the report and data files (created if missing)
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output options"
    )]
    output: PathBuf,

    /// Write the .dat data file
    #[arg(short = 'd', long, help_heading = "Output options")]
    dat: bool,

    /// Write the .csv table (one row per repeat)
    #[arg(long, help_heading = "Output options")]
    csv: bool,

    /// Write the masked sequence file
    #[arg(short = 'm', long, help_heading = "Output options")]
    mask: bool,

    /// Include the flanking sequences around each repeat in the report
    #[arg(short = 'f', long, help_heading = "Output options")]
    flank: bool,

    /// Write the HTML report (repeat tables and alignments)
    #[arg(long, help_heading = "Output options")]
    html: bool,

    /// Keep redundant repeats
    #[arg(long, help_heading = "Output options")]
    no_redundancy: bool,

    /// Print compact single-line records to stdout instead of writing files
    #[arg(short = 'n', long, help_heading = "Output options")]
    ngs: bool,

    /// Suppress progress messages
    #[arg(short = 'q', long, help_heading = "Output options")]
    quiet: bool,

    /// Expected maximum tandem repeat length, in millions of bases
    #[arg(
        short = 'l',
        long = "max-length",
        value_name = "N",
        default_value_t = DEFAULT_MAX_LENGTH,
        value_parser = clap::value_parser!(u32).range(1..=2000),
        help_heading = "Memory options"
    )]
    max_length: u32,
}

/// The search parameters, validated and resolved.
#[derive(Clone, Copy, Debug)]
pub struct SearchParams {
    pub match_weight: u32,
    pub mismatch: u32,
    pub indel: u32,
    pub match_probability_percent: u32,
    pub indel_probability_percent: u32,
    pub minscore: u32,
    pub min_period: u32,
    pub maxperiod: u32,
}

/// The output options.
#[derive(Clone, Debug)]
pub struct OutputOptions {
    /// Directory the report and data files are written to (`--output`).
    pub directory: PathBuf,
    /// Write the `.dat` data file (`--dat`).
    pub data_file: bool,
    /// Write the `.csv` table (`--csv`).
    pub csv_file: bool,
    /// Write the masked sequence file (`--mask`).
    pub masked_file: bool,
    /// Include flanking sequences in the report (`--flank`).
    pub flanking_sequence: bool,
    /// Write the HTML report (`--html`); off by default.
    pub html_output: bool,
    /// Keep redundant repeats (`--no-redundancy`).
    pub no_redundancy: bool,
    /// Print compact records to stdout instead of writing files (`--ngs`).
    pub ngs: bool,
    /// Suppress progress messages (`--quiet`).
    pub quiet: bool,
}

/// Validated command line arguments.
#[derive(Clone, Debug)]
pub struct Cli {
    pub file: PathBuf,
    pub search: SearchParams,
    pub output: OutputOptions,
    /// Expected maximum tandem repeat length in millions of bases
    /// (`--max-length`).
    pub max_length_millions: u32,
}

/// Parses `std::env::args()`.
///
/// Help, version and malformed command lines are handled by `clap` (help and
/// version exit with 0, a command line error with 2).  The cross-field checks
/// report to stderr and exit with 2, so a [`Cli`] is always valid when this
/// returns `Ok`.
pub fn parse() -> Result<Cli, i32> {
    let args = Args::parse();

    /* Period range: an explicit option wins, then the `--str` preset, then
    the built in default. */
    let min_period = args.minperiod.unwrap_or(1);
    let maxperiod = match args.maxperiod {
        Some(value) => value,
        None if args.str => STR_MAX_PERIOD,
        None => DEFAULT_MAXPERIOD,
    };

    if min_period > maxperiod {
        eprintln!(
            "error: invalid period range: --minperiod ({min_period}) must not exceed \
--maxperiod ({maxperiod})\n\nFor more information, try '--help'."
        );
        return Err(2);
    }

    /* Only these two match probabilities have statistical tables. */
    if args.pm != 75 && args.pm != 80 {
        eprintln!(
            "error: invalid value '{pm}' for '--pm <N>': only 75 and 80 are supported \
(no statistical tables exist for other values)\n\nFor more information, try '--help'.",
            pm = args.pm
        );
        return Err(2);
    }

    Ok(Cli {
        file: args.file,
        search: SearchParams {
            match_weight: args.match_weight,
            mismatch: args.mismatch,
            indel: args.indel,
            match_probability_percent: args.pm,
            indel_probability_percent: args.pi,
            minscore: args.minscore,
            min_period,
            maxperiod,
        },
        output: OutputOptions {
            directory: args.output,
            data_file: args.dat,
            csv_file: args.csv,
            masked_file: args.mask,
            flanking_sequence: args.flank,
            html_output: args.html,
            no_redundancy: args.no_redundancy,
            ngs: args.ngs,
            quiet: args.quiet,
        },
        max_length_millions: args.max_length,
    })
}
