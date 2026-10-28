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

//! Command line driver (`trf.c`).
//!
//! This is the former `main.rs` body: it parses the options, validates the
//! search parameters and hands control to [`crate::control`].  It is kept out
//! of `lib.rs` so that the public library surface stays small; the binary is a
//! thin wrapper around [`run_cli`].

use crate::cli::{self, MaxLength};
use crate::control;
use crate::params::{ParamSet, VERSION};

/// Prints the version banner (also used when `-v`/`--version` is given).
pub(crate) fn print_banner() {
    eprint!("\nTandem Repeats Finder, Version {}", VERSION);
    eprint!("\nCopyright (C) Dr. Gary Benson 1999-2012. All rights reserved.\n");
}

/// `GetNamePartAddress()`: file name without the directory part.
fn name_part(name: &str) -> String {
    match name.rfind('/') {
        Some(pos) => name[pos + 1..].to_string(),
        None => name.to_string(),
    }
}

/// `ParseInt()`: `strtol(str, &temp, 0)` plus a full consumption check.
fn parse_int(s: &str) -> Option<i64> {
    let trimmed = s.trim_start();
    if trimmed.is_empty() {
        return None;
    }
    let (neg, rest) = if let Some(r) = trimmed.strip_prefix('-') {
        (true, r)
    } else if let Some(r) = trimmed.strip_prefix('+') {
        (false, r)
    } else {
        (false, trimmed)
    };
    let (radix, digits) =
        if let Some(r) = rest.strip_prefix("0x").or_else(|| rest.strip_prefix("0X")) {
            (16u32, r)
        } else if rest.len() > 1 && rest.starts_with('0') {
            (8u32, &rest[1..])
        } else {
            (10u32, rest)
        };
    if digits.is_empty() {
        if rest == "0" {
            return Some(0);
        }
        return None;
    }
    if !digits.chars().all(|c| c.is_digit(radix)) {
        return None;
    }
    let value = i64::from_str_radix(digits, radix).ok()?;
    Some(if neg { -value } else { value })
}

/// `ParseUInt()`: successful conversion and a non negative value.
fn parse_uint(s: &str, dest: &mut u32) -> bool {
    match parse_int(s) {
        Some(v) => {
            let value = v as i32; /* the reference implementation truncates to int */
            if value >= 0 {
                *dest = value as u32;
                true
            } else {
                false
            }
        }
        None => false,
    }
}

/// Handles `-l`/`-L`/`--maxlength`. Returns an exit code when the run must stop.
fn set_max_length(optarg: &str, maxperiod: &str, paramset: &mut ParamSet) -> Option<i32> {
    if parse_int(optarg).unwrap_or(0) < 1 {
        eprintln!("Error: max TR length must be at least 1 million");
        print_banner();
        return Some(2);
    }

    /* Note: the reference implementation (4.10.0) reads the value from av[8]
    here, i.e. the MaxPeriod argument, instead of optarg.  This behaviour is
    reproduced so that both programs allocate the same amount of memory. */
    if !parse_uint(maxperiod, &mut paramset.max_wrap_length) {
        eprintln!("Error while parsing max TR length (option '-L') value");
        print_banner();
        return Some(1);
    }
    paramset.max_wrap_length = (paramset.max_wrap_length as f64 * 1e6) as u32;
    None
}

/// Runs the command line program and returns its exit code (`main()` of the
/// reference implementation).
pub fn run_cli() -> i32 {
    let cli = match cli::parse() {
        Ok(cli) => cli,
        Err(code) => return code,
    };

    let mut paramset = ParamSet {
        data_file: cli.data_file,
        masked_file: cli.masked_file,
        flanking_sequence: cli.flanking_sequence,
        html_output_off: cli.html_output_off,
        redundancy_removal_off: cli.redundancy_removal_off,
        ngs: cli.ngs,
        ..ParamSet::default()
    };

    match cli.maxlength {
        MaxLength::Absent => {}
        MaxLength::MissingValue => {
            /* getopt reports the missing value and the program keeps going.
            The reference implementation prints `&av[8]` here, i.e. the
            MaxPeriod argument, which is what `cli.maxperiod` holds (its
            default when the option is omitted). */
            eprintln!("{}: option requires an argument -- 'l'", cli.maxperiod);
        }
        MaxLength::Value(value) => {
            if let Some(code) = set_max_length(&value, &cli.maxperiod, &mut paramset) {
                return code;
            }
        }
    }

    /* get input parameters */
    paramset.input_file_name = cli.file.clone();
    paramset.output_prefix = name_part(&cli.file);

    /* validate these parameters */
    if !parse_uint(&cli.match_weight, &mut paramset.match_weight) {
        paramset.end_status =
            Some("Error parsing match parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.mismatch, &mut paramset.mismatch) {
        paramset.end_status =
            Some("Error parsing mismatch parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.indel, &mut paramset.indel) {
        paramset.end_status =
            Some("Error parsing indel parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.pm, &mut paramset.match_probability_percent) {
        paramset.end_status =
            Some("Error parsing PM parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.pi, &mut paramset.indel_probability_percent) {
        paramset.end_status =
            Some("Error parsing PI parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.minscore, &mut paramset.minscore) {
        paramset.end_status =
            Some("Error parsing Minscore parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.maxperiod, &mut paramset.maxperiod)
        || paramset.maxperiod > 2000
        || paramset.maxperiod == 0
    {
        paramset.end_status = Some(
            "Error parsing MaxPeriod parameter. Value must be between 1 and 2000, inclusive."
                .into(),
        );
    }

    /* error if any validation failed */
    if let Some(status) = paramset.end_status.clone() {
        println!("{}\n", status);
        println!(
            "Please run with -h for help, or visit https://tandem.bu.edu/trf/trf.unix.help.html"
        );
        return 1;
    }

    /* paramset.data_file must be set if HTMLoff is set */
    paramset.data_file |= paramset.html_output_off;

    if paramset.ngs {
        paramset.data_file = true;
    } else {
        print_banner();
    }

    if cli.file == "-" {
        if paramset.ngs {
            paramset.use_stdin = true;
            paramset.html_output_off = true;
        } else {
            eprint!("\n\nPlease use -ngs flag if piping input into TRF");
            eprintln!();
            return -1;
        }
    }

    /* call the function that controls execution */
    control::trf_control_routine(&mut paramset);

    if paramset.output_count == 0 && paramset.end_status.is_none() {
        if !paramset.stdout_closed {
            println!("No TRs found. Exiting...");
        }
        return 0;
    } else if let Some(status) = paramset.end_status.clone() {
        if !paramset.stdout_closed {
            println!("Error processing input: {}", status);
        }
        return 1;
    }
    0
}
