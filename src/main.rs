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

//! Tandem Repeats Finder — Rust port of the reference C implementation
//! (version 4.10.0).
//!
//! The command line interface and the control flow mirror `trf.c`; the output
//! is byte for byte identical to the C program.
//!
//! Module map (reference C source -> Rust module):
//!
//! | C source                | Rust module                                   |
//! |-------------------------|-----------------------------------------------|
//! | `trf.c`                 | `cli.rs` (command line), `main.rs` (driver)   |
//! | `trfrun.h`              | `control.rs` (control routine, `TRF()`)       |
//! | `tr30dat.h` + `tr30dat.c` | `state.rs`, `algo.rs`, `scan.rs`            |
//! | `trfclean.h`            | `clean.rs` (redundancy removal, HTML output)  |
//! | FASTA input             | `fasta.rs`                                    |
//! | generated tables        | `tables.rs`                                   |
//!
//! The C program keeps its state in file scope globals; here everything lives in
//! [`state::Trf`], which is created once per FASTA sequence.

mod algo;
mod clean;
mod cli;
mod control;
mod fasta;
mod params;
mod scan;
mod state;
mod tables;

use cli::MaxLength;
use params::{ParamSet, VERSION};

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
    let t = s.trim_start();
    if t.is_empty() {
        return None;
    }
    let (neg, rest) = if let Some(r) = t.strip_prefix('-') {
        (true, r)
    } else if let Some(r) = t.strip_prefix('+') {
        (false, r)
    } else {
        (false, t)
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
    let v = i64::from_str_radix(digits, radix).ok()?;
    Some(if neg { -v } else { v })
}

/// `ParseUInt()`: successful conversion and a non negative value.
fn parse_uint(s: &str, dest: &mut u32) -> bool {
    match parse_int(s) {
        Some(v) => {
            let v = v as i32; /* the reference implementation truncates to int */
            if v >= 0 {
                *dest = v as u32;
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
    if !parse_uint(maxperiod, &mut paramset.maxwraplength) {
        eprintln!("Error while parsing max TR length (option '-L') value");
        print_banner();
        return Some(1);
    }
    paramset.maxwraplength = (paramset.maxwraplength as f64 * 1e6) as u32;
    None
}

fn main() {
    /* The reference program does not ignore SIGPIPE: writing into a closed
    pipe terminates it.  Rust ignores SIGPIPE by default, which would turn the
    same situation into a panic, so restore the default disposition before any
    output is written. */
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    std::process::exit(real_main());
}

fn real_main() -> i32 {
    let cli = match cli::parse() {
        Ok(cli) => cli,
        Err(code) => return code,
    };

    let mut paramset = ParamSet::default();
    paramset.datafile = cli.datafile;
    paramset.maskedfile = cli.maskedfile;
    paramset.flankingsequence = cli.flankingsequence;
    paramset.htmloff = cli.htmloff;
    paramset.redundoff = cli.redundoff;
    paramset.ngs = cli.ngs;

    match cli.maxlength {
        MaxLength::Absent => {}
        MaxLength::MissingValue => {
            /* getopt reports the missing value and the program keeps going.
            `&av[8]` is the program name the reference program reports. */
            eprintln!("{}: option requires an argument -- 'l'", cli.maxperiod);
        }
        MaxLength::Value(value) => {
            if let Some(code) = set_max_length(&value, &cli.maxperiod, &mut paramset) {
                return code;
            }
        }
    }

    /* get input parameters */
    paramset.inputfilename = cli.file.clone();
    paramset.outputprefix = name_part(&cli.file);

    /* validate these parameters */
    if !parse_uint(&cli.match_weight, &mut paramset.match_weight) {
        paramset.endstatus =
            Some("Error parsing match parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.mismatch, &mut paramset.mismatch) {
        paramset.endstatus =
            Some("Error parsing mismatch parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.indel, &mut paramset.indel) {
        paramset.endstatus =
            Some("Error parsing indel parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.pm, &mut paramset.pm) {
        paramset.endstatus =
            Some("Error parsing PM parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.pi, &mut paramset.pi) {
        paramset.endstatus =
            Some("Error parsing PI parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.minscore, &mut paramset.minscore) {
        paramset.endstatus =
            Some("Error parsing Minscore parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&cli.maxperiod, &mut paramset.maxperiod)
        || paramset.maxperiod > 2000
        || paramset.maxperiod == 0
    {
        paramset.endstatus = Some(
            "Error parsing MaxPeriod parameter. Value must be between 1 and 2000, inclusive."
                .into(),
        );
    }

    /* error if any validation failed */
    if let Some(status) = paramset.endstatus.clone() {
        println!("{}\n", status);
        println!(
            "Please run with -h for help, or visit https://tandem.bu.edu/trf/trf.unix.help.html"
        );
        return 1;
    }

    /* paramset.datafile must be set if HTMLoff is set */
    paramset.datafile |= paramset.htmloff;

    if paramset.ngs {
        paramset.datafile = true;
    } else {
        print_banner();
    }

    if cli.file == "-" {
        if paramset.ngs {
            paramset.use_stdin = true;
            paramset.htmloff = true;
        } else {
            eprint!("\n\nPlease use -ngs flag if piping input into TRF");
            eprintln!();
            return -1;
        }
    }

    /* call the function that controls execution */
    control::trf_control_routine(&mut paramset);

    if paramset.outputcount == 0 && paramset.endstatus.is_none() {
        if !paramset.stdout_closed {
            println!("No TRs found. Exiting...");
        }
        return 0;
    } else if let Some(status) = paramset.endstatus.clone() {
        if !paramset.stdout_closed {
            println!("Error processing input: {}", status);
        }
        return 1;
    }
    0
}
