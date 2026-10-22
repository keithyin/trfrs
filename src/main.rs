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
//! | `trf.c`                 | `main.rs` (command line interface)            |
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
mod control;
mod fasta;
mod params;
mod scan;
mod state;
mod tables;

use params::{ParamSet, VERSION};

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

fn usage(program: &str) -> String {
    USAGE_BODY.replacen("%s", program, 1)
}

fn print_banner() {
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

/// Command line options understood by the program.
#[derive(Clone, Copy)]
enum Flag {
    Help,
    Version,
    Dat,
    Mask,
    Flank,
    HtmlOff,
    RedundOff,
    Ngs,
    MaxLength,
}

/// Long option names of `getopt_long_only()`, in table order.
const LONG_OPTIONS: &[(&str, Flag)] = &[
    ("help", Flag::Help),
    ("version", Flag::Version),
    ("dat", Flag::Dat),
    ("mask", Flag::Mask),
    ("flank", Flag::Flank),
    ("html-off", Flag::HtmlOff),
    ("redund-off", Flag::RedundOff),
    ("ngs", Flag::Ngs),
    ("Ngs", Flag::Ngs),
    ("NGS", Flag::Ngs),
    ("maxlength", Flag::MaxLength),
];

/// `getopt_long_only()` lookup: an exact match wins, otherwise a unique
/// prefix is accepted (long option names are case sensitive).
fn long_option(name: &str) -> Option<Flag> {
    if name.is_empty() {
        return None;
    }
    if let Some((_, flag)) = LONG_OPTIONS.iter().find(|(n, _)| *n == name) {
        return Some(*flag);
    }
    let mut matches = LONG_OPTIONS.iter().filter(|(n, _)| n.starts_with(name));
    let first = matches.next()?;
    if matches.next().is_some() {
        return None; /* ambiguous abbreviation */
    }
    Some(first.1)
}

/// Reports an option error the way glibc's `getopt` does.  The reference
/// program passes `&av[8]` to `getopt_long_only()`, so that is what glibc
/// uses as the program name in these messages.
fn getopt_error(prog: &str, message: &str) {
    eprintln!("{}: {}", prog, message);
}

/// Consumes the next argument as an option value (`optarg`).
fn take_value(args: &[String], idx: &mut usize) -> Option<String> {
    if *idx < args.len() {
        let v = args[*idx].clone();
        *idx += 1;
        Some(v)
    } else {
        None
    }
}

/// Handles `-l`/`-L`/`--maxlength`. Returns an exit code when the run must stop.
fn set_max_length(optarg: &str, args: &[String], paramset: &mut ParamSet) -> Option<i32> {
    if parse_int(optarg).unwrap_or(0) < 1 {
        eprintln!("Error: max TR length must be at least 1 million");
        print_banner();
        return Some(2);
    }

    /* Note: the reference implementation (4.10.0) reads the value from av[8]
    here, i.e. the MaxPeriod argument, instead of optarg.  This behaviour is
    reproduced so that both programs allocate the same amount of memory. */
    if !parse_uint(&args[8], &mut paramset.maxwraplength) {
        eprintln!("Error while parsing max TR length (option '-L') value");
        print_banner();
        return Some(1);
    }
    paramset.maxwraplength = (paramset.maxwraplength as f64 * 1e6) as u32;
    None
}

fn main() {
    std::process::exit(real_main());
}

fn real_main() -> i32 {
    let args: Vec<String> = std::env::args().collect();
    let ac = args.len();

    /* Handle a lone -v argument ourselves */
    if ac == 2 && (args[1] == "-v" || args[1] == "-V") {
        print_banner();
        return 0;
    }

    /* Expects exactly 8 non-option arguments */
    if ac < 9 {
        eprint!("{}", usage(&args[0]));
        return 1;
    }

    let mut paramset = ParamSet::default();

    /* Parse command line options; options start at argument 8.  The reference
    program passes `&av[8]` to `getopt_long_only()`, so that is also the
    program name used in getopt's error messages. */
    let prog = args[8].clone();
    let mut idx = 9usize;
    while idx < ac {
        let arg = args[idx].clone();
        idx += 1;
        if arg == "--" {
            break;
        }
        if !arg.starts_with('-') || arg == "-" {
            continue;
        }

        let (body, double_dash) = match arg.strip_prefix("--") {
            Some(stripped) => (stripped, true),
            None => (&arg[1..], false),
        };

        /* getopt_long_only(): try a long option first, then the short options */
        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v)),
            None => (body, None),
        };

        if let Some(flag) = long_option(name) {
            match flag {
                Flag::Help => {
                    print!("{}", usage(&args[0]));
                    return 0;
                }
                Flag::Version => {
                    print_banner();
                    return 0;
                }
                Flag::Dat => paramset.datafile = true,
                Flag::Mask => paramset.maskedfile = true,
                Flag::Flank => paramset.flankingsequence = true,
                Flag::HtmlOff => paramset.htmloff = true,
                Flag::RedundOff => paramset.redundoff = true,
                Flag::Ngs => paramset.ngs = true,
                Flag::MaxLength => {
                    let optarg = match inline {
                        Some(v) => v.to_string(),
                        None => match take_value(&args, &mut idx) {
                            Some(v) => v,
                            None => {
                                getopt_error(
                                    &prog,
                                    &format!("option '--{}' requires an argument", name),
                                );
                                continue;
                            }
                        },
                    };
                    if let Some(code) = set_max_length(&optarg, &args, &mut paramset) {
                        return code;
                    }
                }
            }
            continue;
        }

        if double_dash {
            getopt_error(&prog, &format!("unrecognized option '{}'", arg));
            continue;
        }

        /* short option cluster, e.g. "-dm" or "-l3" */
        let chars: Vec<char> = body.chars().collect();
        if chars.is_empty() {
            continue;
        }
        if !matches!(
            chars[0].to_ascii_lowercase(),
            'u' | 'v' | 'd' | 'm' | 'f' | 'h' | 'r' | 'l'
        ) {
            getopt_error(&prog, &format!("unrecognized option '{}'", arg));
            continue;
        }
        let mut ci = 0usize;
        while ci < chars.len() {
            match chars[ci].to_ascii_lowercase() {
                'u' => {
                    print!("{}", usage(&args[0]));
                    return 0;
                }
                'v' => {
                    print_banner();
                    return 0;
                }
                'd' => paramset.datafile = true,
                'm' => paramset.maskedfile = true,
                'f' => paramset.flankingsequence = true,
                'h' => paramset.htmloff = true,
                'r' => paramset.redundoff = true,
                'l' => {
                    let rest: String = chars[ci + 1..].iter().collect();
                    let optarg = if rest.is_empty() {
                        match take_value(&args, &mut idx) {
                            Some(v) => v,
                            None => {
                                getopt_error(
                                    &prog,
                                    &format!("option requires an argument -- '{}'", chars[ci]),
                                );
                                break;
                            }
                        }
                    } else {
                        rest
                    };
                    if let Some(code) = set_max_length(&optarg, &args, &mut paramset) {
                        return code;
                    }
                    break;
                }
                _ => {
                    /* getopt reports the character as it was given */
                    getopt_error(&prog, &format!("invalid option -- '{}'", chars[ci]));
                }
            }
            ci += 1;
        }
    }

    /* get input parameters */
    paramset.inputfilename = args[1].clone();
    paramset.outputprefix = name_part(&args[1]);

    /* validate these parameters */
    if !parse_uint(&args[2], &mut paramset.match_weight) {
        paramset.endstatus =
            Some("Error parsing match parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&args[3], &mut paramset.mismatch) {
        paramset.endstatus =
            Some("Error parsing mismatch parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&args[4], &mut paramset.indel) {
        paramset.endstatus =
            Some("Error parsing indel parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&args[5], &mut paramset.pm) {
        paramset.endstatus =
            Some("Error parsing PM parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&args[6], &mut paramset.pi) {
        paramset.endstatus =
            Some("Error parsing PI parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&args[7], &mut paramset.minscore) {
        paramset.endstatus =
            Some("Error parsing Minscore parameter. Value must be a positive integer.".into());
    } else if !parse_uint(&args[8], &mut paramset.maxperiod)
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

    if args[1] == "-" {
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
