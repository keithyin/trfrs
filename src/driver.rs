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

//! Command line driver.
//!
//! This is the former `main.rs` body: it turns the validated [`cli::Cli`]
//! into a [`ParamSet`] and hands control to [`crate::control`].  It is kept
//! out of `lib.rs` so that the public library surface stays small; the binary
//! is a thin wrapper around [`run_cli`].
//!
//! Exit code scheme: `0` on success (including a run that finds no repeats),
//! `1` on a run time error, `2` on a command line error.

use std::path::Path;

use crate::cli;
use crate::control;
use crate::params::{ParamSet, VERSION};

/// Prints the version banner (only shown in `--verbose` mode, never in
/// `--ngs` mode).
pub(crate) fn print_banner() {
    eprint!("\nTandem Repeats Finder, Version {}", VERSION);
    eprint!("\nCopyright (C) Dr. Gary Benson 1999-2012. All rights reserved.\n");
}

/// The output file prefix: the input file name without its directory part.
fn name_part(name: &Path) -> String {
    let text = name.to_string_lossy();
    name.file_name()
        .map(|part| part.to_string_lossy().into_owned())
        .unwrap_or_else(|| text.into_owned())
}

/// Runs the command line program and returns its exit code.
pub fn run_cli() -> i32 {
    let cli = match cli::parse() {
        Ok(cli) => cli,
        Err(code) => return code,
    };

    let mut paramset = ParamSet {
        match_weight: cli.search.match_weight,
        mismatch: cli.search.mismatch,
        indel: cli.search.indel,
        minscore: cli.search.minscore,
        min_period: cli.search.min_period,
        maxperiod: cli.search.maxperiod,
        match_probability_percent: cli.search.match_probability_percent,
        indel_probability_percent: cli.search.indel_probability_percent,
        data_file: cli.output.data_file,
        csv_file: cli.output.csv_file,
        masked_file: cli.output.masked_file,
        flanking_sequence: cli.output.flanking_sequence,
        html_output_off: !cli.output.html_output,
        redundancy_removal_off: cli.output.no_redundancy,
        ngs: cli.output.ngs,
        verbose: cli.output.verbose,
        max_wrap_length: (cli.max_length_millions as f64 * 1e6) as u32,
        output_directory: cli.output.directory.to_string_lossy().into_owned(),
        ..ParamSet::default()
    };

    paramset.input_file_name = cli.file.to_string_lossy().into_owned();
    paramset.output_prefix = name_part(&cli.file);

    if paramset.ngs {
        /* the compact records go to stdout; there is nothing else to show */
        paramset.data_file = true;
    } else {
        if cli.output.verbose {
            print_banner();
        }
        /* the report files are opt-in; a run without any of them has no
        visible effect besides the exit code */
        let writes_nothing = !cli.output.html_output
            && !cli.output.data_file
            && !cli.output.csv_file
            && !cli.output.masked_file;
        if writes_nothing {
            eprintln!(
                "note: no output files requested; use --html, -d, --csv or -m to write results"
            );
        }
    }

    if cli.file == Path::new("-") {
        if paramset.ngs {
            paramset.use_stdin = true;
            paramset.html_output_off = true;
        } else {
            eprintln!(
                "error: reading from standard input requires --ngs\n\nFor more information, try '--help'."
            );
            return 2;
        }
    }

    /* call the function that controls execution */
    control::trf_control_routine(&mut paramset);

    if let Some(status) = paramset.end_status.clone() {
        if !paramset.stdout_closed {
            eprintln!("Error processing input: {}", status);
        }
        return 1;
    }
    if paramset.output_count == 0 && !paramset.stdout_closed {
        eprintln!("No TRs found.");
    }
    0
}
