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
//! The crate provides two entry points:
//!
//! * the command line program (`trfrs`), which mirrors `trf.c`, and
//! * a library API for embedding the search into other Rust programs.
//!
//! # Library
//!
//! [`TrfFinder`] is created once from a set of [`Params`] and can then be used
//! for any number of sequences.  [`TrfFinder::find`] runs the search in memory:
//! it writes no files and prints nothing.
//!
//! ```
//! use trfrs::{Params, TrfFinder};
//!
//! fn main() -> Result<(), trfrs::Error> {
//!     let finder = TrfFinder::new(&Params::default())?;
//!     let repeats = finder.find("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")?;
//!     for repeat in &repeats {
//!         println!(
//!             "[{start}, {end}) period={period} pattern={pattern}",
//!             start = repeat.start,
//!             end = repeat.end,
//!             period = repeat.period,
//!             pattern = repeat.pattern,
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! Module map (reference C source -> Rust module):
//!
//! | C source                | Rust module                                   |
//! |-------------------------|-----------------------------------------------|
//! | `trf.c`                 | `cli.rs` (command line), `driver.rs` (driver) |
//! | `trfrun.h`              | `control.rs` (control routine, `TRF()`)       |
//! | `tr30dat.h` + `tr30dat.c` | `state.rs`, `algo.rs`, `scan.rs`            |
//! | `trfclean.h`            | `clean.rs` (redundancy removal, HTML output)  |
//! | FASTA input             | `fasta.rs`                                    |
//! | generated tables        | `tables.rs`                                   |
//!
//! The C program keeps its state in file scope globals; here everything lives in
//! `state::Trf`, which is created once per sequence.

mod algo;
mod clean;
mod cli;
mod control;
mod driver;
mod fasta;
mod params;
mod scan;
mod state;
mod tables;

pub mod api;

pub use api::{Error, Params, Repeat, TrfFinder};
pub use driver::run_cli;
pub use params::VERSION;
