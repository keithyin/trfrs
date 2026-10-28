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

//! Library interface: find tandem repeats in a bare nucleotide sequence.
//!
//! [`TrfFinder`] owns the search parameters and is cheap to reuse:
//!
//! ```
//! use trfrs::{Params, TrfFinder};
//!
//! fn main() -> Result<(), trfrs::Error> {
//!     let finder = TrfFinder::new(&Params::default())?;
//!     let sequence = "TTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAA";
//!     for repeat in finder.find(sequence)? {
//!         println!(
//!             "{}..{} period={} pattern={}",
//!             repeat.start, repeat.end, repeat.period, repeat.pattern
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! The search runs entirely in memory: no files are written and nothing is
//! printed.  Results are the same as the `.dat` records of the command line
//! program, in the same order.
//!
//! # Memory
//!
//! Like the C implementation the scanner keeps a scoring matrix whose size is
//! proportional to the sequence length (about 604 bytes per base, capped by the
//! `max_wrap_length` parameter of the reference implementation).  Sequences of
//! several hundred thousand bases therefore need hundreds of megabytes; callers
//! working with very large sequences should split them into chunks and accept
//! that repeats spanning a chunk boundary are truncated.
//!
//! # Process termination
//!
//! Every error the caller can trigger through [`Params`] or the input sequence
//! is reported as an [`Error`], never by terminating the process.
//!
//! A few paths that the reference C implementation treats as unrecoverable are
//! still reproduced: an internal algorithm invariant violation (for example the
//! alignment bandwidth checks in `algo.rs` or the statistics invariants in
//! `scan.rs`) and a failed large allocation call `exit()` instead of unwinding.
//! Those paths cannot be triggered by valid parameters and normal input, but a
//! caller embedding the library should treat them like an abort.

use std::fmt;

use crate::clean::RepeatRecord;
use crate::fasta::FastaSequence;
use crate::params::ParamSet;
use crate::state::Trf;

/// Search parameters.
///
/// The defaults are the customary human DNA values of the command line program:
/// `2 5 7 80 10 50 2000` (match, mismatch, indel, Pm, Pi, minscore, maxperiod).
///
/// Scaling all weights by the same factor leaves the consensus pattern and the
/// percentage statistics unchanged, but `score` and the `minscore` threshold are
/// expressed in weight units, so scaling them changes which repeats are
/// reported.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params {
    /// Matching weight (`--match`).  Must be at least 1 and fit in a 32-bit
    /// signed integer.
    pub match_weight: u32,
    /// Mismatching penalty (`--mismatch`).  Must be at least 1 and fit in a
    /// 32-bit signed integer.
    pub mismatch: u32,
    /// Insertion/deletion penalty (`--indel`).  Must be at least 1 and fit in a
    /// 32-bit signed integer.
    pub indel: u32,
    /// Matching probability in percent (`--pm`).  Only 75 and 80 are supported,
    /// the reference implementation has no statistical tables for other values.
    pub match_probability_percent: u32,
    /// Insertion/deletion probability in percent (`--pi`).  At most 100.
    pub indel_probability_percent: u32,
    /// Minimum alignment score a repeat must reach to be reported (`--minscore`).
    /// Must be at least 1 and fit in a 32-bit signed integer.
    pub minscore: u32,
    /// Maximum repeat unit length (`--maxperiod`), between 1 and 2000.
    pub maxperiod: u32,
    /// Remove repeats that are redundant with a higher scoring one.  The
    /// command line program does this unless `-r` is given.
    pub remove_redundancy: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            match_weight: 2,
            mismatch: 5,
            indel: 7,
            match_probability_percent: 80,
            indel_probability_percent: 10,
            minscore: 50,
            maxperiod: 2000,
            remove_redundancy: true,
        }
    }
}

/// Errors returned by [`TrfFinder::new`] and [`TrfFinder::find`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// A search parameter is out of range.
    InvalidParam {
        /// Name of the offending field of [`Params`].
        name: &'static str,
        /// Why the value was rejected.
        reason: String,
    },
    /// The sequence contains no letters, so there is nothing to scan.
    EmptySequence,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidParam { name, reason } => {
                write!(f, "invalid parameter `{}`: {}", name, reason)
            }
            Error::EmptySequence => f.write_str("the sequence is empty"),
        }
    }
}

impl std::error::Error for Error {}

/// One reported tandem repeat.
///
/// `start..end` is a 0-based half open interval (`start` inclusive, `end`
/// exclusive), so the repeated region is `&sequence[start..end]`.
///
/// The percentage fields are whole numbers between 0 and 100, computed exactly
/// like the columns of the command line `.dat` output.
#[derive(Clone, Debug, PartialEq)]
pub struct Repeat {
    /// First base of the repeat (0-based, inclusive).
    pub start: usize,
    /// One past the last base of the repeat (0-based, exclusive).
    pub end: usize,
    /// Length of the repeat unit (`period`).
    pub period: usize,
    /// Copy number: repeat length divided by the period.
    pub copies: f32,
    /// Length of the consensus pattern, which differs from `period` when the
    /// repeat contains indels.
    pub consensus_size: usize,
    /// Alignment score of the repeat.
    pub score: i32,
    /// Matching bases in percent (0-100).
    pub percent_matches: u8,
    /// Inserted/deleted bases in percent (0-100).
    pub percent_indels: u8,
    /// Percentage of `A` in the repeat.
    pub percent_a: u8,
    /// Percentage of `C` in the repeat.
    pub percent_c: u8,
    /// Percentage of `G` in the repeat.
    pub percent_g: u8,
    /// Percentage of `T` in the repeat.
    pub percent_t: u8,
    /// Entropy of the consensus pattern (0-2).
    pub entropy: f32,
    /// Consensus pattern of the repeat.
    pub pattern: String,
}

impl Repeat {
    /// Number of bases covered by the repeat.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Always `false`: a reported repeat covers at least one base.
    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    /// The repeated region itself, taken from the sequence that was scanned.
    ///
    /// Returns `None` when `sequence` is shorter than the repeat, e.g. when a
    /// [`Repeat`] is applied to a different sequence.
    pub fn sequence<'a>(&self, sequence: &'a [u8]) -> Option<&'a [u8]> {
        sequence.get(self.start..self.end)
    }

    fn from_record(record: &RepeatRecord) -> Repeat {
        Repeat {
            /* RepeatRecord stores 1-based, inclusive coordinates; the repeat
            spans `first..=last`, i.e. `first - 1 .. last` in 0-based terms. */
            start: (record.first - 1).max(0) as usize,
            end: record.last.max(0) as usize,
            period: record.period.max(0) as usize,
            copies: record.copies,
            consensus_size: record.size.max(0) as usize,
            score: record.score,
            percent_matches: clamp_percent(record.matches),
            percent_indels: clamp_percent(record.indels),
            percent_a: clamp_percent(record.a_count),
            percent_c: clamp_percent(record.c_count),
            percent_g: clamp_percent(record.g_count),
            percent_t: clamp_percent(record.t_count),
            entropy: record.entropy,
            pattern: record.pattern.clone(),
        }
    }
}

#[inline]
fn clamp_percent(value: i32) -> u8 {
    value.clamp(0, 100) as u8
}

/// Finds tandem repeats in nucleotide sequences.
///
/// The parameters are validated once by [`TrfFinder::new`]; the finder is plain
/// data, so it is `Send + Sync` and a single instance can be shared between
/// threads, each call to [`find`](TrfFinder::find) allocating its own scanner
/// state.
///
/// See the [module documentation](self) for an example.
#[derive(Clone, Debug)]
pub struct TrfFinder {
    params: Params,
    paramset: ParamSet,
}

impl TrfFinder {
    /// Validates `params` and creates a finder.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParam`] when a weight or penalty is zero, when a
    /// weight, penalty or `minscore` does not fit in a 32-bit signed integer,
    /// when `maxperiod` is outside `1..=2000`, when `indel_probability_percent`
    /// exceeds 100, or when `match_probability_percent` is neither 75 nor 80.
    pub fn new(params: &Params) -> Result<Self, Error> {
        validate(params)?;
        Ok(TrfFinder {
            params: params.clone(),
            paramset: silent_paramset(params),
        })
    }

    /// The parameters this finder was created with.
    pub fn params(&self) -> &Params {
        &self.params
    }

    /// Scans `sequence` and returns its tandem repeats in the order the
    /// command line program reports them.
    ///
    /// The sequence may be a `&str`, a `String`, a `Vec<u8>` or any other
    /// `AsRef<[u8]>`.  Characters are normalized like the FASTA reader does:
    /// only ASCII letters are kept and lower case is converted to upper case,
    /// so line breaks and whitespace inside the input are ignored.  Bases other
    /// than `A`, `C`, `G` and `T` (e.g. `N`) are scanned as mismatches, exactly
    /// like the command line program.
    ///
    /// # Errors
    ///
    /// Returns [`Error::EmptySequence`] if the normalized sequence contains no
    /// letters at all.
    pub fn find<S: AsRef<[u8]>>(&self, sequence: S) -> Result<Vec<Repeat>, Error> {
        let fasta = FastaSequence::from_bytes(b"query", sequence.as_ref());
        if fasta.length == 0 {
            return Err(Error::EmptySequence);
        }

        let mut trf = Trf::new(self.paramset.clone(), &fasta);
        trf.run();

        Ok(trf.records.iter().map(Repeat::from_record).collect())
    }
}

/// Checks the user facing parameters, mirroring the validation of the command
/// line driver.
fn validate(params: &Params) -> Result<(), Error> {
    let positive = |name: &'static str, value: u32| -> Result<(), Error> {
        if value == 0 {
            Err(Error::InvalidParam {
                name,
                reason: "must be a positive integer".to_string(),
            })
        } else {
            Ok(())
        }
    };

    /* The algorithm stores the weights, penalties and the score threshold as
    32-bit signed integers, so larger values would silently wrap to negative
    numbers instead of producing an error. */
    let fits_i32 = |name: &'static str, value: u32| -> Result<(), Error> {
        if value > i32::MAX as u32 {
            Err(Error::InvalidParam {
                name,
                reason: "must not exceed 2147483647 (used as a 32-bit signed integer)".to_string(),
            })
        } else {
            Ok(())
        }
    };

    positive("match_weight", params.match_weight)?;
    positive("mismatch", params.mismatch)?;
    positive("indel", params.indel)?;
    positive("minscore", params.minscore)?;

    fits_i32("match_weight", params.match_weight)?;
    fits_i32("mismatch", params.mismatch)?;
    fits_i32("indel", params.indel)?;
    fits_i32("minscore", params.minscore)?;

    if params.maxperiod == 0 || params.maxperiod > 2000 {
        return Err(Error::InvalidParam {
            name: "maxperiod",
            reason: "must be between 1 and 2000, inclusive".to_string(),
        });
    }

    if params.match_probability_percent != 75 && params.match_probability_percent != 80 {
        return Err(Error::InvalidParam {
            name: "match_probability_percent",
            reason: "only 75 and 80 are supported (no statistical tables exist for other values)"
                .to_string(),
        });
    }

    if params.indel_probability_percent > 100 {
        return Err(Error::InvalidParam {
            name: "indel_probability_percent",
            reason: "must be between 0 and 100, inclusive".to_string(),
        });
    }

    Ok(())
}

/// Builds the internal parameter set in "silent" mode: no HTML/alignment file,
/// no `.dat`/`.mask` output, no progress messages and no progress bar.
fn silent_paramset(params: &Params) -> ParamSet {
    ParamSet {
        match_weight: params.match_weight,
        mismatch: params.mismatch,
        indel: params.indel,
        minscore: params.minscore,
        maxperiod: params.maxperiod,
        match_probability_percent: params.match_probability_percent,
        indel_probability_percent: params.indel_probability_percent,
        data_file: false,
        masked_file: false,
        flanking_sequence: false,
        html_output_off: true,
        redundancy_removal_off: !params.remove_redundancy,
        ngs: true,
        use_stdin: false,
        percent: -1,
        ..ParamSet::default()
    }
}
