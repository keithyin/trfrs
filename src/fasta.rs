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

//! FASTA input handling, mirroring `LoadSequenceFromFileEugene()`.

use std::io::{BufRead, Read};

use thiserror::Error;

use crate::params::MAX_SEQ_NAME_LEN;

/// Index of a base in the `ACGTcount[letter - 'A']` tables of the reference
/// implementation.
#[inline]
pub fn base_index(letter: u8) -> usize {
    (letter - b'A') as usize
}

/// One FASTA record.
#[derive(Clone, Debug, Default)]
pub struct FastaSequence {
    pub length: usize,
    /// Occurrences of each letter, indexed by `letter - 'A'`.
    pub composition: [i32; 26],
    /// Number of A, C, G and T characters.
    pub nucleotides: i32,
    /// Description text of the header line.  The bytes are kept as they were
    /// read: FASTA files are not guaranteed to be UTF-8 and the reference
    /// implementation copies the header verbatim into its output.
    pub name: Vec<u8>,
    /// Sequence characters, 0-based, upper case.
    pub sequence: Vec<u8>,
}

impl FastaSequence {
    /// Builds a record from a bare sequence, without a FASTA header line.
    ///
    /// The characters are normalized exactly like [`load_sequence`] does: only
    /// ASCII letters are kept (`A`-`Z`), lower case letters are converted to
    /// upper case, and the remaining bytes (line breaks, whitespace, `-`, `.`,
    /// digits, ...) are dropped.
    pub fn from_bytes(name: &[u8], sequence: &[u8]) -> FastaSequence {
        let mut composition = [0i32; 26];
        let mut normalized = Vec::with_capacity(sequence.len());
        for &byte in sequence {
            let upper = byte.to_ascii_uppercase();
            if upper.is_ascii_uppercase() {
                composition[(upper - b'A') as usize] += 1;
                normalized.push(upper);
            }
        }

        let mut name = name.to_vec();
        name.truncate(MAX_SEQ_NAME_LEN - 1);

        let nucleotides = composition[base_index(b'A')]
            + composition[base_index(b'C')]
            + composition[base_index(b'G')]
            + composition[base_index(b'T')];

        FastaSequence {
            length: normalized.len(),
            composition,
            nucleotides,
            name,
            sequence: normalized,
        }
    }
}

/// Error returned by [`load_sequence`].
#[derive(Debug, Error)]
pub enum LoadError {
    /// The file did not start with `>`, or was empty.
    #[error("empty file or bad format")]
    Silent,
    /// An error message that has to be printed through `PrintError()`.
    #[error("{0}")]
    Message(String),
}

/// A byte oriented reader with one byte of push back, mimicking `getc()`/`ungetc()`.
pub struct ByteReader {
    inner: Box<dyn BufRead>,
    pushback: Option<u8>,
}

impl ByteReader {
    pub fn new(inner: Box<dyn BufRead>) -> Self {
        ByteReader {
            inner,
            pushback: None,
        }
    }

    /// Returns the next byte, or `-1` at end of file (same as `getc`).
    pub fn getc(&mut self) -> i32 {
        if let Some(b) = self.pushback.take() {
            return b as i32;
        }
        let mut buf = [0u8; 1];
        match self.inner.read(&mut buf) {
            Ok(1) => buf[0] as i32,
            _ => -1,
        }
    }

    pub fn ungetc(&mut self, b: u8) {
        self.pushback = Some(b);
    }
}

const LF: i32 = 10;
const CR: i32 = 13;
const EOF: i32 = -1;

/// Reads one sequence from the stream.
///
/// The boolean in the returned tuple is `true` when another sequence header was
/// encountered (the reference implementation returns 1 in that case, 0 at end of
/// file).
pub fn load_sequence(reader: &mut ByteReader) -> Result<(FastaSequence, bool), LoadError> {
    let mut seq = FastaSequence::default();

    // read the FASTA '>' symbol
    let mut char_code = reader.getc();
    if char_code < 0 || (char_code as u8) != b'>' {
        return Err(LoadError::Silent);
    }

    // read name and description text
    let mut name = Vec::with_capacity(MAX_SEQ_NAME_LEN);
    while name.len() < MAX_SEQ_NAME_LEN - 1 {
        char_code = reader.getc();
        if char_code == LF || char_code == CR {
            break;
        } else if char_code == EOF {
            return Err(LoadError::Message(
                "FASTA input terminated too early".into(),
            ));
        } else {
            name.push(char_code as u8);
        }
    }
    // if the line was not read completely, flush the rest
    if name.len() == MAX_SEQ_NAME_LEN - 1 {
        char_code = 0;
        while char_code != CR && char_code != LF && char_code != EOF {
            char_code = reader.getc();
        }
    }
    seq.name = name;

    // read sequence characters
    let mut next = -1;
    let mut sequence: Vec<u8> = Vec::new();
    let mut composition = [0i32; 26];
    while char_code != EOF && (char_code as u8) != b'>' {
        char_code = reader.getc();
        let byte = char_code as u8;
        if (b'A'..=b'Z').contains(&byte) {
            composition[(byte - b'A') as usize] += 1;
            sequence.push(byte);
        } else if byte.is_ascii_lowercase() {
            composition[(byte - b'a') as usize] += 1;
            sequence.push(byte.to_ascii_uppercase());
        } else if byte == b'>' {
            next = 1;
            reader.ungetc(b'>');
        } else if char_code == EOF {
            next = 0;
        }
    }

    seq.length = sequence.len();
    seq.sequence = sequence;
    seq.composition = composition;
    seq.nucleotides = seq.composition[base_index(b'A')]
        + seq.composition[base_index(b'C')]
        + seq.composition[base_index(b'G')]
        + seq.composition[base_index(b'T')];

    Ok((seq, next == 1))
}
