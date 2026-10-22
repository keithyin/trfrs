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
    pub name: String,
    /// Sequence characters, 0-based, upper case.
    pub sequence: Vec<u8>,
}

/// Error returned by [`load_sequence`].
pub enum LoadError {
    /// The file did not start with `>`, or was empty.
    Silent,
    /// An error message that has to be printed through `PrintError()`.
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
pub fn load_sequence(rd: &mut ByteReader) -> Result<(FastaSequence, bool), LoadError> {
    let mut seq = FastaSequence::default();

    // read the FASTA '>' symbol
    let mut c = rd.getc();
    if c < 0 || (c as u8) != b'>' {
        return Err(LoadError::Silent);
    }

    // read name and description text
    let mut name = Vec::with_capacity(MAX_SEQ_NAME_LEN);
    while name.len() < MAX_SEQ_NAME_LEN - 1 {
        c = rd.getc();
        if c == LF || c == CR {
            break;
        } else if c == EOF {
            return Err(LoadError::Message(
                "FASTA input terminated too early".into(),
            ));
        } else {
            name.push(c as u8);
        }
    }
    // if the line was not read completely, flush the rest
    if name.len() == MAX_SEQ_NAME_LEN - 1 {
        c = 0;
        while c != CR && c != LF && c != EOF {
            c = rd.getc();
        }
    }
    seq.name = String::from_utf8_lossy(&name).into_owned();

    // read sequence characters
    let mut next = -1;
    let mut sequence: Vec<u8> = Vec::new();
    let mut composition = [0i32; 26];
    while c != EOF && (c as u8) != b'>' {
        c = rd.getc();
        let b = c as u8;
        if (b'A'..=b'Z').contains(&b) {
            composition[(b - b'A') as usize] += 1;
            sequence.push(b);
        } else if b.is_ascii_lowercase() {
            composition[(b - b'a') as usize] += 1;
            sequence.push(b.to_ascii_uppercase());
        } else if b == b'>' {
            next = 1;
            rd.ungetc(b'>');
        } else if c == EOF {
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
