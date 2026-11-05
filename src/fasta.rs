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

//! Sequence input handling, mirroring `LoadSequenceFromFileEugene()`.
//!
//! Besides FASTA, FASTQ files are accepted as an extension of this port.  The
//! format is detected from the first byte of every record (`>` for FASTA, `@`
//! for FASTQ), so a single run can even mix both.  Only the sequence matters:
//! the quality strings of a FASTQ record are skipped.

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

    /// Reads one line, without its `\n` (or `\r\n`) terminator.  `None` means
    /// end of file; a last line that is not terminated is still returned.
    /// Used by the FASTQ reader, where line structure carries meaning.
    pub fn read_line(&mut self) -> Option<Vec<u8>> {
        let mut line = Vec::new();
        loop {
            let c = self.getc();
            if c < 0 {
                return if line.is_empty() { None } else { Some(line) };
            }
            if c == LF {
                if line.last() == Some(&(CR as u8)) {
                    line.pop();
                }
                return Some(line);
            }
            line.push(c as u8);
        }
    }
}

const LF: i32 = 10;
const CR: i32 = 13;
const EOF: i32 = -1;

/// Reads one sequence from the stream.
///
/// FASTA (`>` header) and FASTQ (`@` header) records are both accepted, as
/// detected from the first byte of the record: FASTQ is an extension of this
/// port, the reference implementation only reads FASTA.
///
/// The boolean in the returned tuple is `true` when another sequence header was
/// encountered (the reference implementation returns 1 in that case, 0 at end of
/// file).
pub fn load_sequence(reader: &mut ByteReader) -> Result<(FastaSequence, bool), LoadError> {
    // read the FASTA '>' or FASTQ '@' symbol
    let char_code = reader.getc();
    match char_code {
        EOF => Err(LoadError::Silent),
        _ if (char_code as u8) == b'>' => load_fasta_sequence(reader),
        _ if (char_code as u8) == b'@' => load_fastq_sequence(reader),
        _ => Err(LoadError::Silent),
    }
}

/// Reads one FASTA record; the leading `>` has already been consumed.
fn load_fasta_sequence(reader: &mut ByteReader) -> Result<(FastaSequence, bool), LoadError> {
    let mut seq = FastaSequence::default();

    // read name and description text
    let mut char_code = reader.getc();
    let mut name = Vec::with_capacity(MAX_SEQ_NAME_LEN);
    while name.len() < MAX_SEQ_NAME_LEN - 1 && char_code != LF && char_code != CR {
        if char_code == EOF {
            return Err(LoadError::Message(
                "FASTA input terminated too early".into(),
            ));
        }
        name.push(char_code as u8);
        char_code = reader.getc();
    }
    // if the line was not read completely, flush the rest
    if name.len() == MAX_SEQ_NAME_LEN - 1 {
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

/// The error reported when a FASTQ record ends before its quality strings.
const FASTQ_TRUNCATED: &str = "FASTQ input terminated too early";

/// Reads one FASTQ record; the leading `@` has already been consumed.
///
/// The layout is the usual four part record
///
/// ```text
/// @name description
/// SEQUENCE...
/// +
/// QUALITY...
/// ```
///
/// where sequence and quality may be wrapped over several lines.  Only the
/// sequence is used: the quality strings are counted (one per sequence
/// character) and skipped.  Normalisation is the one of
/// [`FastaSequence::from_bytes`], i.e. the same rules as for FASTA.
fn load_fastq_sequence(reader: &mut ByteReader) -> Result<(FastaSequence, bool), LoadError> {
    // header line: the text after '@', as in a FASTA header
    let name = reader
        .read_line()
        .ok_or_else(|| LoadError::Message(FASTQ_TRUNCATED.into()))?;

    // sequence lines, up to the '+' separator line
    let mut raw: Vec<u8> = Vec::new();
    loop {
        let line = reader
            .read_line()
            .ok_or_else(|| LoadError::Message(FASTQ_TRUNCATED.into()))?;
        if line.first() == Some(&b'+') {
            break;
        }
        raw.extend_from_slice(&line);
        raw.push(b'\n');
    }

    // one quality character per sequence character; wrapped lines are summed
    let length = raw.iter().filter(|b| b.is_ascii_alphabetic()).count();
    let mut remaining = length;
    while remaining > 0 {
        let line = reader
            .read_line()
            .ok_or_else(|| LoadError::Message(FASTQ_TRUNCATED.into()))?;
        remaining = remaining.saturating_sub(line.len());
    }

    // is another record following?  Blank lines between records are tolerated.
    let more = loop {
        match reader.getc() {
            EOF => break false,
            LF | CR => continue,
            c if (c as u8) == b'@' => {
                reader.ungetc(b'@');
                break true;
            }
            _ => return Err(LoadError::Silent),
        }
    };

    Ok((FastaSequence::from_bytes(&name, &raw), more))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn reader(input: &[u8]) -> ByteReader {
        ByteReader::new(Box::new(Cursor::new(input.to_vec())))
    }

    /// Four line record, CRLF line endings, lower case bases.
    #[test]
    fn fastq_record_is_normalized() {
        let mut reader = reader(b"@name with description\r\nacgtACGT\r\n+\r\nIIIIIIII\r\n");
        let (seq, more) = load_sequence(&mut reader).unwrap();
        assert!(!more);
        assert_eq!(seq.name, b"name with description");
        assert_eq!(seq.sequence, b"ACGTACGT");
        assert_eq!(seq.length, 8);
        assert_eq!(seq.nucleotides, 8);
    }

    /// Wrapped sequence and quality lines, and a second record.
    #[test]
    fn fastq_wrapped_lines_and_multiple_records() {
        let mut reader = reader(b"@a\nACGTACGT\nACGT\n+\nIIIIII\nIIIIII\n@b\nTTTT\n+\n!!!!\n");
        let (first, more) = load_sequence(&mut reader).unwrap();
        assert!(more);
        assert_eq!(first.name, b"a");
        assert_eq!(first.sequence, b"ACGTACGTACGT");

        let (second, more) = load_sequence(&mut reader).unwrap();
        assert!(!more);
        assert_eq!(second.name, b"b");
        assert_eq!(second.sequence, b"TTTT");
    }

    /// `@` and `+` are ordinary quality characters, so only the length of the
    /// quality strings can be used to find the next record.
    #[test]
    fn fastq_quality_may_look_like_a_header() {
        let input = b"@q\nACGT\n+\n@+@+\n@r\nTTTT\n+\nIIII\n";
        let mut reader = reader(input);
        let (first, more) = load_sequence(&mut reader).unwrap();
        assert!(more);
        assert_eq!(first.sequence, b"ACGT");
        let (second, more) = load_sequence(&mut reader).unwrap();
        assert!(!more);
        assert_eq!(second.sequence, b"TTTT");
    }

    /// The last quality line may be unterminated.
    #[test]
    fn fastq_last_quality_line_without_newline() {
        let (seq, more) = load_sequence(&mut reader(b"@q\nACGT\n+\nIIII")).unwrap();
        assert!(!more);
        assert_eq!(seq.sequence, b"ACGT");
    }

    /// A record that ends before its quality strings is reported.
    #[test]
    fn truncated_fastq_is_rejected() {
        let err = load_sequence(&mut reader(b"@q\nACGT\n")).unwrap_err();
        assert!(matches!(err, LoadError::Message(_)));
    }

    /// FASTA remains the default and is still detected by its `>` header.
    #[test]
    fn fasta_records_still_load() {
        let mut reader = reader(b">n\nACGT\n>m\nTTTT\n");
        let (first, more) = load_sequence(&mut reader).unwrap();
        assert!(more);
        assert_eq!(first.name, b"n");
        assert_eq!(first.sequence, b"ACGT");
        let (second, more) = load_sequence(&mut reader).unwrap();
        assert!(!more);
        assert_eq!(second.sequence, b"TTTT");
    }
}
