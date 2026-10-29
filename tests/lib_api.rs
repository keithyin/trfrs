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

//! Tests of the library interface.
//!
//! The golden files next to the FASTA fixtures are the `-ngs` data output of
//! the reference program: one `@name` line per sequence followed by one line
//! per repeat, which contains exactly the values the library reports.
//! The library therefore has to reproduce them field by field.

use std::path::{Path, PathBuf};
use std::process::Command;

use trfrs::{Error, Params, Repeat, TrfFinder};

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn read_fixture(name: &str) -> String {
    std::fs::read_to_string(data_dir().join(name)).expect("read fixture")
}

/// Minimal FASTA reader: returns `(header, upper case sequence)` pairs.
fn parse_fasta(text: &str) -> Vec<(String, String)> {
    let mut records: Vec<(String, String)> = Vec::new();
    let mut current: Option<(String, String)> = None;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            current = Some((header.trim().to_string(), String::new()));
        } else if let Some((_, sequence)) = current.as_mut() {
            sequence.push_str(line.trim());
        }
    }
    if let Some(record) = current.take() {
        records.push(record);
    }
    records
}

/// One row of a golden `-ngs` data file.
struct GoldenRecord {
    first: usize,
    last: usize,
    period: usize,
    copies: f32,
    size: usize,
    percent_matches: u8,
    percent_indels: u8,
    score: i32,
    percent_a: u8,
    percent_c: u8,
    percent_g: u8,
    percent_t: u8,
    entropy: f32,
    pattern: String,
    region: String,
}

/// One `@name` group of a golden `-ngs` data file.
struct GoldenSequence {
    name: String,
    records: Vec<GoldenRecord>,
}

fn parse_golden(text: &str) -> Vec<GoldenSequence> {
    let mut sequences: Vec<GoldenSequence> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(name) = line.strip_prefix('@') {
            sequences.push(GoldenSequence {
                name: name.trim().to_string(),
                records: Vec::new(),
            });
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        assert!(fields.len() >= 15, "unexpected golden line: {}", line);
        let sequence = sequences.last_mut().expect("record before '@' header");
        sequence.records.push(GoldenRecord {
            first: fields[0].parse().unwrap(),
            last: fields[1].parse().unwrap(),
            period: fields[2].parse().unwrap(),
            copies: fields[3].parse().unwrap(),
            size: fields[4].parse().unwrap(),
            percent_matches: fields[5].parse().unwrap(),
            percent_indels: fields[6].parse().unwrap(),
            score: fields[7].parse().unwrap(),
            percent_a: fields[8].parse().unwrap(),
            percent_c: fields[9].parse().unwrap(),
            percent_g: fields[10].parse().unwrap(),
            percent_t: fields[11].parse().unwrap(),
            entropy: fields[12].parse().unwrap(),
            pattern: fields[13].to_string(),
            region: fields[14].to_string(),
        });
    }
    sequences
}

fn assert_matches_golden(repeat: &Repeat, golden: &GoldenRecord, sequence: &[u8]) {
    assert_eq!(
        repeat.start,
        golden.first - 1,
        "start of {}",
        golden.pattern
    );
    assert_eq!(repeat.end, golden.last, "end of {}", golden.pattern);
    assert_eq!(repeat.len(), golden.last - golden.first + 1);
    assert_eq!(repeat.period, golden.period, "period of {}", golden.pattern);
    assert_eq!(repeat.consensus_size, golden.size);
    assert_eq!(repeat.score, golden.score);
    assert_eq!(repeat.percent_matches, golden.percent_matches);
    assert_eq!(repeat.percent_indels, golden.percent_indels);
    assert_eq!(repeat.percent_a, golden.percent_a);
    assert_eq!(repeat.percent_c, golden.percent_c);
    assert_eq!(repeat.percent_g, golden.percent_g);
    assert_eq!(repeat.percent_t, golden.percent_t);
    assert_eq!(repeat.pattern, golden.pattern);
    /* the golden files print with `%.1f` and `%.2f`; compare the rendered value
    so that the rounding of the reference output is reproduced exactly */
    assert_eq!(
        format!("{:.1}", repeat.copies),
        format!("{:.1}", golden.copies),
        "copies of {}",
        golden.pattern
    );
    assert_eq!(
        format!("{:.2}", repeat.entropy),
        format!("{:.2}", golden.entropy),
        "entropy of {}",
        golden.pattern
    );

    let region = repeat.sequence(sequence).expect("region inside sequence");
    assert_eq!(
        std::str::from_utf8(region).unwrap(),
        golden.region,
        "region of {}",
        golden.pattern
    );
}

/// Compares one FASTA fixture with its golden `-ngs` output.
fn check_fixture(fasta_file: &str, golden_file: &str) {
    let sequences = parse_fasta(&read_fixture(fasta_file));
    let golden = parse_golden(&read_fixture(golden_file));
    let finder = TrfFinder::new(&Params::default()).expect("valid parameters");

    for expected in &golden {
        let (_, sequence) = sequences
            .iter()
            .find(|(name, _)| name == &expected.name)
            .unwrap_or_else(|| panic!("{} not found in {}", expected.name, fasta_file));

        let repeats = finder
            .find(sequence)
            .unwrap_or_else(|e| panic!("{}: {}", expected.name, e));

        assert_eq!(
            repeats.len(),
            expected.records.len(),
            "number of repeats in {}",
            expected.name
        );
        for (repeat, golden_record) in repeats.iter().zip(&expected.records) {
            assert_matches_golden(repeat, golden_record, sequence.as_bytes());
        }
    }

    /* sequences without repeats must not be reported */
    for (name, sequence) in &sequences {
        if golden.iter().any(|g| &g.name == name) {
            continue;
        }
        assert!(
            finder.find(sequence).unwrap().is_empty(),
            "{} should have no repeats",
            name
        );
    }
}

#[test]
fn golden_single_sequence() {
    check_fixture("small_test.fasta", "small_test.fasta.out");
}

#[test]
fn golden_three_sequences() {
    check_fixture("diff_cn.fasta", "diff_cn.fasta.out");
}

#[test]
fn golden_multi_sequence() {
    check_fixture("multi.fasta", "multi.fasta.out");
}

#[test]
fn synthetic_repeat_region_and_pattern() {
    let finder = TrfFinder::new(&Params::default()).unwrap();
    let sequence = "TTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGGGGG";
    let repeats = finder.find(sequence).unwrap();
    assert_eq!(repeats.len(), 1, "repeats: {:?}", repeats);

    let repeat = &repeats[0];
    assert_eq!(repeat.period, 4);
    assert_eq!(repeat.consensus_size, 4);
    /* the repeat starts one base inside the leading T run, so the consensus
    pattern is the rotated "TACG" */
    assert_eq!(repeat.pattern, "TACG");
    assert_eq!(repeat.start, 5);
    assert_eq!(repeat.end, 46);
    assert_eq!(
        &sequence[repeat.start..repeat.end],
        format!("T{}", "ACGT".repeat(10))
    );
    assert_eq!(repeat.len(), 41);
    assert_eq!(repeat.percent_matches, 100);
}

#[test]
fn input_is_normalized_like_fasta() {
    let finder = TrfFinder::new(&Params::default()).unwrap();
    let sequence = "ACGT".repeat(12);
    let plain = finder.find(&sequence).unwrap();

    /* lower case, line breaks and whitespace are ignored */
    let wrapped_text = format!(
        "{}\n{}\r\n{}\n",
        "acgt".repeat(4),
        "acgt".repeat(4),
        "acgt".repeat(4)
    );
    let wrapped = finder.find(&wrapped_text).unwrap();
    assert_eq!(plain, wrapped);

    /* an owned Vec<u8> is accepted as well */
    let owned: Vec<u8> = sequence.as_bytes().to_vec();
    let bytes = finder.find(&owned).unwrap();
    assert_eq!(plain, bytes);
}

#[test]
fn finder_can_be_reused() {
    let finder = TrfFinder::new(&Params::default()).unwrap();
    let first = finder
        .find("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
        .unwrap();
    let second = finder
        .find("CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT")
        .unwrap();
    let first_again = finder
        .find("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
        .unwrap();

    assert_eq!(first, first_again);
    assert_eq!(first[0].pattern, "ACGT");
    assert_eq!(second[0].pattern, "CAGT");
    assert_eq!(finder.params(), &Params::default());
}

#[test]
fn invalid_parameters_are_rejected() {
    let too_large = i32::MAX as u32 + 1;

    let cases: [(&str, Params); 11] = [
        (
            "match_weight",
            Params {
                match_weight: 0,
                ..Params::default()
            },
        ),
        (
            "mismatch",
            Params {
                mismatch: 0,
                ..Params::default()
            },
        ),
        (
            "indel",
            Params {
                indel: 0,
                ..Params::default()
            },
        ),
        (
            "minscore",
            Params {
                minscore: 0,
                ..Params::default()
            },
        ),
        (
            "maxperiod",
            Params {
                maxperiod: 2001,
                ..Params::default()
            },
        ),
        (
            "match_probability_percent",
            Params {
                match_probability_percent: 70,
                ..Params::default()
            },
        ),
        /* weights and the score threshold are used as 32-bit signed integers
        internally, so larger values must be rejected instead of wrapping */
        (
            "match_weight",
            Params {
                match_weight: u32::MAX,
                ..Params::default()
            },
        ),
        (
            "mismatch",
            Params {
                mismatch: too_large,
                ..Params::default()
            },
        ),
        (
            "minscore",
            Params {
                minscore: too_large,
                ..Params::default()
            },
        ),
        (
            "min_period",
            Params {
                min_period: 0,
                ..Params::default()
            },
        ),
        (
            "min_period",
            Params {
                min_period: 16,
                maxperiod: 15,
                ..Params::default()
            },
        ),
    ];

    for (name, params) in cases {
        match TrfFinder::new(&params) {
            Err(Error::InvalidParam { name: reported, .. }) => assert_eq!(reported, name),
            other => panic!("{} should be rejected, got {:?}", name, other),
        }
    }

    let maxperiod_zero = Params {
        maxperiod: 0,
        ..Params::default()
    };
    assert!(TrfFinder::new(&maxperiod_zero).is_err());

    let pi_too_large = Params {
        indel_probability_percent: 101,
        ..Params::default()
    };
    assert!(TrfFinder::new(&pi_too_large).is_err());

    /* the largest value that still fits into an i32 is accepted */
    let boundary = Params {
        match_weight: i32::MAX as u32,
        mismatch: i32::MAX as u32,
        indel: i32::MAX as u32,
        minscore: i32::MAX as u32,
        ..Params::default()
    };
    assert!(TrfFinder::new(&boundary).is_ok());
}

#[test]
fn very_short_sequences_are_handled() {
    let finder = TrfFinder::new(&Params::default()).unwrap();

    for sequence in ["", "   ", "\n-12.3"] {
        assert_eq!(
            finder.find(sequence).unwrap_err(),
            Error::EmptySequence,
            "{:?} should be empty",
            sequence
        );
    }

    /* shorter than the smallest tuple size (4), nothing can be reported */
    for sequence in ["A", "AC", "ACG", "ACGT", "NNNN", "acgt"] {
        let repeats = finder
            .find(sequence)
            .unwrap_or_else(|e| panic!("{:?}: {}", sequence, e));
        assert!(repeats.is_empty(), "{:?} produced {:?}", sequence, repeats);
    }

    /* a short homopolymer is reported once the score threshold allows it */
    let permissive = TrfFinder::new(&Params {
        minscore: 2,
        maxperiod: 2,
        ..Params::default()
    })
    .unwrap();
    let repeats = permissive.find("AAAAAA").unwrap();
    assert_eq!(repeats.len(), 1, "repeats: {:?}", repeats);
    assert_eq!(repeats[0].pattern, "A");
    assert_eq!(repeats[0].period, 1);
    assert!(repeats[0].end > repeats[0].start);
    assert_eq!(
        repeats[0].sequence(b"AAAAAA").unwrap().len(),
        repeats[0].len()
    );
}

#[test]
fn empty_sequence_is_rejected() {
    let finder = TrfFinder::new(&Params::default()).unwrap();
    assert_eq!(finder.find("").unwrap_err(), Error::EmptySequence);
    assert_eq!(
        finder.find(" \n\r\t-12.3").unwrap_err(),
        Error::EmptySequence
    );
}

/// The STR preset only narrows the period range; it must stay equal to the
/// explicit parameters, and be accepted by the finder.
#[test]
fn str_preset_values() {
    let preset = Params::str_preset();
    assert_eq!(preset.min_period, 1);
    assert_eq!(preset.maxperiod, 15);
    assert_eq!(preset.minscore, 50);
    assert_eq!(
        preset,
        Params {
            min_period: 1,
            maxperiod: 15,
            ..Params::default()
        }
    );
    assert!(TrfFinder::new(&preset).is_ok());
}

/// The default scan finds period 16 in the fixture, the STR preset does not,
/// and the records it keeps are exactly the ones the command line reports.
#[test]
fn str_preset_matches_command_line() {
    let (_, sequence) = parse_fasta(&read_fixture("str_short.fasta"))
        .into_iter()
        .next()
        .expect("fixture sequence");

    let default = TrfFinder::new(&Params::default())
        .unwrap()
        .find(&sequence)
        .unwrap();
    assert!(
        default.iter().any(|repeat| repeat.period == 16),
        "the fixture should contain a period 16 repeat"
    );

    let repeats = TrfFinder::new(&Params::str_preset())
        .unwrap()
        .find(&sequence)
        .unwrap();
    assert!(!repeats.is_empty(), "fixture should contain STRs");
    assert!(repeats
        .iter()
        .all(|repeat| (1..=15).contains(&repeat.period)));

    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .arg(data_dir().join("str_short.fasta"))
        .args(["--str", "-ngs"])
        .output()
        .expect("run trf");
    assert!(output.status.success(), "trf failed: {:?}", output.status);
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let rows: Vec<&str> = stdout
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('@'))
        .collect();

    assert_eq!(rows.len(), repeats.len(), "number of records");
    for (row, repeat) in rows.iter().zip(&repeats) {
        let fields: Vec<&str> = row.split_whitespace().collect();
        assert_eq!(fields[0].parse::<usize>().unwrap(), repeat.start + 1);
        assert_eq!(fields[1].parse::<usize>().unwrap(), repeat.end);
        assert_eq!(fields[2].parse::<usize>().unwrap(), repeat.period);
        assert_eq!(fields[7].parse::<i32>().unwrap(), repeat.score);
        assert_eq!(fields[13], repeat.pattern);
    }
}

/// A lower period bound removes single base repeats, while the default keeps
/// them (the fixture starts with a poly-A run).
#[test]
fn min_period_filters_lower_periods() {
    let (_, sequence) = parse_fasta(&read_fixture("str_short.fasta"))
        .into_iter()
        .next()
        .expect("fixture sequence");

    let default = TrfFinder::new(&Params::default())
        .unwrap()
        .find(&sequence)
        .unwrap();
    assert!(default.iter().any(|repeat| repeat.period == 1));

    let from_two = TrfFinder::new(&Params {
        min_period: 2,
        maxperiod: 15,
        ..Params::default()
    })
    .unwrap()
    .find(&sequence)
    .unwrap();
    assert!(!from_two.is_empty());
    assert!(from_two.iter().all(|repeat| repeat.period >= 2));
}

#[test]
fn redundancy_removal_can_be_disabled() {
    let sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let with_removal = TrfFinder::new(&Params::default()).unwrap();
    let without_removal = TrfFinder::new(&Params {
        remove_redundancy: false,
        ..Params::default()
    })
    .unwrap();

    let kept = with_removal.find(sequence).unwrap();
    let all = without_removal.find(sequence).unwrap();
    assert!(kept.len() <= all.len());
    assert!(!all.is_empty());
}
