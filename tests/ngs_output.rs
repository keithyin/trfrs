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

//! Regression tests against the expected output of the reference C program.
//!
//! The `.out`, `.dat` and `.summary.html` files come from the upstream
//! implementation (for `multi.fasta` they were produced with the reference
//! binary, version 4.10.0).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The search parameters as key-value options; every test below exercises
/// this spelling.
const PARAMS: [&str; 14] = [
    "--match",
    "2",
    "--mismatch",
    "5",
    "--indel",
    "7",
    "--pm",
    "80",
    "--pi",
    "10",
    "--minscore",
    "50",
    "--maxperiod",
    "2000",
];

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("trf-rs-test-{}-{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// Copies `fasta` into a fresh directory and runs the binary there.
fn run_fixture(fasta: &str, args: &[&str], tag: &str) -> (PathBuf, Output) {
    let dir = temp_dir(tag);
    std::fs::copy(data_dir().join(fasta), dir.join(fasta)).expect("copy fixture");
    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .current_dir(&dir)
        .arg(fasta)
        .args(PARAMS)
        .args(args)
        .output()
        .expect("run trf");
    (dir, output)
}

/// Reads a golden output file with LF line endings.
///
/// `.gitattributes` pins `tests/data` to `eol=lf`, but a working tree created
/// before that file existed (or with `core.autocrlf=true`) may still hold CRLF;
/// normalise so the comparison is about content, not checkout settings.  The
/// program itself always writes `\n`, so the produced output is not normalised.
fn golden(name: &str) -> String {
    std::fs::read_to_string(data_dir().join(name))
        .expect("read golden output")
        .replace("\r\n", "\n")
}

fn assert_ngs(fasta: &str, golden_file: &str, tag: &str) {
    let (dir, output) = run_fixture(fasta, &["-ngs"], tag);
    assert!(output.status.success(), "trf failed: {:?}", output.status);
    let actual = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert_eq!(
        actual,
        golden(golden_file),
        "ngs output differs for {}",
        fasta
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn small_test_ngs() {
    assert_ngs("small_test.fasta", "small_test.fasta.out", "small");
}

#[test]
fn diff_cn_ngs() {
    assert_ngs("diff_cn.fasta", "diff_cn.fasta.out", "diffcn");
}

#[test]
fn multi_ngs() {
    assert_ngs("multi.fasta", "multi.fasta.out", "multings");
}

/// A command line without a file name prints the usage text instead of a clap
/// diagnostic.  A lone `-h` only means "suppress html output", so the run
/// reports the usage and exits with 1; `-u` is the help option and exits 0.
#[test]
fn usage_without_file() {
    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .arg("-h")
        .output()
        .expect("run trf");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Please use:"),
        "usage text missing: {stderr}"
    );
    assert!(
        stderr.contains("--maxperiod"),
        "parameter options missing: {stderr}"
    );

    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .arg("-u")
        .output()
        .expect("run trf");
    assert_eq!(output.status.code(), Some(0));
    assert!(String::from_utf8_lossy(&output.stdout).contains("Please use:"));
}

/// The search parameters default to `2 5 7 80 10 50 2000`, so a run that
/// only names the file produces the same output as an explicit one.
#[test]
fn default_parameters() {
    let dir = temp_dir("defaults");
    std::fs::copy(
        data_dir().join("small_test.fasta"),
        dir.join("small_test.fasta"),
    )
    .expect("copy fixture");
    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .current_dir(&dir)
        .args(["small_test.fasta", "-ngs"])
        .output()
        .expect("run trf");
    assert!(output.status.success(), "trf failed: {:?}", output.status);
    let actual = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert_eq!(
        actual,
        golden("small_test.fasta.out"),
        "default parameters differ"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn multi_dat_and_summary() {
    let (dir, output) = run_fixture("multi.fasta", &["-d"], "multidat");
    assert!(output.status.success(), "trf failed: {:?}", output.status);

    let dat = std::fs::read_to_string(dir.join("multi.fasta.2.5.7.80.10.50.2000.dat"))
        .expect("read .dat output");
    assert_eq!(dat, golden("multi.fasta.dat"), ".dat output differs");

    let summary = std::fs::read_to_string(dir.join("multi.fasta.2.5.7.80.10.50.2000.summary.html"))
        .expect("read summary output");
    assert_eq!(
        summary,
        golden("multi.fasta.summary.html"),
        "summary output differs"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// `getopt_long_only()` accepts clustered short options.
#[test]
fn clustered_options() {
    let (plain, _) = run_fixture("multi.fasta", &["-d", "-f"], "plain");
    let (cluster, _) = run_fixture("multi.fasta", &["-df"], "cluster");
    let dat = "multi.fasta.2.5.7.80.10.50.2000.dat";
    let left = std::fs::read(plain.join(dat)).expect("plain .dat");
    let right = std::fs::read(cluster.join(dat)).expect("clustered .dat");
    assert_eq!(left, right, "-d -f and -df must produce the same output");
    let _ = std::fs::remove_dir_all(&plain);
    let _ = std::fs::remove_dir_all(&cluster);
}

/// A missing value makes getopt ignore the option; the run still succeeds.
#[test]
fn missing_maxlength_argument() {
    let (dir, output) = run_fixture("diff_cn.fasta", &["-ngs", "-l"], "missingarg");
    assert!(output.status.success(), "trf should still succeed");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("option requires an argument -- 'l'"),
        "unexpected stderr: {stderr}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `true` when `haystack` contains `needle` (`strstr()` on bytes).
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || haystack.windows(needle.len()).any(|w| w == needle)
}

/// The FASTA header is copied verbatim into every output file.  Header bytes
/// that are not valid UTF-8 must reach the output unchanged: the reference
/// program is byte oriented and `String::from_utf8_lossy` would replace them
/// with U+FFFD.
#[test]
fn non_utf8_header_is_preserved() {
    let dir = temp_dir("badutf8");
    let mut fasta: Vec<u8> = b">bad\xff\xfe name".to_vec();
    fasta.extend_from_slice(b"\nACCCCTCAGGGTACCCCTCAGGGTACCCCTCAGGGTACCCCTCAGGGT\n");
    std::fs::write(dir.join("bad.fasta"), &fasta).expect("write fasta");

    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .current_dir(&dir)
        .arg("bad.fasta")
        .args(PARAMS)
        .args(["-d", "-m"])
        .output()
        .expect("run trf");
    assert!(output.status.success(), "trf failed: {:?}", output.status);

    let stem = "bad.fasta.2.5.7.80.10.50.2000";
    let sequence_line: &[u8] = b"Sequence: bad\xff\xfe name";

    let dat = std::fs::read(dir.join(format!("{stem}.dat"))).expect("read .dat");
    assert!(contains(&dat, sequence_line), ".dat lost the header bytes");

    let html = std::fs::read(dir.join(format!("{stem}.1.html"))).expect("read .html");
    assert!(
        contains(&html, sequence_line),
        ".html lost the header bytes"
    );

    let txt = std::fs::read(dir.join(format!("{stem}.1.txt.html"))).expect("read .txt.html");
    assert!(
        contains(&txt, sequence_line),
        ".txt.html lost the header bytes"
    );

    let mask = std::fs::read(dir.join(format!("{stem}.mask"))).expect("read .mask");
    assert!(
        mask.starts_with(b">bad\xff\xfe name\n"),
        ".mask header was modified"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The reference program does not ignore `SIGPIPE`: a consumer that closes the
/// pipe early terminates it (`-ngs` output is meant to be piped).  Ignoring the
/// signal would turn the same situation into a panic and a bogus exit code.
#[test]
#[cfg(unix)]
fn broken_pipe_terminates_by_signal() {
    use std::os::unix::process::ExitStatusExt;
    use std::process::Stdio;

    let dir = temp_dir("sigpipe");

    /* enough tandem arrays that the `-ngs` output exceeds the pipe buffer */
    let mut state: u32 = 0x2545_f491;
    let mut sequence = String::new();
    for i in 0..800 {
        let period = 21 + i % 40;
        let mut unit = String::new();
        for _ in 0..period {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            unit.push(b"ACGT"[(state % 4) as usize] as char);
        }
        for _ in 0..3 {
            sequence.push_str(&unit);
        }
        sequence.push_str("ACGTACGTACGT");
    }
    std::fs::write(dir.join("big.fasta"), format!(">big\n{sequence}\n")).expect("write fasta");

    let mut child = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .current_dir(&dir)
        .arg("big.fasta")
        .args(PARAMS)
        .arg("-ngs")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn trf");

    /* close the reading end of the pipe */
    drop(child.stdout.take());

    let output = child.wait_with_output().expect("wait for trf");
    assert_eq!(
        output.status.signal(),
        Some(libc::SIGPIPE),
        "expected death by SIGPIPE, got {:?}",
        output.status
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("panicked"), "unexpected panic: {stderr}");

    let _ = std::fs::remove_dir_all(&dir);
}
