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

/// Runs the binary with exactly the given arguments, for options that supply
/// their own parameter defaults (`--str`).
fn run_fixture_args(fasta: &str, args: &[&str], tag: &str) -> (PathBuf, Output) {
    let dir = temp_dir(tag);
    std::fs::copy(data_dir().join(fasta), dir.join(fasta)).expect("copy fixture");
    let output = Command::new(env!("CARGO_BIN_EXE_trfrs"))
        .current_dir(&dir)
        .arg(fasta)
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

/// Version of the reference C program whose output the golden files hold.
const REFERENCE_VERSION: &str = "4.10.0";

/// The `.dat` and `.summary.html` headers carry the program version.  The
/// golden files were produced by the reference program, while trfrs reports
/// its own version from `Cargo.toml`; rewrite that line so the comparison
/// covers the report and not the version number.
fn golden_with_program_version(name: &str) -> String {
    golden(name).replace(
        &format!("Version {}", REFERENCE_VERSION),
        &format!("Version {}", env!("CARGO_PKG_VERSION")),
    )
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
    assert_eq!(
        dat,
        golden_with_program_version("multi.fasta.dat"),
        ".dat output differs"
    );

    let summary = std::fs::read_to_string(dir.join("multi.fasta.2.5.7.80.10.50.2000.summary.html"))
        .expect("read summary output");
    assert_eq!(
        summary,
        golden_with_program_version("multi.fasta.summary.html"),
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

/// Period column (third field) of every record of an `-ngs` output.
fn periods(stdout: &[u8]) -> Vec<usize> {
    String::from_utf8(stdout.to_vec())
        .expect("stdout is utf8")
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('@'))
        .map(|line| {
            line.split_whitespace()
                .nth(2)
                .expect("period field")
                .parse()
                .expect("numeric period")
        })
        .collect()
}

/// `--str` is a shorthand for the period range used by STR detection: it must
/// produce exactly the same bytes as the equivalent explicit parameters (and as
/// `--maxperiod 15`, because the lower bound defaults to 1).
#[test]
fn str_preset_equals_explicit_period_range() {
    let (preset_dir, preset) = run_fixture_args("str_short.fasta", &["--str", "-ngs"], "strpreset");
    let (explicit_dir, explicit) = run_fixture_args(
        "str_short.fasta",
        &["--minperiod", "1", "--maxperiod", "15", "-ngs"],
        "strexplicit",
    );
    let (maxonly_dir, maxonly) = run_fixture_args(
        "str_short.fasta",
        &["--maxperiod", "15", "-ngs"],
        "strmaxonly",
    );
    let (dash_dir, dash) = run_fixture_args("str_short.fasta", &["-str", "-ngs"], "strdash");

    assert!(preset.status.success(), "trf failed: {:?}", preset.status);
    assert_eq!(
        preset.stdout, explicit.stdout,
        "--str and explicit range differ"
    );
    assert_eq!(
        preset.stdout, maxonly.stdout,
        "--str and --maxperiod 15 differ"
    );
    assert_eq!(preset.stdout, dash.stdout, "-str and --str differ");

    for dir in [preset_dir, explicit_dir, maxonly_dir, dash_dir] {
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// The fixture contains period 1, 2, 3, 5, 15 and 16 repeats; STR mode keeps
/// the 1–15 ones (single base repeats included) and drops period 16.
#[test]
fn str_preset_reports_only_short_units() {
    let (all_dir, all) = run_fixture_args("str_short.fasta", &["-ngs"], "strall");
    let (str_dir, str_mode) = run_fixture_args("str_short.fasta", &["--str", "-ngs"], "stronly");

    let mut all_periods = periods(&all.stdout);
    all_periods.sort_unstable();
    assert_eq!(all_periods, vec![1, 2, 3, 5, 15, 16]);

    let mut str_periods = periods(&str_mode.stdout);
    str_periods.sort_unstable();
    assert_eq!(str_periods, vec![1, 2, 3, 5, 15]);
    assert!(str_periods.iter().all(|period| (1..=15).contains(period)));

    let _ = std::fs::remove_dir_all(&all_dir);
    let _ = std::fs::remove_dir_all(&str_dir);
}

/// Explicit options win over the preset, and `--minperiod` removes the lower
/// bound of the reported period range.
#[test]
fn period_range_options_take_precedence() {
    let (expanded_dir, expanded) = run_fixture_args(
        "str_short.fasta",
        &["--str", "--maxperiod", "2000", "-ngs"],
        "strovr",
    );
    let mut expanded_periods = periods(&expanded.stdout);
    expanded_periods.sort_unstable();
    assert_eq!(expanded_periods, vec![1, 2, 3, 5, 15, 16]);

    let (lower_dir, lower) = run_fixture_args(
        "str_short.fasta",
        &["--minperiod", "2", "--maxperiod", "15", "-ngs"],
        "strlower",
    );
    let mut lower_periods = periods(&lower.stdout);
    lower_periods.sort_unstable();
    assert!(!lower_periods.contains(&1), "period 1 should be filtered");
    assert!(lower_periods.contains(&2), "period 2 should be kept");

    let _ = std::fs::remove_dir_all(&expanded_dir);
    let _ = std::fs::remove_dir_all(&lower_dir);
}

/// The default output is unchanged: `--minperiod 1` is the default and must
/// reproduce the reference output byte for byte.
#[test]
fn default_minperiod_keeps_output() {
    let (dir, output) = run_fixture(
        "small_test.fasta",
        &["-ngs", "--minperiod", "1"],
        "mindefault",
    );
    assert!(output.status.success(), "trf failed: {:?}", output.status);
    assert_eq!(
        String::from_utf8(output.stdout).expect("stdout is utf8"),
        golden("small_test.fasta.out"),
        "default minperiod changed the output"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// A zero lower bound, a lower bound above the upper bound and a non numeric
/// value are rejected with the usual driver message and exit code 1.
#[test]
fn invalid_minperiod_is_rejected() {
    let cases: [(&[&str], &str); 3] = [
        (&["--minperiod", "0", "-ngs"], "positive integer"),
        (
            &["--minperiod", "16", "--maxperiod", "15", "-ngs"],
            "must not exceed MaxPeriod",
        ),
        (&["--minperiod", "abc", "-ngs"], "positive integer"),
    ];

    for (index, (args, expected)) in cases.iter().enumerate() {
        let (dir, output) = run_fixture_args("str_short.fasta", args, &format!("badmin{index}"));
        assert_eq!(output.status.code(), Some(1), "args: {args:?}");
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            stdout.contains(expected),
            "missing {expected:?} for {args:?}: {stdout}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
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
