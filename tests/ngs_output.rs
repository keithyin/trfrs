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

const PARAMS: [&str; 7] = ["2", "5", "7", "80", "10", "50", "2000"];

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
    let output = Command::new(env!("CARGO_BIN_EXE_trf"))
        .current_dir(&dir)
        .arg(fasta)
        .args(PARAMS)
        .args(args)
        .output()
        .expect("run trf");
    (dir, output)
}

fn golden(name: &str) -> String {
    std::fs::read_to_string(data_dir().join(name)).expect("read golden output")
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
