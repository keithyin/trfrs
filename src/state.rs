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

//! Core state of the tandem repeat scan (`tr30dat.h` in the reference
//! implementation).
//!
//! The C program keeps all of this in file scope globals.  Here the state lives
//! in a single [`Trf`] value that is created per FASTA sequence.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::fasta::FastaSequence;
use crate::params::ParamSet;

/// Maximum width of the alignment band, `MAX_BANDWIDTH`.
pub const MAX_BANDWIDTH: usize = 150;
/// Row stride of the `S` scoring matrix.
pub const SCORE_MATRIX_STRIDE: usize = MAX_BANDWIDTH + 1;
/// `MAXDISTANCECONSTANT` / `MAX_PATTERN_SIZE_CONSTANT`.
pub const MAX_PATTERN_SIZE_CONSTANT: usize = 2000;
/// Distances up to this size use the full (non banded) wraparound DP.
pub const SMALL_DISTANCE: i32 = 20;
/// Minimum band radius, `MIN_BAND_RADIUS`.
pub const MIN_BAND_RADIUS: i32 = 6;
/// Number of consecutive matches needed to recenter the band.
pub const RECENTER_CRITERION: i32 = 3;
/// Tag separation for the linked list of active distances.
pub const TAG_SEPARATION: i32 = 50;
/// Number of best periods kept for a repeat.
pub const NUMBER_OF_PERIODS: usize = 5;
/// Number of best periods actually tested.
pub const NUMBER_OF_PERIODS_TO_TEST: usize = 3;
/// Number of best periods copied into `sort_multiples`.
pub const NUMBER_OF_PERIODS_INTO_SORTMULTIPLES: usize = 5;
/// Maximum number of different tuple sizes.
pub const MAX_TUPLE_SIZES: usize = 10;
/// Alignment traceback flavor.
pub const LOCAL: u8 = 1;
/// Alignment traceback flavor.
pub const GLOBAL: u8 = 0;

pub const DASH: u8 = b'-';

/// `four_to_the[]` from the reference implementation.
pub const FOUR_TO_THE: [i32; 11] = [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576];

/// Formats a float like glibc's `printf("%*.*f")`, including the `-nan` spelling
/// that C produces for 0/0.
pub fn format_c_float(value: f64, width: usize, precision: usize) -> String {
    if value.is_nan() {
        let body = if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        };
        return format!("{:>width$}", body, width = width);
    }
    format!(
        "{:>width$.precision$}",
        value,
        width = width,
        precision = precision
    )
}

#[derive(Clone, Copy, Default)]
pub struct DistanceEntry {
    pub location: i32,
    pub size: i32,
}

pub struct DistanceList {
    pub k_run_sums_criteria: i32,
    pub waiting_time_criteria: i32,
    pub min_distance_range: i32,
    pub max_distance_range: i32,
    pub entry_count: i32,
    pub match_count: i32,
    pub low_index: i32,
    pub high_index: i32,
    pub linked: bool,
    pub link_down: i32,
    pub link_up: i32,
    pub entries: Vec<DistanceEntry>,
}

impl DistanceList {
    fn new(capacity: usize) -> Self {
        DistanceList {
            k_run_sums_criteria: 0,
            waiting_time_criteria: 0,
            min_distance_range: 0,
            max_distance_range: 0,
            entry_count: 0,
            match_count: 0,
            low_index: 0,
            high_index: 0,
            linked: false,
            link_down: 0,
            link_up: 0,
            entries: vec![DistanceEntry::default(); capacity],
        }
    }
}

#[derive(Clone, Copy, Default)]
pub struct DistanceSeenEntry {
    /// Location in the sequence where the distance was last seen.
    pub location: i32,
    pub end: i32,
    pub score: i32,
}

#[derive(Clone, Copy, Default)]
pub struct BestPeriodElement {
    pub index_high: i32,
    pub index_low: i32,
    pub best_1: i32,
    pub best_2: i32,
    pub best_3: i32,
    pub best_4: i32,
    pub best_5: i32,
}

#[derive(Clone, Copy, Default)]
pub struct HistoryEntry {
    pub location: i32,
    pub previous: i32,
    pub code: i32,
}

/// Working copy of the consensus pattern plus the nucleotide counters.
pub struct Consensus {
    pub pattern: Vec<u8>,
    pub a: Vec<i32>,
    pub c: Vec<i32>,
    pub g: Vec<i32>,
    pub t: Vec<i32>,
    pub dash: Vec<i32>,
    pub insert: Vec<i32>,
    pub letters: Vec<i32>,
    pub total: Vec<i32>,
}

impl Default for Consensus {
    fn default() -> Self {
        let capacity = 2 * (MAX_PATTERN_SIZE_CONSTANT + 1);
        Consensus {
            pattern: vec![DASH; capacity],
            a: vec![0; capacity],
            c: vec![0; capacity],
            g: vec![0; capacity],
            t: vec![0; capacity],
            dash: vec![0; capacity],
            insert: vec![0; capacity],
            letters: vec![0; capacity],
            total: vec![0; capacity],
        }
    }
}

impl Consensus {
    pub fn clear(&mut self) {
        let last_index = 2 * MAX_PATTERN_SIZE_CONSTANT;
        for index in 0..=last_index {
            self.a[index] = 0;
            self.c[index] = 0;
            self.g[index] = 0;
            self.t[index] = 0;
            self.dash[index] = 0;
            self.insert[index] = 0;
            self.letters[index] = 0;
            self.total[index] = 0;
            self.pattern[index] = DASH;
        }
    }
}

/// Traceback alignment between the sequence and the current pattern (`EC`).
pub struct AlignPair {
    pub length: i32,
    pub score: i32,
    pub sequence_text: Vec<u8>,
    pub pattern_text: Vec<u8>,
    pub sequence_index: Vec<i32>,
    pub pattern_index: Vec<i32>,
}

impl AlignPair {
    fn new(size: usize) -> Self {
        AlignPair {
            length: 0,
            score: 0,
            sequence_text: vec![0; size + 2],
            pattern_text: vec![0; size + 2],
            sequence_index: vec![0; size + 2],
            pattern_index: vec![0; size + 2],
        }
    }
}

/// The complete algorithm state for one sequence.
pub struct Trf {
    pub params: ParamSet,
    /// The input record (name, composition and raw sequence).
    pub fasta: FastaSequence,
    /// Sequence characters; index 0 is a sentinel, `sequence[1..=length]` are bases.
    pub sequence: Vec<u8>,
    pub length: i32,
    /// `max_wrap_length`: number of rows of `S`.
    pub band_width: usize,

    /// Wraparound dynamic programming scoring matrix (`S`).
    pub score_matrix: Vec<i32>,
    pub band_up: [i32; MAX_BANDWIDTH + 1],
    pub band_diag: [i32; MAX_BANDWIDTH + 1],
    /// Pattern currently being aligned against the sequence (`EC`).
    pub current_pattern: Vec<u8>,
    pub consensus: Consensus,
    pub align_pair: AlignPair,

    /// Similarity matrix, `SM[256*a+b]`.
    pub similarity_matrix: Vec<i32>,
    /// Maps each base to its nucleotide index (A/C/G/T -> 0-3), `Index[c]`.
    pub nucleotide_index: [i32; 256],

    pub alpha: i32,
    #[allow(dead_code)]
    pub beta: i32,
    pub delta: i32,
    pub match_probability_percent: i32,
    pub indel_probability_percent: i32,
    pub indel_probability: f64,

    pub minscore: i32,
    pub maxperiod: i32,
    pub max_distance: i32,
    #[allow(dead_code)]
    pub maxpatternsize: i32,
    pub max_real_col: i32,
    pub report_min: i32,
    pub heading: i32,
    pub class_length: i32,
    pub consensus_class_length: i32,
    pub period: i32,
    pub minsize: i32,

    /// Counter of processed rows, a `double` in the reference implementation.
    pub rows: f64,
    pub total_characters: f64,
    pub wrap_end: i32,
    pub max_real_row: i32,
    pub max_row: i32,
    pub max_col: i32,
    pub max_score: i32,

    pub tags: Vec<i32>,
    pub top_tag: i32,
    pub distances: Vec<DistanceList>,
    pub distance_seen_entries: Vec<DistanceSeenEntry>,
    pub best_period_list: VecDeque<BestPeriodElement>,

    pub criteria_counts: Vec<i32>,
    pub consensus_counts: Vec<i32>,
    pub output_size_counts: Vec<i32>,
    pub cell_counts: Vec<f64>,
    pub distance_statistics: Vec<i32>,
    pub sort_multiples: Vec<i32>,

    pub num_tuple_sizes: i32,
    pub tuple_sizes: [i32; MAX_TUPLE_SIZES + 1],
    pub tuple_max_distances: [i32; MAX_TUPLE_SIZES + 1],
    pub tuple_codes: [i32; MAX_TUPLE_SIZES + 1],
    pub tuple_hashes: Vec<Vec<i32>>,
    pub histories: Vec<Vec<HistoryEntry>>,
    pub history_sizes: [i32; MAX_TUPLE_SIZES + 1],
    pub next_free_history_indices: [i32; MAX_TUPLE_SIZES + 1],

    pub band_centers: Vec<i32>,

    pub copynumber: f64,
    pub wraparound_dp_count: f64,
    pub output_count: f64,
    #[allow(dead_code)]
    pub try_waiting_time_count: f64,
    #[allow(dead_code)]
    pub fail_waiting_time_count: f64,
    pub cell_total: f64,
    pub wasted_total: f64,

    pub counter_in_sequence: i32,
    pub print_flanking: bool,
    /// `ready` flag of `SetProgressBar()`.
    pub progress_ready: bool,

    /// Header lines used by the HTML/data/masked output.  `header_sequence` holds
    /// raw bytes because it contains the (possibly non UTF-8) sequence name.
    pub header_sequence: Vec<u8>,
    pub header_parameters: String,
    pub header_length: String,

    /// Alignment file currently being written (`Fptxt`).
    pub alignment_file: Option<BufWriter<File>>,

    /// `GlobalIndexList`: repeats found for the current sequence.
    pub records: Vec<crate::clean::RepeatRecord>,
}

impl Trf {
    pub fn new(params: ParamSet, seq: &FastaSequence) -> Trf {
        let length = seq.length as i32;

        let mut sequence = Vec::with_capacity(seq.length + 2);
        sequence.push(0u8);
        sequence.extend_from_slice(&seq.sequence);
        sequence.push(0u8);

        let band_width = std::cmp::min(params.max_wrap_length as usize, seq.length);

        let alpha = params.match_weight as i32;
        let beta = -(params.mismatch as i32);
        let delta = -(params.indel as i32);

        let mut max_distance = params.maxperiod as i32;
        if max_distance < 500 {
            max_distance = 500;
        }
        max_distance = std::cmp::min(max_distance, (length as f64 * 0.6) as i32);
        max_distance = std::cmp::max(max_distance, 200);

        // init_sm() / init_index()
        let mut similarity_matrix = vec![beta; 256 * 256];
        similarity_matrix[b'A' as usize * 256 + b'A' as usize] = alpha;
        similarity_matrix[b'C' as usize * 256 + b'C' as usize] = alpha;
        similarity_matrix[b'G' as usize * 256 + b'G' as usize] = alpha;
        similarity_matrix[b'T' as usize * 256 + b'T' as usize] = alpha;

        let mut nucleotide_index = [0i32; 256];
        nucleotide_index[b'A' as usize] = 0;
        nucleotide_index[b'C' as usize] = 1;
        nucleotide_index[b'G' as usize] = 2;
        nucleotide_index[b'T' as usize] = 3;

        // new_distancelist()
        let mut distances = Vec::with_capacity(max_distance as usize + 1);
        for distance_value in 0..=max_distance {
            let capacity = std::cmp::max(distance_value, 20) as usize + 1;
            distances.push(DistanceList::new(capacity));
        }
        // clear_distancelist()
        for distance_value in 1..=max_distance {
            let list = &mut distances[distance_value as usize];
            list.low_index = 0;
            list.high_index = std::cmp::max(distance_value, 20);
            list.entry_count = 0;
            list.match_count = 0;
        }
        // init_links()
        distances[0].link_up = max_distance + 1;

        let top_tag = max_distance / TAG_SEPARATION;

        let mut trf = Trf {
            alpha,
            beta,
            delta,
            match_probability_percent: params.match_probability_percent as i32,
            indel_probability_percent: params.indel_probability_percent as i32,
            indel_probability: 0.0,
            minscore: params.minscore as i32,
            maxperiod: params.maxperiod as i32,
            max_distance,
            maxpatternsize: max_distance,
            max_real_col: 0,
            report_min: 0,
            heading: 0,
            class_length: 0,
            consensus_class_length: 0,
            period: 0,
            minsize: 1,
            rows: 0.0,
            total_characters: 0.0,
            wrap_end: 0,
            max_real_row: 0,
            max_row: 0,
            max_col: 0,
            max_score: 0,
            tags: vec![0; (max_distance / TAG_SEPARATION) as usize + 1],
            top_tag,
            distances,
            distance_seen_entries: vec![
                DistanceSeenEntry::default();
                MAX_PATTERN_SIZE_CONSTANT + 1
            ],
            best_period_list: VecDeque::new(),
            criteria_counts: vec![0; 2 * (max_distance as usize + 1)],
            consensus_counts: vec![0; 2 * (max_distance as usize + 1)],
            output_size_counts: vec![0; 2 * (max_distance as usize + 1)],
            cell_counts: vec![0.0; 2 * (max_distance as usize + 1)],
            distance_statistics: vec![0; 4 * max_distance as usize],
            sort_multiples: vec![0; max_distance as usize + 1],
            num_tuple_sizes: 0,
            tuple_sizes: [0; MAX_TUPLE_SIZES + 1],
            tuple_max_distances: [0; MAX_TUPLE_SIZES + 1],
            tuple_codes: [0; MAX_TUPLE_SIZES + 1],
            tuple_hashes: Vec::new(),
            histories: Vec::new(),
            history_sizes: [0; MAX_TUPLE_SIZES + 1],
            next_free_history_indices: [0; MAX_TUPLE_SIZES + 1],
            band_centers: vec![0; band_width + 1],
            copynumber: 0.0,
            wraparound_dp_count: 0.0,
            output_count: 0.0,
            try_waiting_time_count: 0.0,
            fail_waiting_time_count: 0.0,
            cell_total: 0.0,
            wasted_total: 0.0,
            counter_in_sequence: 0,
            print_flanking: false,
            progress_ready: false,
            header_sequence: Vec::new(),
            header_parameters: String::new(),
            header_length: String::new(),
            alignment_file: None,
            records: Vec::new(),
            params,
            fasta: seq.clone(),
            sequence,
            length,
            band_width,
            score_matrix: vec![0; (band_width + 1) * SCORE_MATRIX_STRIDE],
            band_up: [0; MAX_BANDWIDTH + 1],
            band_diag: [0; MAX_BANDWIDTH + 1],
            current_pattern: vec![0; 2 * (MAX_PATTERN_SIZE_CONSTANT + 1)],
            consensus: Consensus::default(),
            align_pair: AlignPair::new(2 * band_width),
            similarity_matrix,
            nucleotide_index,
        };
        trf.score_matrix[0] = 1; // S[0][0]=1
        trf
    }

    #[inline]
    pub fn similarity_score(&self, first_base: u8, second_base: u8) -> i32 {
        self.similarity_matrix[256 * first_base as usize + second_base as usize]
    }

    #[inline]
    pub fn score(&self, row: usize, col: usize) -> i32 {
        self.score_matrix[row * SCORE_MATRIX_STRIDE + col]
    }

    #[inline]
    pub fn set_score(&mut self, row: usize, col: usize, value: i32) {
        self.score_matrix[row * SCORE_MATRIX_STRIDE + col] = value;
    }

    /// 95% range of the random walk caused by indels over `distance` bases.
    pub fn indel_random_walk_range(&self, distance: i32) -> i32 {
        (2.3 * (self.indel_probability * distance as f64).sqrt()).floor() as i32
    }

    /// `trf_message()`: write to the alignment file unless HTML output is off.
    pub fn write_to_alignment_file(&mut self, text: &str) {
        self.write_bytes_to_alignment_file(text.as_bytes());
    }

    /// `trf_message()` for text that is not valid UTF-8 (the sequence name).
    pub fn write_bytes_to_alignment_file(&mut self, text: &[u8]) {
        if !self.params.html_output_off {
            if let Some(f) = self.alignment_file.as_mut() {
                let _ = f.write_all(text);
            }
        }
    }

    /// Writes the `fill_align_pair()` macro.
    #[inline]
    pub fn fill_align_pair(
        &mut self,
        sequence_char: u8,
        pattern_char: u8,
        column: i32,
        sequence_position: i32,
        pattern_position: i32,
    ) {
        let column = column as usize;
        self.align_pair.sequence_text[column] = sequence_char;
        self.align_pair.pattern_text[column] = pattern_char;
        self.align_pair.sequence_index[column] = sequence_position;
        self.align_pair.pattern_index[column] = pattern_position;
    }

    pub fn param_string(&self) -> String {
        self.params.param_string()
    }

    /// `PrintProgress()`: progress messages go to stderr unless `-ngs` is used.
    pub fn progress(&self, text: &str) {
        if !self.params.ngs {
            eprint!("\n{}", text);
        }
    }

    /// `PrintError()`.
    pub fn print_error(&self, text: &str) {
        eprintln!("Error: {}", text);
    }

    /// Terminates the program the way C's `exit(3)` does: buffered output is
    /// flushed before the process ends.  `std::process::exit` alone would drop
    /// the pending HTML written through `BufWriter`.
    pub fn fatal(&mut self, code: i32) -> ! {
        if let Some(mut f) = self.alignment_file.take() {
            let _ = f.flush();
        }
        let _ = std::io::stdout().flush();
        std::process::exit(code)
    }

    /// `SetProgressBar()` for the console build: one dot per 5 percent.
    pub fn set_progress_bar(&mut self) {
        if self.params.percent < 0 {
            if !self.progress_ready {
                return;
            }
            self.progress_ready = false;
            println!();
            return;
        }
        if !self.progress_ready {
            self.progress_ready = true;
            println!();
        }
        if self.params.percent > 100 {
            self.params.percent = 100;
        }
        if self.params.percent % 5 == 0 {
            print!(".");
            let _ = std::io::stdout().flush();
        }
    }
}
