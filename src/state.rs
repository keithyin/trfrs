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
    /// `new_DistanceList()` + `clear_DistanceList()`: the ring holds
    /// `max(distance, 20) + 1` entries.  `high_index` starts on the last slot so
    /// the first recorded match lands in slot 0; `low_index` starts at slot 0.
    fn new(distance: i32) -> Self {
        let window_size = distance.max(20) + 1;
        DistanceList {
            k_run_sums_criteria: 0,
            waiting_time_criteria: 0,
            min_distance_range: 0,
            max_distance_range: 0,
            entry_count: 0,
            match_count: 0,
            low_index: 0,
            high_index: window_size - 1,
            linked: false,
            link_down: 0,
            link_up: 0,
            entries: vec![DistanceEntry::default(); window_size as usize],
        }
    }

    /// Drops entries that have left the trailing window ending at `location`
    /// and advances `low_index` around the ring.
    ///
    /// This is the eviction loop the reference implementation repeats at the
    /// start of `add_tuple_match_to_Distance_entry()` and
    /// `no_matches_so_unlink_Distance()`; the statement order mirrors it.  The
    /// modulus is the allocated ring capacity, `max(distance, 20) + 1`.
    #[inline]
    pub fn evict_outside_window(&mut self, location: i32, distance: i32) {
        let window_size = self.entries.len() as i32;
        let window_left_end = std::cmp::min(location - distance, location - 20) + 1;
        while self.entry_count != 0 {
            let lo = self.low_index as usize;
            if self.entries[lo].location >= window_left_end {
                break;
            }
            let entry_size = self.entries[lo].size;
            self.entry_count -= 1;
            self.match_count -= entry_size;
            self.low_index += 1;
            self.low_index %= window_size;
        }
    }

    /// `add_tuple_match_to_Distance_entry()`: records one tuple occurrence
    /// ending at `location`.
    ///
    /// Expired entries are dropped first.  If the newest surviving entry ends
    /// at `location - 1`, the occurrence extends it by one tuple; otherwise a
    /// new entry holding `tuple_size` tuples is written at `high_index`,
    /// wrapping around the ring.
    #[inline]
    pub fn record_match(&mut self, location: i32, tuple_size: i32, distance: i32) {
        self.evict_outside_window(location, distance);

        let hi = self.high_index as usize;
        if self.entry_count != 0 && self.entries[hi].location == location - 1 {
            /* adjacent tuple: extend the last entry */
            self.entries[hi].location += 1;
            self.entries[hi].size += 1;
            self.match_count += 1;
        } else {
            self.entry_count += 1;
            self.high_index += 1;
            self.high_index %= self.entries.len() as i32;
            let hi = self.high_index as usize;
            self.entries[hi].location = location;
            self.entries[hi].size = tuple_size;
            self.match_count += tuple_size;
        }
    }

    /// Start of the oldest match kept by this list:
    /// `entries[low_index].location - entries[low_index].size + tuple_size`.
    /// Only meaningful while `entry_count != 0`.
    #[inline]
    pub fn first_match_location(&self, tuple_size: i32) -> i32 {
        let entry = self.entries[self.low_index as usize];
        entry.location - entry.size + tuple_size
    }
}

#[derive(Clone, Copy, Default)]
pub struct DistanceSeenEntry {
    /// Location in the sequence where the distance was last seen.  Written by
    /// [`DistanceSeen::record`] but never read, exactly as in the reference
    /// implementation; kept to mirror the C `distance_seen` struct.
    #[allow(dead_code)]
    pub location: i32,
    /// Row reached by the most recent DP run for this distance.
    pub end: i32,
    /// Score of the most recent DP run.  Written but never read, as in C.
    #[allow(dead_code)]
    pub score: i32,
}

/// `distance_seen[]` from the reference implementation: for each distance, the
/// outcome of the last wraparound DP run, used to skip candidates it covered.
pub struct DistanceSeen {
    /// One slot per distance, indexed exactly as the C array (`0..=2000`).
    pub entries: Vec<DistanceSeenEntry>,
}

impl DistanceSeen {
    /// Allocates `distance_seen[MAX_PATTERN_SIZE_CONSTANT + 1]`, zeroed as in
    /// the reference implementation.  Private: only [`Trf::new`] creates it.
    fn new() -> Self {
        DistanceSeen {
            entries: vec![DistanceSeenEntry::default(); MAX_PATTERN_SIZE_CONSTANT + 1],
        }
    }

    /// `add_to_distance_seen()`: records the result of the DP run for
    /// `distance`.  Argument order follows the C function.
    #[inline]
    pub fn record(&mut self, location: i32, distance: i32, end: i32, score: i32) {
        let entry = &mut self.entries[distance as usize];
        entry.location = location;
        entry.end = end;
        entry.score = score;
    }

    /// `search_for_distance_match_in_distance_seen()`: true when the last run
    /// for `distance` already reached `start` (`end >= start`).
    #[inline]
    pub fn covers(&self, distance: i32, start: i32) -> bool {
        self.entries[distance as usize].end >= start
    }
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

/// The rolling tuple-hash machinery for one tuple size: `tuple_hash[k]`,
/// `history[k]`, `history_size[k]` and `next_free_history_index[k]` in the
/// reference implementation.
///
/// `hash[code]` points into `entries`; slot/location `0` is the null sentinel
/// (`location == 0` marks a free slot, `previous == 0` ends a code's chain).
/// When [`Self::push`] recycles a slot it also clears the hash entry that still
/// points at it, so the two arrays are only consistent when mutated together.
#[derive(Default)]
pub struct TupleHistory {
    /// `tuple_hash[k]`: slot of the newest occurrence of each `4^tuple_size` code.
    pub hash: Vec<i32>,
    /// `history[k]`: ring of occurrences; slot 0 is never written.
    pub entries: Vec<HistoryEntry>,
    /// `history_size[k]`: `entries.len()` as an `i32`, mirroring the C comparison.
    pub size: i32,
    /// `next_free_history_index[k]`: slot written by the next `push`, wrapping to 1.
    pub next_free: i32,
}

impl TupleHistory {
    /// Capacities exactly as allocated by the reference `scan_sequence`:
    /// `4^tuple_size` hash entries and `2 * (max_distance + 1) + 2` slots.
    /// Private: only [`Trf::init_tuple_tables`] (same module) creates these.
    fn new(tuple_size: i32, max_distance: i32) -> Self {
        let size = 2 * (max_distance + 1) + 2;
        TupleHistory {
            hash: vec![0; FOUR_TO_THE[tuple_size as usize] as usize],
            entries: vec![HistoryEntry::default(); size as usize],
            size,
            next_free: 1,
        }
    }

    /// Records the occurrence of `code` ending at `position`, evicting the slot
    /// about to be recycled when it is still its code's chain head.
    ///
    /// Returns `(write_slot, previous_head)`: the slot just written (the walk's
    /// `prev_node`) and the head to start the distance walk from.  The statement
    /// order mirrors `scan_sequence`; in particular the eviction check inspects
    /// `next_free_slot`, not `write_slot`.
    #[inline]
    pub fn push(&mut self, position: i32, code: i32) -> (i32, i32) {
        let previous_head = self.hash[code as usize];
        let write_slot = self.next_free;
        let mut next_free_slot = write_slot + 1;
        if next_free_slot == self.size {
            next_free_slot = 1;
        }

        // HistoryEntry is Copy, so this value read holds no borrow across the
        // hash update below.
        let evicted = self.entries[next_free_slot as usize];
        if evicted.location != 0 && next_free_slot == self.hash[evicted.code as usize] {
            self.hash[evicted.code as usize] = 0;
        }
        self.next_free = next_free_slot;

        self.hash[code as usize] = write_slot;
        self.entries[write_slot as usize] = HistoryEntry {
            location: position,
            previous: previous_head,
            code,
        };

        (write_slot, previous_head)
    }
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
    /// `distance_seen[]`: per-distance DP coverage bookkeeping.
    pub distance_seen: DistanceSeen,
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
    /// Per tuple size `k`: hash table and rolling history.  Index 0 is an unused
    /// sentinel so `k` stays 1-based; live entries are `1..=num_tuple_sizes`.
    pub tuple_histories: Vec<TupleHistory>,

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

        // new_distancelist() + clear_distancelist()
        let mut distances = Vec::with_capacity(max_distance as usize + 1);
        for distance_value in 0..=max_distance {
            distances.push(DistanceList::new(distance_value));
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
            distance_seen: DistanceSeen::new(),
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
            tuple_histories: Vec::new(),
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

    /// Allocates [`Self::tuple_histories`] from the configured tuple sizes.
    /// Index 0 is left as an inert sentinel so every tuple loop keeps the
    /// reference implementation's 1-based indexing.
    pub fn init_tuple_tables(&mut self) {
        let num_tuple_sizes = self.num_tuple_sizes as usize;
        self.tuple_histories = Vec::with_capacity(num_tuple_sizes + 1);
        self.tuple_histories.push(TupleHistory::default());
        for tuple_index in 1..=num_tuple_sizes {
            self.tuple_histories.push(TupleHistory::new(
                self.tuple_sizes[tuple_index],
                self.tuple_max_distances[tuple_index],
            ));
        }
        debug_assert_eq!(self.tuple_histories.len(), num_tuple_sizes + 1);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_history_capacities_match_reference() {
        let history = TupleHistory::new(4, 29);
        assert_eq!(history.hash.len(), FOUR_TO_THE[4] as usize);
        assert_eq!(history.size, 2 * (29 + 1) + 2);
        assert_eq!(history.entries.len(), history.size as usize);
        assert_eq!(history.next_free, 1);
        assert!(history.entries.iter().all(|entry| entry.location == 0));
    }

    #[test]
    fn push_links_new_entry_to_previous_head() {
        let mut history = TupleHistory::new(4, 3);
        assert_eq!(history.push(10, 2), (1, 0));
        assert_eq!(history.hash[2], 1);
        assert_eq!(history.push(20, 2), (2, 1));
        assert_eq!(history.hash[2], 2);
        assert_eq!(history.entries[2].location, 20);
        assert_eq!(history.entries[2].previous, 1);
        assert_eq!(history.entries[2].code, 2);
    }

    #[test]
    fn push_wraps_and_clears_evicted_hash_entries() {
        // size = 2 * (1 + 1) + 2 = 6, usable slots 1..=5.
        let mut history = TupleHistory::new(1, 1);
        for (position, code) in (1..=4).zip(0..=3) {
            history.push(position, code);
        }
        // Writing slot 5 inspects slot 1, which is still code 0's chain head.
        assert_eq!(history.push(5, 0), (5, 1));
        assert_eq!(history.hash[0], 5);
        // The next write recycles slot 1 and clears slot 2's stale head entry.
        assert_eq!(history.push(6, 2), (1, 3));
        assert_eq!(history.hash[2], 1);
        assert_eq!(history.hash[1], 0);
        assert_eq!(history.hash[0], 5);
    }

    #[test]
    fn default_tuple_history_is_an_inert_sentinel() {
        let history = TupleHistory::default();
        assert!(history.hash.is_empty());
        assert!(history.entries.is_empty());
        assert_eq!(history.size, 0);
        assert_eq!(history.next_free, 0);
    }

    #[test]
    fn init_tuple_tables_keeps_index_zero_as_sentinel() {
        let fasta = FastaSequence::from_bytes(b"t", b"ACGTACGT");
        let mut trf = Trf::new(ParamSet::default(), &fasta);
        trf.num_tuple_sizes = 3;
        trf.tuple_sizes[1] = 4;
        trf.tuple_sizes[2] = 5;
        trf.tuple_sizes[3] = 7;
        trf.tuple_max_distances[1] = 29;
        trf.tuple_max_distances[2] = 159;
        trf.tuple_max_distances[3] = trf.max_distance;

        trf.init_tuple_tables();

        assert_eq!(trf.tuple_histories.len(), 4);
        assert!(trf.tuple_histories[0].entries.is_empty());
        assert_eq!(trf.tuple_histories[1].hash.len(), FOUR_TO_THE[4] as usize);
        assert_eq!(trf.tuple_histories[2].size, 2 * (159 + 1) + 2);
        assert_eq!(trf.tuple_histories[3].size, 2 * (trf.max_distance + 1) + 2);
        assert!(trf.tuple_histories[1..]
            .iter()
            .all(|history| history.next_free == 1));
    }

    #[test]
    fn distance_list_new_initializes_ring_cursors() {
        let list = DistanceList::new(4);
        assert_eq!(list.entries.len(), 21);
        assert_eq!(list.low_index, 0);
        assert_eq!(list.high_index, 20);
        assert_eq!(list.entry_count, 0);
        assert_eq!(list.match_count, 0);

        let list = DistanceList::new(30);
        assert_eq!(list.entries.len(), 31);
        assert_eq!(list.high_index, 30);
    }

    #[test]
    fn record_match_places_first_entry_in_wrapped_slot() {
        let mut list = DistanceList::new(30);
        list.record_match(100, 5, 30);
        assert_eq!(list.entry_count, 1);
        assert_eq!(list.match_count, 5);
        assert_eq!(list.low_index, 0);
        assert_eq!(list.high_index, 0);
        assert_eq!((list.entries[0].location, list.entries[0].size), (100, 5));
    }

    #[test]
    fn record_match_extends_adjacent_entry() {
        let mut list = DistanceList::new(30);
        list.record_match(100, 5, 30);
        list.record_match(101, 5, 30);
        list.record_match(102, 5, 30);
        assert_eq!(list.entry_count, 1);
        assert_eq!(list.match_count, 7);
        assert_eq!(list.high_index, 0);
        assert_eq!((list.entries[0].location, list.entries[0].size), (102, 7));
    }

    #[test]
    fn record_match_starts_new_entry_after_gap() {
        let mut list = DistanceList::new(30);
        list.record_match(100, 5, 30);
        list.record_match(120, 4, 30);
        assert_eq!(list.entry_count, 2);
        assert_eq!(list.match_count, 9);
        assert_eq!(list.low_index, 0);
        assert_eq!(list.high_index, 1);
        assert_eq!((list.entries[0].location, list.entries[0].size), (100, 5));
        assert_eq!((list.entries[1].location, list.entries[1].size), (120, 4));
    }

    #[test]
    fn record_match_evicts_entries_outside_window() {
        let mut list = DistanceList::new(30);
        list.record_match(100, 5, 30);
        list.record_match(120, 4, 30);
        list.record_match(140, 3, 30); // window_left_end = 111 drops location 100
        assert_eq!(list.entry_count, 2);
        assert_eq!(list.match_count, 7);
        assert_eq!(list.low_index, 1);
        assert_eq!(list.high_index, 2);
        assert_eq!((list.entries[2].location, list.entries[2].size), (140, 3));
        // Eviction only moves the counters and cursors; it never clears the slot.
        assert_eq!((list.entries[0].location, list.entries[0].size), (100, 5));
    }

    #[test]
    fn record_match_evicts_to_empty_then_writes() {
        let mut list = DistanceList::new(30);
        list.record_match(100, 5, 30);
        list.record_match(200, 4, 30); // window_left_end = 171 drops the only entry
        assert_eq!(list.entry_count, 1);
        assert_eq!(list.match_count, 4);
        assert_eq!(list.low_index, 1);
        assert_eq!(list.high_index, 1);
        assert_eq!((list.entries[1].location, list.entries[1].size), (200, 4));
        assert_eq!((list.entries[0].location, list.entries[0].size), (100, 5));
    }

    #[test]
    fn short_distances_use_twenty_position_window() {
        let mut list = DistanceList::new(5);
        list.record_match(100, 2, 5);
        list.record_match(105, 2, 5);
        list.record_match(120, 1, 5); // window_left_end = 101: 100 goes, 105 stays
        assert_eq!(list.entry_count, 2);
        assert_eq!(list.match_count, 3);
        assert_eq!(list.low_index, 1);
        assert_eq!(list.high_index, 2);
        assert_eq!((list.entries[1].location, list.entries[1].size), (105, 2));
        assert_eq!((list.entries[2].location, list.entries[2].size), (120, 1));
    }

    #[test]
    fn record_match_wraps_high_index() {
        let mut list = DistanceList::new(5); // 21 slots; the first write lands in slot 0
        for step in 0..22 {
            list.record_match(100 + 2 * step, 2, 5); // gap 2: never adjacent
        }
        assert_eq!(list.high_index, 0); // slot 0 again after 22 inserts
        assert_eq!(list.low_index, 12);
        assert_eq!(list.entry_count, 10);
        assert_eq!(list.match_count, 20);
        assert_eq!((list.entries[0].location, list.entries[0].size), (142, 2));

        list.record_match(143, 2, 5); // adjacent to the entry in wrapped slot 0
        assert_eq!(list.high_index, 0);
        assert_eq!(list.low_index, 12);
        assert_eq!(list.entry_count, 10);
        assert_eq!(list.match_count, 21);
        assert_eq!((list.entries[0].location, list.entries[0].size), (143, 3));
    }

    #[test]
    fn first_match_location_returns_oldest_entry_start() {
        let mut list = DistanceList::new(30);
        list.record_match(100, 5, 30);
        list.record_match(120, 4, 30);
        assert_eq!(list.first_match_location(5), 100);
        assert_eq!(list.first_match_location(4), 99);

        list.record_match(140, 3, 30);
        assert_eq!(list.low_index, 1);
        assert_eq!(list.first_match_location(3), 119);
    }

    #[test]
    fn distance_seen_new_allocates_one_zeroed_slot_per_distance() {
        let seen = DistanceSeen::new();
        assert_eq!(seen.entries.len(), MAX_PATTERN_SIZE_CONSTANT + 1);
        assert!(seen
            .entries
            .iter()
            .all(|entry| entry.location == 0 && entry.end == 0 && entry.score == 0));
    }

    #[test]
    fn record_writes_all_fields_at_distance_index() {
        let mut seen = DistanceSeen::new();
        seen.record(10, 3, 42, 7);
        assert_eq!(
            (
                seen.entries[3].location,
                seen.entries[3].end,
                seen.entries[3].score
            ),
            (10, 42, 7)
        );
        assert_eq!((seen.entries[2].end, seen.entries[4].end), (0, 0));
    }

    #[test]
    fn record_overwrites_the_entry_for_the_same_distance() {
        let mut seen = DistanceSeen::new();
        seen.record(10, 3, 42, 7);
        seen.record(20, 3, 60, 9);
        assert_eq!(
            (
                seen.entries[3].location,
                seen.entries[3].end,
                seen.entries[3].score
            ),
            (20, 60, 9)
        );
    }

    #[test]
    fn covers_compares_end_with_start() {
        let mut seen = DistanceSeen::new();
        assert!(!seen.covers(5, 100));
        seen.record(50, 5, 100, 3);
        assert!(seen.covers(5, 100)); // equality counts
        assert!(seen.covers(5, 99));
        assert!(!seen.covers(5, 101));
        assert!(!seen.covers(4, 100)); // slots are independent
    }

    #[test]
    fn zeroed_slot_covers_only_non_positive_starts() {
        // `end` starts at 0, so `end >= start` holds for start <= 0; candidate
        // positions are >= 1, so the sentinel is inert in the scan.
        let seen = DistanceSeen::new();
        assert!(seen.covers(7, 0));
        assert!(!seen.covers(7, 1));
    }

    #[test]
    fn record_accepts_max_pattern_size_index() {
        let mut seen = DistanceSeen::new();
        let max = MAX_PATTERN_SIZE_CONSTANT as i32;
        seen.record(500, max, 520, 11);
        assert_eq!(seen.entries[MAX_PATTERN_SIZE_CONSTANT].end, 520);
        assert!(seen.covers(max, 520));
    }
}
