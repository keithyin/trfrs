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

/// Maximum width of the alignment band, `MAXBANDWIDTH`.
pub const MAXBANDWIDTH: usize = 150;
/// Row stride of the `S` scoring matrix.
pub const STRIDE: usize = MAXBANDWIDTH + 1;
/// `MAXDISTANCECONSTANT` / `MAXPATTERNSIZECONSTANT`.
pub const MAXPATTERNSIZECONSTANT: usize = 2000;
/// Distances up to this size use the full (non banded) wraparound DP.
pub const SMALLDISTANCE: i32 = 20;
/// Minimum band radius, `MINBANDRADIUS`.
pub const MINBANDRADIUS: i32 = 6;
/// Number of consecutive matches needed to recenter the band.
pub const RECENTERCRITERION: i32 = 3;
/// Tag separation for the linked list of active distances.
pub const TAGSEP: i32 = 50;
/// Number of best periods kept for a repeat.
pub const NUMBER_OF_PERIODS: usize = 5;
/// Number of best periods actually tested.
pub const NUMBER_OF_PERIODS_TO_TEST: usize = 3;
/// Number of best periods copied into `sort_multiples`.
pub const NUMBER_OF_PERIODS_INTO_SORTMULTIPLES: usize = 5;
/// Maximum number of different tuple sizes.
pub const MAXTUPLESIZES: usize = 10;
/// Alignment traceback flavor.
pub const LOCAL: u8 = 1;
/// Alignment traceback flavor.
pub const GLOBAL: u8 = 0;

pub const DASH: u8 = b'-';

/// `four_to_the[]` from the reference implementation.
pub const FOUR_TO_THE: [i32; 11] = [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576];

/// Formats a float like glibc's `printf("%*.*f")`, including the `-nan` spelling
/// that C produces for 0/0.
pub fn cfloat(value: f64, width: usize, precision: usize) -> String {
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
    pub lo_d_range: i32,
    pub hi_d_range: i32,
    pub numentries: i32,
    pub nummatches: i32,
    pub lowindex: i32,
    pub highindex: i32,
    pub linked: bool,
    pub linkdown: i32,
    pub linkup: i32,
    pub entry: Vec<DistanceEntry>,
}

impl DistanceList {
    fn new(capacity: usize) -> Self {
        DistanceList {
            k_run_sums_criteria: 0,
            waiting_time_criteria: 0,
            lo_d_range: 0,
            hi_d_range: 0,
            numentries: 0,
            nummatches: 0,
            lowindex: 0,
            highindex: 0,
            linked: false,
            linkdown: 0,
            linkup: 0,
            entry: vec![DistanceEntry::default(); capacity],
        }
    }
}

#[derive(Clone, Copy, Default)]
pub struct DistanceseenArrayElement {
    pub index: i32,
    pub end: i32,
    pub score: i32,
}

#[derive(Clone, Copy, Default)]
pub struct BestPeriodElement {
    pub indexhigh: i32,
    pub indexlow: i32,
    pub best1: i32,
    pub best2: i32,
    pub best3: i32,
    pub best4: i32,
    pub best5: i32,
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
        let n = 2 * (MAXPATTERNSIZECONSTANT + 1);
        Consensus {
            pattern: vec![DASH; n],
            a: vec![0; n],
            c: vec![0; n],
            g: vec![0; n],
            t: vec![0; n],
            dash: vec![0; n],
            insert: vec![0; n],
            letters: vec![0; n],
            total: vec![0; n],
        }
    }
}

impl Consensus {
    pub fn clear(&mut self) {
        let n = 2 * MAXPATTERNSIZECONSTANT;
        for i in 0..=n {
            self.a[i] = 0;
            self.c[i] = 0;
            self.g[i] = 0;
            self.t[i] = 0;
            self.dash[i] = 0;
            self.insert[i] = 0;
            self.letters[i] = 0;
            self.total[i] = 0;
            self.pattern[i] = DASH;
        }
    }
}

/// Traceback alignment between the sequence and the current pattern (`EC`).
pub struct AlignPair {
    pub length: i32,
    pub score: i32,
    pub textprime: Vec<u8>,
    pub textsecnd: Vec<u8>,
    pub indexprime: Vec<i32>,
    pub indexsecnd: Vec<i32>,
}

impl AlignPair {
    fn new(size: usize) -> Self {
        AlignPair {
            length: 0,
            score: 0,
            textprime: vec![0; size + 2],
            textsecnd: vec![0; size + 2],
            indexprime: vec![0; size + 2],
            indexsecnd: vec![0; size + 2],
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
    /// `maxwraplength`: number of rows of `S`.
    pub bw: usize,

    /// Wraparound dynamic programming scoring matrix (`S`).
    pub s: Vec<i32>,
    pub up: [i32; MAXBANDWIDTH + 1],
    pub diag: [i32; MAXBANDWIDTH + 1],
    /// Pattern currently being aligned against the sequence (`EC`).
    pub ec: Vec<u8>,
    pub consensus: Consensus,
    pub align_pair: AlignPair,

    /// Similarity matrix, `SM[256*a+b]`.
    pub sm: Vec<i32>,
    /// Index of the four nucleotides, `Index[c]`.
    pub index: [i32; 256],

    pub alpha: i32,
    #[allow(dead_code)]
    pub beta: i32,
    pub delta: i32,
    pub pm: i32,
    pub pi: i32,
    pub pindel: f64,

    pub minscore: i32,
    pub maxperiod: i32,
    pub maxdistance: i32,
    #[allow(dead_code)]
    pub maxpatternsize: i32,
    pub maxrealcol: i32,
    pub reportmin: i32,
    pub heading: i32,
    pub classlength: i32,
    pub cons_classlength: i32,
    pub period: i32,
    pub minsize: i32,

    /// Counter of processed rows, a `double` in the reference implementation.
    pub rows: f64,
    pub totalcharacters: f64,
    pub wrapend: i32,
    pub maxrealrow: i32,
    pub maxrow: i32,
    pub maxcol: i32,
    pub maxscore: i32,

    pub tag: Vec<i32>,
    pub top_tag: i32,
    pub distance: Vec<DistanceList>,
    pub distanceseenarray: Vec<DistanceseenArrayElement>,
    pub best_period_list: VecDeque<BestPeriodElement>,

    pub criteria_count: Vec<i32>,
    pub consensus_count: Vec<i32>,
    pub outputsize_count: Vec<i32>,
    pub cell_count: Vec<f64>,
    pub statistics_distance: Vec<i32>,
    pub sort_multiples: Vec<i32>,

    pub nts: i32,
    pub tuple_size: [i32; MAXTUPLESIZES + 1],
    pub tuple_max_distance: [i32; MAXTUPLESIZES + 1],
    pub tuple_code: [i32; MAXTUPLESIZES + 1],
    pub tuple_hash: Vec<Vec<i32>>,
    pub history: Vec<Vec<HistoryEntry>>,
    pub history_size: [i32; MAXTUPLESIZES + 1],
    pub next_free_history_index: [i32; MAXTUPLESIZES + 1],

    pub band_center: Vec<i32>,

    pub copynumber: f64,
    pub wdp_count: f64,
    pub output_count: f64,
    #[allow(dead_code)]
    pub try_waiting_time_count: f64,
    #[allow(dead_code)]
    pub fail_waiting_time_count: f64,
    pub cell_total: f64,
    pub wasted_total: f64,

    pub counter_in_seq: i32,
    pub print_flanking: bool,
    /// `ready` flag of `SetProgressBar()`.
    pub progress_ready: bool,

    /// Header lines used by the HTML/data/masked output.  `hsequence` holds
    /// raw bytes because it contains the (possibly non UTF-8) sequence name.
    pub hsequence: Vec<u8>,
    pub hparameters: String,
    pub hlength: String,

    /// Alignment file currently being written (`Fptxt`).
    pub fptxt: Option<BufWriter<File>>,

    /// `GlobalIndexList`: repeats found for the current sequence.
    pub records: Vec<crate::clean::IlRecord>,
}

impl Trf {
    pub fn new(params: ParamSet, seq: &FastaSequence) -> Trf {
        let length = seq.length as i32;

        let mut sequence = Vec::with_capacity(seq.length + 2);
        sequence.push(0u8);
        sequence.extend_from_slice(&seq.sequence);
        sequence.push(0u8);

        let bw = std::cmp::min(params.maxwraplength as usize, seq.length);

        let alpha = params.match_weight as i32;
        let beta = -(params.mismatch as i32);
        let delta = -(params.indel as i32);

        let mut maxdistance = params.maxperiod as i32;
        if maxdistance < 500 {
            maxdistance = 500;
        }
        maxdistance = std::cmp::min(maxdistance, (length as f64 * 0.6) as i32);
        maxdistance = std::cmp::max(maxdistance, 200);

        // init_sm() / init_index()
        let mut sm = vec![beta; 256 * 256];
        sm[b'A' as usize * 256 + b'A' as usize] = alpha;
        sm[b'C' as usize * 256 + b'C' as usize] = alpha;
        sm[b'G' as usize * 256 + b'G' as usize] = alpha;
        sm[b'T' as usize * 256 + b'T' as usize] = alpha;

        let mut index = [0i32; 256];
        index[b'A' as usize] = 0;
        index[b'C' as usize] = 1;
        index[b'G' as usize] = 2;
        index[b'T' as usize] = 3;

        // new_distancelist()
        let mut distance = Vec::with_capacity(maxdistance as usize + 1);
        for g in 0..=maxdistance {
            let capacity = std::cmp::max(g, 20) as usize + 1;
            distance.push(DistanceList::new(capacity));
        }
        // clear_distancelist()
        for g in 1..=maxdistance {
            let d = &mut distance[g as usize];
            d.lowindex = 0;
            d.highindex = std::cmp::max(g, 20);
            d.numentries = 0;
            d.nummatches = 0;
        }
        // init_links()
        distance[0].linkup = maxdistance + 1;

        let top_tag = maxdistance / TAGSEP;

        let mut trf = Trf {
            alpha,
            beta,
            delta,
            pm: params.pm as i32,
            pi: params.pi as i32,
            pindel: 0.0,
            minscore: params.minscore as i32,
            maxperiod: params.maxperiod as i32,
            maxdistance,
            maxpatternsize: maxdistance,
            maxrealcol: 0,
            reportmin: 0,
            heading: 0,
            classlength: 0,
            cons_classlength: 0,
            period: 0,
            minsize: 1,
            rows: 0.0,
            totalcharacters: 0.0,
            wrapend: 0,
            maxrealrow: 0,
            maxrow: 0,
            maxcol: 0,
            maxscore: 0,
            tag: vec![0; (maxdistance / TAGSEP) as usize + 1],
            top_tag,
            distance,
            distanceseenarray: vec![
                DistanceseenArrayElement::default();
                MAXPATTERNSIZECONSTANT + 1
            ],
            best_period_list: VecDeque::new(),
            criteria_count: vec![0; 2 * (maxdistance as usize + 1)],
            consensus_count: vec![0; 2 * (maxdistance as usize + 1)],
            outputsize_count: vec![0; 2 * (maxdistance as usize + 1)],
            cell_count: vec![0.0; 2 * (maxdistance as usize + 1)],
            statistics_distance: vec![0; 4 * maxdistance as usize],
            sort_multiples: vec![0; maxdistance as usize + 1],
            nts: 0,
            tuple_size: [0; MAXTUPLESIZES + 1],
            tuple_max_distance: [0; MAXTUPLESIZES + 1],
            tuple_code: [0; MAXTUPLESIZES + 1],
            tuple_hash: Vec::new(),
            history: Vec::new(),
            history_size: [0; MAXTUPLESIZES + 1],
            next_free_history_index: [0; MAXTUPLESIZES + 1],
            band_center: vec![0; bw + 1],
            copynumber: 0.0,
            wdp_count: 0.0,
            output_count: 0.0,
            try_waiting_time_count: 0.0,
            fail_waiting_time_count: 0.0,
            cell_total: 0.0,
            wasted_total: 0.0,
            counter_in_seq: 0,
            print_flanking: false,
            progress_ready: false,
            hsequence: Vec::new(),
            hparameters: String::new(),
            hlength: String::new(),
            fptxt: None,
            records: Vec::new(),
            params,
            fasta: seq.clone(),
            sequence,
            length,
            bw,
            s: vec![0; (bw + 1) * STRIDE],
            up: [0; MAXBANDWIDTH + 1],
            diag: [0; MAXBANDWIDTH + 1],
            ec: vec![0; 2 * (MAXPATTERNSIZECONSTANT + 1)],
            consensus: Consensus::default(),
            align_pair: AlignPair::new(2 * bw),
            sm,
            index,
        };
        trf.s[0] = 1; // S[0][0]=1
        trf
    }

    #[inline]
    pub fn mm(&self, a: u8, b: u8) -> i32 {
        self.sm[256 * a as usize + b as usize]
    }

    #[inline]
    pub fn get_s(&self, r: usize, i: usize) -> i32 {
        self.s[r * STRIDE + i]
    }

    #[inline]
    pub fn put_s(&mut self, r: usize, i: usize, v: i32) {
        self.s[r * STRIDE + i] = v;
    }

    /// `d_range(d)`: 95% range of the random walk caused by indels.
    pub fn d_range(&self, d: i32) -> i32 {
        (2.3 * (self.pindel * d as f64).sqrt()).floor() as i32
    }

    /// `trf_message()`: write to the alignment file unless HTML output is off.
    pub fn msg(&mut self, text: &str) {
        self.msg_bytes(text.as_bytes());
    }

    /// `trf_message()` for text that is not valid UTF-8 (the sequence name).
    pub fn msg_bytes(&mut self, text: &[u8]) {
        if !self.params.htmloff {
            if let Some(f) = self.fptxt.as_mut() {
                let _ = f.write_all(text);
            }
        }
    }

    /// Writes the `fill_align_pair()` macro.
    #[inline]
    pub fn fill_align_pair(&mut self, c1: u8, c2: u8, l: i32, i: i32, j: i32) {
        let l = l as usize;
        self.align_pair.textprime[l] = c1;
        self.align_pair.textsecnd[l] = c2;
        self.align_pair.indexprime[l] = i;
        self.align_pair.indexsecnd[l] = j;
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
        if let Some(mut f) = self.fptxt.take() {
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
