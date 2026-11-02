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

//! Wraparound dynamic programming and traceback.
//!
//! Direct translation of `wraparound_dp()`, `narrow_band_wraparound_dp()`,
//! `get_pair_alignment_with_copynumber()` and
//! `get_narrow_band_pair_alignment_with_copynumber()` from the reference
//! implementation.  [`Trf`] keeps thin wrappers with the original signatures.

use crate::state::{
    AlignPair, Trf, DASH, GLOBAL, LOCAL, MAX_BANDWIDTH, MAX_PATTERN_SIZE_CONSTANT,
    SCORE_MATRIX_STRIDE,
};

#[inline]
fn max2(a: i32, b: i32) -> i32 {
    if a >= b {
        a
    } else {
        b
    }
}

#[inline]
fn max3(a: i32, b: i32, c: i32) -> i32 {
    if a >= b {
        if a >= c {
            a
        } else {
            c
        }
    } else if b >= c {
        b
    } else {
        c
    }
}

#[inline]
fn max4(a: i32, b: i32, c: i32, d: i32) -> i32 {
    if a >= b {
        if a >= c {
            if a >= d {
                a
            } else {
                d
            }
        } else if c >= d {
            c
        } else {
            d
        }
    } else if b >= c {
        if b >= d {
            b
        } else {
            d
        }
    } else if c >= d {
        c
    } else {
        d
    }
}

/// Column bookkeeping shared by the two halves of the banded alignment.
#[derive(Default)]
struct TraceState {
    max_score: i32,
    max_real_row: i32,
    max_row: i32,
    max_col: i32,
    min_real_row: i32,
    min_col: i32,
    min_col_band_center: i32,
    min_col_position: i32,
    max_row_score: i32,
    match_at_max_col: i32,
    matches_in_diagonal: i32,
    end_of_trace: bool,
}

/// Read-only inputs borrowed from [`Trf`] for one DP/traceback call.
pub struct DpContext<'a> {
    sequence: &'a [u8],
    length: i32,
    alpha: i32,
    delta: i32,
    similarity_matrix: &'a [i32],
    consensus_pattern: &'a [u8],
}

impl DpContext<'_> {
    #[inline]
    fn similarity_score(&self, first_base: u8, second_base: u8) -> i32 {
        self.similarity_matrix[256 * first_base as usize + second_base as usize]
    }
}

/// Wraparound DP workspace and results (`S`, band caches, `EC` and `max_*`).
pub struct WraparoundDp {
    score_matrix: Vec<i32>,
    band_up: [i32; MAX_BANDWIDTH + 1],
    band_diag: [i32; MAX_BANDWIDTH + 1],
    band_centers: Vec<i32>,
    max_real_col: i32,
    max_row: i32,
    max_col: i32,
    band_width: usize,
    /// Write-only in the reference implementation; kept for C parity.
    #[allow(dead_code)]
    wrap_end: i32,

    // Accessed directly by scan.rs, control.rs and algo.rs.
    pub current_pattern: Vec<u8>,
    pub max_real_row: i32,
    pub max_score: i32,
    /// Rows processed by the current DP run (C's `rows` counter).
    pub rows: f64,
    /// Report threshold; always 0, kept to mirror the C structure.
    pub report_min: i32,
}

impl WraparoundDp {
    pub fn new(band_width: usize) -> Self {
        let mut dp = WraparoundDp {
            score_matrix: vec![0; (band_width + 1) * SCORE_MATRIX_STRIDE],
            band_up: [0; MAX_BANDWIDTH + 1],
            band_diag: [0; MAX_BANDWIDTH + 1],
            band_centers: vec![0; band_width + 1],
            max_real_col: 0,
            max_row: 0,
            max_col: 0,
            band_width,
            wrap_end: 0,
            current_pattern: vec![0; 2 * (MAX_PATTERN_SIZE_CONSTANT + 1)],
            max_real_row: 0,
            max_score: 0,
            rows: 0.0,
            report_min: 0,
        };
        dp.score_matrix[0] = 1; // S[0][0]=1
        dp
    }

    #[inline]
    fn score(&self, row: usize, col: usize) -> i32 {
        self.score_matrix[row * SCORE_MATRIX_STRIDE + col]
    }

    #[inline]
    fn set_score(&mut self, row: usize, col: usize, value: i32) {
        self.score_matrix[row * SCORE_MATRIX_STRIDE + col] = value;
    }
}

impl WraparoundDp {
    /// Mirrors the `test_trace_and_backwards_maxscore`,
    /// `test_maxrowscore_with_match` and `test_maxrowscore_without_match` macros.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn backward_trace_cell(
        &mut self,
        ctx: &DpContext,
        trace: &mut TraceState,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        col: i32,
        pattern_col: i32,
        cell_score: &mut i32,
        left_score: &mut i32,
        with_match: Option<(i32, i32)>,
    ) {
        if real_row <= start - size && *cell_score == 0 {
            *cell_score = -1000;
            self.set_score(row as usize, col as usize, *cell_score);
            *left_score = -1000;
        } else {
            trace.end_of_trace = false;
        }
        if *cell_score >= trace.max_score {
            trace.max_score = *cell_score;
            trace.min_real_row = real_row;
            trace.min_col = pattern_col;
            trace.min_col_band_center = self.band_centers[row as usize];
            trace.min_col_position = col;
        }
        if *cell_score > trace.max_row_score {
            trace.max_row_score = *cell_score;
            let set = match with_match {
                Some((diag_val, match_score)) => {
                    *cell_score == diag_val && match_score == ctx.alpha
                }
                None => false,
            };
            if set {
                trace.match_at_max_col = pattern_col;
            } else {
                trace.match_at_max_col = -2;
            }
        }
    }

    /// Forward half of the banded alignment.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn forward_trace_cell(
        &mut self,
        ctx: &DpContext,
        trace: &mut TraceState,
        real_row: i32,
        start: i32,
        row: i32,
        col: i32,
        pattern_col: i32,
        cell_score: &mut i32,
        left_score: &mut i32,
        with_match: Option<(i32, i32)>,
    ) {
        if real_row >= start && *cell_score == 0 {
            *cell_score = -1000;
            self.set_score(row as usize, col as usize, *cell_score);
            *left_score = -1000;
        } else {
            trace.end_of_trace = false;
        }
        if *cell_score >= trace.max_score {
            trace.max_score = *cell_score;
            trace.max_real_row = real_row;
            trace.max_row = row;
            trace.max_col = pattern_col;
        }
        if *cell_score > trace.max_row_score {
            trace.max_row_score = *cell_score;
            let set = match with_match {
                Some((diag_val, match_score)) => {
                    *cell_score == diag_val && match_score == ctx.alpha
                }
                None => false,
            };
            if set {
                trace.match_at_max_col = pattern_col;
            } else {
                trace.match_at_max_col = -2;
            }
        }
    }

    /// `narrow_band_wraparound_dp()`: banded wraparound dynamic programming.
    #[allow(clippy::too_many_arguments)]
    pub fn narrow_band_wraparound_dp(
        &mut self,
        ctx: &DpContext,
        start: i32,
        size: i32,
        band_radius: i32,
        band_radius_forward: i32,
        with_consensus: bool,
        recenter_criterion: i32,
    ) {
        debug_assert!(2 * band_radius < MAX_BANDWIDTH as i32);
        debug_assert!(2 * band_radius_forward < MAX_BANDWIDTH as i32);
        let delta = ctx.delta;

        // fill EC
        self.fill_current_pattern(ctx, start, size, with_consensus);

        /* ---------------- backward wraparound DP ---------------- */
        let mut trace = TraceState {
            max_score: 0,
            match_at_max_col: -2,
            max_row_score: 0,
            ..Default::default()
        };
        let mut real_row = start + 1;
        let mut row = self.band_width as i32;
        self.band_centers[row as usize] = 0;
        trace.matches_in_diagonal = 0;

        for col in 0..=band_radius {
            let cell_score = delta * (band_radius - col);
            self.set_score(row as usize, col as usize, cell_score);
            self.band_diag[col as usize] = cell_score;
            self.band_up[col as usize] = cell_score + delta;
        }
        for col in (band_radius + 1)..=(2 * band_radius) {
            self.set_score(row as usize, col as usize, -1000);
            self.band_diag[col as usize] = -1000;
            self.band_up[col as usize] = -1000 + delta;
        }

        trace.end_of_trace = false;
        while !trace.end_of_trace && real_row > 1 && row > 0 {
            row -= 1;
            real_row -= 1;
            self.rows += 1.0;
            self.banded_backward_row(
                ctx,
                &mut trace,
                real_row,
                start,
                size,
                row,
                band_radius,
                recenter_criterion,
            );
        }

        /* ---------------- forward wraparound DP ---------------- */
        let mut row = 0i32;
        let mut real_row = trace.min_real_row - 1;

        self.band_centers[0] = (trace.min_col_band_center - 1 + size) % size;
        trace.matches_in_diagonal = 0;
        trace.match_at_max_col = -2;
        trace.max_score = 0;

        let zero_column = trace.min_col_position - band_radius;

        /* widen the narrow band in the forward direction */
        let first_end = std::cmp::max(0, band_radius_forward + zero_column);
        for col in 0..first_end {
            self.set_score(0, col as usize, -1000);
            self.band_diag[col as usize] = -1000;
            self.band_up[col as usize] = -1000 + delta;
        }
        for col in first_end..=(2 * band_radius_forward) {
            let cell_score = delta * (col - (band_radius_forward + zero_column));
            self.set_score(0, col as usize, cell_score);
            self.band_diag[col as usize] = cell_score;
            self.band_up[col as usize] = cell_score + delta;
        }

        trace.end_of_trace = false;
        while !trace.end_of_trace && real_row < ctx.length && row < self.band_width as i32 {
            row += 1;
            real_row += 1;
            self.rows += 1.0;
            self.banded_forward_row(
                ctx,
                &mut trace,
                real_row,
                start,
                size,
                row,
                band_radius_forward,
                recenter_criterion,
            );
        }

        /* store last position scanned with this pattern */
        self.wrap_end = real_row;
        /* test for report */
        if trace.max_score >= self.report_min {
            self.max_real_row = trace.max_real_row;
            self.max_row = trace.max_row;
            self.max_col = trace.max_col;
            self.max_score = trace.max_score;
        }
    }

    /// Fills `current_pattern` (the `EC` working pattern) either with the
    /// consensus or with the sequence window ending at `start`.
    fn fill_current_pattern(
        &mut self,
        ctx: &DpContext,
        start: i32,
        size: i32,
        with_consensus: bool,
    ) {
        if with_consensus {
            for pattern_pos in 0..size {
                self.current_pattern[pattern_pos as usize] =
                    ctx.consensus_pattern[pattern_pos as usize];
            }
        } else {
            for pattern_pos in 0..size {
                self.current_pattern[pattern_pos as usize] =
                    ctx.sequence[(start - size + pattern_pos + 1) as usize];
            }
        }
    }

    /// Copies the finished DP row into the `band_diag` / `band_up` caches so
    /// that the next row can use them as its diagonal / vertical predecessors.
    fn store_banded_row(&mut self, ctx: &DpContext, row: i32, band_radius: i32) {
        let delta = ctx.delta;
        for col in 0..=(2 * band_radius) {
            let cell_score = self.score(row as usize, col as usize);
            self.band_diag[col as usize] = cell_score;
            self.band_up[col as usize] = cell_score + delta;
        }
    }

    /// One row of the backward half: recenter the band, fill the band cells
    /// and update the trace bookkeeping.
    #[allow(clippy::too_many_arguments)]
    fn banded_backward_row(
        &mut self,
        ctx: &DpContext,
        trace: &mut TraceState,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        band_radius: i32,
        recenter_criterion: i32,
    ) {
        trace.max_row_score = -1;
        let last_match_at_max_col = trace.match_at_max_col;
        trace.end_of_trace = true;
        let current_base = ctx.sequence[real_row as usize];

        if trace.matches_in_diagonal >= recenter_criterion {
            self.band_centers[row as usize] = (trace.match_at_max_col - 1 + size) % size;
        } else {
            self.band_centers[row as usize] =
                (self.band_centers[(row + 1) as usize] - 1 + size) % size;
        }

        let mut band_shift =
            (self.band_centers[row as usize] - self.band_centers[(row + 1) as usize] + size) % size;
        if size - band_shift <= band_shift {
            band_shift = -(size - band_shift);
        }
        let pattern_col = (self.band_centers[row as usize] + band_radius) % size;

        if band_shift <= -1 {
            /* band shifts left */
            self.banded_backward_row_shifts_left(
                ctx,
                trace,
                current_base,
                real_row,
                start,
                size,
                row,
                band_radius,
                -band_shift,
                pattern_col,
            );
        } else {
            /* band shifts right */
            self.banded_backward_row_shifts_right(
                ctx,
                trace,
                current_base,
                real_row,
                start,
                size,
                row,
                band_radius,
                band_shift,
                pattern_col,
            );
        }

        self.store_banded_row(ctx, row, band_radius);

        if (trace.match_at_max_col - last_match_at_max_col + size) % size == size - 1 {
            trace.matches_in_diagonal += 1;
        } else {
            trace.matches_in_diagonal = 0;
        }
    }

    /// Backward row cells when the band shifts left.  `band_shift` is the
    /// (positive) magnitude of the shift, `pattern_col` the pattern position
    /// of the rightmost band column.
    #[allow(clippy::too_many_arguments)]
    fn banded_backward_row_shifts_left(
        &mut self,
        ctx: &DpContext,
        mut trace: &mut TraceState,
        current_base: u8,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        band_radius: i32,
        band_shift: i32,
        mut pattern_col: i32,
    ) {
        let delta = ctx.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = 2 * band_radius - band_shift + 1;
        let mut up_index = 2 * band_radius - band_shift;
        let mut col = 2 * band_radius;
        while col >= band_shift {
            let match_score =
                ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                size,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                Some((diag_score, match_score)),
            );
            col -= 1;
            diag_index -= 1;
            up_index -= 1;
            pattern_col = (pattern_col - 1 + size) % size;
        }
        /* col == band_shift-1: diagonal and left only */
        let match_score =
            ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
        self.band_diag[0] += match_score;
        let diag_score = self.band_diag[0];
        let mut cell_score = max3(0, diag_score, left_score);
        self.set_score(row as usize, (band_shift - 1) as usize, cell_score);
        left_score = cell_score + delta;
        self.backward_trace_cell(
            ctx,
            &mut trace,
            real_row,
            start,
            size,
            row,
            band_shift - 1,
            pattern_col,
            &mut cell_score,
            &mut left_score,
            Some((diag_score, match_score)),
        );
        pattern_col = (pattern_col - 1 + size) % size;
        let mut col = band_shift - 2;
        while col >= 0 {
            let mut cell_score = max2(0, left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                size,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                None,
            );
            col -= 1;
            pattern_col = (pattern_col - 1 + size) % size;
        }
    }

    /// Backward row cells when the band shifts right (or stays).  `band_shift`
    /// is the signed shift (>= 0 here), `pattern_col` the pattern position of
    /// the rightmost band column.
    #[allow(clippy::too_many_arguments)]
    fn banded_backward_row_shifts_right(
        &mut self,
        ctx: &DpContext,
        mut trace: &mut TraceState,
        current_base: u8,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        band_radius: i32,
        band_shift: i32,
        mut pattern_col: i32,
    ) {
        let delta = ctx.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = 2 * band_radius;
        let mut up_index = 2 * band_radius;
        let mut col = 2 * band_radius;
        while col >= 2 * band_radius - band_shift + 1 {
            let mut cell_score = max2(0, left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                size,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                None,
            );
            col -= 1;
            pattern_col = (pattern_col - 1 + size) % size;
        }
        /* col == 2*band_radius-band_shift: up and left only */
        let mut cell_score = max3(0, self.band_up[up_index as usize], left_score);
        self.set_score(row as usize, col as usize, cell_score);
        left_score = cell_score + delta;
        self.backward_trace_cell(
            ctx,
            &mut trace,
            real_row,
            start,
            size,
            row,
            col,
            pattern_col,
            &mut cell_score,
            &mut left_score,
            None,
        );
        up_index -= 1;
        col -= 1;
        pattern_col = (pattern_col - 1 + size) % size;
        while col >= 0 {
            let match_score =
                ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                size,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                Some((diag_score, match_score)),
            );
            col -= 1;
            diag_index -= 1;
            up_index -= 1;
            pattern_col = (pattern_col - 1 + size) % size;
        }
    }

    /// One row of the forward half: recenter the band, fill the band cells
    /// and update the trace bookkeeping.
    #[allow(clippy::too_many_arguments)]
    fn banded_forward_row(
        &mut self,
        ctx: &DpContext,
        trace: &mut TraceState,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        band_radius: i32,
        recenter_criterion: i32,
    ) {
        trace.end_of_trace = true;
        trace.max_row_score = -1;
        let last_match_at_max_col = trace.match_at_max_col;
        let current_base = ctx.sequence[real_row as usize];

        if trace.matches_in_diagonal >= recenter_criterion {
            self.band_centers[row as usize] = (trace.match_at_max_col + 1) % size;
        } else {
            self.band_centers[row as usize] = (self.band_centers[(row - 1) as usize] + 1) % size;
        }

        let mut band_shift =
            (self.band_centers[row as usize] - self.band_centers[(row - 1) as usize] + size) % size;
        if size - band_shift <= band_shift {
            band_shift = -(size - band_shift);
        }
        let pattern_col = (self.band_centers[row as usize] - band_radius + size) % size;

        if band_shift >= 1 {
            /* band shifts right */
            self.banded_forward_row_shifts_right(
                ctx,
                trace,
                current_base,
                real_row,
                start,
                size,
                row,
                band_radius,
                band_shift,
                pattern_col,
            );
        } else {
            /* band shifts left */
            self.banded_forward_row_shifts_left(
                ctx,
                trace,
                current_base,
                real_row,
                start,
                size,
                row,
                band_radius,
                -band_shift,
                pattern_col,
            );
        }

        self.store_banded_row(ctx, row, band_radius);

        if (trace.match_at_max_col - last_match_at_max_col + size) % size == 1 {
            trace.matches_in_diagonal += 1;
        } else {
            trace.matches_in_diagonal = 0;
        }
    }

    /// Forward row cells when the band shifts right.  `band_shift` is the
    /// signed shift (>= 1 here), `pattern_col` the pattern position of the
    /// leftmost band column.
    #[allow(clippy::too_many_arguments)]
    fn banded_forward_row_shifts_right(
        &mut self,
        ctx: &DpContext,
        mut trace: &mut TraceState,
        current_base: u8,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        band_radius: i32,
        band_shift: i32,
        mut pattern_col: i32,
    ) {
        let delta = ctx.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = band_shift - 1;
        let mut up_index = band_shift;
        let mut col = 0;
        while col <= 2 * band_radius - band_shift {
            let match_score =
                ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                Some((diag_score, match_score)),
            );
            col += 1;
            diag_index += 1;
            up_index += 1;
            pattern_col = (pattern_col + 1) % size;
        }
        /* col == 2*band_radius-band_shift+1: diagonal and left only */
        let match_score =
            ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
        self.band_diag[diag_index as usize] += match_score;
        let diag_score = self.band_diag[diag_index as usize];
        let mut cell_score = max3(0, diag_score, left_score);
        self.set_score(row as usize, col as usize, cell_score);
        left_score = cell_score + delta;
        self.forward_trace_cell(
            ctx,
            &mut trace,
            real_row,
            start,
            row,
            col,
            pattern_col,
            &mut cell_score,
            &mut left_score,
            Some((diag_score, match_score)),
        );
        col += 1;
        pattern_col = (pattern_col + 1) % size;
        while col <= 2 * band_radius {
            let mut cell_score = max2(0, left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                None,
            );
            col += 1;
            pattern_col = (pattern_col + 1) % size;
        }
    }

    /// Forward row cells when the band shifts left (or stays).  `band_shift`
    /// is the (positive) magnitude of the shift, `pattern_col` the pattern
    /// position of the leftmost band column.
    #[allow(clippy::too_many_arguments)]
    fn banded_forward_row_shifts_left(
        &mut self,
        ctx: &DpContext,
        mut trace: &mut TraceState,
        current_base: u8,
        real_row: i32,
        start: i32,
        size: i32,
        row: i32,
        band_radius: i32,
        band_shift: i32,
        mut pattern_col: i32,
    ) {
        let delta = ctx.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = 0;
        let mut up_index = 0;
        let mut col = 0;
        while col <= band_shift - 1 {
            let mut cell_score = max2(0, left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                None,
            );
            col += 1;
            pattern_col = (pattern_col + 1) % size;
        }
        /* col == band_shift: up and left only */
        let mut cell_score = max3(0, self.band_up[up_index as usize], left_score);
        self.set_score(row as usize, col as usize, cell_score);
        left_score = cell_score + delta;
        self.forward_trace_cell(
            ctx,
            &mut trace,
            real_row,
            start,
            row,
            col,
            pattern_col,
            &mut cell_score,
            &mut left_score,
            None,
        );
        col += 1;
        up_index += 1;
        pattern_col = (pattern_col + 1) % size;
        while col <= 2 * band_radius {
            let match_score =
                ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
                ctx,
                &mut trace,
                real_row,
                start,
                row,
                col,
                pattern_col,
                &mut cell_score,
                &mut left_score,
                Some((diag_score, match_score)),
            );
            col += 1;
            diag_index += 1;
            up_index += 1;
            pattern_col = (pattern_col + 1) % size;
        }
    }

    /// `wraparound_dp()`: full (non banded) wraparound dynamic programming.
    pub fn wraparound_dp(
        &mut self,
        ctx: &DpContext,
        start: i32,
        size: i32,
        consensus_present: bool,
    ) {
        let delta = ctx.delta;
        let size_minus_one = size - 1;
        let size_minus_two = size_minus_one - 1;

        // fill EC
        self.fill_current_pattern(ctx, start, size, consensus_present);

        /* ---------------- backward wraparound DP ---------------- */
        let mut max_score: i32 = 0;
        let mut real_row = start + 1;
        let mut row = self.band_width as i32;
        let mut min_row: i32 = 0;
        let mut max_row: i32 = 0;
        let mut max_col: i32 = 0;
        let mut max_real_row: i32 = 0;

        for pattern_col in 0..size {
            self.set_score(row as usize, pattern_col as usize, max_score);
            self.band_diag[pattern_col as usize] = max_score;
            self.band_up[pattern_col as usize] = max_score + delta;
        }

        let mut end_of_trace = false;
        while !end_of_trace && real_row > 1 && row > 0 {
            row -= 1;
            real_row -= 1;
            self.rows += 1.0;
            let current_base = ctx.sequence[real_row as usize];

            // first pass
            let mut left_score = delta;
            let mut pattern_col = size_minus_one;
            while pattern_col >= 0 {
                self.band_diag[pattern_col as usize] +=
                    ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
                let cell_score = max4(
                    0,
                    self.band_diag[pattern_col as usize],
                    self.band_up[pattern_col as usize],
                    left_score,
                );
                self.set_score(row as usize, pattern_col as usize, cell_score);
                left_score = cell_score + delta;
                pattern_col -= 1;
            }

            end_of_trace = true;
            // second pass
            let mut pattern_col = size_minus_one;
            while pattern_col >= 0 {
                let mut cell_score = max4(
                    0,
                    self.band_diag[pattern_col as usize],
                    self.band_up[pattern_col as usize],
                    left_score,
                );
                left_score = cell_score + delta;
                self.set_score(row as usize, pattern_col as usize, cell_score);
                self.band_up[pattern_col as usize] = left_score;
                if real_row <= start - std::cmp::max(size, 20) && cell_score == 0 {
                    cell_score = -1000;
                    self.set_score(row as usize, pattern_col as usize, cell_score);
                    self.band_up[pattern_col as usize] = -1000;
                    left_score = -1000;
                } else {
                    // real_row >  start - std::cmp::max(size, 20)  || cell_score != 0
                    end_of_trace = false;
                }
                if cell_score >= max_score {
                    max_score = cell_score;
                    min_row = real_row;
                }
                pattern_col -= 1;
            }

            let mut pattern_col = size_minus_two;
            while pattern_col >= 0 {
                self.band_diag[pattern_col as usize] =
                    self.score(row as usize, (pattern_col + 1) as usize);
                pattern_col -= 1;
            }
            self.band_diag[size_minus_one as usize] = self.score(row as usize, 0);
        }

        /***********************************/
        let mut row = 0i32;
        let mut real_row = start - 1;

        for pattern_col in 0..size {
            self.set_score(0, pattern_col as usize, max_score);
            self.band_diag[pattern_col as usize] = max_score;
            self.band_up[pattern_col as usize] = max_score + delta;
        }

        end_of_trace = false;
        while !end_of_trace && real_row < ctx.length && row < self.band_width as i32 {
            row += 1;
            real_row += 1;
            self.rows += 1.0;
            let current_base = ctx.sequence[real_row as usize];

            // first pass
            let mut left_score = delta;
            let mut pattern_col = 0;
            while pattern_col < size {
                self.band_diag[pattern_col as usize] +=
                    ctx.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
                let cell_score = max4(
                    0,
                    self.band_diag[pattern_col as usize],
                    self.band_up[pattern_col as usize],
                    left_score,
                );
                self.set_score(row as usize, pattern_col as usize, cell_score);
                left_score = cell_score + delta;
                pattern_col += 1;
            }

            end_of_trace = true;
            // second pass
            let mut pattern_col = 0;
            while pattern_col < size {
                let mut cell_score = max4(
                    0,
                    self.band_diag[pattern_col as usize],
                    self.band_up[pattern_col as usize],
                    left_score,
                );
                left_score = cell_score + delta;
                self.set_score(row as usize, pattern_col as usize, cell_score);
                self.band_up[pattern_col as usize] = left_score;
                if real_row >= start && cell_score == 0 {
                    cell_score = -1000;
                    self.set_score(row as usize, pattern_col as usize, cell_score);
                    self.band_up[pattern_col as usize] = -1000;
                    left_score = -1000;
                } else {
                    end_of_trace = false;
                }
                if cell_score >= max_score {
                    max_score = cell_score;
                    max_row = real_row;
                }
                pattern_col += 1;
            }

            let mut pattern_col = 0;
            while pattern_col < size {
                self.band_diag[(pattern_col + 1) as usize] =
                    self.score(row as usize, pattern_col as usize);
                pattern_col += 1;
            }
            self.band_diag[0] = self.band_diag[size as usize];
        }

        /* store last position scanned with this pattern */
        self.wrap_end = real_row;

        /* redo for an accurate alignment if score is big enough to report */
        self.report_min = 0;
        if max_score >= self.report_min {
            let mut row = 0i32;
            let start_plus_one = start + 1;
            let mut real_row = min_row - size - 1;
            if real_row < 0 {
                real_row = 0;
            }
            max_score = 0;

            for pattern_col in 0..size {
                self.set_score(0, pattern_col as usize, max_score);
                self.band_diag[pattern_col as usize] = max_score;
                self.band_up[pattern_col as usize] = max_score + delta;
            }

            end_of_trace = false;
            while !end_of_trace && real_row < ctx.length && row < self.band_width as i32 {
                row += 1;
                real_row += 1;
                self.rows += 1.0;
                let current_base = ctx.sequence[real_row as usize];

                let mut left_score = delta;
                let mut pattern_col = 0;
                while pattern_col < size {
                    self.band_diag[pattern_col as usize] += ctx
                        .similarity_score(current_base, self.current_pattern[pattern_col as usize]);
                    let cell_score = max4(
                        0,
                        self.band_diag[pattern_col as usize],
                        self.band_up[pattern_col as usize],
                        left_score,
                    );
                    self.set_score(row as usize, pattern_col as usize, cell_score);
                    left_score = cell_score + delta;
                    pattern_col += 1;
                }

                end_of_trace = true;
                let mut pattern_col = 0;
                while pattern_col < size {
                    let mut cell_score = max4(
                        0,
                        self.band_diag[pattern_col as usize],
                        self.band_up[pattern_col as usize],
                        left_score,
                    );
                    left_score = cell_score + delta;
                    self.set_score(row as usize, pattern_col as usize, cell_score);
                    self.band_up[pattern_col as usize] = left_score;
                    if real_row >= start_plus_one && cell_score == 0 {
                        cell_score = -1000;
                        self.set_score(row as usize, pattern_col as usize, cell_score);
                        self.band_up[pattern_col as usize] = -1000;
                        left_score = -1000;
                    } else {
                        end_of_trace = false;
                    }
                    if cell_score > max_score {
                        max_score = cell_score;
                        max_real_row = real_row;
                        max_row = row;
                        max_col = pattern_col;
                    }
                    pattern_col += 1;
                }

                let mut pattern_col = 0;
                while pattern_col < size {
                    self.band_diag[(pattern_col + 1) as usize] =
                        self.score(row as usize, pattern_col as usize);
                    pattern_col += 1;
                }
                self.band_diag[0] = self.band_diag[size as usize];
            }
        }

        if max_score >= self.report_min {
            self.max_real_row = max_real_row;
            self.max_row = max_row;
            self.max_col = max_col;
            self.max_score = max_score;
        }
    }

    #[inline]
    fn score_at(&self, row: i32, col: i32) -> i32 {
        if row < 0 || col < 0 {
            return 0;
        }
        self.score(row as usize, col as usize)
    }

    /// `get_pair_alignment_with_copynumber()`: traceback for the small distances.
    pub fn get_pair_alignment_with_copynumber(
        &mut self,
        ctx: &DpContext,
        align_pair: &mut AlignPair,
        copynumber: &mut f64,
        size: i32,
    ) -> Vec<String> {
        let mut diagnostics = Vec::new();
        let mut real_row = self.max_real_row;
        let mut row = self.max_row;
        let mut pattern_col = self.max_col;
        let size_minus_one = size - 1;
        let fullcopy = (pattern_col + 1) % size;

        align_pair.score = self.score_at(row, pattern_col);
        let mut length = 0i32;
        *copynumber = 0.0;

        loop {
            let cell_score = self.score_at(row, pattern_col);
            let mut done = false;
            if cell_score <= 0 {
                done = true;
            } else if cell_score
                == self.score_at(row - 1, (pattern_col + size_minus_one) % size)
                    + ctx.similarity_score(
                        ctx.sequence[real_row as usize],
                        self.current_pattern[pattern_col as usize],
                    )
            {
                length += 1;
                if pattern_col == fullcopy {
                    *copynumber += 1.0;
                }
                align_pair.fill(
                    ctx.sequence[real_row as usize],
                    self.current_pattern[pattern_col as usize],
                    length,
                    real_row,
                    pattern_col,
                );
                real_row -= 1;
                row -= 1;
                pattern_col = (pattern_col + size_minus_one) % size;
            } else if cell_score == self.score_at(row - 1, pattern_col) + ctx.delta {
                length += 1;
                align_pair.fill(
                    ctx.sequence[real_row as usize],
                    DASH,
                    length,
                    real_row,
                    (pattern_col + 1) % size,
                );
                real_row -= 1;
                row -= 1;
            } else if cell_score
                == self.score_at(row, (pattern_col + size_minus_one) % size) + ctx.delta
            {
                length += 1;
                if pattern_col == fullcopy {
                    *copynumber += 1.0;
                }
                align_pair.fill(
                    DASH,
                    self.current_pattern[pattern_col as usize],
                    length,
                    real_row + 1,
                    pattern_col,
                );
                pattern_col = (pattern_col + size_minus_one) % size;
            } else {
                diagnostics
                    .push("\nget_pair_alignment_with_copynumber: error in trace back".to_string());
                diagnostics.push(format!("\nrow={}  column={}", real_row, pattern_col));
                diagnostics.push(format!(
                    "\nS={}  Sleft={}  Sup={}  Sdiag={}  match={}",
                    cell_score,
                    self.score_at(row, (pattern_col + size_minus_one) % size),
                    self.score_at(row - 1, pattern_col),
                    self.score_at(row - 1, (pattern_col + size_minus_one) % size),
                    ctx.similarity_score(
                        ctx.sequence[real_row as usize],
                        self.current_pattern[pattern_col as usize]
                    )
                ));
                done = true;
            }

            if done {
                align_pair.length = length;
                if self.max_col >= pattern_col {
                    *copynumber += (self.max_col - pattern_col) as f64 / size as f64;
                } else {
                    *copynumber += (self.max_col + size - pattern_col) as f64 / size as f64;
                }
                return diagnostics;
            }
        }
    }

    /// `get_narrow_band_pair_alignment_with_copynumber()`.
    pub fn get_narrow_band_pair_alignment_with_copynumber(
        &mut self,
        ctx: &DpContext,
        align_pair: &mut AlignPair,
        copynumber: &mut f64,
        size: i32,
        band_radius: i32,
        option: u8,
    ) -> Vec<String> {
        let mut diagnostics = Vec::new();
        let band_radius = band_radius;
        let mut real_row = self.max_real_row;
        let mut row = self.max_row;
        let mut pattern_col = self.max_col;
        let mut band_shift = (self.max_col - self.band_centers[row as usize] + size) % size;
        if size - band_shift <= band_shift {
            band_shift = -(size - band_shift);
        }
        let mut col = band_radius + band_shift;
        let fullcopy = (self.max_col + 1) % size;

        align_pair.score = self.score_at(row, col);
        let mut length = 0i32;
        *copynumber = 0.0;

        loop {
            /* stop at zeros or -1000 for local; stop at row=0 for global */
            let mut legitimate_zero = false;
            let mut cell_score = self.score_at(row, col);
            if option == LOCAL && cell_score == 0 && row != 0 {
                /* check to see if the LOCAL zero score is a legitimate zero */
                let mut column_shift =
                    (pattern_col - self.band_centers[(row - 1) as usize] + size) % size;
                if size - column_shift <= column_shift {
                    column_shift = -(size - column_shift);
                }
                let up_col = band_radius + column_shift;
                let mut center_shift = (self.band_centers[row as usize]
                    - self.band_centers[(row - 1) as usize]
                    + size)
                    % size;
                if size - center_shift <= center_shift {
                    center_shift = -(size - center_shift);
                }
                let match_score = ctx.similarity_score(
                    ctx.sequence[real_row as usize],
                    self.current_pattern[pattern_col as usize],
                );
                if center_shift >= 1 {
                    /* band shifts right */
                    if col <= 2 * band_radius - center_shift {
                        if cell_score == self.score_at(row - 1, up_col - 1) + match_score
                            || cell_score == self.score_at(row - 1, up_col) + ctx.delta
                            || cell_score == self.score_at(row, col - 1) + ctx.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if col == 2 * band_radius - center_shift + 1 {
                        if cell_score == self.score_at(row - 1, up_col - 1) + match_score
                            || cell_score == self.score_at(row, col - 1) + ctx.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if cell_score == self.score_at(row, col - 1) + ctx.delta {
                        legitimate_zero = true;
                    }
                } else {
                    /* band shifts left */
                    let center_shift = -center_shift;
                    if col <= center_shift - 1 {
                        if cell_score == self.score_at(row, col - 1) + ctx.delta {
                            legitimate_zero = true;
                        }
                    } else if col == center_shift {
                        if cell_score == self.score_at(row - 1, up_col) + ctx.delta
                            || cell_score == self.score_at(row, col - 1) + ctx.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if cell_score == self.score_at(row - 1, up_col - 1) + match_score
                        || cell_score == self.score_at(row - 1, up_col) + ctx.delta
                        || cell_score == self.score_at(row, col - 1) + ctx.delta
                    {
                        legitimate_zero = true;
                    }
                }
            }

            /* Stop at zeros or -1000 for local, at row=0 for global.  The
            reference implementation also tests `GLOBAL && r == 0 && c == -1`,
            which is already implied by the leading `row == 0`. */
            let stop = row == 0 || (option == LOCAL && cell_score <= 0 && !legitimate_zero);
            if stop {
                align_pair.length = length;
                if self.max_col >= pattern_col {
                    *copynumber += (self.max_col - pattern_col) as f64 / size as f64;
                } else {
                    *copynumber += (self.max_col + size - pattern_col) as f64 / size as f64;
                }
                return diagnostics;
            }

            let mut column_shift =
                (pattern_col - self.band_centers[(row - 1) as usize] + size) % size;
            if size - column_shift <= column_shift {
                column_shift = -(size - column_shift);
            }
            let up_col = band_radius + column_shift;

            /* try the three predecessors: diagonal, up, left */
            let mut advanced = false;
            if !(option == GLOBAL && pattern_col == -1) {
                let mut center_shift = (self.band_centers[row as usize]
                    - self.band_centers[(row - 1) as usize]
                    + size)
                    % size;
                if size - center_shift <= center_shift {
                    center_shift = -(size - center_shift);
                }
                let shifts_right = center_shift >= 1;
                let shift_magnitude = if shifts_right {
                    center_shift
                } else {
                    -center_shift
                };

                let try_diag = if shifts_right {
                    col <= 2 * band_radius - shift_magnitude
                        || col == 2 * band_radius - shift_magnitude + 1
                } else {
                    col > shift_magnitude
                };
                let try_up = if shifts_right {
                    col <= 2 * band_radius - shift_magnitude
                } else {
                    col >= shift_magnitude
                };
                /* the left branch is always tested */
                let try_left = true;

                cell_score = self.score_at(row, col);
                if try_diag
                    && cell_score
                        == self.score_at(row - 1, up_col - 1)
                            + ctx.similarity_score(
                                ctx.sequence[real_row as usize],
                                self.current_pattern[pattern_col as usize],
                            )
                {
                    length += 1;
                    if pattern_col == fullcopy {
                        *copynumber += 1.0;
                    }
                    if option == LOCAL {
                        align_pair.fill(
                            ctx.sequence[real_row as usize],
                            self.current_pattern[pattern_col as usize],
                            length,
                            real_row,
                            pattern_col,
                        );
                        pattern_col = (pattern_col - 1 + size) % size;
                    } else {
                        align_pair.fill(
                            ctx.sequence[real_row as usize],
                            self.current_pattern[pattern_col as usize],
                            length,
                            real_row,
                            pattern_col + self.max_real_col - self.max_col,
                        );
                        pattern_col -= 1;
                    }
                    real_row -= 1;
                    row -= 1;
                    col = up_col - 1;
                    advanced = true;
                }

                if !advanced && try_up && cell_score == self.score_at(row - 1, up_col) + ctx.delta {
                    length += 1;
                    if option == LOCAL {
                        align_pair.fill(
                            ctx.sequence[real_row as usize],
                            DASH,
                            length,
                            real_row,
                            (pattern_col + 1) % size,
                        );
                    } else {
                        align_pair.fill(
                            ctx.sequence[real_row as usize],
                            DASH,
                            length,
                            real_row,
                            pattern_col + 1 + self.max_real_col - self.max_col,
                        );
                    }
                    real_row -= 1;
                    row -= 1;
                    col = up_col;
                    advanced = true;
                }

                if !advanced && try_left {
                    if col == 0 {
                        diagnostics.push(
                            "\nget_pair_alignment_with_copynumber: error in trace back".to_string(),
                        );
                        diagnostics
                            .push("\nattempted to compute left branch when i==0".to_string());
                        diagnostics.push(format!("\nS[{}][{}]={}", row, col, cell_score));
                        break;
                    } else if cell_score == self.score_at(row, col - 1) + ctx.delta {
                        length += 1;
                        if pattern_col == fullcopy {
                            *copynumber += 1.0;
                        }
                        if option == LOCAL {
                            align_pair.fill(
                                DASH,
                                self.current_pattern[pattern_col as usize],
                                length,
                                real_row + 1,
                                pattern_col,
                            );
                            pattern_col = (pattern_col - 1 + size) % size;
                        } else {
                            align_pair.fill(
                                DASH,
                                self.current_pattern[pattern_col as usize],
                                length,
                                real_row + 1,
                                pattern_col + self.max_real_col - self.max_col,
                            );
                            pattern_col -= 1;
                        }
                        col -= 1;
                        advanced = true;
                    }
                }
            } else {
                /* global alignment past the left edge: the only option is up */
                cell_score = self.score_at(row, col);
                if cell_score == self.score_at(row - 1, up_col) + ctx.delta {
                    length += 1;
                    align_pair.fill(
                        ctx.sequence[real_row as usize],
                        DASH,
                        length,
                        real_row,
                        (pattern_col + 1) % size,
                    );
                    real_row -= 1;
                    row -= 1;
                    col = up_col;
                    advanced = true;
                }
            }

            if !advanced {
                diagnostics
                    .push("\nget_pair_alignment_with_copynumber: error in trace back".to_string());
                diagnostics.push(format!(
                    "\nS: row={}  column={}  upi={}  Sequence: realrow={}  EC: realcol={}",
                    row, col, up_col, real_row, pattern_col
                ));
                break;
            }
        }
        diagnostics
    }
}

impl Trf {
    pub fn wraparound_dp(&mut self, start: i32, size: i32, consensus_present: bool) {
        let ctx = DpContext {
            sequence: &self.sequence,
            length: self.length,
            alpha: self.alpha,
            delta: self.delta,
            similarity_matrix: &self.similarity_matrix,
            consensus_pattern: &self.consensus.pattern,
        };
        self.dp.wraparound_dp(&ctx, start, size, consensus_present);
    }

    pub fn narrow_band_wraparound_dp(
        &mut self,
        start: i32,
        size: i32,
        band_radius: i32,
        band_radius_forward: i32,
        with_consensus: bool,
        recenter_criterion: i32,
    ) {
        if (MAX_BANDWIDTH as i32) < 2 * band_radius + 1 {
            let text = format!(
                "\nIn narrow_band_wraparound_dp, MAXBANDWIDTH: {} exceeded by 2*w+1: {}\n",
                MAX_BANDWIDTH,
                2 * band_radius + 1
            );
            self.write_to_alignment_file(&text);
            self.fatal(-1);
        }
        if (MAX_BANDWIDTH as i32) < 2 * band_radius_forward + 1 {
            let text = format!(
                "\nIn narrow_band_wraparound_dp, MAXBANDWIDTH: {} exceeded by bandradiusforward 2*w+1: {}\n",
                MAX_BANDWIDTH,
                2 * band_radius_forward + 1
            );
            self.write_to_alignment_file(&text);
            self.fatal(-1);
        }

        let ctx = DpContext {
            sequence: &self.sequence,
            length: self.length,
            alpha: self.alpha,
            delta: self.delta,
            similarity_matrix: &self.similarity_matrix,
            consensus_pattern: &self.consensus.pattern,
        };
        self.dp.narrow_band_wraparound_dp(
            &ctx,
            start,
            size,
            band_radius,
            band_radius_forward,
            with_consensus,
            recenter_criterion,
        );
    }

    pub fn get_pair_alignment_with_copynumber(&mut self, size: i32) {
        let ctx = DpContext {
            sequence: &self.sequence,
            length: self.length,
            alpha: self.alpha,
            delta: self.delta,
            similarity_matrix: &self.similarity_matrix,
            consensus_pattern: &self.consensus.pattern,
        };
        let diagnostics = self.dp.get_pair_alignment_with_copynumber(
            &ctx,
            &mut self.align_pair,
            &mut self.copynumber,
            size,
        );
        for text in diagnostics {
            self.write_to_alignment_file(&text);
        }
    }

    pub fn get_narrow_band_pair_alignment_with_copynumber(
        &mut self,
        size: i32,
        band_radius: i32,
        option: u8,
    ) {
        let ctx = DpContext {
            sequence: &self.sequence,
            length: self.length,
            alpha: self.alpha,
            delta: self.delta,
            similarity_matrix: &self.similarity_matrix,
            consensus_pattern: &self.consensus.pattern,
        };
        let diagnostics = self.dp.get_narrow_band_pair_alignment_with_copynumber(
            &ctx,
            &mut self.align_pair,
            &mut self.copynumber,
            size,
            band_radius,
            option,
        );
        for text in diagnostics {
            self.write_to_alignment_file(&text);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::FastaSequence;
    use crate::params::ParamSet;
    use crate::state::{MIN_BAND_RADIUS, RECENTER_CRITERION};

    fn scoring_trf(sequence: &[u8]) -> Trf {
        let fasta = FastaSequence::from_bytes(b"test", sequence);
        let params = ParamSet {
            match_weight: 2,
            mismatch: 5,
            indel: 7,
            ..Default::default()
        };
        Trf::new(params, &fasta)
    }

    #[test]
    fn new_allocates_workspace_and_sets_origin() {
        let dp = WraparoundDp::new(7);
        assert_eq!(dp.score_matrix.len(), 8 * SCORE_MATRIX_STRIDE);
        assert_eq!(dp.score(0, 0), 1);
        assert_eq!(dp.band_centers.len(), 8);
        assert_eq!(
            dp.current_pattern.len(),
            2 * (MAX_PATTERN_SIZE_CONSTANT + 1)
        );
        assert_eq!(dp.rows, 0.0);
        assert_eq!(dp.max_score, 0);
        assert_eq!(dp.max_real_row, 0);
    }

    #[test]
    fn dp_context_reads_similarity_matrix_row_major() {
        let mut matrix = vec![0i32; 256 * 256];
        matrix[256 * b'A' as usize + b'C' as usize] = 42;
        let ctx = DpContext {
            sequence: &[0],
            length: 0,
            alpha: 0,
            delta: 0,
            similarity_matrix: &matrix,
            consensus_pattern: &[],
        };
        assert_eq!(ctx.similarity_score(b'A', b'C'), 42);
        assert_eq!(ctx.similarity_score(b'C', b'A'), 0);
    }

    #[test]
    fn wraparound_dp_and_traceback_smoke() {
        let mut trf = scoring_trf(&[b'A'; 30]);
        trf.wraparound_dp(15, 5, false);
        assert!(trf.dp.max_score > 0);
        assert!(trf.dp.rows > 0.0);
        assert!(trf.dp.wrap_end <= trf.length);

        let ctx = DpContext {
            sequence: &trf.sequence,
            length: trf.length,
            alpha: trf.alpha,
            delta: trf.delta,
            similarity_matrix: &trf.similarity_matrix,
            consensus_pattern: &trf.consensus.pattern,
        };
        let diagnostics = trf.dp.get_pair_alignment_with_copynumber(
            &ctx,
            &mut trf.align_pair,
            &mut trf.copynumber,
            5,
        );
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        assert!(trf.align_pair.length > 0);
        assert!(trf.copynumber > 0.0);
    }

    #[test]
    fn narrow_band_dp_and_traceback_smoke() {
        let mut trf = scoring_trf(&[b'A'; 200]);
        let size = 25;
        let band_radius = std::cmp::max(MIN_BAND_RADIUS, trf.indel_random_walk_range(size));
        let band_radius_forward = std::cmp::min(2 * band_radius, size / 3);
        trf.narrow_band_wraparound_dp(
            100,
            size,
            band_radius,
            band_radius_forward,
            false,
            RECENTER_CRITERION,
        );
        assert!(trf.dp.max_score > 0);
        assert!(trf.dp.rows > 0.0);

        let trace_band_radius = std::cmp::min(2 * band_radius, size / 3);
        let ctx = DpContext {
            sequence: &trf.sequence,
            length: trf.length,
            alpha: trf.alpha,
            delta: trf.delta,
            similarity_matrix: &trf.similarity_matrix,
            consensus_pattern: &trf.consensus.pattern,
        };
        let diagnostics = trf.dp.get_narrow_band_pair_alignment_with_copynumber(
            &ctx,
            &mut trf.align_pair,
            &mut trf.copynumber,
            size,
            trace_band_radius,
            crate::state::LOCAL,
        );
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        assert!(trf.align_pair.length > 0);
        assert!(trf.copynumber > 0.0);
    }
}
