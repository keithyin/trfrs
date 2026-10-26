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

//! Wraparound dynamic programming, traceback and consensus computation.
//!
//! This is a direct translation of `narrow_band_wraparound_dp()`, `wraparound_dp()`,
//! `get_pair_alignment_with_copynumber()`,
//! `get_narrow_band_pair_alignment_with_copynumber()`, `get_consensus()` and the
//! alignment printing routines from the reference implementation.

use crate::state::{Trf, DASH, GLOBAL, LOCAL, MAX_BANDWIDTH};

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

impl Trf {
    /// Mirrors the `test_trace_and_backwards_maxscore`,
    /// `test_maxrowscore_with_match` and `test_maxrowscore_without_match` macros.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn backward_trace_cell(
        &mut self,
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
            trace.min_col_band_center = self.band_center[row as usize];
            trace.min_col_position = col;
        }
        if *cell_score > trace.max_row_score {
            trace.max_row_score = *cell_score;
            let set = match with_match {
                Some((diag_val, match_score)) => *cell_score == diag_val && match_score == self.alpha,
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
                Some((diag_val, match_score)) => *cell_score == diag_val && match_score == self.alpha,
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
    pub fn narrow_band_wraparound_dp(
        &mut self,
        start: i32,
        size: i32,
        band_radius: i32,
        band_radius_forward: i32,
        with_consensus: bool,
        recenter_criterion: i32,
    ) {
        let delta = self.delta;
        if (MAX_BANDWIDTH as i32) < 2 * band_radius + 1 {
            let text = format!(
                "\nIn narrow_band_wraparound_dp, MAXBANDWIDTH: {} exceeded by 2*w+1: {}\n",
                MAX_BANDWIDTH,
                2 * band_radius + 1
            );
            self.write_to_alignment_file(&text);
            self.fatal(-1);
        }

        // fill EC
        self.fill_current_pattern(start, size, with_consensus);

        /* ---------------- backward wraparound DP ---------------- */
        let mut trace = TraceState {
            max_score: 0,
            match_at_max_col: -2,
            max_row_score: 0,
            ..Default::default()
        };
        let mut real_row = start + 1;
        let mut row = self.band_width as i32;
        self.band_center[row as usize] = 0;
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

        self.band_center[0] = (trace.min_col_band_center - 1 + size) % size;
        trace.matches_in_diagonal = 0;
        trace.match_at_max_col = -2;
        trace.max_score = 0;

        let zero_column = trace.min_col_position - band_radius;

        /* widen the narrow band in the forward direction */
        if (MAX_BANDWIDTH as i32) < 2 * band_radius_forward + 1 {
            let text = format!(
                "\nIn narrow_band_wraparound_dp, MAXBANDWIDTH: {} exceeded by bandradiusforward 2*w+1: {}\n",
                MAX_BANDWIDTH,
                2 * band_radius_forward + 1
            );
            self.write_to_alignment_file(&text);
            self.fatal(-1);
        }

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
        while !trace.end_of_trace && real_row < self.length && row < self.band_width as i32 {
            row += 1;
            real_row += 1;
            self.rows += 1.0;
            self.banded_forward_row(
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
    fn fill_current_pattern(&mut self, start: i32, size: i32, with_consensus: bool) {
        if with_consensus {
            for pattern_pos in 0..size {
                self.current_pattern[pattern_pos as usize] =
                    self.consensus.pattern[pattern_pos as usize];
            }
        } else {
            for pattern_pos in 0..size {
                self.current_pattern[pattern_pos as usize] =
                    self.sequence[(start - size + pattern_pos + 1) as usize];
            }
        }
    }

    /// Copies the finished DP row into the `band_diag` / `band_up` caches so
    /// that the next row can use them as its diagonal / vertical predecessors.
    fn store_banded_row(&mut self, row: i32, band_radius: i32) {
        let delta = self.delta;
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
        let current_base = self.sequence[real_row as usize];

        if trace.matches_in_diagonal >= recenter_criterion {
            self.band_center[row as usize] = (trace.match_at_max_col - 1 + size) % size;
        } else {
            self.band_center[row as usize] =
                (self.band_center[(row + 1) as usize] - 1 + size) % size;
        }

        let mut band_shift =
            (self.band_center[row as usize] - self.band_center[(row + 1) as usize] + size) % size;
        if size - band_shift <= band_shift {
            band_shift = -(size - band_shift);
        }
        let pattern_col = (self.band_center[row as usize] + band_radius) % size;

        if band_shift <= -1 {
            /* band shifts left */
            self.banded_backward_row_shifts_left(
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

        self.store_banded_row(row, band_radius);

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
        let delta = self.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = 2 * band_radius - band_shift + 1;
        let mut up_index = 2 * band_radius - band_shift;
        let mut col = 2 * band_radius;
        while col >= band_shift {
            let match_score = self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
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
        let match_score = self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
        self.band_diag[0] += match_score;
        let diag_score = self.band_diag[0];
        let mut cell_score = max3(0, diag_score, left_score);
        self.set_score(row as usize, (band_shift - 1) as usize, cell_score);
        left_score = cell_score + delta;
        self.backward_trace_cell(
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
                &mut trace, real_row, start, size, row, col, pattern_col, &mut cell_score, &mut left_score, None,
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
        let delta = self.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = 2 * band_radius;
        let mut up_index = 2 * band_radius;
        let mut col = 2 * band_radius;
        while col >= 2 * band_radius - band_shift + 1 {
            let mut cell_score = max2(0, left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
                &mut trace, real_row, start, size, row, col, pattern_col, &mut cell_score, &mut left_score, None,
            );
            col -= 1;
            pattern_col = (pattern_col - 1 + size) % size;
        }
        /* col == 2*band_radius-band_shift: up and left only */
        let mut cell_score = max3(0, self.band_up[up_index as usize], left_score);
        self.set_score(row as usize, col as usize, cell_score);
        left_score = cell_score + delta;
        self.backward_trace_cell(
            &mut trace, real_row, start, size, row, col, pattern_col, &mut cell_score, &mut left_score, None,
        );
        up_index -= 1;
        col -= 1;
        pattern_col = (pattern_col - 1 + size) % size;
        while col >= 0 {
            let match_score = self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.backward_trace_cell(
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
        let current_base = self.sequence[real_row as usize];

        if trace.matches_in_diagonal >= recenter_criterion {
            self.band_center[row as usize] = (trace.match_at_max_col + 1) % size;
        } else {
            self.band_center[row as usize] = (self.band_center[(row - 1) as usize] + 1) % size;
        }

        let mut band_shift =
            (self.band_center[row as usize] - self.band_center[(row - 1) as usize] + size) % size;
        if size - band_shift <= band_shift {
            band_shift = -(size - band_shift);
        }
        let pattern_col = (self.band_center[row as usize] - band_radius + size) % size;

        if band_shift >= 1 {
            /* band shifts right */
            self.banded_forward_row_shifts_right(
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

        self.store_banded_row(row, band_radius);

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
        let delta = self.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = band_shift - 1;
        let mut up_index = band_shift;
        let mut col = 0;
        while col <= 2 * band_radius - band_shift {
            let match_score = self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
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
        let match_score = self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
        self.band_diag[diag_index as usize] += match_score;
        let diag_score = self.band_diag[diag_index as usize];
        let mut cell_score = max3(0, diag_score, left_score);
        self.set_score(row as usize, col as usize, cell_score);
        left_score = cell_score + delta;
        self.forward_trace_cell(
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
                &mut trace, real_row, start, row, col, pattern_col, &mut cell_score, &mut left_score, None,
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
        let delta = self.delta;
        let mut left_score: i32 = -1000;
        let mut diag_index = 0;
        let mut up_index = 0;
        let mut col = 0;
        while col <= band_shift - 1 {
            let mut cell_score = max2(0, left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
                &mut trace, real_row, start, row, col, pattern_col, &mut cell_score, &mut left_score, None,
            );
            col += 1;
            pattern_col = (pattern_col + 1) % size;
        }
        /* col == band_shift: up and left only */
        let mut cell_score = max3(0, self.band_up[up_index as usize], left_score);
        self.set_score(row as usize, col as usize, cell_score);
        left_score = cell_score + delta;
        self.forward_trace_cell(&mut trace, real_row, start, row, col, pattern_col, &mut cell_score, &mut left_score, None);
        col += 1;
        up_index += 1;
        pattern_col = (pattern_col + 1) % size;
        while col <= 2 * band_radius {
            let match_score = self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
            self.band_diag[diag_index as usize] += match_score;
            let diag_score = self.band_diag[diag_index as usize];
            let mut cell_score = max4(0, diag_score, self.band_up[up_index as usize], left_score);
            self.set_score(row as usize, col as usize, cell_score);
            left_score = cell_score + delta;
            self.forward_trace_cell(
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
    pub fn wraparound_dp(&mut self, start: i32, size: i32, consensus_present: bool) {
        let delta = self.delta;
        let size_minus_one = size - 1;
        let size_minus_two = size_minus_one - 1;

        // fill EC
        self.fill_current_pattern(start, size, consensus_present);

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
            let current_base = self.sequence[real_row as usize];

            // first pass
            let mut left_score = delta;
            let mut pattern_col = size_minus_one;
            while pattern_col >= 0 {
                self.band_diag[pattern_col as usize] += self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
                let cell_score = max4(0, self.band_diag[pattern_col as usize], self.band_up[pattern_col as usize], left_score);
                self.set_score(row as usize, pattern_col as usize, cell_score);
                left_score = cell_score + delta;
                pattern_col -= 1;
            }

            end_of_trace = true;
            // second pass
            let mut pattern_col = size_minus_one;
            while pattern_col >= 0 {
                let mut cell_score = max4(0, self.band_diag[pattern_col as usize], self.band_up[pattern_col as usize], left_score);
                left_score = cell_score + delta;
                self.set_score(row as usize, pattern_col as usize, cell_score);
                self.band_up[pattern_col as usize] = left_score;
                if real_row <= start - std::cmp::max(size, 20) && cell_score == 0 {
                    cell_score = -1000;
                    self.set_score(row as usize, pattern_col as usize, cell_score);
                    self.band_up[pattern_col as usize] = -1000;
                    left_score = -1000;
                } else {
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
                self.band_diag[pattern_col as usize] = self.score(row as usize, (pattern_col + 1) as usize);
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
        while !end_of_trace && real_row < self.length && row < self.band_width as i32 {
            row += 1;
            real_row += 1;
            self.rows += 1.0;
            let current_base = self.sequence[real_row as usize];

            // first pass
            let mut left_score = delta;
            let mut pattern_col = 0;
            while pattern_col < size {
                self.band_diag[pattern_col as usize] += self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
                let cell_score = max4(0, self.band_diag[pattern_col as usize], self.band_up[pattern_col as usize], left_score);
                self.set_score(row as usize, pattern_col as usize, cell_score);
                left_score = cell_score + delta;
                pattern_col += 1;
            }

            end_of_trace = true;
            // second pass
            let mut pattern_col = 0;
            while pattern_col < size {
                let mut cell_score = max4(0, self.band_diag[pattern_col as usize], self.band_up[pattern_col as usize], left_score);
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
                self.band_diag[(pattern_col + 1) as usize] = self.score(row as usize, pattern_col as usize);
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
            while !end_of_trace && real_row < self.length && row < self.band_width as i32 {
                row += 1;
                real_row += 1;
                self.rows += 1.0;
                let current_base = self.sequence[real_row as usize];

                let mut left_score = delta;
                let mut pattern_col = 0;
                while pattern_col < size {
                    self.band_diag[pattern_col as usize] += self.similarity_score(current_base, self.current_pattern[pattern_col as usize]);
                    let cell_score = max4(0, self.band_diag[pattern_col as usize], self.band_up[pattern_col as usize], left_score);
                    self.set_score(row as usize, pattern_col as usize, cell_score);
                    left_score = cell_score + delta;
                    pattern_col += 1;
                }

                end_of_trace = true;
                let mut pattern_col = 0;
                while pattern_col < size {
                    let mut cell_score = max4(0, self.band_diag[pattern_col as usize], self.band_up[pattern_col as usize], left_score);
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
                    self.band_diag[(pattern_col + 1) as usize] = self.score(row as usize, pattern_col as usize);
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
    pub fn get_pair_alignment_with_copynumber(&mut self, size: i32) {
        let mut real_row = self.max_real_row;
        let mut row = self.max_row;
        let mut pattern_col = self.max_col;
        let size_minus_one = size - 1;
        let fullcopy = (pattern_col + 1) % size;

        self.align_pair.score = self.score_at(row, pattern_col);
        let mut length = 0i32;
        self.copynumber = 0.0;

        loop {
            let cell_score = self.score_at(row, pattern_col);
            let mut done = false;
            if cell_score <= 0 {
                done = true;
            } else if cell_score
                == self.score_at(row - 1, (pattern_col + size_minus_one) % size)
                    + self.similarity_score(self.sequence[real_row as usize], self.current_pattern[pattern_col as usize])
            {
                length += 1;
                if pattern_col == fullcopy {
                    self.copynumber += 1.0;
                }
                self.fill_align_pair(self.sequence[real_row as usize], self.current_pattern[pattern_col as usize], length, real_row, pattern_col);
                real_row -= 1;
                row -= 1;
                pattern_col = (pattern_col + size_minus_one) % size;
            } else if cell_score == self.score_at(row - 1, pattern_col) + self.delta {
                length += 1;
                self.fill_align_pair(self.sequence[real_row as usize], DASH, length, real_row, (pattern_col + 1) % size);
                real_row -= 1;
                row -= 1;
            } else if cell_score == self.score_at(row, (pattern_col + size_minus_one) % size) + self.delta {
                length += 1;
                if pattern_col == fullcopy {
                    self.copynumber += 1.0;
                }
                self.fill_align_pair(DASH, self.current_pattern[pattern_col as usize], length, real_row + 1, pattern_col);
                pattern_col = (pattern_col + size_minus_one) % size;
            } else {
                self.write_to_alignment_file("\nget_pair_alignment_with_copynumber: error in trace back");
                let text = format!("\nrow={}  column={}", real_row, pattern_col);
                self.write_to_alignment_file(&text);
                let text = format!(
                    "\nS={}  Sleft={}  Sup={}  Sdiag={}  match={}",
                    cell_score,
                    self.score_at(row, (pattern_col + size_minus_one) % size),
                    self.score_at(row - 1, pattern_col),
                    self.score_at(row - 1, (pattern_col + size_minus_one) % size),
                    self.similarity_score(self.sequence[real_row as usize], self.current_pattern[pattern_col as usize])
                );
                self.write_to_alignment_file(&text);
                done = true;
            }

            if done {
                self.align_pair.length = length;
                if self.max_col >= pattern_col {
                    self.copynumber += (self.max_col - pattern_col) as f64 / size as f64;
                } else {
                    self.copynumber += (self.max_col + size - pattern_col) as f64 / size as f64;
                }
                return;
            }
        }
    }

    /// `get_narrow_band_pair_alignment_with_copynumber()`.
    pub fn get_narrow_band_pair_alignment_with_copynumber(
        &mut self,
        size: i32,
        band_radius: i32,
        option: u8,
    ) {
        let band_radius = band_radius;
        let mut real_row = self.max_real_row;
        let mut row = self.max_row;
        let mut pattern_col = self.max_col;
        let mut band_shift = (self.max_col - self.band_center[row as usize] + size) % size;
        if size - band_shift <= band_shift {
            band_shift = -(size - band_shift);
        }
        let mut col = band_radius + band_shift;
        let fullcopy = (self.max_col + 1) % size;

        self.align_pair.score = self.score_at(row, col);
        let mut length = 0i32;
        self.copynumber = 0.0;

        loop {
            /* stop at zeros or -1000 for local; stop at row=0 for global */
            let mut legitimate_zero = false;
            let mut cell_score = self.score_at(row, col);
            if option == LOCAL && cell_score == 0 && row != 0 {
                /* check to see if the LOCAL zero score is a legitimate zero */
                let mut column_shift = (pattern_col - self.band_center[(row - 1) as usize] + size) % size;
                if size - column_shift <= column_shift {
                    column_shift = -(size - column_shift);
                }
                let up_col = band_radius + column_shift;
                let mut center_shift = (self.band_center[row as usize] - self.band_center[(row - 1) as usize]
                    + size)
                    % size;
                if size - center_shift <= center_shift {
                    center_shift = -(size - center_shift);
                }
                let match_score = self.similarity_score(self.sequence[real_row as usize], self.current_pattern[pattern_col as usize]);
                if center_shift >= 1 {
                    /* band shifts right */
                    if col <= 2 * band_radius - center_shift {
                        if cell_score == self.score_at(row - 1, up_col - 1) + match_score
                            || cell_score == self.score_at(row - 1, up_col) + self.delta
                            || cell_score == self.score_at(row, col - 1) + self.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if col == 2 * band_radius - center_shift + 1 {
                        if cell_score == self.score_at(row - 1, up_col - 1) + match_score
                            || cell_score == self.score_at(row, col - 1) + self.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if cell_score == self.score_at(row, col - 1) + self.delta {
                        legitimate_zero = true;
                    }
                } else {
                    /* band shifts left */
                    let center_shift = -center_shift;
                    if col <= center_shift - 1 {
                        if cell_score == self.score_at(row, col - 1) + self.delta {
                            legitimate_zero = true;
                        }
                    } else if col == center_shift {
                        if cell_score == self.score_at(row - 1, up_col) + self.delta
                            || cell_score == self.score_at(row, col - 1) + self.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if cell_score == self.score_at(row - 1, up_col - 1) + match_score
                        || cell_score == self.score_at(row - 1, up_col) + self.delta
                        || cell_score == self.score_at(row, col - 1) + self.delta
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
                self.align_pair.length = length;
                if self.max_col >= pattern_col {
                    self.copynumber += (self.max_col - pattern_col) as f64 / size as f64;
                } else {
                    self.copynumber += (self.max_col + size - pattern_col) as f64 / size as f64;
                }
                return;
            }

            let mut column_shift = (pattern_col - self.band_center[(row - 1) as usize] + size) % size;
            if size - column_shift <= column_shift {
                column_shift = -(size - column_shift);
            }
            let up_col = band_radius + column_shift;

            /* try the three predecessors: diagonal, up, left */
            let mut advanced = false;
            if !(option == GLOBAL && pattern_col == -1) {
                let mut center_shift = (self.band_center[row as usize] - self.band_center[(row - 1) as usize]
                    + size)
                    % size;
                if size - center_shift <= center_shift {
                    center_shift = -(size - center_shift);
                }
                let shifts_right = center_shift >= 1;
                let shift_magnitude = if shifts_right { center_shift } else { -center_shift };

                let try_diag = if shifts_right {
                    col <= 2 * band_radius - shift_magnitude || col == 2 * band_radius - shift_magnitude + 1
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
                            + self.similarity_score(self.sequence[real_row as usize], self.current_pattern[pattern_col as usize])
                {
                    length += 1;
                    if pattern_col == fullcopy {
                        self.copynumber += 1.0;
                    }
                    if option == LOCAL {
                        self.fill_align_pair(
                            self.sequence[real_row as usize],
                            self.current_pattern[pattern_col as usize],
                            length,
                            real_row,
                            pattern_col,
                        );
                        pattern_col = (pattern_col - 1 + size) % size;
                    } else {
                        self.fill_align_pair(
                            self.sequence[real_row as usize],
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

                if !advanced && try_up && cell_score == self.score_at(row - 1, up_col) + self.delta {
                    length += 1;
                    if option == LOCAL {
                        self.fill_align_pair(
                            self.sequence[real_row as usize],
                            DASH,
                            length,
                            real_row,
                            (pattern_col + 1) % size,
                        );
                    } else {
                        self.fill_align_pair(
                            self.sequence[real_row as usize],
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
                        self.write_to_alignment_file("\nget_pair_alignment_with_copynumber: error in trace back");
                        self.write_to_alignment_file("\nattempted to compute left branch when i==0");
                        let text = format!("\nS[{}][{}]={}", row, col, cell_score);
                        self.write_to_alignment_file(&text);
                        break;
                    } else if cell_score == self.score_at(row, col - 1) + self.delta {
                        length += 1;
                        if pattern_col == fullcopy {
                            self.copynumber += 1.0;
                        }
                        if option == LOCAL {
                            self.fill_align_pair(DASH, self.current_pattern[pattern_col as usize], length, real_row + 1, pattern_col);
                            pattern_col = (pattern_col - 1 + size) % size;
                        } else {
                            self.fill_align_pair(
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
                if cell_score == self.score_at(row - 1, up_col) + self.delta {
                    length += 1;
                    self.fill_align_pair(
                        self.sequence[real_row as usize],
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
                self.write_to_alignment_file("\nget_pair_alignment_with_copynumber: error in trace back");
                let text = format!(
                    "\nS: row={}  column={}  upi={}  Sequence: realrow={}  EC: realcol={}",
                    row, col, up_col, real_row, pattern_col
                );
                self.write_to_alignment_file(&text);
                break;
            }
        }
    }

    /// `reverse()`: reverses the alignment in `AlignPair`.
    pub fn reverse(&mut self) {
        let length = self.align_pair.length;
        let mut front = 1;
        while front <= length / 2 {
            let back = length - front + 1;
            self.align_pair.sequence_text.swap(front as usize, back as usize);
            self.align_pair.pattern_text.swap(front as usize, back as usize);
            self.align_pair.sequence_index.swap(front as usize, back as usize);
            self.align_pair.pattern_index.swap(front as usize, back as usize);
            front += 1;
        }
    }

    /// `shift_pattern_indices()`: makes the pattern indices run from 0 to size-1.
    pub fn shift_pattern_indices(&mut self, pattern_size: i32) {
        let downshift = self.align_pair.pattern_index[1];
        let upshift = pattern_size - downshift;
        for column in 1..=self.align_pair.length {
            if self.align_pair.pattern_index[column as usize] >= downshift {
                self.align_pair.pattern_index[column as usize] -= downshift;
            } else {
                self.align_pair.pattern_index[column as usize] += upshift;
            }
        }
    }

    /// `print_alignment()`.
    pub fn print_alignment(&mut self, pattern_width: i32) {
        self.shift_pattern_indices(pattern_width);
        let print_width = 75i32;
        let mut out = String::new();

        if print_width <= 0 {
            out.push_str("Error, pwidth<=0, can't report alignments.\n");
            self.write_to_alignment_file(&out);
            return;
        }

        let align_length = self.align_pair.length;
        if self.align_pair.sequence_index[1] != 1 {
            let mut display_start = self.align_pair.sequence_index[1] - 10;
            if display_start < 1 {
                display_start = 1;
            }
            let mut col = display_start;
            out.push_str(&format!("  {:9} ", col));
            for _ in 1..=10 {
                out.push(self.sequence[col as usize] as char);
                col += 1;
                if col == self.align_pair.sequence_index[1] {
                    break;
                }
            }
            out.push_str("\n\n");
        }

        let mut block_start = 1i32;
        loop {
            let mut col = block_start;
            out.push_str("            ");
            let mut first = true;
            let mut segment_count = 0;
            for _ in 0..print_width - 10 {
                if col == align_length + 1 {
                    break;
                }
                if !first
                    && self.align_pair.pattern_index[col as usize] == 0
                    && self.align_pair.pattern_index[(col - 1) as usize] != 0
                {
                    if pattern_width > 6 || print_width - 10 - segment_count < 2 * pattern_width {
                        break;
                    }
                    out.push(' ');
                    segment_count += 1;
                }
                first = false;
                if self.align_pair.sequence_text[col as usize] != self.align_pair.pattern_text[col as usize]
                    && self.align_pair.sequence_text[col as usize] != b'-'
                    && self.align_pair.pattern_text[col as usize] != b'-'
                {
                    out.push('*');
                } else {
                    out.push(' ');
                }
                col += 1;
                segment_count += 1;
            }
            out.push('\n');

            col = block_start;
            out.push_str(&format!("  {:9} ", self.align_pair.sequence_index[col as usize]));
            first = true;
            segment_count = 0;
            for _ in 0..print_width - 10 {
                if col == align_length + 1 {
                    break;
                }
                if !first
                    && self.align_pair.pattern_index[col as usize] == 0
                    && self.align_pair.pattern_index[(col - 1) as usize] != 0
                {
                    if pattern_width > 6 || print_width - 10 - segment_count < 2 * pattern_width {
                        break;
                    }
                    out.push(' ');
                    segment_count += 1;
                }
                first = false;
                out.push(self.align_pair.sequence_text[col as usize] as char);
                col += 1;
                segment_count += 1;
            }
            out.push('\n');

            col = block_start;
            out.push_str(&format!(
                "  {:9} ",
                self.align_pair.pattern_index[col as usize] + 1
            ));
            first = true;
            segment_count = 0;
            for _ in 0..print_width - 10 {
                if col == align_length + 1 {
                    break;
                }
                if !first
                    && self.align_pair.pattern_index[col as usize] == 0
                    && self.align_pair.pattern_index[(col - 1) as usize] != 0
                {
                    if pattern_width > 6 || print_width - 10 - segment_count < 2 * pattern_width {
                        break;
                    }
                    out.push(' ');
                    segment_count += 1;
                }
                first = false;
                out.push(self.align_pair.pattern_text[col as usize] as char);
                col += 1;
                segment_count += 1;
            }
            out.push_str("\n\n");
            block_start = col;
            if col == align_length + 1 {
                break;
            }
        }

        if self.align_pair.sequence_index[align_length as usize] != self.length {
            let display_end =
                std::cmp::min(self.align_pair.sequence_index[align_length as usize] + 10, self.length);
            let mut col = self.align_pair.sequence_index[align_length as usize] + 1;
            out.push_str(&format!("  {:9} ", col));
            for _ in 1..=10 {
                out.push(self.sequence[col as usize] as char);
                col += 1;
                if col > display_end {
                    break;
                }
            }
            out.push_str("\n\n");
        }

        self.write_to_alignment_file(&out);
    }

    /// `print_alignment_headings()`.
    pub fn print_alignment_headings(&mut self, consensus_length: i32) {
        self.heading = 1;
        let mut out = String::new();
        out.push_str(&format!(
            "\n\n<A NAME=\"{0}--{1},{2},{3:3.1},{4},{5}\">",
            self.align_pair.sequence_index[1],
            self.align_pair.sequence_index[self.align_pair.length as usize],
            self.period,
            self.copynumber,
            consensus_length,
            self.output_count as i32
        ));
        out.push_str("</A>");
        out.push_str("<A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n");
        out.push_str(&format!(
            "    Indices: {}--{}",
            self.align_pair.sequence_index[1],
            self.align_pair.sequence_index[self.align_pair.length as usize]
        ));
        out.push_str(&format!("  Score: {}", self.max_score));
        out.push_str(&format!(
            "\n    Period size: {}  Copynumber: {:3.1}  Consensus size: {}\n\n",
            self.period, self.copynumber, consensus_length
        ));
        self.write_to_alignment_file(&out);
    }

    /// `print_flanking_sequence()`.
    pub fn print_flanking_sequence(&mut self, flank_length: i32) {
        let print_width = 75i32;
        let mut out = String::new();
        let mut flank_start = self.align_pair.sequence_index[1] - flank_length;
        if flank_start < 1 {
            flank_start = 1;
        }
        let mut flank_end =
            self.align_pair.sequence_index[self.align_pair.length as usize] + flank_length;
        if flank_end > self.length {
            flank_end = self.length;
        }

        if flank_start == self.align_pair.sequence_index[1] {
            out.push_str("\nLeft flanking sequence: None");
        } else {
            out.push_str(&format!(
                "\nLeft flanking sequence: Indices {} -- {}\n",
                flank_start,
                self.align_pair.sequence_index[1] - 1
            ));
            let repeat_start = self.align_pair.sequence_index[1];
            let mut position = flank_start;
            loop {
                for _ in 1..=print_width - 10 {
                    out.push(self.sequence[position as usize] as char);
                    position += 1;
                    if position >= repeat_start {
                        break;
                    }
                }
                out.push('\n');
                if position >= repeat_start {
                    break;
                }
            }
        }

        let last = self.align_pair.sequence_index[self.align_pair.length as usize];
        if flank_end == last {
            out.push_str("\n\nRight flanking sequence: None");
        } else {
            out.push_str(&format!(
                "\n\nRight flanking sequence: Indices {} -- {}\n",
                last + 1,
                flank_end
            ));
            let mut position = last + 1;
            loop {
                for _ in 1..=print_width - 10 {
                    out.push(self.sequence[position as usize] as char);
                    position += 1;
                    if position > flank_end {
                        break;
                    }
                }
                out.push('\n');
                if position > flank_end {
                    break;
                }
            }
        }
        out.push_str("\n\n");
        self.write_to_alignment_file(&out);
    }

    /// `printECtoAlignments()`.
    pub fn print_consensus_to_alignment_file(&mut self, start: i32, width: i32) {
        let mut out = String::new();
        let mut chars = 0;
        out.push_str(&format!("\nConsensus pattern ({} bp):   ", width));
        for counter in start..width {
            if chars % 65 == 0 {
                out.push('\n');
            }
            out.push(self.current_pattern[counter as usize] as char);
            chars += 1;
        }
        for counter in 0..start {
            if chars % 65 == 0 {
                out.push('\n');
            }
            out.push(self.current_pattern[counter as usize] as char);
            chars += 1;
        }
        out.push('\n');
        self.write_to_alignment_file(&out);
    }

    /// `printECtoBuffer()`.
    pub fn consensus_pattern_to_string(&self, start: i32, width: i32) -> String {
        let mut out = String::new();
        for counter in start..width {
            out.push(self.current_pattern[counter as usize] as char);
        }
        for counter in 0..start {
            out.push(self.current_pattern[counter as usize] as char);
        }
        out
    }

    /// `printECtoData()`.
    #[allow(dead_code)]
    pub fn print_consensus_to_data(&mut self, start: i32, width: i32) {
        let mut out = String::new();
        for counter in start..width {
            out.push(self.current_pattern[counter as usize] as char);
        }
        for counter in 0..start {
            out.push(self.current_pattern[counter as usize] as char);
        }
        out.push('\n');
        self.write_to_alignment_file(&out);
    }

    /// `get_consensus()`.
    #[allow(unused_assignments)]
    pub fn get_consensus(&mut self, pattern_length: i32) {
        self.consensus.clear();

        let mut previous_pattern_pos: i32 = -1;
        let mut align_col: i32 = 1;
        while align_col <= self.align_pair.length {
            let pattern_position = self.align_pair.pattern_index[align_col as usize];
            if pattern_position != previous_pattern_pos {
                match self.align_pair.sequence_text[align_col as usize] {
                    b'A' => self.consensus.a[2 * pattern_position as usize + 1] += 1,
                    b'C' => self.consensus.c[2 * pattern_position as usize + 1] += 1,
                    b'G' => self.consensus.g[2 * pattern_position as usize + 1] += 1,
                    b'T' => self.consensus.t[2 * pattern_position as usize + 1] += 1,
                    DASH => self.consensus.dash[2 * pattern_position as usize + 1] += 1,
                    _ => {}
                }
                if previous_pattern_pos != -1 {
                    if pattern_position == pattern_length - 1 {
                        self.consensus.total[0] += 1;
                    } else {
                        self.consensus.total[2 * pattern_position as usize + 2] += 1;
                    }
                }
                previous_pattern_pos = pattern_position;
                align_col += 1;
            } else {
                self.consensus.insert[2 * pattern_position as usize] += 1;
                while align_col <= self.align_pair.length
                    && self.align_pair.pattern_index[align_col as usize] == previous_pattern_pos
                {
                    let pattern_position = self.align_pair.pattern_index[align_col as usize] as usize;
                    match self.align_pair.sequence_text[align_col as usize] {
                        b'A' => self.consensus.a[2 * pattern_position] += 1,
                        b'C' => self.consensus.c[2 * pattern_position] += 1,
                        b'G' => self.consensus.g[2 * pattern_position] += 1,
                        b'T' => self.consensus.t[2 * pattern_position] += 1,
                        DASH => self.consensus.dash[2 * pattern_position] += 1,
                        _ => {}
                    }
                    align_col += 1;
                }
            }
        }

        for pattern_pos in (1..=2 * pattern_length).step_by(2) {
            let mut max_count = self.consensus.dash[pattern_pos as usize];
            let mut max_char = DASH;
            if max_count < self.consensus.a[pattern_pos as usize] {
                max_count = self.consensus.a[pattern_pos as usize];
                max_char = b'A';
            }
            if max_count < self.consensus.c[pattern_pos as usize] {
                max_count = self.consensus.c[pattern_pos as usize];
                max_char = b'C';
            }
            if max_count < self.consensus.g[pattern_pos as usize] {
                max_count = self.consensus.g[pattern_pos as usize];
                max_char = b'G';
            }
            if max_count < self.consensus.t[pattern_pos as usize] {
                max_count = self.consensus.t[pattern_pos as usize];
                max_char = b'T';
            }
            self.consensus.pattern[pattern_pos as usize] = max_char;
        }
        for pattern_pos in (0..=2 * pattern_length).step_by(2) {
            let total = self.consensus.total[pattern_pos as usize];
            let insert = self.consensus.insert[pattern_pos as usize];
            if total != 0 && (insert as f32 / total as f32) >= 0.5 {
                let mut max_count = self.consensus.a[pattern_pos as usize];
                let mut max_char = b'A';
                if max_count < self.consensus.c[pattern_pos as usize] {
                    max_count = self.consensus.c[pattern_pos as usize];
                    max_char = b'C';
                }
                if max_count < self.consensus.g[pattern_pos as usize] {
                    max_count = self.consensus.g[pattern_pos as usize];
                    max_char = b'G';
                }
                if max_count < self.consensus.t[pattern_pos as usize] {
                    max_count = self.consensus.t[pattern_pos as usize];
                    max_char = b'T';
                }
                self.consensus.pattern[pattern_pos as usize] = max_char;
            } else {
                self.consensus.pattern[pattern_pos as usize] = DASH;
            }
        }

        let mut write_index: i32 = 0;
        for pattern_pos in 0..=2 * pattern_length {
            if self.consensus.pattern[pattern_pos as usize] != DASH {
                self.consensus.pattern[write_index as usize] = self.consensus.pattern[pattern_pos as usize];
                write_index += 1;
            }
        }
        self.consensus_class_length = write_index;
    }
}
