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
//! This is a direct translation of `narrowbandwrap()`, `newwrap()`,
//! `get_pair_alignment_with_copynumber()`,
//! `get_narrowband_pair_alignment_with_copynumber()`, `get_consensus()` and the
//! alignment printing routines from the reference implementation.

use crate::state::{Trf, DASH, GLOBAL, LOCAL, MAXBANDWIDTH};

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
    maxscore: i32,
    maxrealrow: i32,
    maxrow: i32,
    maxcol: i32,
    minrealrow: i32,
    mincol: i32,
    mincolbandcenter: i32,
    mincolposition: i32,
    maxrowscore: i32,
    matchatmax_col: i32,
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
        ts: &mut TraceState,
        realr: i32,
        start: i32,
        size: i32,
        r: i32,
        i: i32,
        c: i32,
        v: &mut i32,
        pleft: &mut i32,
        with_match: Option<(i32, i32)>,
    ) {
        if realr <= start - size && *v == 0 {
            *v = -1000;
            self.put_s(r as usize, i as usize, *v);
            *pleft = -1000;
        } else {
            ts.end_of_trace = false;
        }
        if *v >= ts.maxscore {
            ts.maxscore = *v;
            ts.minrealrow = realr;
            ts.mincol = c;
            ts.mincolbandcenter = self.band_center[r as usize];
            ts.mincolposition = i;
        }
        if *v > ts.maxrowscore {
            ts.maxrowscore = *v;
            let set = match with_match {
                Some((diag_val, m)) => *v == diag_val && m == self.alpha,
                None => false,
            };
            if set {
                ts.matchatmax_col = c;
            } else {
                ts.matchatmax_col = -2;
            }
        }
    }

    /// Forward half of the banded alignment.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn forward_trace_cell(
        &mut self,
        ts: &mut TraceState,
        realr: i32,
        start: i32,
        r: i32,
        i: i32,
        c: i32,
        v: &mut i32,
        pleft: &mut i32,
        with_match: Option<(i32, i32)>,
    ) {
        if realr >= start && *v == 0 {
            *v = -1000;
            self.put_s(r as usize, i as usize, *v);
            *pleft = -1000;
        } else {
            ts.end_of_trace = false;
        }
        if *v >= ts.maxscore {
            ts.maxscore = *v;
            ts.maxrealrow = realr;
            ts.maxrow = r;
            ts.maxcol = c;
        }
        if *v > ts.maxrowscore {
            ts.maxrowscore = *v;
            let set = match with_match {
                Some((diag_val, m)) => *v == diag_val && m == self.alpha,
                None => false,
            };
            if set {
                ts.matchatmax_col = c;
            } else {
                ts.matchatmax_col = -2;
            }
        }
    }

    /// `narrowbandwrap()`: banded wraparound dynamic programming.
    pub fn narrowbandwrap(
        &mut self,
        start: i32,
        size: i32,
        bandradius: i32,
        bandradiusforward: i32,
        with_consensus: bool,
        tuplesize: i32,
    ) {
        let delta = self.delta;
        let mut w = bandradius;
        if (MAXBANDWIDTH as i32) < 2 * w + 1 {
            let text = format!(
                "\nIn narrowbandwrap, MAXBANDWIDTH: {} exceeded by 2*w+1: {}\n",
                MAXBANDWIDTH,
                2 * w + 1
            );
            self.msg(&text);
            self.fatal(-1);
        }

        // fill EC
        if with_consensus {
            for g in 0..size {
                self.ec[g as usize] = self.consensus.pattern[g as usize];
            }
        } else {
            for g in 0..size {
                self.ec[g as usize] = self.sequence[(start - size + g + 1) as usize];
            }
        }

        /* ---------------- backward wdp ---------------- */
        let mut ts = TraceState {
            maxscore: 0,
            matchatmax_col: -2,
            maxrowscore: 0,
            ..Default::default()
        };
        let mut realr = start + 1;
        let mut r = self.bw as i32;
        self.band_center[r as usize] = 0;
        ts.matches_in_diagonal = 0;

        for i in 0..=w {
            let v = delta * (w - i);
            self.put_s(r as usize, i as usize, v);
            self.diag[i as usize] = v;
            self.up[i as usize] = v + delta;
        }
        for i in (w + 1)..=(2 * w) {
            self.put_s(r as usize, i as usize, -1000);
            self.diag[i as usize] = -1000;
            self.up[i as usize] = -1000 + delta;
        }

        ts.end_of_trace = false;
        while !ts.end_of_trace && realr > 1 && r > 0 {
            r -= 1;
            realr -= 1;
            self.rows += 1.0;
            ts.maxrowscore = -1;
            let lastmatchatmax_col = ts.matchatmax_col;
            ts.end_of_trace = true;
            let currchar = self.sequence[realr as usize];
            let mut pleft: i32 = -1000;

            if ts.matches_in_diagonal >= tuplesize {
                self.band_center[r as usize] = (ts.matchatmax_col - 1 + size) % size;
            } else {
                self.band_center[r as usize] =
                    (self.band_center[(r + 1) as usize] - 1 + size) % size;
            }

            let mut k =
                (self.band_center[r as usize] - self.band_center[(r + 1) as usize] + size) % size;
            if size - k <= k {
                k = -(size - k);
            }
            let mut c = (self.band_center[r as usize] + w) % size;

            if k <= -1 {
                /* band shifts left */
                k = -k;
                let mut pd = 2 * w - k + 1;
                let mut pu = 2 * w - k;
                let mut i = 2 * w;
                while i >= k {
                    let m = self.mm(currchar, self.ec[c as usize]);
                    self.diag[pd as usize] += m;
                    let dv = self.diag[pd as usize];
                    let mut v = max4(0, dv, self.up[pu as usize], pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.backward_trace_cell(
                        &mut ts,
                        realr,
                        start,
                        size,
                        r,
                        i,
                        c,
                        &mut v,
                        &mut pleft,
                        Some((dv, m)),
                    );
                    i -= 1;
                    pd -= 1;
                    pu -= 1;
                    c = (c - 1 + size) % size;
                }
                /* i == k-1: diagonal and left only */
                let m = self.mm(currchar, self.ec[c as usize]);
                self.diag[0] += m;
                let dv = self.diag[0];
                let mut v = max3(0, dv, pleft);
                self.put_s(r as usize, (k - 1) as usize, v);
                pleft = v + delta;
                self.backward_trace_cell(
                    &mut ts,
                    realr,
                    start,
                    size,
                    r,
                    k - 1,
                    c,
                    &mut v,
                    &mut pleft,
                    Some((dv, m)),
                );
                c = (c - 1 + size) % size;
                let mut i = k - 2;
                while i >= 0 {
                    let mut v = max2(0, pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.backward_trace_cell(
                        &mut ts, realr, start, size, r, i, c, &mut v, &mut pleft, None,
                    );
                    i -= 1;
                    c = (c - 1 + size) % size;
                }
            } else {
                /* band shifts right */
                let mut pd = 2 * w;
                let mut pu = 2 * w;
                let mut i = 2 * w;
                while i >= 2 * w - k + 1 {
                    let mut v = max2(0, pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.backward_trace_cell(
                        &mut ts, realr, start, size, r, i, c, &mut v, &mut pleft, None,
                    );
                    i -= 1;
                    c = (c - 1 + size) % size;
                }
                /* i == 2*w-k: up and left only */
                let mut v = max3(0, self.up[pu as usize], pleft);
                self.put_s(r as usize, i as usize, v);
                pleft = v + delta;
                self.backward_trace_cell(
                    &mut ts, realr, start, size, r, i, c, &mut v, &mut pleft, None,
                );
                pu -= 1;
                i -= 1;
                c = (c - 1 + size) % size;
                while i >= 0 {
                    let m = self.mm(currchar, self.ec[c as usize]);
                    self.diag[pd as usize] += m;
                    let dv = self.diag[pd as usize];
                    let mut v = max4(0, dv, self.up[pu as usize], pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.backward_trace_cell(
                        &mut ts,
                        realr,
                        start,
                        size,
                        r,
                        i,
                        c,
                        &mut v,
                        &mut pleft,
                        Some((dv, m)),
                    );
                    i -= 1;
                    pd -= 1;
                    pu -= 1;
                    c = (c - 1 + size) % size;
                }
            }

            for i in 0..=(2 * w) {
                let v = self.get_s(r as usize, i as usize);
                self.diag[i as usize] = v;
                self.up[i as usize] = v + delta;
            }

            if (ts.matchatmax_col - lastmatchatmax_col + size) % size == size - 1 {
                ts.matches_in_diagonal += 1;
            } else {
                ts.matches_in_diagonal = 0;
            }
        }

        /* ---------------- forward wdp ---------------- */
        let mut r = 0i32;
        let mut realr = ts.minrealrow - 1;

        self.band_center[0] = (ts.mincolbandcenter - 1 + size) % size;
        ts.matches_in_diagonal = 0;
        ts.matchatmax_col = -2;
        ts.maxscore = 0;

        let zeroat = ts.mincolposition - w;

        /* widen the narrow band in the forward direction */
        w = bandradiusforward;
        if (MAXBANDWIDTH as i32) < 2 * w + 1 {
            let text = format!(
                "\nIn narrowbandwrap, MAXBANDWIDTH: {} exceeded by bandradiusforward 2*w+1: {}\n",
                MAXBANDWIDTH,
                2 * w + 1
            );
            self.msg(&text);
            self.fatal(-1);
        }

        let first_end = std::cmp::max(0, w + zeroat);
        for i in 0..first_end {
            self.put_s(0, i as usize, -1000);
            self.diag[i as usize] = -1000;
            self.up[i as usize] = -1000 + delta;
        }
        for i in first_end..=(2 * w) {
            let v = delta * (i - (w + zeroat));
            self.put_s(0, i as usize, v);
            self.diag[i as usize] = v;
            self.up[i as usize] = v + delta;
        }

        ts.end_of_trace = false;
        while !ts.end_of_trace && realr < self.length && r < self.bw as i32 {
            r += 1;
            realr += 1;
            self.rows += 1.0;
            ts.end_of_trace = true;
            ts.maxrowscore = -1;
            let lastmatchatmax_col = ts.matchatmax_col;
            let currchar = self.sequence[realr as usize];
            let mut pleft: i32 = -1000;

            if ts.matches_in_diagonal >= tuplesize {
                self.band_center[r as usize] = (ts.matchatmax_col + 1) % size;
            } else {
                self.band_center[r as usize] = (self.band_center[(r - 1) as usize] + 1) % size;
            }

            let mut k =
                (self.band_center[r as usize] - self.band_center[(r - 1) as usize] + size) % size;
            if size - k <= k {
                k = -(size - k);
            }
            let mut c = (self.band_center[r as usize] - w + size) % size;

            if k >= 1 {
                /* band shifts right */
                let mut pd = k - 1;
                let mut pu = k;
                let mut i = 0;
                while i <= 2 * w - k {
                    let m = self.mm(currchar, self.ec[c as usize]);
                    self.diag[pd as usize] += m;
                    let dv = self.diag[pd as usize];
                    let mut v = max4(0, dv, self.up[pu as usize], pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.forward_trace_cell(
                        &mut ts,
                        realr,
                        start,
                        r,
                        i,
                        c,
                        &mut v,
                        &mut pleft,
                        Some((dv, m)),
                    );
                    i += 1;
                    pd += 1;
                    pu += 1;
                    c = (c + 1) % size;
                }
                /* i == 2*w-k+1: diagonal and left only */
                let m = self.mm(currchar, self.ec[c as usize]);
                self.diag[pd as usize] += m;
                let dv = self.diag[pd as usize];
                let mut v = max3(0, dv, pleft);
                self.put_s(r as usize, i as usize, v);
                pleft = v + delta;
                self.forward_trace_cell(
                    &mut ts,
                    realr,
                    start,
                    r,
                    i,
                    c,
                    &mut v,
                    &mut pleft,
                    Some((dv, m)),
                );
                i += 1;
                c = (c + 1) % size;
                while i <= 2 * w {
                    let mut v = max2(0, pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.forward_trace_cell(
                        &mut ts, realr, start, r, i, c, &mut v, &mut pleft, None,
                    );
                    i += 1;
                    c = (c + 1) % size;
                }
            } else {
                /* band shifts left */
                k = -k;
                let mut pd = 0;
                let mut pu = 0;
                let mut i = 0;
                while i <= k - 1 {
                    let mut v = max2(0, pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.forward_trace_cell(
                        &mut ts, realr, start, r, i, c, &mut v, &mut pleft, None,
                    );
                    i += 1;
                    c = (c + 1) % size;
                }
                /* i == k: up and left only */
                let mut v = max3(0, self.up[pu as usize], pleft);
                self.put_s(r as usize, i as usize, v);
                pleft = v + delta;
                self.forward_trace_cell(&mut ts, realr, start, r, i, c, &mut v, &mut pleft, None);
                i += 1;
                pu += 1;
                c = (c + 1) % size;
                while i <= 2 * w {
                    let m = self.mm(currchar, self.ec[c as usize]);
                    self.diag[pd as usize] += m;
                    let dv = self.diag[pd as usize];
                    let mut v = max4(0, dv, self.up[pu as usize], pleft);
                    self.put_s(r as usize, i as usize, v);
                    pleft = v + delta;
                    self.forward_trace_cell(
                        &mut ts,
                        realr,
                        start,
                        r,
                        i,
                        c,
                        &mut v,
                        &mut pleft,
                        Some((dv, m)),
                    );
                    i += 1;
                    pd += 1;
                    pu += 1;
                    c = (c + 1) % size;
                }
            }

            for i in 0..=(2 * w) {
                let v = self.get_s(r as usize, i as usize);
                self.diag[i as usize] = v;
                self.up[i as usize] = v + delta;
            }

            if (ts.matchatmax_col - lastmatchatmax_col + size) % size == 1 {
                ts.matches_in_diagonal += 1;
            } else {
                ts.matches_in_diagonal = 0;
            }
        }

        /* store last position scanned with this pattern */
        self.wrapend = realr;
        /* test for report */
        if ts.maxscore >= self.reportmin {
            self.maxrealrow = ts.maxrealrow;
            self.maxrow = ts.maxrow;
            self.maxcol = ts.maxcol;
            self.maxscore = ts.maxscore;
        }
    }

    /// `newwrap()`: full (non banded) wraparound dynamic programming.
    pub fn newwrap(&mut self, start: i32, size: i32, consensus_present: bool) {
        let delta = self.delta;
        let adjlength = size - 1;
        let adjmone = adjlength - 1;

        // fill EC
        if consensus_present {
            for g in 0..size {
                self.ec[g as usize] = self.consensus.pattern[g as usize];
            }
        } else {
            for g in 0..size {
                self.ec[g as usize] = self.sequence[(start - size + g + 1) as usize];
            }
        }

        /* ---------------- backward wdp ---------------- */
        let mut maxscore: i32 = 0;
        let mut realr = start + 1;
        let mut r = self.bw as i32;
        let mut minrow: i32 = 0;
        let mut maxrow: i32 = 0;
        let mut maxcol: i32 = 0;
        let mut maxrealrow: i32 = 0;

        for c in 0..size {
            self.put_s(r as usize, c as usize, maxscore);
            self.diag[c as usize] = maxscore;
            self.up[c as usize] = maxscore + delta;
        }

        let mut end_of_trace = false;
        while !end_of_trace && realr > 1 && r > 0 {
            r -= 1;
            realr -= 1;
            self.rows += 1.0;
            let currchar = self.sequence[realr as usize];

            // first pass
            let mut pleft = delta;
            let mut c = adjlength;
            while c >= 0 {
                self.diag[c as usize] += self.mm(currchar, self.ec[c as usize]);
                let v = max4(0, self.diag[c as usize], self.up[c as usize], pleft);
                self.put_s(r as usize, c as usize, v);
                pleft = v + delta;
                c -= 1;
            }

            end_of_trace = true;
            // second pass
            let mut c = adjlength;
            while c >= 0 {
                let mut v = max4(0, self.diag[c as usize], self.up[c as usize], pleft);
                pleft = v + delta;
                self.put_s(r as usize, c as usize, v);
                self.up[c as usize] = pleft;
                if realr <= start - std::cmp::max(size, 20) && v == 0 {
                    v = -1000;
                    self.put_s(r as usize, c as usize, v);
                    self.up[c as usize] = -1000;
                    pleft = -1000;
                } else {
                    end_of_trace = false;
                }
                if v >= maxscore {
                    maxscore = v;
                    minrow = realr;
                }
                c -= 1;
            }

            let mut c = adjmone;
            while c >= 0 {
                self.diag[c as usize] = self.get_s(r as usize, (c + 1) as usize);
                c -= 1;
            }
            self.diag[adjlength as usize] = self.get_s(r as usize, 0);
        }

        /***********************************/
        let mut r = 0i32;
        let mut realr = start - 1;

        for c in 0..size {
            self.put_s(0, c as usize, maxscore);
            self.diag[c as usize] = maxscore;
            self.up[c as usize] = maxscore + delta;
        }

        end_of_trace = false;
        while !end_of_trace && realr < self.length && r < self.bw as i32 {
            r += 1;
            realr += 1;
            self.rows += 1.0;
            let currchar = self.sequence[realr as usize];

            // first pass
            let mut pleft = delta;
            let mut c = 0;
            while c < size {
                self.diag[c as usize] += self.mm(currchar, self.ec[c as usize]);
                let v = max4(0, self.diag[c as usize], self.up[c as usize], pleft);
                self.put_s(r as usize, c as usize, v);
                pleft = v + delta;
                c += 1;
            }

            end_of_trace = true;
            // second pass
            let mut c = 0;
            while c < size {
                let mut v = max4(0, self.diag[c as usize], self.up[c as usize], pleft);
                pleft = v + delta;
                self.put_s(r as usize, c as usize, v);
                self.up[c as usize] = pleft;
                if realr >= start && v == 0 {
                    v = -1000;
                    self.put_s(r as usize, c as usize, v);
                    self.up[c as usize] = -1000;
                    pleft = -1000;
                } else {
                    end_of_trace = false;
                }
                if v >= maxscore {
                    maxscore = v;
                    maxrow = realr;
                }
                c += 1;
            }

            let mut c = 0;
            while c < size {
                self.diag[(c + 1) as usize] = self.get_s(r as usize, c as usize);
                c += 1;
            }
            self.diag[0] = self.diag[size as usize];
        }

        /* store last position scanned with this pattern */
        self.wrapend = realr;

        /* redo for an accurate alignment if score is big enough to report */
        self.reportmin = 0;
        if maxscore >= self.reportmin {
            let mut r = 0i32;
            let modstart = start + 1;
            let mut realr = minrow - size - 1;
            if realr < 0 {
                realr = 0;
            }
            maxscore = 0;

            for c in 0..size {
                self.put_s(0, c as usize, maxscore);
                self.diag[c as usize] = maxscore;
                self.up[c as usize] = maxscore + delta;
            }

            end_of_trace = false;
            while !end_of_trace && realr < self.length && r < self.bw as i32 {
                r += 1;
                realr += 1;
                self.rows += 1.0;
                let currchar = self.sequence[realr as usize];

                let mut pleft = delta;
                let mut c = 0;
                while c < size {
                    self.diag[c as usize] += self.mm(currchar, self.ec[c as usize]);
                    let v = max4(0, self.diag[c as usize], self.up[c as usize], pleft);
                    self.put_s(r as usize, c as usize, v);
                    pleft = v + delta;
                    c += 1;
                }

                end_of_trace = true;
                let mut c = 0;
                while c < size {
                    let mut v = max4(0, self.diag[c as usize], self.up[c as usize], pleft);
                    pleft = v + delta;
                    self.put_s(r as usize, c as usize, v);
                    self.up[c as usize] = pleft;
                    if realr >= modstart && v == 0 {
                        v = -1000;
                        self.put_s(r as usize, c as usize, v);
                        self.up[c as usize] = -1000;
                        pleft = -1000;
                    } else {
                        end_of_trace = false;
                    }
                    if v > maxscore {
                        maxscore = v;
                        maxrealrow = realr;
                        maxrow = r;
                        maxcol = c;
                    }
                    c += 1;
                }

                let mut c = 0;
                while c < size {
                    self.diag[(c + 1) as usize] = self.get_s(r as usize, c as usize);
                    c += 1;
                }
                self.diag[0] = self.diag[size as usize];
            }
        }

        if maxscore >= self.reportmin {
            self.maxrealrow = maxrealrow;
            self.maxrow = maxrow;
            self.maxcol = maxcol;
            self.maxscore = maxscore;
        }
    }

    #[inline]
    fn s_at(&self, r: i32, i: i32) -> i32 {
        if r < 0 || i < 0 {
            return 0;
        }
        self.get_s(r as usize, i as usize)
    }

    /// `get_pair_alignment_with_copynumber()`: traceback for the small distances.
    pub fn get_pair_alignment_with_copynumber(&mut self, size: i32) {
        let mut i = self.maxrealrow;
        let mut si = self.maxrow;
        let mut j = self.maxcol;
        let adjlength = size - 1;
        let fullcopy = (j + 1) % size;

        self.align_pair.score = self.s_at(si, j);
        let mut length = 0i32;
        self.copynumber = 0.0;

        loop {
            let sval = self.s_at(si, j);
            let mut done = false;
            if sval <= 0 {
                done = true;
            } else if sval
                == self.s_at(si - 1, (j + adjlength) % size)
                    + self.mm(self.sequence[i as usize], self.ec[j as usize])
            {
                length += 1;
                if j == fullcopy {
                    self.copynumber += 1.0;
                }
                self.fill_align_pair(self.sequence[i as usize], self.ec[j as usize], length, i, j);
                i -= 1;
                si -= 1;
                j = (j + adjlength) % size;
            } else if sval == self.s_at(si - 1, j) + self.delta {
                length += 1;
                self.fill_align_pair(self.sequence[i as usize], DASH, length, i, (j + 1) % size);
                i -= 1;
                si -= 1;
            } else if sval == self.s_at(si, (j + adjlength) % size) + self.delta {
                length += 1;
                if j == fullcopy {
                    self.copynumber += 1.0;
                }
                self.fill_align_pair(DASH, self.ec[j as usize], length, i + 1, j);
                j = (j + adjlength) % size;
            } else {
                self.msg("\nget_pair_alignment_with_copynumber: error in trace back");
                let text = format!("\nrow={}  column={}", i, j);
                self.msg(&text);
                let text = format!(
                    "\nS={}  Sleft={}  Sup={}  Sdiag={}  match={}",
                    sval,
                    self.s_at(si, (j + adjlength) % size),
                    self.s_at(si - 1, j),
                    self.s_at(si - 1, (j + adjlength) % size),
                    self.mm(self.sequence[i as usize], self.ec[j as usize])
                );
                self.msg(&text);
                done = true;
            }

            if done {
                self.align_pair.length = length;
                if self.maxcol >= j {
                    self.copynumber += (self.maxcol - j) as f64 / size as f64;
                } else {
                    self.copynumber += (self.maxcol + size - j) as f64 / size as f64;
                }
                return;
            }
        }
    }

    /// `get_narrowband_pair_alignment_with_copynumber()`.
    pub fn get_narrowband_pair_alignment_with_copynumber(
        &mut self,
        size: i32,
        bandradius: i32,
        option: u8,
    ) {
        let w = bandradius;
        let mut realr = self.maxrealrow;
        let mut r = self.maxrow;
        let mut c = self.maxcol;
        let mut k = (self.maxcol - self.band_center[r as usize] + size) % size;
        if size - k <= k {
            k = -(size - k);
        }
        let mut i = w + k;
        let fullcopy = (self.maxcol + 1) % size;

        self.align_pair.score = self.s_at(r, i);
        let mut length = 0i32;
        self.copynumber = 0.0;

        loop {
            /* stop at zeros or -1000 for local; stop at r=0 for global */
            let mut legitimate_zero = false;
            let mut sval = self.s_at(r, i);
            if option == LOCAL && sval == 0 && r != 0 {
                /* check to see if the LOCAL zero score is a legitimate zero */
                let mut kk = (c - self.band_center[(r - 1) as usize] + size) % size;
                if size - kk <= kk {
                    kk = -(size - kk);
                }
                let upi = w + kk;
                let mut k2 = (self.band_center[r as usize] - self.band_center[(r - 1) as usize]
                    + size)
                    % size;
                if size - k2 <= k2 {
                    k2 = -(size - k2);
                }
                let m = self.mm(self.sequence[realr as usize], self.ec[c as usize]);
                if k2 >= 1 {
                    /* band shifts right */
                    if i <= 2 * w - k2 {
                        if sval == self.s_at(r - 1, upi - 1) + m
                            || sval == self.s_at(r - 1, upi) + self.delta
                            || sval == self.s_at(r, i - 1) + self.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if i == 2 * w - k2 + 1 {
                        if sval == self.s_at(r - 1, upi - 1) + m
                            || sval == self.s_at(r, i - 1) + self.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if sval == self.s_at(r, i - 1) + self.delta {
                        legitimate_zero = true;
                    }
                } else {
                    /* band shifts left */
                    let k2 = -k2;
                    if i <= k2 - 1 {
                        if sval == self.s_at(r, i - 1) + self.delta {
                            legitimate_zero = true;
                        }
                    } else if i == k2 {
                        if sval == self.s_at(r - 1, upi) + self.delta
                            || sval == self.s_at(r, i - 1) + self.delta
                        {
                            legitimate_zero = true;
                        }
                    } else if sval == self.s_at(r - 1, upi - 1) + m
                        || sval == self.s_at(r - 1, upi) + self.delta
                        || sval == self.s_at(r, i - 1) + self.delta
                    {
                        legitimate_zero = true;
                    }
                }
            }

            /* Stop at zeros or -1000 for local, at r=0 for global.  The
            reference implementation also tests `GLOBAL && r == 0 && c == -1`,
            which is already implied by the leading `r == 0`. */
            let stop = r == 0 || (option == LOCAL && sval <= 0 && !legitimate_zero);
            if stop {
                self.align_pair.length = length;
                if self.maxcol >= c {
                    self.copynumber += (self.maxcol - c) as f64 / size as f64;
                } else {
                    self.copynumber += (self.maxcol + size - c) as f64 / size as f64;
                }
                return;
            }

            let mut kk = (c - self.band_center[(r - 1) as usize] + size) % size;
            if size - kk <= kk {
                kk = -(size - kk);
            }
            let upi = w + kk;

            /* try the three predecessors: diagonal, up, left */
            let mut advanced = false;
            if !(option == GLOBAL && c == -1) {
                let mut k2 = (self.band_center[r as usize] - self.band_center[(r - 1) as usize]
                    + size)
                    % size;
                if size - k2 <= k2 {
                    k2 = -(size - k2);
                }
                let shifts_right = k2 >= 1;
                let kabs = if shifts_right { k2 } else { -k2 };

                let try_diag = if shifts_right {
                    i <= 2 * w - kabs || i == 2 * w - kabs + 1
                } else {
                    i > kabs
                };
                let try_up = if shifts_right {
                    i <= 2 * w - kabs
                } else {
                    i >= kabs
                };
                /* the left branch is always tested */
                let try_left = true;

                sval = self.s_at(r, i);
                if try_diag
                    && sval
                        == self.s_at(r - 1, upi - 1)
                            + self.mm(self.sequence[realr as usize], self.ec[c as usize])
                {
                    length += 1;
                    if c == fullcopy {
                        self.copynumber += 1.0;
                    }
                    if option == LOCAL {
                        self.fill_align_pair(
                            self.sequence[realr as usize],
                            self.ec[c as usize],
                            length,
                            realr,
                            c,
                        );
                        c = (c - 1 + size) % size;
                    } else {
                        self.fill_align_pair(
                            self.sequence[realr as usize],
                            self.ec[c as usize],
                            length,
                            realr,
                            c + self.maxrealcol - self.maxcol,
                        );
                        c -= 1;
                    }
                    realr -= 1;
                    r -= 1;
                    i = upi - 1;
                    advanced = true;
                }

                if !advanced && try_up && sval == self.s_at(r - 1, upi) + self.delta {
                    length += 1;
                    if option == LOCAL {
                        self.fill_align_pair(
                            self.sequence[realr as usize],
                            DASH,
                            length,
                            realr,
                            (c + 1) % size,
                        );
                    } else {
                        self.fill_align_pair(
                            self.sequence[realr as usize],
                            DASH,
                            length,
                            realr,
                            c + 1 + self.maxrealcol - self.maxcol,
                        );
                    }
                    realr -= 1;
                    r -= 1;
                    i = upi;
                    advanced = true;
                }

                if !advanced && try_left {
                    if i == 0 {
                        self.msg("\nget_pair_alignment_with_copynumber: error in trace back");
                        self.msg("\nattempted to compute left branch when i==0");
                        let text = format!("\nS[{}][{}]={}", r, i, sval);
                        self.msg(&text);
                        break;
                    } else if sval == self.s_at(r, i - 1) + self.delta {
                        length += 1;
                        if c == fullcopy {
                            self.copynumber += 1.0;
                        }
                        if option == LOCAL {
                            self.fill_align_pair(DASH, self.ec[c as usize], length, realr + 1, c);
                            c = (c - 1 + size) % size;
                        } else {
                            self.fill_align_pair(
                                DASH,
                                self.ec[c as usize],
                                length,
                                realr + 1,
                                c + self.maxrealcol - self.maxcol,
                            );
                            c -= 1;
                        }
                        i -= 1;
                        advanced = true;
                    }
                }
            } else {
                /* global alignment past the left edge: the only option is up */
                sval = self.s_at(r, i);
                if sval == self.s_at(r - 1, upi) + self.delta {
                    length += 1;
                    self.fill_align_pair(
                        self.sequence[realr as usize],
                        DASH,
                        length,
                        realr,
                        (c + 1) % size,
                    );
                    realr -= 1;
                    r -= 1;
                    i = upi;
                    advanced = true;
                }
            }

            if !advanced {
                self.msg("\nget_pair_alignment_with_copynumber: error in trace back");
                let text = format!(
                    "\nS: row={}  column={}  upi={}  Sequence: realrow={}  EC: realcol={}",
                    r, i, upi, realr, c
                );
                self.msg(&text);
                break;
            }
        }
    }

    /// `reverse()`: reverses the alignment in `AlignPair`.
    pub fn reverse(&mut self) {
        let ml = self.align_pair.length;
        let mut j = 1;
        while j <= ml / 2 {
            let k = ml - j + 1;
            self.align_pair.textprime.swap(j as usize, k as usize);
            self.align_pair.textsecnd.swap(j as usize, k as usize);
            self.align_pair.indexprime.swap(j as usize, k as usize);
            self.align_pair.indexsecnd.swap(j as usize, k as usize);
            j += 1;
        }
    }

    /// `shift_pattern_indices()`: makes the pattern indices run from 0 to size-1.
    pub fn shift_pattern_indices(&mut self, patternsize: i32) {
        let downshift = self.align_pair.indexsecnd[1];
        let upshift = patternsize - downshift;
        for l in 1..=self.align_pair.length {
            if self.align_pair.indexsecnd[l as usize] >= downshift {
                self.align_pair.indexsecnd[l as usize] -= downshift;
            } else {
                self.align_pair.indexsecnd[l as usize] += upshift;
            }
        }
    }

    /// `alt3_print_alignment()`.
    pub fn alt3_print_alignment(&mut self, patternwidth: i32) {
        self.shift_pattern_indices(patternwidth);
        let pwidth = 75i32;
        let mut out = String::new();

        if pwidth <= 0 {
            out.push_str("Error, pwidth<=0, can't report alignments.\n");
            self.msg(&out);
            return;
        }

        let alen = self.align_pair.length;
        if self.align_pair.indexprime[1] != 1 {
            let mut m = self.align_pair.indexprime[1] - 10;
            if m < 1 {
                m = 1;
            }
            let mut j = m;
            out.push_str(&format!("  {:9} ", j));
            for _ in 1..=10 {
                out.push(self.sequence[j as usize] as char);
                j += 1;
                if j == self.align_pair.indexprime[1] {
                    break;
                }
            }
            out.push_str("\n\n");
        }

        let mut g = 1i32;
        loop {
            let mut j = g;
            out.push_str("            ");
            let mut first = true;
            let mut h = 0;
            for _ in 0..pwidth - 10 {
                if j == alen + 1 {
                    break;
                }
                if !first
                    && self.align_pair.indexsecnd[j as usize] == 0
                    && self.align_pair.indexsecnd[(j - 1) as usize] != 0
                {
                    if patternwidth > 6 || pwidth - 10 - h < 2 * patternwidth {
                        break;
                    }
                    out.push(' ');
                    h += 1;
                }
                first = false;
                if self.align_pair.textprime[j as usize] != self.align_pair.textsecnd[j as usize]
                    && self.align_pair.textprime[j as usize] != b'-'
                    && self.align_pair.textsecnd[j as usize] != b'-'
                {
                    out.push('*');
                } else {
                    out.push(' ');
                }
                j += 1;
                h += 1;
            }
            out.push('\n');

            j = g;
            out.push_str(&format!("  {:9} ", self.align_pair.indexprime[j as usize]));
            first = true;
            h = 0;
            for _ in 0..pwidth - 10 {
                if j == alen + 1 {
                    break;
                }
                if !first
                    && self.align_pair.indexsecnd[j as usize] == 0
                    && self.align_pair.indexsecnd[(j - 1) as usize] != 0
                {
                    if patternwidth > 6 || pwidth - 10 - h < 2 * patternwidth {
                        break;
                    }
                    out.push(' ');
                    h += 1;
                }
                first = false;
                out.push(self.align_pair.textprime[j as usize] as char);
                j += 1;
                h += 1;
            }
            out.push('\n');

            j = g;
            out.push_str(&format!(
                "  {:9} ",
                self.align_pair.indexsecnd[j as usize] + 1
            ));
            first = true;
            h = 0;
            for _ in 0..pwidth - 10 {
                if j == alen + 1 {
                    break;
                }
                if !first
                    && self.align_pair.indexsecnd[j as usize] == 0
                    && self.align_pair.indexsecnd[(j - 1) as usize] != 0
                {
                    if patternwidth > 6 || pwidth - 10 - h < 2 * patternwidth {
                        break;
                    }
                    out.push(' ');
                    h += 1;
                }
                first = false;
                out.push(self.align_pair.textsecnd[j as usize] as char);
                j += 1;
                h += 1;
            }
            out.push_str("\n\n");
            g = j;
            if j == alen + 1 {
                break;
            }
        }

        if self.align_pair.indexprime[alen as usize] != self.length {
            let m = std::cmp::min(self.align_pair.indexprime[alen as usize] + 10, self.length);
            let mut j = self.align_pair.indexprime[alen as usize] + 1;
            out.push_str(&format!("  {:9} ", j));
            for _ in 1..=10 {
                out.push(self.sequence[j as usize] as char);
                j += 1;
                if j > m {
                    break;
                }
            }
            out.push_str("\n\n");
        }

        self.msg(&out);
    }

    /// `print_alignment_headings()`.
    pub fn print_alignment_headings(&mut self, consensuslength: i32) {
        self.heading = 1;
        let mut out = String::new();
        out.push_str(&format!(
            "\n\n<A NAME=\"{0}--{1},{2},{3:3.1},{4},{5}\">",
            self.align_pair.indexprime[1],
            self.align_pair.indexprime[self.align_pair.length as usize],
            self.period,
            self.copynumber,
            consensuslength,
            self.output_count as i32
        ));
        out.push_str("</A>");
        out.push_str("<A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n");
        out.push_str(&format!(
            "    Indices: {}--{}",
            self.align_pair.indexprime[1],
            self.align_pair.indexprime[self.align_pair.length as usize]
        ));
        out.push_str(&format!("  Score: {}", self.maxscore));
        out.push_str(&format!(
            "\n    Period size: {}  Copynumber: {:3.1}  Consensus size: {}\n\n",
            self.period, self.copynumber, consensuslength
        ));
        self.msg(&out);
    }

    /// `print_flanking_sequence()`.
    pub fn print_flanking_sequence(&mut self, flank_length: i32) {
        let pwidth = 75i32;
        let mut out = String::new();
        let mut m = self.align_pair.indexprime[1] - flank_length;
        if m < 1 {
            m = 1;
        }
        let mut n = self.align_pair.indexprime[self.align_pair.length as usize] + flank_length;
        if n > self.length {
            n = self.length;
        }

        if m == self.align_pair.indexprime[1] {
            out.push_str("\nLeft flanking sequence: None");
        } else {
            out.push_str(&format!(
                "\nLeft flanking sequence: Indices {} -- {}\n",
                m,
                self.align_pair.indexprime[1] - 1
            ));
            let k = self.align_pair.indexprime[1];
            let mut j = m;
            loop {
                for _ in 1..=pwidth - 10 {
                    out.push(self.sequence[j as usize] as char);
                    j += 1;
                    if j >= k {
                        break;
                    }
                }
                out.push('\n');
                if j >= k {
                    break;
                }
            }
        }

        let last = self.align_pair.indexprime[self.align_pair.length as usize];
        if n == last {
            out.push_str("\n\nRight flanking sequence: None");
        } else {
            out.push_str(&format!(
                "\n\nRight flanking sequence: Indices {} -- {}\n",
                last + 1,
                n
            ));
            let mut j = last + 1;
            loop {
                for _ in 1..=pwidth - 10 {
                    out.push(self.sequence[j as usize] as char);
                    j += 1;
                    if j > n {
                        break;
                    }
                }
                out.push('\n');
                if j > n {
                    break;
                }
            }
        }
        out.push_str("\n\n");
        self.msg(&out);
    }

    /// `printECtoAlignments()`.
    pub fn print_ec_to_alignments(&mut self, start: i32, width: i32) {
        let mut out = String::new();
        let mut chars = 0;
        out.push_str(&format!("\nConsensus pattern ({} bp):   ", width));
        for counter in start..width {
            if chars % 65 == 0 {
                out.push('\n');
            }
            out.push(self.ec[counter as usize] as char);
            chars += 1;
        }
        for counter in 0..start {
            if chars % 65 == 0 {
                out.push('\n');
            }
            out.push(self.ec[counter as usize] as char);
            chars += 1;
        }
        out.push('\n');
        self.msg(&out);
    }

    /// `printECtoBuffer()`.
    pub fn print_ec_to_buffer(&self, start: i32, width: i32) -> String {
        let mut out = String::new();
        for counter in start..width {
            out.push(self.ec[counter as usize] as char);
        }
        for counter in 0..start {
            out.push(self.ec[counter as usize] as char);
        }
        out
    }

    /// `printECtoData()`.
    #[allow(dead_code)]
    pub fn print_ec_to_data(&mut self, start: i32, width: i32) {
        let mut out = String::new();
        for counter in start..width {
            out.push(self.ec[counter as usize] as char);
        }
        for counter in 0..start {
            out.push(self.ec[counter as usize] as char);
        }
        out.push('\n');
        self.msg(&out);
    }

    /// `get_consensus()`.
    #[allow(unused_assignments)]
    pub fn get_consensus(&mut self, patternlength: i32) {
        self.consensus.clear();

        let mut lastindex: i32 = -1;
        let mut i: i32 = 1;
        while i <= self.align_pair.length {
            let secnd = self.align_pair.indexsecnd[i as usize];
            if secnd != lastindex {
                match self.align_pair.textprime[i as usize] {
                    b'A' => self.consensus.a[2 * secnd as usize + 1] += 1,
                    b'C' => self.consensus.c[2 * secnd as usize + 1] += 1,
                    b'G' => self.consensus.g[2 * secnd as usize + 1] += 1,
                    b'T' => self.consensus.t[2 * secnd as usize + 1] += 1,
                    DASH => self.consensus.dash[2 * secnd as usize + 1] += 1,
                    _ => {}
                }
                if lastindex != -1 {
                    if secnd == patternlength - 1 {
                        self.consensus.total[0] += 1;
                    } else {
                        self.consensus.total[2 * secnd as usize + 2] += 1;
                    }
                }
                lastindex = secnd;
                i += 1;
            } else {
                self.consensus.insert[2 * secnd as usize] += 1;
                while i <= self.align_pair.length
                    && self.align_pair.indexsecnd[i as usize] == lastindex
                {
                    let secnd = self.align_pair.indexsecnd[i as usize] as usize;
                    match self.align_pair.textprime[i as usize] {
                        b'A' => self.consensus.a[2 * secnd] += 1,
                        b'C' => self.consensus.c[2 * secnd] += 1,
                        b'G' => self.consensus.g[2 * secnd] += 1,
                        b'T' => self.consensus.t[2 * secnd] += 1,
                        DASH => self.consensus.dash[2 * secnd] += 1,
                        _ => {}
                    }
                    i += 1;
                }
            }
        }

        for i in (1..=2 * patternlength).step_by(2) {
            let mut max = self.consensus.dash[i as usize];
            let mut maxchar = DASH;
            if max < self.consensus.a[i as usize] {
                max = self.consensus.a[i as usize];
                maxchar = b'A';
            }
            if max < self.consensus.c[i as usize] {
                max = self.consensus.c[i as usize];
                maxchar = b'C';
            }
            if max < self.consensus.g[i as usize] {
                max = self.consensus.g[i as usize];
                maxchar = b'G';
            }
            if max < self.consensus.t[i as usize] {
                max = self.consensus.t[i as usize];
                maxchar = b'T';
            }
            self.consensus.pattern[i as usize] = maxchar;
        }
        for i in (0..=2 * patternlength).step_by(2) {
            let total = self.consensus.total[i as usize];
            let insert = self.consensus.insert[i as usize];
            if total != 0 && (insert as f32 / total as f32) >= 0.5 {
                let mut max = self.consensus.a[i as usize];
                let mut maxchar = b'A';
                if max < self.consensus.c[i as usize] {
                    max = self.consensus.c[i as usize];
                    maxchar = b'C';
                }
                if max < self.consensus.g[i as usize] {
                    max = self.consensus.g[i as usize];
                    maxchar = b'G';
                }
                if max < self.consensus.t[i as usize] {
                    max = self.consensus.t[i as usize];
                    maxchar = b'T';
                }
                self.consensus.pattern[i as usize] = maxchar;
            } else {
                self.consensus.pattern[i as usize] = DASH;
            }
        }

        let mut j: i32 = 0;
        for i in 0..=2 * patternlength {
            if self.consensus.pattern[i as usize] != DASH {
                self.consensus.pattern[j as usize] = self.consensus.pattern[i as usize];
                j += 1;
            }
        }
        self.cons_classlength = j;
    }
}
