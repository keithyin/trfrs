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

//! Distance bookkeeping, statistical criteria, the main scan loop and the
//! per-repeat statistics (`tr30dat.c`).

use std::collections::VecDeque;

use crate::clean::IlRecord;
use crate::fasta::base_index;
use crate::state::{
    BestPeriodElement, HistoryEntry, Trf, DASH, FOUR_TO_THE, MAXPATTERNSIZECONSTANT, MINBANDRADIUS,
    NUMBER_OF_PERIODS, NUMBER_OF_PERIODS_INTO_SORTMULTIPLES, NUMBER_OF_PERIODS_TO_TEST,
    RECENTERCRITERION, SMALLDISTANCE, TAGSEP,
};
use crate::tables::{SUMDATA75, SUMDATA80, WAITDATA75, WAITDATA80};

#[inline]
fn max2(a: i32, b: i32) -> i32 {
    if a >= b {
        a
    } else {
        b
    }
}

#[inline]
fn min2(a: i32, b: i32) -> i32 {
    if a <= b {
        a
    } else {
        b
    }
}

#[inline]
fn is_acgt(c: u8) -> bool {
    matches!(c, b'a' | b'c' | b'g' | b't' | b'A' | b'C' | b'G' | b'T')
}

/// `GetTopPeriods()`.
fn get_top_periods(
    index: &[i32; 256],
    pattern: &[u8],
    length: i32,
    maxdistance: i32,
    toparray: &mut [i32],
) -> i32 {
    let mut counts = vec![0f64; length as usize];
    let mut history = vec![0i32; length as usize];
    let mut heads = [-1i32; 16];

    let mut i = 0i32;
    let end = length - 2;
    while i <= end {
        let tupid =
            index[pattern[i as usize] as usize] * 4 + index[pattern[(i + 1) as usize] as usize];
        history[i as usize] = heads[tupid as usize];
        heads[tupid as usize] = i;

        let mut dist = 0;
        let mut curr = i;
        while history[curr as usize] != -1 && dist < MAXPATTERNSIZECONSTANT as i32 * 3 {
            dist = i - history[curr as usize];
            counts[dist as usize] += 1.0;
            curr = history[curr as usize];
        }
        i += 1;
    }

    let mut xysum = 0f64;
    let mut xsum = 0f64;
    let mut ysum = 0f64;
    let mut x2sum = 0f64;
    let end = length - 2;
    for i in 1..=end {
        xysum += i as f64 * counts[i as usize];
        xsum += i as f64;
        ysum += counts[i as usize];
        x2sum += (i * i) as f64;
    }
    let n = end as f64;
    let s = (n * xysum - xsum * ysum) / (n * x2sum - xsum * xsum);

    for i in 1..=end {
        counts[i as usize] -= i as f64 * s;
    }

    let mut end = length - 2;
    if end > maxdistance {
        /* accepts smaller multiples if the best ones are too large */
        end = maxdistance;
    }
    for t in 0..NUMBER_OF_PERIODS {
        let mut topind = 0;
        let mut topval = 0.0;
        for i in 1..=end {
            if counts[i as usize] > topval {
                topind = i;
                topval = counts[i as usize];
            }
        }
        toparray[t] = topind;
        counts[topind as usize] = 0.0;
    }

    0
}

impl Trf {
    /* ------------------------------------------------------------------ */
    /* distance lists                                                      */
    /* ------------------------------------------------------------------ */

    /// `add_tuple_match_to_Distance_entry()`.
    pub fn add_tuple_match_to_distance_entry(&mut self, location: i32, size: i32, d: i32) {
        let windowsize = max2(d, 20) + 1;
        let list = &mut self.distance[d as usize];

        if list.numentries != 0 {
            let windowleftend = min2(location - d, location - 20) + 1;
            loop {
                if list.numentries == 0 {
                    break;
                }
                let lo = list.lowindex as usize;
                if list.entry[lo].location < windowleftend {
                    let sz = list.entry[lo].size;
                    list.numentries -= 1;
                    list.nummatches -= sz;
                    list.lowindex += 1;
                    list.lowindex %= windowsize;
                } else {
                    break;
                }
            }
        }

        let hi = list.highindex as usize;
        if list.numentries != 0 && list.entry[hi].location == location - 1 {
            /* adjacent tuple: extend the last entry */
            list.entry[hi].location += 1;
            list.entry[hi].size += 1;
            list.nummatches += 1;
        } else {
            list.numentries += 1;
            list.highindex += 1;
            list.highindex %= windowsize;
            let hi = list.highindex as usize;
            list.entry[hi].location = location;
            list.entry[hi].size = size;
            list.nummatches += size;
        }
    }

    /// `link_Distance_window()`.
    pub fn link_distance_window(&mut self, d: i32) {
        /* the reference implementation applies ceil() to an integer division */
        let t = d / TAGSEP;
        let mut f;
        if self.tag[t as usize] < d {
            f = self.tag[t as usize];
            let mut tt = t;
            while tt <= self.top_tag && self.tag[tt as usize] < d {
                self.tag[tt as usize] = d;
                tt += 1;
            }
        } else if self.tag[t as usize] > d {
            f = self.tag[t as usize];
            while f > d {
                f = self.distance[f as usize].linkdown;
            }
            if f == d {
                let text = format!("\nTag error following links.  f==d={}", d);
                self.msg(&text);
            }
        } else {
            let text = format!("\nTag error Tag[{}]={}", t, d);
            self.msg(&text);
            self.fatal(-2);
        }

        /* link in d */
        self.distance[d as usize].linkdown = f;
        let h = self.distance[f as usize].linkup;
        self.distance[d as usize].linkup = h;
        self.distance[f as usize].linkup = d;
        if h <= self.maxdistance {
            self.distance[h as usize].linkdown = d;
        }
        self.distance[d as usize].linked = true;
    }

    /// `untag_Distance_window()`.
    pub fn untag_distance_window(&mut self, d: i32, linkdown: i32) {
        let t = d / TAGSEP;
        if self.tag[t as usize] != d {
            return;
        }
        let mut tt = t;
        while tt <= self.top_tag && self.tag[tt as usize] == d {
            self.tag[tt as usize] = linkdown;
            tt += 1;
        }
    }

    /// `no_matches_so_unlink_Distance()`; returns 1 when the distance was unlinked.
    pub fn no_matches_so_unlink_distance(&mut self, d: i32, location: i32) -> i32 {
        let windowsize = max2(d, 20) + 1;
        let mut numentries = self.distance[d as usize].numentries;
        if numentries != 0 {
            let windowleftend = min2(location - d, location - 20) + 1;
            loop {
                if numentries == 0 {
                    break;
                }
                let list = &mut self.distance[d as usize];
                let lo = list.lowindex as usize;
                if list.entry[lo].location < windowleftend {
                    let sz = list.entry[lo].size;
                    list.numentries -= 1;
                    list.nummatches -= sz;
                    list.lowindex += 1;
                    list.lowindex %= windowsize;
                    numentries = list.numentries;
                } else {
                    break;
                }
            }
        }
        if self.distance[d as usize].numentries == 0 {
            let g = self.distance[d as usize].linkdown;
            let h = self.distance[d as usize].linkup;

            self.distance[g as usize].linkup = h;
            if h <= self.maxdistance {
                self.distance[h as usize].linkdown = g;
            }
            self.distance[d as usize].linked = false;

            self.untag_distance_window(d, g);
            return 1;
        }
        0
    }

    /* ------------------------------------------------------------------ */
    /* distances seen                                                      */
    /* ------------------------------------------------------------------ */

    pub fn add_to_distanceseenarray(&mut self, location: i32, distance: i32, end: i32, score: i32) {
        let ptr = &mut self.distanceseenarray[distance as usize];
        ptr.index = location;
        ptr.end = end;
        ptr.score = score;
    }

    pub fn search_for_distance_match_in_distanceseenarray(
        &self,
        distance: i32,
        start: i32,
    ) -> bool {
        self.distanceseenarray[distance as usize].end >= start
    }

    /* ------------------------------------------------------------------ */
    /* best period list                                                    */
    /* ------------------------------------------------------------------ */

    pub fn init_bestperiodlist(&mut self) {
        self.best_period_list = VecDeque::new();
    }

    /// `add_to_bestperiodlist()`.
    pub fn add_to_bestperiodlist(&mut self, d: i32) {
        if d != 1 {
            let e = BestPeriodElement {
                indexlow: self.align_pair.indexprime[self.align_pair.length as usize],
                indexhigh: self.align_pair.indexprime[1],
                best1: self.sort_multiples[0],
                best2: self.sort_multiples[1],
                best3: self.sort_multiples[2],
                best4: self.sort_multiples[3],
                best5: self.sort_multiples[4],
            };
            self.best_period_list.push_front(e);
        }
    }

    /// `adjust_bestperiod_entry()`.
    pub fn adjust_bestperiod_entry(&mut self, d: i32) {
        if d != 1 {
            if let Some(entry) = self.best_period_list.front_mut() {
                if entry.indexhigh > self.align_pair.indexprime[1] {
                    entry.indexhigh = self.align_pair.indexprime[1];
                }
            }
        }
    }

    /// `search_for_range_in_bestperiodlist()`.
    pub fn search_for_range_in_bestperiodlist(&mut self, start: i32, distance: i32) -> bool {
        let mut range_covered = false;
        let mut i = 0usize;
        while i < self.best_period_list.len() {
            let e = self.best_period_list[i];
            if e.indexhigh < start - 2 * self.maxdistance {
                self.best_period_list.remove(i);
            } else {
                if e.indexlow
                    <= start - 2 * distance
                        + 1
                        + self.distance[distance as usize].waiting_time_criteria
                    && e.indexhigh >= start
                {
                    range_covered = true;
                    if e.best1 == distance
                        || e.best2 == distance
                        || e.best3 == distance
                        || e.best4 == distance
                        || e.best5 == distance
                    {
                        return true;
                    }
                }
                i += 1;
            }
        }
        !range_covered
    }

    /* ------------------------------------------------------------------ */
    /* criteria                                                            */
    /* ------------------------------------------------------------------ */

    /// `multiples_criteria_4()`.
    pub fn multiples_criteria_4(&mut self, found_d: i32) -> bool {
        let lowerindex = self.align_pair.indexprime[self.align_pair.length as usize];
        let upperindex = self.align_pair.indexprime[1];
        let length = upperindex - lowerindex + 1;

        if found_d == 1 {
            let mut comps = [0i32; 4];
            let mut total = 0;
            let mut maxind = 0usize;
            for g in 0..length {
                let b = self.sequence[(lowerindex + g) as usize];
                comps[self.index[b as usize] as usize] += 1;
            }
            for g in 0..4 {
                total += comps[g];
                if comps[g] > comps[maxind] {
                    maxind = g;
                }
            }
            let percentmatch = comps[maxind] as f32 * 100.0 / total as f32;
            return percentmatch >= 80.0;
        }

        let pattern = &self.sequence[lowerindex as usize..(lowerindex + length) as usize];
        let mut topperiods = [0i32; NUMBER_OF_PERIODS];
        if get_top_periods(
            &self.index,
            pattern,
            length,
            self.maxdistance,
            &mut topperiods,
        ) != 0
        {
            eprint!("\nUnable to allocate counts array in GetTopPeriods()!!!");
            self.fatal(-1);
        }

        for g in 0..NUMBER_OF_PERIODS_INTO_SORTMULTIPLES {
            self.sort_multiples[g] = topperiods[g];
        }
        for g in 0..NUMBER_OF_PERIODS_TO_TEST {
            if found_d == topperiods[g] {
                return true;
            }
        }
        false
    }

    /// `new_meet_criteria_3()`.
    pub fn new_meet_criteria_3(&mut self, d: i32, location: i32, tuplesize: i32) -> bool {
        let info = &self.distance[d as usize];
        let min_krun_matches = info.k_run_sums_criteria;
        let max_first_match_location = max2(0, location - max2(d, 20)) + info.waiting_time_criteria;
        let low_end_of_range = info.lo_d_range;
        let high_end_of_range = info.hi_d_range;
        let main_d_matches = info.nummatches;
        let range_d_min_for_waiting_time_test =
            (0.35 * min2(min_krun_matches, main_d_matches) as f64) as i32;
        let main_d_first_match_location = info.entry[info.lowindex as usize].location
            - info.entry[info.lowindex as usize].size
            + tuplesize;

        let mut waiting_time_ok = false;
        let mut waiting_time_d = 0;
        if main_d_first_match_location <= max_first_match_location {
            waiting_time_ok = true;
            waiting_time_d = d;
        }

        if main_d_matches >= min_krun_matches && waiting_time_ok {
            return true;
        }

        /* look in the low range to see if d has the most matches */
        let mut hipointer;
        let mut lopointer = d;
        let mut t = self.distance[d as usize].linkdown;
        let mut m = main_d_matches;
        let mut d_still_best = true;
        while t >= low_end_of_range && d_still_best {
            let s = self.distance[t as usize].linkdown;
            if self.no_matches_so_unlink_distance(t, location) == 0 {
                let range_d_matches = self.distance[t as usize].nummatches;
                if range_d_matches > main_d_matches {
                    d_still_best = false;
                } else {
                    lopointer = t;
                    m += range_d_matches;
                    if !waiting_time_ok && range_d_matches >= range_d_min_for_waiting_time_test {
                        let info = &self.distance[t as usize];
                        let loc = info.entry[info.lowindex as usize].location
                            - info.entry[info.lowindex as usize].size
                            + tuplesize;
                        if loc <= max_first_match_location {
                            waiting_time_ok = true;
                            waiting_time_d = t;
                        }
                    }
                }
            }
            t = s;
        }

        if !d_still_best {
            return false;
        }

        /* now check that d is the best in the upper range */
        t = self.distance[d as usize].linkup;
        while t <= high_end_of_range && d_still_best {
            let s = self.distance[t as usize].linkup;
            if self.no_matches_so_unlink_distance(t, location) == 0
                && self.distance[t as usize].nummatches > main_d_matches
            {
                d_still_best = false;
            }
            t = s;
        }

        if !d_still_best {
            return false;
        }

        if m >= min_krun_matches && waiting_time_ok {
            return true;
        }

        /* lower range didn't work, now test higher ranges */
        let d_range = d - low_end_of_range + 1;
        t = self.distance[d as usize].linkup;
        while t <= high_end_of_range {
            let s = self.distance[t as usize].linkup;
            hipointer = t;
            let range_d_matches = self.distance[t as usize].nummatches;
            m += range_d_matches;
            while lopointer < hipointer - d_range + 1 {
                m -= self.distance[lopointer as usize].nummatches;
                if m <= 0 {
                    self.msg("\n*** error in meet crieria, m<=0");
                    self.fatal(-16);
                }
                if waiting_time_ok && lopointer == waiting_time_d {
                    waiting_time_ok = false;
                    waiting_time_d = 0;
                }
                lopointer = self.distance[lopointer as usize].linkup;
            }
            if range_d_matches >= range_d_min_for_waiting_time_test {
                let info = &self.distance[t as usize];
                let loc = info.entry[info.lowindex as usize].location
                    - info.entry[info.lowindex as usize].size
                    + tuplesize;
                if loc <= max_first_match_location {
                    waiting_time_ok = true;
                    waiting_time_d = t;
                }
            }
            if m >= min_krun_matches && waiting_time_ok {
                return true;
            }
            t = s;
        }

        false
    }

    /* ------------------------------------------------------------------ */
    /* statistics                                                          */
    /* ------------------------------------------------------------------ */

    /// `init_and_fill_coin_toss_stats2000_with_4tuplesizes()`.
    pub fn init_and_fill_coin_toss_stats(&mut self) {
        let text = format!(
            "\nPmatch={:3.2},Pindel={:3.2}",
            self.pm as f32 / 100.0,
            self.pi as f32 / 100.0
        );
        self.msg(&text);

        self.pindel = (self.pi as f32 / 100.0f32) as f64;

        for g in 1..=self.maxdistance {
            if g <= SMALLDISTANCE {
                self.distance[g as usize].lo_d_range = g;
                self.distance[g as usize].hi_d_range = g;
            } else {
                self.distance[g as usize].lo_d_range = max2(g - self.d_range(g), 1);
                self.distance[g as usize].hi_d_range = min2(g + self.d_range(g), self.maxdistance);
            }
        }

        let waitdata: &[i32; 2001];
        let sumdata: &[i32; 2001];
        if self.pm == 80 {
            self.nts = 3;
            self.tuple_size[0] = 0;
            self.tuple_size[1] = 4;
            self.tuple_size[2] = 5;
            self.tuple_size[3] = 7;
            self.msg("\ntuple sizes 0,4,5,7");
            self.tuple_max_distance[0] = 0;
            self.tuple_max_distance[1] = 29;
            self.tuple_max_distance[2] = 159;
            self.tuple_max_distance[3] = self.maxdistance;
            let text = format!("\ntuple distances 0, 29, 159, {}", self.maxdistance);
            self.msg(&text);
            waitdata = &WAITDATA80;
            sumdata = &SUMDATA80;
        } else if self.pm == 75 {
            self.nts = 4;
            self.tuple_size[0] = 0;
            self.tuple_size[1] = 3;
            self.tuple_size[2] = 4;
            self.tuple_size[3] = 5;
            self.tuple_size[4] = 7;
            self.msg("\ntuple sizes 0,3,4,5,7");
            self.tuple_max_distance[0] = 0;
            self.tuple_max_distance[1] = 29;
            self.tuple_max_distance[2] = 43;
            self.tuple_max_distance[3] = 159;
            self.tuple_max_distance[4] = self.maxdistance;
            let text = format!("\ntuple distances 0, 29, 43, 159, {}", self.maxdistance);
            self.msg(&text);
            waitdata = &WAITDATA75;
            sumdata = &SUMDATA75;
        } else {
            let text = format!("\nNo wait table file for PM={}", self.pm);
            self.msg(&text);
            eprint!("\nNo wait table file for PM={}", self.pm);
            self.fatal(-13);
        }

        for d in 1..=self.maxdistance {
            let idx = std::cmp::min(MAXPATTERNSIZECONSTANT as i32, d) as usize;
            self.distance[d as usize].waiting_time_criteria = waitdata[idx];
        }
        for d in 1..=self.maxdistance {
            let idx = std::cmp::min(MAXPATTERNSIZECONSTANT as i32, d) as usize;
            self.distance[d as usize].k_run_sums_criteria = sumdata[idx];
        }
    }

    /// Handles one candidate distance found while scanning tuples.
    fn process_candidate(&mut self, i: i32, d: i32, tuplesize: i32) {
        self.add_tuple_match_to_distance_entry(i, tuplesize, d);

        let found = self.search_for_distance_match_in_distanceseenarray(d, i);
        if found {
            return;
        }

        if !self.distance[d as usize].linked {
            self.link_distance_window(d);
        }

        if !(self.new_meet_criteria_3(d, i, tuplesize)
            && (d <= 250 || self.search_for_range_in_bestperiodlist(i, d)))
        {
            return;
        }

        self.wdp_count += 1.0;
        self.criteria_count[d as usize] += 1;
        self.rows = 0.0;
        if d <= SMALLDISTANCE {
            self.newwrap(i, d, false);
            self.cell_count[d as usize] += self.rows * d as f64;
            self.get_pair_alignment_with_copynumber(d);
        } else {
            let band = max2(MINBANDRADIUS, self.d_range(d));
            let bandforward = min2(2 * band, d / 3);
            self.narrowbandwrap(i, d, band, bandforward, false, RECENTERCRITERION);
            self.cell_count[d as usize] +=
                self.rows * (2 * max2(MINBANDRADIUS, self.d_range(d)) + 1) as f64;
            let traceband = min2(2 * max2(MINBANDRADIUS, self.d_range(d)), d / 3);
            self.get_narrowband_pair_alignment_with_copynumber(d, traceband, crate::state::LOCAL);
        }
        self.add_to_distanceseenarray(i, d, self.maxrealrow, self.maxscore);

        /* number of copies required, less restrictive for larger patterns */
        let copynumber = self.copynumber;
        if (d <= 50 && copynumber < 1.9)
            || (d > 50 && d <= 100 && copynumber < 1.9 - 0.002 * (d - 50) as f64)
            || (d > 100 && copynumber < 1.8)
        {
            return;
        }

        let pass_multiples_test = self.multiples_criteria_4(d);
        self.add_to_bestperiodlist(d);
        if !pass_multiples_test {
            return;
        }

        /* get consensus and realign */
        self.classlength = d;
        self.get_consensus(d);
        if self.cons_classlength != self.classlength {
            self.classlength = self.cons_classlength;
        }
        self.consensus_count[self.classlength as usize] += 1;
        self.rows = 0.0;
        if self.classlength <= SMALLDISTANCE {
            self.newwrap(i, self.classlength, true);
            self.cell_count[self.classlength as usize] += self.rows * d as f64;
            self.get_pair_alignment_with_copynumber(self.classlength);
        } else {
            let cl = self.classlength;
            let band = max2(MINBANDRADIUS, self.d_range(cl));
            let bandforward = min2(2 * band, cl / 3);
            self.narrowbandwrap(i, cl, band, bandforward, true, RECENTERCRITERION);
            self.cell_count[cl as usize] +=
                self.rows * (2 * max2(MINBANDRADIUS, self.d_range(cl)) + 1) as f64;
            let traceband = min2(2 * max2(MINBANDRADIUS, self.d_range(cl)), cl / 3);
            self.get_narrowband_pair_alignment_with_copynumber(cl, traceband, crate::state::LOCAL);
        }
        self.add_to_distanceseenarray(i, d, self.maxrealrow, self.maxscore);
        self.adjust_bestperiod_entry(d);

        let cl = self.classlength;
        let copynumber = self.copynumber;
        if (cl <= 50 && copynumber < 1.9)
            || (cl > 50 && cl <= 100 && copynumber < 1.9 - 0.002 * (d - 50) as f64)
            || (cl > 100 && copynumber < 1.8)
        {
            return;
        }

        if cl >= self.minsize && self.align_pair.score >= self.minscore {
            self.output_count += 1.0;
            self.outputsize_count[cl as usize] += 1;
            let text = format!("\nFound at i:{} original size:{} final size:{}", i, d, cl);
            self.msg(&text);
            self.get_statistics(cl);
        }
    }

    /// `newtupbo()`: the main function of the algorithm.
    pub fn newtupbo(&mut self) {
        self.params.percent = 0;
        if !self.params.ngs {
            self.set_progress_bar();
        }

        self.output_count = 0.0;

        let mintuplesize = self.tuple_size[1];
        let maxtuplesize = self.tuple_size[self.nts as usize];

        let nts = self.nts as usize;
        self.tuple_hash = vec![Vec::new(); nts + 1];
        self.history = vec![Vec::new(); nts + 1];
        for g in 1..=nts {
            let codes = FOUR_TO_THE[self.tuple_size[g] as usize] as usize;
            let hsize = 2 * (self.tuple_max_distance[g] + 1) + 2;
            self.tuple_hash[g] = vec![0i32; codes];
            self.history_size[g] = hsize;
            self.history[g] = vec![HistoryEntry::default(); hsize as usize];
            self.next_free_history_index[g] = 1;
        }

        self.sort_multiples = vec![0i32; self.maxdistance as usize + 1];

        let mut build_entire_code = true;
        let onepercent = self.length / 100;
        let mut percentincrease = 0;
        let mut progbarpos = 0;

        let mut i: i32 = 0;
        let mut badcharindex: i32 = 0;
        let mut code: i32 = 0;

        while i <= self.length {
            percentincrease += 1;
            if percentincrease == onepercent {
                percentincrease = 0;
                progbarpos += 1;
                self.params.percent = progbarpos;
                if !self.params.ngs {
                    self.set_progress_bar();
                }
            }

            if i == 0 || !is_acgt(self.sequence[i as usize]) {
                badcharindex = i;
                build_entire_code = true;
                /* find first good string of mintuplesize characters */
                let mut g = 0;
                while g < mintuplesize && i < self.length {
                    i += 1;
                    if !is_acgt(self.sequence[i as usize]) {
                        badcharindex = i;
                        g = 0;
                    } else {
                        g += 1;
                    }
                }
                if g < mintuplesize {
                    break;
                }
            }

            if build_entire_code {
                code = 0;
                let mut g = badcharindex + 1;
                while g <= i {
                    code = code * 4 + self.index[self.sequence[g as usize] as usize];
                    g += 1;
                }
                if i - badcharindex >= maxtuplesize {
                    build_entire_code = false;
                }
            } else {
                code = (code % FOUR_TO_THE[(self.tuple_size[self.nts as usize] - 1) as usize]) * 4
                    + self.index[self.sequence[i as usize] as usize];
            }

            self.tuple_code[self.nts as usize] = code;
            for h in (1..self.nts as usize).rev() {
                self.tuple_code[h] = code % FOUR_TO_THE[self.tuple_size[h] as usize];
            }

            /* process index i using all the tuple sizes */
            let mut g = 1usize;
            while g <= self.nts as usize && i - badcharindex >= self.tuple_size[g] {
                let y = self.tuple_hash[g][self.tuple_code[g] as usize];
                let h = self.next_free_history_index[g];
                let mut j = h + 1;
                if j == self.history_size[g] {
                    j = 1;
                }

                if self.history[g][j as usize].location != 0
                    && j == self.tuple_hash[g][self.history[g][j as usize].code as usize]
                {
                    let c = self.history[g][j as usize].code as usize;
                    self.tuple_hash[g][c] = 0;
                }
                self.next_free_history_index[g] = j;

                self.tuple_hash[g][self.tuple_code[g] as usize] = h;
                self.history[g][h as usize].location = i;
                self.history[g][h as usize].previous = y;
                self.history[g][h as usize].code = self.tuple_code[g];

                let mut yy = h;
                let mut y = y;
                while y != 0 {
                    let d = i - self.history[g][y as usize].location;
                    if d > self.tuple_max_distance[g] {
                        self.history[g][yy as usize].previous = 0;
                        y = 0;
                    } else {
                        yy = y;
                        y = self.history[g][y as usize].previous;
                        if d > self.tuple_max_distance[g - 1] {
                            let tuplesize = self.tuple_size[g];
                            self.process_candidate(i, d, tuplesize);
                        }
                    }
                }
                g += 1;
            }

            i += 1;
        }

        /* close progress indicator */
        self.params.percent = -1;
        if !self.params.ngs {
            self.set_progress_bar();
        }
    }

    /// `get_statistics()`: prints the alignment and records the repeat.
    pub fn get_statistics(&mut self, consensussize: i32) {
        let size = consensussize;

        for g in 1..=self.maxdistance + self.d_range(self.maxdistance) {
            self.statistics_distance[g as usize] = 0;
        }

        let mut matchcount = 0i32;
        let mut mismatch = 0i32;
        let mut indel = 0i32;

        let alen = self.align_pair.length;

        let mut lp: i32 = 1;
        while self.align_pair.textsecnd[lp as usize] == DASH {
            lp += 1;
            if lp > alen {
                self.msg("\nError in statistics.");
                self.msg("\nInitial left pointer exceeds AlignPair.length while");
                self.msg("\nlooking for first non -");
                self.fatal(-10);
            }
        }
        let mut rp = lp + 1;
        if rp > alen {
            self.msg("\nError in statistics.");
            self.msg("\nInitial right pointer exceeds AlignPair.length");
            self.fatal(-10);
        }
        while self.align_pair.indexsecnd[rp as usize] != self.align_pair.indexsecnd[lp as usize] {
            rp += 1;
            if rp > alen {
                self.msg("\nError in statistics.");
                self.msg(
                    "\nInitial right pointer exceeds AlignPair.length while\nlooking for AlignPair.indexsecnd[lp]",
                );
                self.fatal(-10);
            }
        }
        while self.align_pair.textsecnd[rp as usize] == DASH {
            rp += 1;
            if rp > alen {
                self.msg("\nError in statistics.");
                self.msg("\nInitial right pointer exceeds AlignPair.length while looking for first non -");
                self.fatal(-10);
            }
        }
        if self.align_pair.indexsecnd[lp as usize] != self.align_pair.indexsecnd[rp as usize] {
            self.msg("\nError in statistics.");
            self.msg("\nInitial left pointer index not the same as\ninitial right pointer index");
            self.fatal(-10);
        }

        let mut mindistance = size;
        let mut maxdistance = size;

        while rp <= alen && lp < rp {
            let tp_lp = self.align_pair.textprime[lp as usize];
            let tp_rp = self.align_pair.textprime[rp as usize];
            let ts_lp = self.align_pair.textsecnd[lp as usize];
            let ts_rp = self.align_pair.textsecnd[rp as usize];

            if ts_lp != DASH && ts_rp != DASH {
                if tp_lp != DASH && tp_rp != DASH {
                    if self.align_pair.indexsecnd[lp as usize]
                        != self.align_pair.indexsecnd[rp as usize]
                    {
                        self.msg("\nError in statistics.");
                        self.msg("\nLeft pointer index not the same as right pointer index");
                        self.msg("\non match or mismatch when nothing is dash.");
                        let text = format!(
                            "\nAlignPair.indexsecnd[{}]: {}, AlignPair.indexprime[{}]: {}\nAlignPair.indexsecnd[{}]: {}, AlignPair.indexprime[{}]: {}",
                            lp,
                            self.align_pair.indexsecnd[lp as usize],
                            lp,
                            self.align_pair.indexprime[lp as usize],
                            rp,
                            self.align_pair.indexsecnd[rp as usize],
                            rp,
                            self.align_pair.indexprime[rp as usize]
                        );
                        self.msg(&text);
                        self.fatal(-10);
                    }
                    if tp_lp == tp_rp {
                        matchcount += 1;
                        let d = (self.align_pair.indexprime[rp as usize]
                            - self.align_pair.indexprime[lp as usize])
                            .abs();
                        if d < 4 * self.maxdistance {
                            self.statistics_distance[d as usize] += 1;
                            if d < mindistance {
                                mindistance = d;
                            }
                            if d > maxdistance {
                                maxdistance = d;
                            }
                        }
                    } else {
                        mismatch += 1;
                    }
                    lp += 1;
                    rp += 1;
                } else if tp_lp == DASH && tp_rp == DASH {
                    /* do nothing */
                    if self.align_pair.indexsecnd[lp as usize]
                        != self.align_pair.indexsecnd[rp as usize]
                    {
                        self.msg("\nError in statistics.");
                        self.msg("\nLeft pointer index not the same as right pointer index");
                        let text = format!("\non do nothing.\nlp: {},  rp: {}", lp, rp);
                        self.msg(&text);
                        self.fatal(-10);
                    }
                    lp += 1;
                    rp += 1;
                } else {
                    if self.align_pair.indexsecnd[lp as usize]
                        != self.align_pair.indexsecnd[rp as usize]
                    {
                        self.msg("\nError in statistics.");
                        self.msg("\nLeft pointer index not the same as right pointer index");
                        self.msg("\non indel caused by dash in textprime.");
                        let text = format!(
                            "\nAlignPair.indexsecnd[{}]: {}, AlignPair.indexprime[{}]: {}\nAlignPair.indexsecnd[{}]: {}, AlignPair.indexprime[{}]: {}",
                            lp,
                            self.align_pair.indexsecnd[lp as usize],
                            lp,
                            self.align_pair.indexprime[lp as usize],
                            rp,
                            self.align_pair.indexsecnd[rp as usize],
                            rp,
                            self.align_pair.indexprime[rp as usize]
                        );
                        self.msg(&text);
                        self.fatal(-10);
                    }
                    indel += 1;
                    lp += 1;
                    rp += 1;
                }
            } else if ts_lp == DASH && ts_rp == DASH {
                if self.align_pair.indexsecnd[lp as usize]
                    != self.align_pair.indexsecnd[rp as usize]
                {
                    self.msg("\nError in statistics.");
                    self.msg("\nLeft pointer index not the same as right pointer index");
                    self.msg("\non match or mismatch when both secnd are dash.");
                    let text = format!("\nlp: {},  rp: {}", lp, rp);
                    self.msg(&text);
                    self.fatal(-10);
                }
                if tp_lp == tp_rp {
                    matchcount += 1;
                    let d = (self.align_pair.indexprime[rp as usize]
                        - self.align_pair.indexprime[lp as usize])
                        .abs();
                    if d < 4 * self.maxdistance {
                        self.statistics_distance[d as usize] += 1;
                        if d < mindistance {
                            mindistance = d;
                        }
                        if d > maxdistance {
                            maxdistance = d;
                        }
                    }
                } else {
                    mismatch += 1;
                }
                lp += 1;
                rp += 1;
            } else if ts_lp == DASH {
                indel += 1;
                lp += 1;
            } else if ts_rp == DASH {
                indel += 1;
                rp += 1;
            }
        }
        if lp >= rp {
            self.msg("\nError in statistics.");
            self.msg("\nLeft pointer >= right pointer");
            self.fatal(-10);
        }

        let x = matchcount + mismatch + indel;

        let mut best_match_distance = 0;
        let mut best_match_count = 0;
        for g in mindistance..=maxdistance {
            if self.statistics_distance[g as usize] != 0
                && self.statistics_distance[g as usize] > best_match_count
            {
                best_match_count = self.statistics_distance[g as usize];
                best_match_distance = g;
            }
        }

        self.reverse();
        self.period = best_match_distance;

        if !self.params.htmloff {
            let cl = self.classlength;
            self.print_alignment_headings(cl);
        }

        /* save the starting position of consensus in EC */
        let start_ecpos = self.align_pair.indexsecnd[1];

        if !self.params.htmloff {
            let cl = self.classlength;
            self.alt3_print_alignment(cl);
        }

        self.reverse();

        if !self.params.htmloff {
            let mut out = String::new();
            out.push_str("\nStatistics");
            out.push_str(&format!(
                "\nMatches: {},  Mismatches: {}, Indels: {}",
                matchcount, mismatch, indel
            ));
            out.push_str(&format!(
                "\n        {:0.2}            {:0.2}        {:0.2}",
                matchcount as f32 / x as f32,
                mismatch as f32 / x as f32,
                indel as f32 / x as f32
            ));
            out.push('\n');
            out.push_str("\nMatches are distributed among these distances:");
            for g in mindistance..=maxdistance {
                if self.statistics_distance[g as usize] != 0 {
                    out.push_str(&format!(
                        "\n {:3}  {:3}  {:0.2}",
                        g,
                        self.statistics_distance[g as usize],
                        self.statistics_distance[g as usize] as f32 / matchcount as f32
                    ));
                }
            }
            self.msg(&out);
        }

        let mut acgtcount = [0i32; 26];
        let mut count = 0i32;
        for i in 1..=alen {
            let b = self.align_pair.textprime[i as usize];
            if b != DASH && (b'A'..=b'Z').contains(&b) {
                acgtcount[(b - b'A') as usize] += 1;
                count += 1;
            }
        }

        if count != 0 {
            let text = format!(
                "\n\nACGTcount: A:{:3.2}, C:{:3.2}, G:{:3.2}, T:{:3.2}",
                acgtcount[base_index(b'A')] as f64 / count as f64,
                acgtcount[base_index(b'C')] as f64 / count as f64,
                acgtcount[base_index(b'G')] as f64 / count as f64,
                acgtcount[base_index(b'T')] as f64 / count as f64
            );
            self.msg(&text);
        } else {
            self.msg("\nError in statistics: ACGTcount=0");
        }

        let mut diversity = [
            acgtcount[base_index(b'A')] as f64 / count as f64,
            acgtcount[base_index(b'C')] as f64 / count as f64,
            acgtcount[base_index(b'G')] as f64 / count as f64,
            acgtcount[base_index(b'T')] as f64 / count as f64,
        ];
        for e in (0..=2).rev() {
            for f in 0..=e {
                if diversity[f] < diversity[f + 1] {
                    diversity.swap(f, f + 1);
                }
            }
        }

        let mut entropy = 0.0f64;
        for d in diversity.iter() {
            if *d != 0.0 {
                entropy += d * (d.ln() / std::f64::consts::LN_2);
            }
        }
        if entropy < 0.0 {
            entropy = -entropy;
        }

        if !self.params.htmloff {
            self.msg("\n");
            self.msg("\n");
        }

        if !self.params.htmloff {
            self.print_ec_to_alignments(start_ecpos, consensussize);

            if self.print_flanking {
                self.reverse();
                let fl = self.params.flankinglength as i32;
                self.print_flanking_sequence(fl);
                self.reverse();
            }
        }

        /* add the record */
        self.counter_in_seq += 1;
        let reference = format!(
            "{}--{},{},{:3.1},{},{}",
            self.align_pair.indexprime[self.align_pair.length as usize],
            self.align_pair.indexprime[1],
            best_match_distance,
            self.copynumber,
            self.classlength,
            self.output_count as i32
        );
        let pattern = self.print_ec_to_buffer(start_ecpos, consensussize);
        let record = IlRecord {
            count: self.counter_in_seq,
            reference,
            first: self.align_pair.indexprime[self.align_pair.length as usize],
            last: self.align_pair.indexprime[1],
            period: best_match_distance,
            copies: self.copynumber as f32,
            size: self.classlength,
            matches: (100.0 * matchcount as f32 / x as f32) as i32,
            indels: (100.0 * indel as f32 / x as f32) as i32,
            score: self.align_pair.score,
            acount: (100.0 * acgtcount[0] as f64 / count as f64) as i32,
            ccount: (100.0 * acgtcount[2] as f64 / count as f64) as i32,
            gcount: (100.0 * acgtcount[6] as f64 / count as f64) as i32,
            tcount: (100.0 * acgtcount[19] as f64 / count as f64) as i32,
            entropy: entropy as f32,
            pattern,
        };
        self.records.push(record);
    }
}
