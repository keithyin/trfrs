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

use crate::clean::RepeatRecord;
use crate::fasta::base_index;
use crate::state::{
    BestPeriodElement, HistoryEntry, Trf, DASH, FOUR_TO_THE, MAX_PATTERN_SIZE_CONSTANT, MIN_BAND_RADIUS,
    NUMBER_OF_PERIODS, NUMBER_OF_PERIODS_INTO_SORTMULTIPLES, NUMBER_OF_PERIODS_TO_TEST,
    RECENTER_CRITERION, SMALL_DISTANCE, TAG_SEPARATION,
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
fn is_acgt(base: u8) -> bool {
    matches!(base, b'a' | b'c' | b'g' | b't' | b'A' | b'C' | b'G' | b'T')
}

/// `GetTopPeriods()`.
fn get_top_periods(
    nucleotide_index: &[i32; 256],
    pattern: &[u8],
    length: i32,
    max_distance: i32,
    top_periods: &mut [i32],
) -> i32 {
    let mut counts = vec![0f64; length as usize];
    let mut history = vec![0i32; length as usize];
    let mut heads = [-1i32; 16];

    let mut position = 0i32;
    let end = length - 2;
    while position <= end {
        let tuple_id =
            nucleotide_index[pattern[position as usize] as usize] * 4 + nucleotide_index[pattern[(position + 1) as usize] as usize];
        history[position as usize] = heads[tuple_id as usize];
        heads[tuple_id as usize] = position;

        let mut distance = 0;
        let mut current = position;
        while history[current as usize] != -1 && distance < MAX_PATTERN_SIZE_CONSTANT as i32 * 3 {
            distance = position - history[current as usize];
            counts[distance as usize] += 1.0;
            current = history[current as usize];
        }
        position += 1;
    }

    let mut xy_sum = 0f64;
    let mut x_sum = 0f64;
    let mut y_sum = 0f64;
    let mut x2_sum = 0f64;
    let end = length - 2;
    for position in 1..=end {
        xy_sum += position as f64 * counts[position as usize];
        x_sum += position as f64;
        y_sum += counts[position as usize];
        x2_sum += (position * position) as f64;
    }
    let point_count = end as f64;
    let slope = (point_count * xy_sum - x_sum * y_sum) / (point_count * x2_sum - x_sum * x_sum);

    for position in 1..=end {
        counts[position as usize] -= position as f64 * slope;
    }

    let mut end = length - 2;
    if end > max_distance {
        /* accepts smaller multiples if the best ones are too large */
        end = max_distance;
    }
    for slot in 0..NUMBER_OF_PERIODS {
        let mut top_index = 0;
        let mut top_value = 0.0;
        for position in 1..=end {
            if counts[position as usize] > top_value {
                top_index = position;
                top_value = counts[position as usize];
            }
        }
        top_periods[slot] = top_index;
        counts[top_index as usize] = 0.0;
    }

    0
}

impl Trf {
    /* ------------------------------------------------------------------ */
    /* distance lists                                                      */
    /* ------------------------------------------------------------------ */

    /// `add_tuple_match_to_Distance_entry()`.
    pub fn add_tuple_match_to_distance_entry(
        &mut self,
        location: i32,
        tuple_size: i32,
        distance: i32,
    ) {
        let window_size = max2(distance, 20) + 1;
        let list = &mut self.distance[distance as usize];

        if list.entry_count != 0 {
            let window_left_end = min2(location - distance, location - 20) + 1;
            loop {
                if list.entry_count == 0 {
                    break;
                }
                let lo = list.low_index as usize;
                if list.entries[lo].location < window_left_end {
                    let entry_size = list.entries[lo].size;
                    list.entry_count -= 1;
                    list.match_count -= entry_size;
                    list.low_index += 1;
                    list.low_index %= window_size;
                } else {
                    break;
                }
            }
        }

        let hi = list.high_index as usize;
        if list.entry_count != 0 && list.entries[hi].location == location - 1 {
            /* adjacent tuple: extend the last entry */
            list.entries[hi].location += 1;
            list.entries[hi].size += 1;
            list.match_count += 1;
        } else {
            list.entry_count += 1;
            list.high_index += 1;
            list.high_index %= window_size;
            let hi = list.high_index as usize;
            list.entries[hi].location = location;
            list.entries[hi].size = tuple_size;
            list.match_count += tuple_size;
        }
    }

    /// `link_Distance_window()`.
    pub fn link_distance_window(&mut self, distance: i32) {
        /* the reference implementation applies ceil() to an integer division */
        let tag_bucket = distance / TAG_SEPARATION;
        let mut lower_distance;
        if self.tag[tag_bucket as usize] < distance {
            lower_distance = self.tag[tag_bucket as usize];
            let mut tag_bucket_cursor = tag_bucket;
            while tag_bucket_cursor <= self.top_tag && self.tag[tag_bucket_cursor as usize] < distance {
                self.tag[tag_bucket_cursor as usize] = distance;
                tag_bucket_cursor += 1;
            }
        } else if self.tag[tag_bucket as usize] > distance {
            lower_distance = self.tag[tag_bucket as usize];
            while lower_distance > distance {
                lower_distance = self.distance[lower_distance as usize].link_down;
            }
            if lower_distance == distance {
                let text = format!("\nTag error following links.  f==d={}", distance);
                self.write_to_alignment_file(&text);
            }
        } else {
            let text = format!("\nTag error Tag[{}]={}", tag_bucket, distance);
            self.write_to_alignment_file(&text);
            self.fatal(-2);
        }

        /* link in distance */
        self.distance[distance as usize].link_down = lower_distance;
        let higher_distance = self.distance[lower_distance as usize].link_up;
        self.distance[distance as usize].link_up = higher_distance;
        self.distance[lower_distance as usize].link_up = distance;
        if higher_distance <= self.max_distance {
            self.distance[higher_distance as usize].link_down = distance;
        }
        self.distance[distance as usize].linked = true;
    }

    /// `untag_Distance_window()`.
    pub fn untag_distance_window(&mut self, distance: i32, link_down: i32) {
        let tag_bucket = distance / TAG_SEPARATION;
        if self.tag[tag_bucket as usize] != distance {
            return;
        }
        let mut tag_bucket_cursor = tag_bucket;
        while tag_bucket_cursor <= self.top_tag && self.tag[tag_bucket_cursor as usize] == distance {
            self.tag[tag_bucket_cursor as usize] = link_down;
            tag_bucket_cursor += 1;
        }
    }

    /// `no_matches_so_unlink_Distance()`; returns 1 when the distance was unlinked.
    pub fn no_matches_so_unlink_distance(&mut self, distance: i32, location: i32) -> i32 {
        let window_size = max2(distance, 20) + 1;
        let mut entry_count = self.distance[distance as usize].entry_count;
        if entry_count != 0 {
            let window_left_end = min2(location - distance, location - 20) + 1;
            loop {
                if entry_count == 0 {
                    break;
                }
                let list = &mut self.distance[distance as usize];
                let lo = list.low_index as usize;
                if list.entries[lo].location < window_left_end {
                    let entry_size = list.entries[lo].size;
                    list.entry_count -= 1;
                    list.match_count -= entry_size;
                    list.low_index += 1;
                    list.low_index %= window_size;
                    entry_count = list.entry_count;
                } else {
                    break;
                }
            }
        }
        if self.distance[distance as usize].entry_count == 0 {
            let lower_distance = self.distance[distance as usize].link_down;
            let higher_distance = self.distance[distance as usize].link_up;

            self.distance[lower_distance as usize].link_up = higher_distance;
            if higher_distance <= self.max_distance {
                self.distance[higher_distance as usize].link_down = lower_distance;
            }
            self.distance[distance as usize].linked = false;

            self.untag_distance_window(distance, lower_distance);
            return 1;
        }
        0
    }

    /* ------------------------------------------------------------------ */
    /* distances seen                                                      */
    /* ------------------------------------------------------------------ */

    pub fn add_to_distance_seen(&mut self, location: i32, distance: i32, end: i32, score: i32) {
        let ptr = &mut self.distance_seen[distance as usize];
        ptr.location = location;
        ptr.end = end;
        ptr.score = score;
    }

    pub fn search_for_distance_match_in_distance_seen(
        &self,
        distance: i32,
        start: i32,
    ) -> bool {
        self.distance_seen[distance as usize].end >= start
    }

    /* ------------------------------------------------------------------ */
    /* best period list                                                    */
    /* ------------------------------------------------------------------ */

    pub fn init_best_period_list(&mut self) {
        self.best_period_list = VecDeque::new();
    }

    /// `add_to_best_period_list()`.
    pub fn add_to_best_period_list(&mut self, distance: i32) {
        if distance != 1 {
            let entry = BestPeriodElement {
                index_low: self.align_pair.sequence_index[self.align_pair.length as usize],
                index_high: self.align_pair.sequence_index[1],
                best_1: self.sort_multiples[0],
                best_2: self.sort_multiples[1],
                best_3: self.sort_multiples[2],
                best_4: self.sort_multiples[3],
                best_5: self.sort_multiples[4],
            };
            self.best_period_list.push_front(entry);
        }
    }

    /// `adjust_best_period_entry()`.
    pub fn adjust_best_period_entry(&mut self, distance: i32) {
        if distance != 1 {
            if let Some(entry) = self.best_period_list.front_mut() {
                if entry.index_high > self.align_pair.sequence_index[1] {
                    entry.index_high = self.align_pair.sequence_index[1];
                }
            }
        }
    }

    /// `search_for_range_in_best_period_list()`.
    pub fn search_for_range_in_best_period_list(&mut self, start: i32, distance: i32) -> bool {
        let mut range_covered = false;
        let mut entry_index = 0usize;
        while entry_index < self.best_period_list.len() {
            let entry = self.best_period_list[entry_index];
            if entry.index_high < start - 2 * self.max_distance {
                self.best_period_list.remove(entry_index);
            } else {
                if entry.index_low
                    <= start - 2 * distance
                        + 1
                        + self.distance[distance as usize].waiting_time_criteria
                    && entry.index_high >= start
                {
                    range_covered = true;
                    if entry.best_1 == distance
                        || entry.best_2 == distance
                        || entry.best_3 == distance
                        || entry.best_4 == distance
                        || entry.best_5 == distance
                    {
                        return true;
                    }
                }
                entry_index += 1;
            }
        }
        !range_covered
    }

    /* ------------------------------------------------------------------ */
    /* criteria                                                            */
    /* ------------------------------------------------------------------ */

    /// `passes_top_period_test()`.
    pub fn passes_top_period_test(&mut self, found_distance: i32) -> bool {
        let lowerindex = self.align_pair.sequence_index[self.align_pair.length as usize];
        let upperindex = self.align_pair.sequence_index[1];
        let length = upperindex - lowerindex + 1;

        if found_distance == 1 {
            let mut comps = [0i32; 4];
            let mut total = 0;
            let mut maxind = 0usize;
            for offset in 0..length {
                let base = self.sequence[(lowerindex + offset) as usize];
                comps[self.nucleotide_index[base as usize] as usize] += 1;
            }
            for base in 0..4 {
                total += comps[base];
                if comps[base] > comps[maxind] {
                    maxind = base;
                }
            }
            let percent_match = comps[maxind] as f32 * 100.0 / total as f32;
            return percent_match >= 80.0;
        }

        let pattern = &self.sequence[lowerindex as usize..(lowerindex + length) as usize];
        let mut top_periods = [0i32; NUMBER_OF_PERIODS];
        if get_top_periods(
            &self.nucleotide_index,
            pattern,
            length,
            self.max_distance,
            &mut top_periods,
        ) != 0
        {
            eprint!("\nUnable to allocate counts array in GetTopPeriods()!!!");
            self.fatal(-1);
        }

        for slot in 0..NUMBER_OF_PERIODS_INTO_SORTMULTIPLES {
            self.sort_multiples[slot] = top_periods[slot];
        }
        for slot in 0..NUMBER_OF_PERIODS_TO_TEST {
            if found_distance == top_periods[slot] {
                return true;
            }
        }
        false
    }

    /// `distance_meets_criteria()`.
    pub fn distance_meets_criteria(&mut self, distance: i32, location: i32, tuple_size: i32) -> bool {
        let info = &self.distance[distance as usize];
        let min_krun_matches = info.k_run_sums_criteria;
        let max_first_match_location = max2(0, location - max2(distance, 20)) + info.waiting_time_criteria;
        let low_end_of_range = info.min_distance_range;
        let high_end_of_range = info.max_distance_range;
        let main_d_matches = info.match_count;
        let range_d_min_for_waiting_time_test =
            (0.35 * min2(min_krun_matches, main_d_matches) as f64) as i32;
        let main_d_first_match_location = info.entries[info.low_index as usize].location
            - info.entries[info.low_index as usize].size
            + tuple_size;

        let mut waiting_time_ok = false;
        let mut waiting_time_d = 0;
        if main_d_first_match_location <= max_first_match_location {
            waiting_time_ok = true;
            waiting_time_d = distance;
        }

        if main_d_matches >= min_krun_matches && waiting_time_ok {
            return true;
        }

        /* look in the low range to see if distance has the most matches */
        let mut high_pointer;
        let mut low_pointer = distance;
        let mut candidate = self.distance[distance as usize].link_down;
        let mut accumulated_matches = main_d_matches;
        let mut d_still_best = true;
        while candidate >= low_end_of_range && d_still_best {
            let next_candidate = self.distance[candidate as usize].link_down;
            if self.no_matches_so_unlink_distance(candidate, location) == 0 {
                let range_d_matches = self.distance[candidate as usize].match_count;
                if range_d_matches > main_d_matches {
                    d_still_best = false;
                } else {
                    low_pointer = candidate;
                    accumulated_matches += range_d_matches;
                    if !waiting_time_ok && range_d_matches >= range_d_min_for_waiting_time_test {
                        let info = &self.distance[candidate as usize];
                        let loc = info.entries[info.low_index as usize].location
                            - info.entries[info.low_index as usize].size
                            + tuple_size;
                        if loc <= max_first_match_location {
                            waiting_time_ok = true;
                            waiting_time_d = candidate;
                        }
                    }
                }
            }
            candidate = next_candidate;
        }

        if !d_still_best {
            return false;
        }

        /* now check that distance is the best in the upper range */
        candidate = self.distance[distance as usize].link_up;
        while candidate <= high_end_of_range && d_still_best {
            let next_candidate = self.distance[candidate as usize].link_up;
            if self.no_matches_so_unlink_distance(candidate, location) == 0
                && self.distance[candidate as usize].match_count > main_d_matches
            {
                d_still_best = false;
            }
            candidate = next_candidate;
        }

        if !d_still_best {
            return false;
        }

        if accumulated_matches >= min_krun_matches && waiting_time_ok {
            return true;
        }

        /* lower range didn'candidate work, now test higher ranges */
        let indel_random_walk_range = distance - low_end_of_range + 1;
        candidate = self.distance[distance as usize].link_up;
        while candidate <= high_end_of_range {
            let next_candidate = self.distance[candidate as usize].link_up;
            high_pointer = candidate;
            let range_d_matches = self.distance[candidate as usize].match_count;
            accumulated_matches += range_d_matches;
            while low_pointer < high_pointer - indel_random_walk_range + 1 {
                accumulated_matches -= self.distance[low_pointer as usize].match_count;
                if accumulated_matches <= 0 {
                    self.write_to_alignment_file("\n*** error in meet crieria, m<=0");
                    self.fatal(-16);
                }
                if waiting_time_ok && low_pointer == waiting_time_d {
                    waiting_time_ok = false;
                    waiting_time_d = 0;
                }
                low_pointer = self.distance[low_pointer as usize].link_up;
            }
            if range_d_matches >= range_d_min_for_waiting_time_test {
                let info = &self.distance[candidate as usize];
                let loc = info.entries[info.low_index as usize].location
                    - info.entries[info.low_index as usize].size
                    + tuple_size;
                if loc <= max_first_match_location {
                    waiting_time_ok = true;
                    waiting_time_d = candidate;
                }
            }
            if accumulated_matches >= min_krun_matches && waiting_time_ok {
                return true;
            }
            candidate = next_candidate;
        }

        false
    }

    /* ------------------------------------------------------------------ */
    /* statistics                                                          */
    /* ------------------------------------------------------------------ */

    /// `init_and_fill_coin_toss_stats2000_with_4tuplesizes()`.
    pub fn init_criteria_tables(&mut self) {
        let text = format!(
            "\nPmatch={:3.2},Pindel={:3.2}",
            self.match_probability_percent as f32 / 100.0,
            self.indel_probability_percent as f32 / 100.0
        );
        self.write_to_alignment_file(&text);

        self.indel_probability = (self.indel_probability_percent as f32 / 100.0f32) as f64;

        for distance in 1..=self.max_distance {
            if distance <= SMALL_DISTANCE {
                self.distance[distance as usize].min_distance_range = distance;
                self.distance[distance as usize].max_distance_range = distance;
            } else {
                self.distance[distance as usize].min_distance_range = max2(distance - self.indel_random_walk_range(distance), 1);
                self.distance[distance as usize].max_distance_range = min2(distance + self.indel_random_walk_range(distance), self.max_distance);
            }
        }

        let waitdata: &[i32; 2001];
        let sumdata: &[i32; 2001];
        if self.match_probability_percent == 80 {
            self.num_tuple_sizes = 3;
            self.tuple_size[0] = 0;
            self.tuple_size[1] = 4;
            self.tuple_size[2] = 5;
            self.tuple_size[3] = 7;
            self.write_to_alignment_file("\ntuple sizes 0,4,5,7");
            self.tuple_max_distance[0] = 0;
            self.tuple_max_distance[1] = 29;
            self.tuple_max_distance[2] = 159;
            self.tuple_max_distance[3] = self.max_distance;
            let text = format!("\ntuple distances 0, 29, 159, {}", self.max_distance);
            self.write_to_alignment_file(&text);
            waitdata = &WAITDATA80;
            sumdata = &SUMDATA80;
        } else if self.match_probability_percent == 75 {
            self.num_tuple_sizes = 4;
            self.tuple_size[0] = 0;
            self.tuple_size[1] = 3;
            self.tuple_size[2] = 4;
            self.tuple_size[3] = 5;
            self.tuple_size[4] = 7;
            self.write_to_alignment_file("\ntuple sizes 0,3,4,5,7");
            self.tuple_max_distance[0] = 0;
            self.tuple_max_distance[1] = 29;
            self.tuple_max_distance[2] = 43;
            self.tuple_max_distance[3] = 159;
            self.tuple_max_distance[4] = self.max_distance;
            let text = format!("\ntuple distances 0, 29, 43, 159, {}", self.max_distance);
            self.write_to_alignment_file(&text);
            waitdata = &WAITDATA75;
            sumdata = &SUMDATA75;
        } else {
            let text = format!("\nNo wait table file for PM={}", self.match_probability_percent);
            self.write_to_alignment_file(&text);
            eprint!("\nNo wait table file for PM={}", self.match_probability_percent);
            self.fatal(-13);
        }

        for distance in 1..=self.max_distance {
            let idx = std::cmp::min(MAX_PATTERN_SIZE_CONSTANT as i32, distance) as usize;
            self.distance[distance as usize].waiting_time_criteria = waitdata[idx];
        }
        for distance in 1..=self.max_distance {
            let idx = std::cmp::min(MAX_PATTERN_SIZE_CONSTANT as i32, distance) as usize;
            self.distance[distance as usize].k_run_sums_criteria = sumdata[idx];
        }
    }

    /// Handles one candidate distance found while scanning tuples.
    fn process_candidate(&mut self, position: i32, distance: i32, tuple_size: i32) {
        self.add_tuple_match_to_distance_entry(position, tuple_size, distance);

        let found = self.search_for_distance_match_in_distance_seen(distance, position);
        if found {
            return;
        }

        if !self.distance[distance as usize].linked {
            self.link_distance_window(distance);
        }

        if !(self.distance_meets_criteria(distance, position, tuple_size)
            && (distance <= 250 || self.search_for_range_in_best_period_list(position, distance)))
        {
            return;
        }

        self.wraparound_dp_count += 1.0;
        self.criteria_count[distance as usize] += 1;
        self.rows = 0.0;
        if distance <= SMALL_DISTANCE {
            self.wraparound_dp(position, distance, false);
            self.cell_count[distance as usize] += self.rows * distance as f64;
            self.get_pair_alignment_with_copynumber(distance);
        } else {
            let band_radius = max2(MIN_BAND_RADIUS, self.indel_random_walk_range(distance));
            let band_radius_forward = min2(2 * band_radius, distance / 3);
            self.narrow_band_wraparound_dp(position, distance, band_radius, band_radius_forward, false, RECENTER_CRITERION);
            self.cell_count[distance as usize] +=
                self.rows * (2 * max2(MIN_BAND_RADIUS, self.indel_random_walk_range(distance)) + 1) as f64;
            let trace_band_radius = min2(2 * max2(MIN_BAND_RADIUS, self.indel_random_walk_range(distance)), distance / 3);
            self.get_narrow_band_pair_alignment_with_copynumber(distance, trace_band_radius, crate::state::LOCAL);
        }
        self.add_to_distance_seen(position, distance, self.max_real_row, self.max_score);

        /* number of copies required, less restrictive for larger patterns */
        let copynumber = self.copynumber;
        if (distance <= 50 && copynumber < 1.9)
            || (distance > 50 && distance <= 100 && copynumber < 1.9 - 0.002 * (distance - 50) as f64)
            || (distance > 100 && copynumber < 1.8)
        {
            return;
        }

        let pass_multiples_test = self.passes_top_period_test(distance);
        self.add_to_best_period_list(distance);
        if !pass_multiples_test {
            return;
        }

        /* get consensus and realign */
        self.class_length = distance;
        self.get_consensus(distance);
        if self.consensus_class_length != self.class_length {
            self.class_length = self.consensus_class_length;
        }
        self.consensus_count[self.class_length as usize] += 1;
        self.rows = 0.0;
        if self.class_length <= SMALL_DISTANCE {
            self.wraparound_dp(position, self.class_length, true);
            self.cell_count[self.class_length as usize] += self.rows * distance as f64;
            self.get_pair_alignment_with_copynumber(self.class_length);
        } else {
            let class_length = self.class_length;
            let band_radius = max2(MIN_BAND_RADIUS, self.indel_random_walk_range(class_length));
            let band_radius_forward = min2(2 * band_radius, class_length / 3);
            self.narrow_band_wraparound_dp(position, class_length, band_radius, band_radius_forward, true, RECENTER_CRITERION);
            self.cell_count[class_length as usize] +=
                self.rows * (2 * max2(MIN_BAND_RADIUS, self.indel_random_walk_range(class_length)) + 1) as f64;
            let trace_band_radius = min2(2 * max2(MIN_BAND_RADIUS, self.indel_random_walk_range(class_length)), class_length / 3);
            self.get_narrow_band_pair_alignment_with_copynumber(class_length, trace_band_radius, crate::state::LOCAL);
        }
        self.add_to_distance_seen(position, distance, self.max_real_row, self.max_score);
        self.adjust_best_period_entry(distance);

        let class_length = self.class_length;
        let copynumber = self.copynumber;
        if (class_length <= 50 && copynumber < 1.9)
            || (class_length > 50 && class_length <= 100 && copynumber < 1.9 - 0.002 * (distance - 50) as f64)
            || (class_length > 100 && copynumber < 1.8)
        {
            return;
        }

        if class_length >= self.minsize && self.align_pair.score >= self.minscore {
            self.output_count += 1.0;
            self.output_size_count[class_length as usize] += 1;
            let text = format!("\nFound at i:{} original size:{} final size:{}", position, distance, class_length);
            self.write_to_alignment_file(&text);
            self.report_repeat(class_length);
        }
    }

    /// `scan_sequence()`: the main function of the algorithm.
    pub fn scan_sequence(&mut self) {
        self.params.percent = 0;
        if !self.params.ngs {
            self.set_progress_bar();
        }

        self.output_count = 0.0;

        let min_tuple_size = self.tuple_size[1];
        let max_tuple_size = self.tuple_size[self.num_tuple_sizes as usize];

        let num_tuple_sizes = self.num_tuple_sizes as usize;
        self.tuple_hash = vec![Vec::new(); num_tuple_sizes + 1];
        self.history = vec![Vec::new(); num_tuple_sizes + 1];
        for tuple_index in 1..=num_tuple_sizes {
            let num_codes = FOUR_TO_THE[self.tuple_size[tuple_index] as usize] as usize;
            let history_capacity = 2 * (self.tuple_max_distance[tuple_index] + 1) + 2;
            self.tuple_hash[tuple_index] = vec![0i32; num_codes];
            self.history_size[tuple_index] = history_capacity;
            self.history[tuple_index] = vec![HistoryEntry::default(); history_capacity as usize];
            self.next_free_history_index[tuple_index] = 1;
        }

        self.sort_multiples = vec![0i32; self.max_distance as usize + 1];

        let mut build_entire_code = true;
        let one_percent = self.length / 100;
        let mut percent_increase = 0;
        let mut progress_position = 0;

        let mut position: i32 = 0;
        let mut bad_char_index: i32 = 0;
        let mut rolling_code: i32 = 0;

        while position <= self.length {
            percent_increase += 1;
            if percent_increase == one_percent {
                percent_increase = 0;
                progress_position += 1;
                self.params.percent = progress_position;
                if !self.params.ngs {
                    self.set_progress_bar();
                }
            }

            if position == 0 || !is_acgt(self.sequence[position as usize]) {
                bad_char_index = position;
                build_entire_code = true;
                /* find first good string of min_tuple_size characters */
                let mut good_run_length = 0;
                while good_run_length < min_tuple_size && position < self.length {
                    position += 1;
                    if !is_acgt(self.sequence[position as usize]) {
                        bad_char_index = position;
                        good_run_length = 0;
                    } else {
                        good_run_length += 1;
                    }
                }
                if good_run_length < min_tuple_size {
                    break;
                }
            }

            if build_entire_code {
                rolling_code = 0;
                let mut code_pos = bad_char_index + 1;
                while code_pos <= position {
                    rolling_code = rolling_code * 4 + self.nucleotide_index[self.sequence[code_pos as usize] as usize];
                    code_pos += 1;
                }
                if position - bad_char_index >= max_tuple_size {
                    build_entire_code = false;
                }
            } else {
                rolling_code = (rolling_code % FOUR_TO_THE[(self.tuple_size[self.num_tuple_sizes as usize] - 1) as usize]) * 4
                    + self.nucleotide_index[self.sequence[position as usize] as usize];
            }

            self.tuple_code[self.num_tuple_sizes as usize] = rolling_code;
            for smaller_tuple_index in (1..self.num_tuple_sizes as usize).rev() {
                self.tuple_code[smaller_tuple_index] = rolling_code % FOUR_TO_THE[self.tuple_size[smaller_tuple_index] as usize];
            }

            /* process this position using all the tuple sizes */
            let mut tuple_index = 1usize;
            while tuple_index <= self.num_tuple_sizes as usize && position - bad_char_index >= self.tuple_size[tuple_index] {
                let node = self.tuple_hash[tuple_index][self.tuple_code[tuple_index] as usize];
                let write_slot = self.next_free_history_index[tuple_index];
                let mut next_free_slot = write_slot + 1;
                if next_free_slot == self.history_size[tuple_index] {
                    next_free_slot = 1;
                }

                if self.history[tuple_index][next_free_slot as usize].location != 0
                    && next_free_slot == self.tuple_hash[tuple_index][self.history[tuple_index][next_free_slot as usize].code as usize]
                {
                    let evicted_code = self.history[tuple_index][next_free_slot as usize].code as usize;
                    self.tuple_hash[tuple_index][evicted_code] = 0;
                }
                self.next_free_history_index[tuple_index] = next_free_slot;

                self.tuple_hash[tuple_index][self.tuple_code[tuple_index] as usize] = write_slot;
                self.history[tuple_index][write_slot as usize].location = position;
                self.history[tuple_index][write_slot as usize].previous = node;
                self.history[tuple_index][write_slot as usize].code = self.tuple_code[tuple_index];

                let mut prev_node = write_slot;
                let mut node = node;
                while node != 0 {
                    let distance = position - self.history[tuple_index][node as usize].location;
                    if distance > self.tuple_max_distance[tuple_index] {
                        self.history[tuple_index][prev_node as usize].previous = 0;
                        node = 0;
                    } else {
                        prev_node = node;
                        node = self.history[tuple_index][node as usize].previous;
                        if distance > self.tuple_max_distance[tuple_index - 1] {
                            let tuple_size = self.tuple_size[tuple_index];
                            self.process_candidate(position, distance, tuple_size);
                        }
                    }
                }
                tuple_index += 1;
            }

            position += 1;
        }

        /* close progress indicator */
        self.params.percent = -1;
        if !self.params.ngs {
            self.set_progress_bar();
        }
    }

    /// `get_statistics()` in the reference implementation: prints the alignment
    /// to the output file and appends one [`RepeatRecord`] to `self.records`.
    pub fn report_repeat(&mut self, consensus_size: i32) {
        for distance in 1..=self.max_distance + self.indel_random_walk_range(self.max_distance) {
            self.statistics_distance[distance as usize] = 0;
        }

        let mut match_count = 0i32;
        let mut mismatch = 0i32;
        let mut indel = 0i32;

        let align_length = self.align_pair.length;

        let mut left_col: i32 = 1;
        while self.align_pair.pattern_text[left_col as usize] == DASH {
            left_col += 1;
            if left_col > align_length {
                self.write_to_alignment_file("\nError in statistics.");
                self.write_to_alignment_file("\nInitial left pointer exceeds AlignPair.length while");
                self.write_to_alignment_file("\nlooking for first non -");
                self.fatal(-10);
            }
        }
        let mut right_col = left_col + 1;
        if right_col > align_length {
            self.write_to_alignment_file("\nError in statistics.");
            self.write_to_alignment_file("\nInitial right pointer exceeds AlignPair.length");
            self.fatal(-10);
        }
        while self.align_pair.pattern_index[right_col as usize] != self.align_pair.pattern_index[left_col as usize] {
            right_col += 1;
            if right_col > align_length {
                self.write_to_alignment_file("\nError in statistics.");
                self.write_to_alignment_file(
                    "\nInitial right pointer exceeds AlignPair.length while\nlooking for AlignPair.pattern_index[lp]",
                );
                self.fatal(-10);
            }
        }
        while self.align_pair.pattern_text[right_col as usize] == DASH {
            right_col += 1;
            if right_col > align_length {
                self.write_to_alignment_file("\nError in statistics.");
                self.write_to_alignment_file("\nInitial right pointer exceeds AlignPair.length while looking for first non -");
                self.fatal(-10);
            }
        }
        if self.align_pair.pattern_index[left_col as usize] != self.align_pair.pattern_index[right_col as usize] {
            self.write_to_alignment_file("\nError in statistics.");
            self.write_to_alignment_file("\nInitial left pointer index not the same as\ninitial right pointer index");
            self.fatal(-10);
        }

        let mut min_distance = consensus_size;
        let mut max_distance = consensus_size;

        while right_col <= align_length && left_col < right_col {
            let left_seq_char = self.align_pair.sequence_text[left_col as usize];
            let right_seq_char = self.align_pair.sequence_text[right_col as usize];
            let left_pat_char = self.align_pair.pattern_text[left_col as usize];
            let right_pat_char = self.align_pair.pattern_text[right_col as usize];

            if left_pat_char != DASH && right_pat_char != DASH {
                if left_seq_char != DASH && right_seq_char != DASH {
                    if self.align_pair.pattern_index[left_col as usize]
                        != self.align_pair.pattern_index[right_col as usize]
                    {
                        self.write_to_alignment_file("\nError in statistics.");
                        self.write_to_alignment_file("\nLeft pointer index not the same as right pointer index");
                        self.write_to_alignment_file("\non match or mismatch when nothing is dash.");
                        let text = format!(
                            "\nAlignPair.pattern_index[{}]: {}, AlignPair.sequence_index[{}]: {}\nAlignPair.pattern_index[{}]: {}, AlignPair.sequence_index[{}]: {}",
                            left_col,
                            self.align_pair.pattern_index[left_col as usize],
                            left_col,
                            self.align_pair.sequence_index[left_col as usize],
                            right_col,
                            self.align_pair.pattern_index[right_col as usize],
                            right_col,
                            self.align_pair.sequence_index[right_col as usize]
                        );
                        self.write_to_alignment_file(&text);
                        self.fatal(-10);
                    }
                    if left_seq_char == right_seq_char {
                        match_count += 1;
                        let distance = (self.align_pair.sequence_index[right_col as usize]
                            - self.align_pair.sequence_index[left_col as usize])
                            .abs();
                        if distance < 4 * self.max_distance {
                            self.statistics_distance[distance as usize] += 1;
                            if distance < min_distance {
                                min_distance = distance;
                            }
                            if distance > max_distance {
                                max_distance = distance;
                            }
                        }
                    } else {
                        mismatch += 1;
                    }
                    left_col += 1;
                    right_col += 1;
                } else if left_seq_char == DASH && right_seq_char == DASH {
                    /* do nothing */
                    if self.align_pair.pattern_index[left_col as usize]
                        != self.align_pair.pattern_index[right_col as usize]
                    {
                        self.write_to_alignment_file("\nError in statistics.");
                        self.write_to_alignment_file("\nLeft pointer index not the same as right pointer index");
                        let text = format!("\non do nothing.\nlp: {},  rp: {}", left_col, right_col);
                        self.write_to_alignment_file(&text);
                        self.fatal(-10);
                    }
                    left_col += 1;
                    right_col += 1;
                } else {
                    if self.align_pair.pattern_index[left_col as usize]
                        != self.align_pair.pattern_index[right_col as usize]
                    {
                        self.write_to_alignment_file("\nError in statistics.");
                        self.write_to_alignment_file("\nLeft pointer index not the same as right pointer index");
                        self.write_to_alignment_file("\non indel caused by dash in sequence_text.");
                        let text = format!(
                            "\nAlignPair.pattern_index[{}]: {}, AlignPair.sequence_index[{}]: {}\nAlignPair.pattern_index[{}]: {}, AlignPair.sequence_index[{}]: {}",
                            left_col,
                            self.align_pair.pattern_index[left_col as usize],
                            left_col,
                            self.align_pair.sequence_index[left_col as usize],
                            right_col,
                            self.align_pair.pattern_index[right_col as usize],
                            right_col,
                            self.align_pair.sequence_index[right_col as usize]
                        );
                        self.write_to_alignment_file(&text);
                        self.fatal(-10);
                    }
                    indel += 1;
                    left_col += 1;
                    right_col += 1;
                }
            } else if left_pat_char == DASH && right_pat_char == DASH {
                if self.align_pair.pattern_index[left_col as usize]
                    != self.align_pair.pattern_index[right_col as usize]
                {
                    self.write_to_alignment_file("\nError in statistics.");
                    self.write_to_alignment_file("\nLeft pointer index not the same as right pointer index");
                    self.write_to_alignment_file("\non match or mismatch when both secnd are dash.");
                    let text = format!("\nlp: {},  rp: {}", left_col, right_col);
                    self.write_to_alignment_file(&text);
                    self.fatal(-10);
                }
                if left_seq_char == right_seq_char {
                    match_count += 1;
                    let distance = (self.align_pair.sequence_index[right_col as usize]
                        - self.align_pair.sequence_index[left_col as usize])
                        .abs();
                    if distance < 4 * self.max_distance {
                        self.statistics_distance[distance as usize] += 1;
                        if distance < min_distance {
                            min_distance = distance;
                        }
                        if distance > max_distance {
                            max_distance = distance;
                        }
                    }
                } else {
                    mismatch += 1;
                }
                left_col += 1;
                right_col += 1;
            } else if left_pat_char == DASH {
                indel += 1;
                left_col += 1;
            } else if right_pat_char == DASH {
                indel += 1;
                right_col += 1;
            }
        }
        if left_col >= right_col {
            self.write_to_alignment_file("\nError in statistics.");
            self.write_to_alignment_file("\nLeft pointer >= right pointer");
            self.fatal(-10);
        }

        let total_compared = match_count + mismatch + indel;

        let mut best_match_distance = 0;
        let mut best_match_count = 0;
        for distance in min_distance..=max_distance {
            if self.statistics_distance[distance as usize] != 0
                && self.statistics_distance[distance as usize] > best_match_count
            {
                best_match_count = self.statistics_distance[distance as usize];
                best_match_distance = distance;
            }
        }

        self.reverse();
        self.period = best_match_distance;

        if !self.params.html_output_off {
            let class_length = self.class_length;
            self.print_alignment_headings(class_length);
        }

        /* save the starting position of consensus in EC */
        let consensus_start_pos = self.align_pair.pattern_index[1];

        if !self.params.html_output_off {
            let class_length = self.class_length;
            self.print_alignment(class_length);
        }

        self.reverse();

        if !self.params.html_output_off {
            let mut out = String::new();
            out.push_str("\nStatistics");
            out.push_str(&format!(
                "\nMatches: {},  Mismatches: {}, Indels: {}",
                match_count, mismatch, indel
            ));
            out.push_str(&format!(
                "\n        {:0.2}            {:0.2}        {:0.2}",
                match_count as f32 / total_compared as f32,
                mismatch as f32 / total_compared as f32,
                indel as f32 / total_compared as f32
            ));
            out.push('\n');
            out.push_str("\nMatches are distributed among these distances:");
            for distance in min_distance..=max_distance {
                if self.statistics_distance[distance as usize] != 0 {
                    out.push_str(&format!(
                        "\n {:3}  {:3}  {:0.2}",
                        distance,
                        self.statistics_distance[distance as usize],
                        self.statistics_distance[distance as usize] as f32 / match_count as f32
                    ));
                }
            }
            self.write_to_alignment_file(&out);
        }

        let mut acgt_count = [0i32; 26];
        let mut count = 0i32;
        for align_col in 1..=align_length {
            let letter = self.align_pair.sequence_text[align_col as usize];
            if letter != DASH && (b'A'..=b'Z').contains(&letter) {
                acgt_count[(letter - b'A') as usize] += 1;
                count += 1;
            }
        }

        if count != 0 {
            let text = format!(
                "\n\nACGTcount: A:{:3.2}, C:{:3.2}, G:{:3.2}, T:{:3.2}",
                acgt_count[base_index(b'A')] as f64 / count as f64,
                acgt_count[base_index(b'C')] as f64 / count as f64,
                acgt_count[base_index(b'G')] as f64 / count as f64,
                acgt_count[base_index(b'T')] as f64 / count as f64
            );
            self.write_to_alignment_file(&text);
        } else {
            self.write_to_alignment_file("\nError in statistics: ACGTcount=0");
        }

        let mut diversity = [
            acgt_count[base_index(b'A')] as f64 / count as f64,
            acgt_count[base_index(b'C')] as f64 / count as f64,
            acgt_count[base_index(b'G')] as f64 / count as f64,
            acgt_count[base_index(b'T')] as f64 / count as f64,
        ];
        for upper in (0..=2).rev() {
            for idx in 0..=upper {
                if diversity[idx] < diversity[idx + 1] {
                    diversity.swap(idx, idx + 1);
                }
            }
        }

        let mut entropy = 0.0f64;
        for frequency in diversity.iter() {
            if *frequency != 0.0 {
                entropy += frequency * (frequency.ln() / std::f64::consts::LN_2);
            }
        }
        if entropy < 0.0 {
            entropy = -entropy;
        }

        if !self.params.html_output_off {
            self.write_to_alignment_file("\n");
            self.write_to_alignment_file("\n");
        }

        if !self.params.html_output_off {
            self.print_consensus_to_alignment_file(consensus_start_pos, consensus_size);

            if self.print_flanking {
                self.reverse();
                let flank_length = self.params.flankinglength as i32;
                self.print_flanking_sequence(flank_length);
                self.reverse();
            }
        }

        /* add the record */
        self.counter_in_sequence += 1;
        let reference = format!(
            "{}--{},{},{:3.1},{},{}",
            self.align_pair.sequence_index[self.align_pair.length as usize],
            self.align_pair.sequence_index[1],
            best_match_distance,
            self.copynumber,
            self.class_length,
            self.output_count as i32
        );
        let pattern = self.consensus_pattern_to_string(consensus_start_pos, consensus_size);
        let record = RepeatRecord {
            count: self.counter_in_sequence,
            reference,
            first: self.align_pair.sequence_index[self.align_pair.length as usize],
            last: self.align_pair.sequence_index[1],
            period: best_match_distance,
            copies: self.copynumber as f32,
            size: self.class_length,
            matches: (100.0 * match_count as f32 / total_compared as f32) as i32,
            indels: (100.0 * indel as f32 / total_compared as f32) as i32,
            score: self.align_pair.score,
            a_count: (100.0 * acgt_count[0] as f64 / count as f64) as i32,
            c_count: (100.0 * acgt_count[2] as f64 / count as f64) as i32,
            g_count: (100.0 * acgt_count[6] as f64 / count as f64) as i32,
            t_count: (100.0 * acgt_count[19] as f64 / count as f64) as i32,
            entropy: entropy as f32,
            pattern,
        };
        self.records.push(record);
    }
}
