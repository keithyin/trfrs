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

//! Consensus computation and alignment printing.
//!
//! This is a direct translation of `get_consensus()` and the alignment
//! printing routines from the reference implementation.

use crate::state::{Trf, DASH};

impl Trf {
    /// `reverse()`: reverses the alignment in `AlignPair`.
    pub fn reverse(&mut self) {
        let length = self.align_pair.length;
        let mut front = 1;
        while front <= length / 2 {
            let back = length - front + 1;
            self.align_pair
                .sequence_text
                .swap(front as usize, back as usize);
            self.align_pair
                .pattern_text
                .swap(front as usize, back as usize);
            self.align_pair
                .sequence_index
                .swap(front as usize, back as usize);
            self.align_pair
                .pattern_index
                .swap(front as usize, back as usize);
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
                if self.align_pair.sequence_text[col as usize]
                    != self.align_pair.pattern_text[col as usize]
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
            out.push_str(&format!(
                "  {:9} ",
                self.align_pair.sequence_index[col as usize]
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
            let display_end = std::cmp::min(
                self.align_pair.sequence_index[align_length as usize] + 10,
                self.length,
            );
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
        out.push_str(&format!("  Score: {}", self.dp.max_score));
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
            out.push(self.dp.current_pattern[counter as usize] as char);
            chars += 1;
        }
        for counter in 0..start {
            if chars % 65 == 0 {
                out.push('\n');
            }
            out.push(self.dp.current_pattern[counter as usize] as char);
            chars += 1;
        }
        out.push('\n');
        self.write_to_alignment_file(&out);
    }

    /// `printECtoBuffer()`.
    pub fn consensus_pattern_to_string(&self, start: i32, width: i32) -> String {
        let mut out = String::new();
        for counter in start..width {
            out.push(self.dp.current_pattern[counter as usize] as char);
        }
        for counter in 0..start {
            out.push(self.dp.current_pattern[counter as usize] as char);
        }
        out
    }

    /// `printECtoData()`.
    #[allow(dead_code)]
    pub fn print_consensus_to_data(&mut self, start: i32, width: i32) {
        let mut out = String::new();
        for counter in start..width {
            out.push(self.dp.current_pattern[counter as usize] as char);
        }
        for counter in 0..start {
            out.push(self.dp.current_pattern[counter as usize] as char);
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
                    let pattern_position =
                        self.align_pair.pattern_index[align_col as usize] as usize;
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
                self.consensus.pattern[write_index as usize] =
                    self.consensus.pattern[pattern_pos as usize];
                write_index += 1;
            }
        }
        self.consensus_class_length = write_index;
    }
}
