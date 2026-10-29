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

//! Redundancy removal and output generation (`trfclean.h`).

use std::fs;

use crate::params::VERSION;
use crate::state::Trf;

/// Maximum number of repeats listed in one table / alignment file.
pub const MAX_REPEATS_PER_TABLE: usize = 120;

/// One reported tandem repeat (`struct index_list` / `IL`).
#[derive(Clone, Debug)]
pub struct RepeatRecord {
    /// Order in which the repeat was found.
    pub count: i32,
    /// Anchor used to link the summary table with the alignment file.
    pub reference: String,
    pub first: i32,
    pub last: i32,
    pub period: i32,
    pub copies: f32,
    pub size: i32,
    pub matches: i32,
    pub indels: i32,
    pub score: i32,
    pub a_count: i32,
    pub c_count: i32,
    pub g_count: i32,
    pub t_count: i32,
    pub entropy: f32,
    pub pattern: String,
}

fn write_file(path: &str, content: &[u8]) {
    if let Err(e) = fs::write(path, content) {
        eprintln!("Error: unable to write {}: {}", path, e);
        std::process::exit(-1);
    }
}

/// `strstr()`: whether `haystack` contains `needle`.
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || haystack.windows(needle.len()).any(|w| w == needle)
}

/// `MakeFileName()`: replace the number in front of the extension.
pub fn make_file_name(oldname: &str, tag: i32) -> String {
    let bytes = oldname.as_bytes();
    let mut numindex = 0usize;
    for (i, b) in bytes.iter().enumerate() {
        if b.is_ascii_digit() {
            numindex = i;
        }
    }
    let oldext = if numindex + 1 <= bytes.len() {
        &oldname[numindex + 1..]
    } else {
        ""
    };
    let mut newname = String::new();
    newname.push_str(&oldname[..numindex + 1]);
    newname.push_str(&format!(".{}{}", tag, oldext));
    newname
}

/// `IntervalOverlap()`.
fn interval_overlap(i: &RepeatRecord, j: &RepeatRecord) -> i32 {
    let beg = if i.first > j.first { i.first } else { j.first };
    let end = if i.last < j.last { i.last } else { j.last };
    let overlap = end - beg + 1;
    if overlap > 0 {
        overlap
    } else {
        0
    }
}

/// `IsRedundant()`: returns true if `i` is redundant with respect to `j`.
fn is_redundant(i: &RepeatRecord, j: &RepeatRecord) -> bool {
    if i.period > j.period && i.period % j.period == 0 && (i.score as f64) <= 1.1 * j.score as f64 {
        return true;
    }
    if i.period == j.period && i.score <= j.score {
        return true;
    }
    false
}

impl Trf {
    /* ------------------------------------------------------------------ */
    /* list handling                                                       */
    /* ------------------------------------------------------------------ */

    /// `RemoveBySize()`, extended with a lower period bound (`--minperiod`).
    ///
    /// Runs before [`Trf::remove_redundancy`]: a repeat that is a multiple of a
    /// reported period (for example a 30 bp repeat in a 15 bp STR search) has to
    /// be dropped first, otherwise the shorter repeat would be discarded as
    /// redundant with a period that is not reported at all.
    fn remove_by_period(&mut self, min_period: i32, max_period: i32) {
        self.records
            .retain(|r| r.period >= min_period && r.period <= max_period);
    }

    /// `SortByIndex()`.
    fn sort_by_index(&mut self) {
        self.records.sort_by_key(|r| r.first);
    }

    /// `SortByCount()`.
    fn sort_by_count(&mut self) {
        self.records.sort_by_key(|r| r.count);
    }

    /// `RemoveRedundancy()`.
    fn remove_redundancy(&mut self) {
        if self.records.len() < 2 {
            return;
        }
        let mut record_index = 0usize;
        while record_index < self.records.len() {
            let mut other_index = record_index + 1;
            loop {
                if other_index >= self.records.len() {
                    break;
                }
                let iinterval =
                    self.records[record_index].last - self.records[record_index].first + 1;
                let jinterval =
                    self.records[other_index].last - self.records[other_index].first + 1;
                let overlap =
                    interval_overlap(&self.records[record_index], &self.records[other_index]);
                if overlap == 0 {
                    break;
                }

                if !(overlap as f64 / (iinterval as f64) < 0.9)
                    && is_redundant(&self.records[record_index], &self.records[other_index])
                {
                    self.records.remove(record_index);
                    if record_index >= self.records.len() {
                        break;
                    }
                    other_index = record_index + 1;
                    continue;
                }

                if !(overlap as f64 / (jinterval as f64) < 0.9)
                    && is_redundant(&self.records[other_index], &self.records[record_index])
                {
                    self.records.remove(other_index);
                    continue;
                }

                other_index += 1;
            }
            record_index += 1;
        }
    }

    /* ------------------------------------------------------------------ */
    /* alignment file post processing                                      */
    /* ------------------------------------------------------------------ */

    /// `CleanAlignments()`: drop the alignment blocks of removed repeats.
    fn clean_alignments(&mut self, alignmentfile: &str) {
        let content = match fs::read(alignmentfile) {
            Ok(c) => c,
            Err(_) => {
                self.print_error(
                    "Unable to open alignment file for reading in CleanAlignments routine!",
                );
                std::process::exit(-1);
            }
        };
        let lines: Vec<&[u8]> = content.split_inclusive(|&b| b == b'\n').collect();
        let mut out: Vec<u8> = Vec::new();
        let mut moving = true;
        let mut curr = 0usize;
        let mut line_index = 0usize;
        while line_index < lines.len() {
            let line = lines[line_index];
            if line.starts_with(b"Fo") {
                let second_line: &[u8] = if line_index + 1 < lines.len() {
                    lines[line_index + 1]
                } else {
                    b""
                };
                let third_line: &[u8] = if line_index + 2 < lines.len() {
                    lines[line_index + 2]
                } else {
                    b""
                };
                if curr < self.records.len()
                    && contains(third_line, self.records[curr].reference.as_bytes())
                {
                    out.extend_from_slice(line);
                    out.extend_from_slice(second_line);
                    out.extend_from_slice(third_line);
                    curr += 1;
                    moving = true;
                    line_index += 3;
                    continue;
                }
                moving = false;
                line_index += 3;
                continue;
            }
            if line.starts_with(b"Do") {
                let second_line: &[u8] = if line_index + 1 < lines.len() {
                    lines[line_index + 1]
                } else {
                    b""
                };
                out.extend_from_slice(line);
                out.extend_from_slice(second_line);
                out.push(b'\n');
                break;
            }
            if moving {
                out.extend_from_slice(line);
            }
            line_index += 1;
        }

        write_file(alignmentfile, &out);
    }

    /// `BreakAlignments()`: split the alignment file into pages.
    fn break_alignments(&mut self, alignmentfile: &str) {
        let alignments = self.records.len();
        let mut nfiles = alignments / MAX_REPEATS_PER_TABLE;
        if alignments % MAX_REPEATS_PER_TABLE > 0 {
            nfiles += 1;
        }
        if nfiles == 0 {
            nfiles = 1;
        }

        if nfiles == 1 {
            let outfile = make_file_name(alignmentfile, 1);
            let _ = fs::remove_file(&outfile);
            if let Err(e) = fs::rename(alignmentfile, &outfile) {
                eprintln!("Error: unable to rename {}: {}", alignmentfile, e);
                std::process::exit(-1);
            }
            return;
        }

        let content = match fs::read(alignmentfile) {
            Ok(c) => c,
            Err(_) => {
                self.print_error(
                    "Unable to open alignment file for reading in BreakAlignments routine!",
                );
                std::process::exit(-1);
            }
        };
        let lines: Vec<&[u8]> = content.split_inclusive(|&b| b == b'\n').collect();

        /* get heading lines */
        let mut headcnt = 0usize;
        for line_index in 0..30 {
            if line_index >= lines.len() {
                break;
            }
            if lines[line_index].starts_with(b"F") {
                break;
            }
            headcnt += 1;
        }

        let mut idx = headcnt;
        for file_no in 1..=nfiles {
            let outfile = make_file_name(alignmentfile, file_no as i32);
            let mut out: Vec<u8> = Vec::new();
            for line in lines.iter().take(headcnt) {
                out.extend_from_slice(line);
            }
            out.extend_from_slice(format!("File {} of {}\n\n", file_no, nfiles).as_bytes());

            for _ in 0..MAX_REPEATS_PER_TABLE {
                if idx >= lines.len() {
                    break;
                }
                out.extend_from_slice(lines[idx]);
                idx += 1;
                let mut saw_done = false;
                while idx < lines.len() {
                    let line = lines[idx];
                    if line.starts_with(b"F") || line.starts_with(b"D") {
                        saw_done = line.starts_with(b"D");
                        break;
                    }
                    out.extend_from_slice(line);
                    idx += 1;
                }
                if idx >= lines.len() || saw_done {
                    break;
                }
            }

            out.extend_from_slice(b"\nDone.\n</PRE></BODY></HTML>\n");
            write_file(&outfile, &out);
        }

        let _ = fs::remove_file(alignmentfile);
    }

    /// `OutputHeading()`.
    fn output_heading(&self, tablefile: &str, alignmentfile: &str) -> Vec<u8> {
        let mut out = String::new();
        out.push_str(&format!(
            "<HTML><HEAD><TITLE>{}</TITLE><BASE TARGET=\"{}\"></HEAD><BODY bgcolor=\"#FBF8BC\"><BR><PRE>Tandem Repeats Finder Program written by:</PRE><PRE><CENTER>Gary Benson<BR>Program in Bioinformatics<BR>Boston University<BR>Version {}<BR></CENTER>",
            tablefile, alignmentfile, VERSION
        ));
        out.push_str("\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");
        let mut bytes = out.into_bytes();
        bytes.push(b'\n');
        bytes.extend_from_slice(&self.header_sequence);
        bytes.extend_from_slice(self.header_parameters.as_bytes());
        bytes.extend_from_slice(self.header_length.as_bytes());
        bytes.extend_from_slice(b"</PRE>\n");
        bytes
    }

    /// `OutputHTML()`.
    fn output_html(&self, tablefile: &str, alignmentfile: &str) {
        let alignments = self.records.len();
        let mut nfiles = alignments / MAX_REPEATS_PER_TABLE;
        if alignments % MAX_REPEATS_PER_TABLE > 0 {
            nfiles += 1;
        }
        if nfiles == 0 {
            nfiles = 1;
        }

        const HEADER_ROW: &str = "<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n";

        let mut curr = 0usize;
        for page_number in 1..=nfiles {
            let outfile = make_file_name(tablefile, page_number as i32);
            let linkfile = make_file_name(alignmentfile, page_number as i32);
            /* the heading contains the raw sequence name, so it is kept as
            bytes and prepended to the (ASCII) body */
            let heading = self.output_heading(&outfile, alignmentfile);
            let mut out = String::new();

            out.push_str("\n<P><PRE>Tables:   ");
            for link_page in 1..=nfiles {
                if link_page != page_number {
                    let namebuffer = make_file_name(tablefile, link_page as i32);
                    out.push_str(&format!(
                        "<A HREF=\"{}\" target=\"_self\">{}</A>   ",
                        namebuffer, link_page
                    ));
                } else {
                    out.push_str(&format!("{}   ", link_page));
                }
                if link_page % 16 == 0 && link_page < nfiles {
                    out.push_str("\n          ");
                }
            }

            out.push_str(&format!(
                "\n\nThis is table  {}  of  {}  ( {} repeats found )\n</PRE>",
                page_number, nfiles, alignments
            ));
            out.push_str("<PRE>\nClick on indices to view alignment\n");
            out.push_str("</PRE><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");
            out.push_str("<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");

            let mut record_in_page = 0usize;
            while record_in_page < MAX_REPEATS_PER_TABLE && curr < self.records.len() {
                if record_in_page % 22 == 0 {
                    out.push_str(HEADER_ROW);
                }
                let record = &self.records[curr];
                out.push_str(&format!(
                    "<TR><TD><CENTER><A HREF=\"{}#{}\">{}--{}</A></CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{:1.1}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{:1.2}</CENTER></TD></TR>\n",
                    linkfile,
                    record.reference,
                    record.first,
                    record.last,
                    record.period,
                    record.copies,
                    record.size,
                    record.matches,
                    record.indels,
                    record.score,
                    record.a_count,
                    record.c_count,
                    record.g_count,
                    record.t_count,
                    record.entropy
                ));
                curr += 1;
                record_in_page += 1;
            }

            if self.records.is_empty() {
                out.push_str(HEADER_ROW);
            }

            out.push_str("\n</TABLE>\n");

            if self.records.is_empty() {
                out.push_str("\nNo Repeats Found!<BR>");
            }

            out.push_str("\n<P><PRE>Tables:   ");
            for link_page in 1..=nfiles {
                if link_page != page_number {
                    let namebuffer = make_file_name(tablefile, link_page as i32);
                    out.push_str(&format!(
                        "<A HREF=\"{}\" target=\"_self\">{}</A>   ",
                        namebuffer, link_page
                    ));
                } else {
                    out.push_str(&format!("{}   ", link_page));
                }
                if link_page % 16 == 0 && link_page < nfiles {
                    out.push_str("\n          ");
                }
            }
            out.push_str("\n</PRE>");

            if page_number == nfiles {
                out.push_str("<P>The End!\n");
            }

            out.push_str("\n</BODY></HTML>\n");

            let mut page = heading;
            page.extend_from_slice(out.as_bytes());
            write_file(&outfile, &page);
        }
    }

    /// `MakeMaskedFile()`: FASTA with all repeat positions replaced by `N`.
    pub fn make_masked_file(&mut self, masked: bool, maskfile: &str) {
        if !masked {
            return;
        }

        for record in &self.records {
            for position in record.first..=record.last {
                self.sequence[position as usize] = b'N';
            }
        }

        /* `&header_sequence[10]` skips the "Sequence: " prefix */
        let header: &[u8] = if self.header_sequence.len() > 10 {
            &self.header_sequence[10..]
        } else {
            b""
        };
        let mut out: Vec<u8> = Vec::with_capacity(self.sequence.len() + 16);
        out.push(b'>');
        out.extend_from_slice(header);

        let mut printcr = 0;
        let mut count = 1usize;
        while self.sequence[count] != 0 {
            out.push(self.sequence[count]);
            printcr += 1;
            if printcr >= 60 {
                printcr = 0;
                out.push(b'\n');
            }
            count += 1;
        }
        out.push(b'\n');
        out.push(b'\n');

        write_file(maskfile, &out);
    }

    /// `TRFClean()`: removes redundancy and writes the HTML/data output.
    ///
    /// The reported period range (`min_period` ..= `maxperiod`) is taken from
    /// the parameters, so the command line and the library apply exactly the
    /// same filter.
    pub fn trf_clean(&mut self, tablefile: &str, alignmentfile: &str) {
        let min_period = self.params.min_period as i32;
        let max_period = self.maxperiod;
        self.remove_by_period(min_period, max_period);
        self.sort_by_index();
        if !self.params.redundancy_removal_off {
            self.remove_redundancy();
        }
        self.sort_by_count();

        if !self.params.html_output_off {
            self.clean_alignments(alignmentfile);
            self.break_alignments(alignmentfile);
            self.output_html(tablefile, alignmentfile);
        }

        self.params.output_count = self.records.len() as i32;
    }
}
