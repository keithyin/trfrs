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
pub const EO_MAX_TBL: usize = 120;

/// One reported tandem repeat (`struct index_list` / `IL`).
#[derive(Clone, Debug)]
pub struct IlRecord {
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
    pub acount: i32,
    pub ccount: i32,
    pub gcount: i32,
    pub tcount: i32,
    pub entropy: f32,
    pub pattern: String,
}

fn write_file(path: &str, content: &str) {
    if let Err(e) = fs::write(path, content) {
        eprintln!("Error: unable to write {}: {}", path, e);
        std::process::exit(-1);
    }
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
fn interval_overlap(i: &IlRecord, j: &IlRecord) -> i32 {
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
fn is_redundant(i: &IlRecord, j: &IlRecord) -> bool {
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

    /// `RemoveBySize()`.
    fn remove_by_size(&mut self, maxsize: i32) {
        self.records.retain(|r| r.period <= maxsize);
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
        let mut i = 0usize;
        while i < self.records.len() {
            let mut j = i + 1;
            loop {
                if j >= self.records.len() {
                    break;
                }
                let iinterval = self.records[i].last - self.records[i].first + 1;
                let jinterval = self.records[j].last - self.records[j].first + 1;
                let overlap = interval_overlap(&self.records[i], &self.records[j]);
                if overlap == 0 {
                    break;
                }

                if !(overlap as f64 / (iinterval as f64) < 0.9)
                    && is_redundant(&self.records[i], &self.records[j])
                {
                    self.records.remove(i);
                    if i >= self.records.len() {
                        break;
                    }
                    j = i + 1;
                    continue;
                }

                if !(overlap as f64 / (jinterval as f64) < 0.9)
                    && is_redundant(&self.records[j], &self.records[i])
                {
                    self.records.remove(j);
                    continue;
                }

                j += 1;
            }
            i += 1;
        }
    }

    /* ------------------------------------------------------------------ */
    /* alignment file post processing                                      */
    /* ------------------------------------------------------------------ */

    /// `CleanAlignments()`: drop the alignment blocks of removed repeats.
    fn clean_alignments(&mut self, alignmentfile: &str) {
        let content = match fs::read_to_string(alignmentfile) {
            Ok(c) => c,
            Err(_) => {
                self.print_error(
                    "Unable to open alignment file for reading in CleanAlignments routine!",
                );
                std::process::exit(-1);
            }
        };
        let lines: Vec<&str> = content.split_inclusive('\n').collect();
        let mut out = String::new();
        let mut moving = true;
        let mut curr = 0usize;
        let mut i = 0usize;
        while i < lines.len() {
            let line = lines[i];
            if line.starts_with("Fo") {
                let l2 = if i + 1 < lines.len() {
                    lines[i + 1]
                } else {
                    ""
                };
                let l3 = if i + 2 < lines.len() {
                    lines[i + 2]
                } else {
                    ""
                };
                if curr < self.records.len() && l3.contains(&self.records[curr].reference) {
                    out.push_str(line);
                    out.push_str(l2);
                    out.push_str(l3);
                    curr += 1;
                    moving = true;
                    i += 3;
                    continue;
                }
                moving = false;
                i += 3;
                continue;
            }
            if line.starts_with("Do") {
                let l2 = if i + 1 < lines.len() {
                    lines[i + 1]
                } else {
                    ""
                };
                out.push_str(line);
                out.push_str(l2);
                out.push('\n');
                break;
            }
            if moving {
                out.push_str(line);
            }
            i += 1;
        }

        write_file(alignmentfile, &out);
    }

    /// `BreakAlignments()`: split the alignment file into pages.
    fn break_alignments(&mut self, alignmentfile: &str) {
        let alignments = self.records.len();
        let mut nfiles = alignments / EO_MAX_TBL;
        if alignments % EO_MAX_TBL > 0 {
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

        let content = match fs::read_to_string(alignmentfile) {
            Ok(c) => c,
            Err(_) => {
                self.print_error(
                    "Unable to open alignment file for reading in BreakAlignments routine!",
                );
                std::process::exit(-1);
            }
        };
        let lines: Vec<&str> = content.split_inclusive('\n').collect();

        /* get heading lines */
        let mut headcnt = 0usize;
        for i in 0..30 {
            if i >= lines.len() {
                break;
            }
            if lines[i].starts_with('F') {
                break;
            }
            headcnt += 1;
        }

        let mut idx = headcnt;
        for file_no in 1..=nfiles {
            let outfile = make_file_name(alignmentfile, file_no as i32);
            let mut out = String::new();
            for line in lines.iter().take(headcnt) {
                out.push_str(line);
            }
            out.push_str(&format!("File {} of {}\n\n", file_no, nfiles));

            for _ in 0..EO_MAX_TBL {
                if idx >= lines.len() {
                    break;
                }
                out.push_str(lines[idx]);
                idx += 1;
                let mut saw_done = false;
                while idx < lines.len() {
                    let l = lines[idx];
                    if l.starts_with('F') || l.starts_with('D') {
                        saw_done = l.starts_with('D');
                        break;
                    }
                    out.push_str(l);
                    idx += 1;
                }
                if idx >= lines.len() || saw_done {
                    break;
                }
            }

            out.push_str("\nDone.\n</PRE></BODY></HTML>\n");
            write_file(&outfile, &out);
        }

        let _ = fs::remove_file(alignmentfile);
    }

    /// `OutputHeading()`.
    fn output_heading(&self, tablefile: &str, alignmentfile: &str) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "<HTML><HEAD><TITLE>{}</TITLE><BASE TARGET=\"{}\"></HEAD><BODY bgcolor=\"#FBF8BC\"><BR><PRE>Tandem Repeats Finder Program written by:</PRE><PRE><CENTER>Gary Benson<BR>Program in Bioinformatics<BR>Boston University<BR>Version {}<BR></CENTER>",
            tablefile, alignmentfile, VERSION
        ));
        out.push_str("\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");
        out.push_str(&format!("\n{}", self.hsequence));
        out.push_str(&self.hparameters);
        out.push_str(&format!("{}</PRE>\n", self.hlength));
        out
    }

    /// `OutputHTML()`.
    fn output_html(&self, tablefile: &str, alignmentfile: &str) {
        let alignments = self.records.len();
        let mut nfiles = alignments / EO_MAX_TBL;
        if alignments % EO_MAX_TBL > 0 {
            nfiles += 1;
        }
        if nfiles == 0 {
            nfiles = 1;
        }

        const HEADER_ROW: &str = "<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n";

        let mut curr = 0usize;
        for i in 1..=nfiles {
            let outfile = make_file_name(tablefile, i as i32);
            let linkfile = make_file_name(alignmentfile, i as i32);
            let mut out = String::new();
            out.push_str(&self.output_heading(&outfile, alignmentfile));

            out.push_str("\n<P><PRE>Tables:   ");
            for j in 1..=nfiles {
                if j != i {
                    let namebuffer = make_file_name(tablefile, j as i32);
                    out.push_str(&format!(
                        "<A HREF=\"{}\" target=\"_self\">{}</A>   ",
                        namebuffer, j
                    ));
                } else {
                    out.push_str(&format!("{}   ", j));
                }
                if j % 16 == 0 && j < nfiles {
                    out.push_str("\n          ");
                }
            }

            out.push_str(&format!(
                "\n\nThis is table  {}  of  {}  ( {} repeats found )\n</PRE>",
                i, nfiles, alignments
            ));
            out.push_str("<PRE>\nClick on indices to view alignment\n");
            out.push_str("</PRE><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");
            out.push_str("<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");

            let mut j = 0usize;
            while j < EO_MAX_TBL && curr < self.records.len() {
                if j % 22 == 0 {
                    out.push_str(HEADER_ROW);
                }
                let r = &self.records[curr];
                out.push_str(&format!(
                    "<TR><TD><CENTER><A HREF=\"{}#{}\">{}--{}</A></CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{:1.1}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{}</CENTER></TD><TD><CENTER>{:1.2}</CENTER></TD></TR>\n",
                    linkfile,
                    r.reference,
                    r.first,
                    r.last,
                    r.period,
                    r.copies,
                    r.size,
                    r.matches,
                    r.indels,
                    r.score,
                    r.acount,
                    r.ccount,
                    r.gcount,
                    r.tcount,
                    r.entropy
                ));
                curr += 1;
                j += 1;
            }

            if self.records.is_empty() {
                out.push_str(HEADER_ROW);
            }

            out.push_str("\n</TABLE>\n");

            if self.records.is_empty() {
                out.push_str("\nNo Repeats Found!<BR>");
            }

            out.push_str("\n<P><PRE>Tables:   ");
            for j in 1..=nfiles {
                if j != i {
                    let namebuffer = make_file_name(tablefile, j as i32);
                    out.push_str(&format!(
                        "<A HREF=\"{}\" target=\"_self\">{}</A>   ",
                        namebuffer, j
                    ));
                } else {
                    out.push_str(&format!("{}   ", j));
                }
                if j % 16 == 0 && j < nfiles {
                    out.push_str("\n          ");
                }
            }
            out.push_str("\n</PRE>");

            if i == nfiles {
                out.push_str("<P>The End!\n");
            }

            out.push_str("\n</BODY></HTML>\n");
            write_file(&outfile, &out);
        }
    }

    /// `MakeMaskedFile()`: FASTA with all repeat positions replaced by `N`.
    pub fn make_masked_file(&mut self, masked: bool, maskfile: &str) {
        if !masked {
            return;
        }

        for r in &self.records {
            for m in r.first..=r.last {
                self.sequence[m as usize] = b'N';
            }
        }

        let header = if self.hsequence.len() > 10 {
            &self.hsequence[10..]
        } else {
            ""
        };
        let mut out = String::new();
        out.push('>');
        out.push_str(header);

        let mut printcr = 0;
        let mut count = 1usize;
        while self.sequence[count] != 0 {
            out.push(self.sequence[count] as char);
            printcr += 1;
            if printcr >= 60 {
                printcr = 0;
                out.push('\n');
            }
            count += 1;
        }
        out.push('\n');
        out.push('\n');

        write_file(maskfile, &out);
    }

    /// `TRFClean()`: removes redundancy and writes the HTML/data output.
    pub fn trf_clean(&mut self, tablefile: &str, alignmentfile: &str, maxsize: i32) {
        self.remove_by_size(maxsize);
        self.sort_by_index();
        if !self.params.redundoff {
            self.remove_redundancy();
        }
        self.sort_by_count();

        if !self.params.htmloff {
            self.clean_alignments(alignmentfile);
            self.break_alignments(alignmentfile);
            self.output_html(tablefile, alignmentfile);
        }

        self.params.outputcount = self.records.len() as i32;
    }
}
