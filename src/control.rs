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

//! Control routine and per sequence driver (`trfrun.h`).

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

use crate::clean::make_file_name;
use crate::fasta::{load_sequence, ByteReader, FastaSequence, LoadError};
use crate::params::{ParamSet, VERSION};
use crate::state::{cfloat, Trf, SMALLDISTANCE};

impl Trf {
    /// `TRF()`: runs the algorithm on a single sequence and writes the output.
    pub fn run(&mut self) {
        let paramstring = self.param_string();
        let htmlstring = format!("{}.{}.html", self.params.outputprefix, paramstring);
        let txtstring = format!("{}.{}.txt.html", self.params.outputprefix, paramstring);
        let _datstring = format!("{}.{}.dat", self.params.outputprefix, paramstring);
        let _maskstring = format!("{}.{}.mask", self.params.outputprefix, paramstring);

        self.init_bestperiodlist();
        self.print_flanking = self.params.flankingsequence;
        self.progress("Allocating Memory...");

        self.reportmin = 0;
        self.rows = 0.0;
        self.totalcharacters = 0.0;

        /* start the alignment (txt) file */
        if !self.params.htmloff {
            match File::create(&txtstring) {
                Ok(f) => self.fptxt = Some(BufWriter::new(f)),
                Err(_) => {
                    self.print_error("Unable to open alignment file for writing in TRF() routine!");
                    std::process::exit(-1);
                }
            }
            let mut head = String::new();
            head.push_str("<HTML>");
            head.push_str("<HEAD>");
            head.push_str(&format!("<TITLE>{}</TITLE>", txtstring));
            head.push_str("</HEAD>");
            head.push_str("<BODY bgcolor=\"#FBF8BC\">");
            head.push_str("<PRE>");
            head.push_str("\nTandem Repeats Finder Program written by:");
            head.push_str("\n\n                 Gary Benson");
            head.push_str("\n      Program in Bioinformatics");
            head.push_str("\n          Boston University");
            head.push_str(&format!("\n\nVersion {}", VERSION));
            head.push_str(&format!(
                "\n\nSequence: {}\n\nParameters: {} {} {} {} {} {} {}\n",
                self.fasta.name,
                self.params.match_weight,
                self.params.mismatch,
                self.params.indel,
                self.params.pm,
                self.params.pi,
                self.params.minscore,
                self.params.maxperiod
            ));
            self.msg(&head);
        }

        self.progress("Initializing data structures...");

        self.progress("Computing TR Model Statistics...");
        self.init_and_fill_coin_toss_stats();

        if !self.params.htmloff {
            let text = format!("\n\nLength: {}", self.length);
            self.msg(&text);
            let text = format!(
                "\nACGTcount: A:{}, C:{}, G:{}, T:{}\n\n",
                cfloat(self.fasta.composition[0] as f64 / self.length as f64, 3, 2),
                cfloat(self.fasta.composition[2] as f64 / self.length as f64, 3, 2),
                cfloat(self.fasta.composition[6] as f64 / self.length as f64, 3, 2),
                cfloat(self.fasta.composition[19] as f64 / self.length as f64, 3, 2)
            );
            self.msg(&text);
            if self.fasta.length as i32 - self.fasta.nucleotides > 0 {
                let text = format!(
                    "Warning! {} characters in sequence are not A, C, G, or T\n\n",
                    self.fasta.length as i32 - self.fasta.nucleotides
                );
                self.msg(&text);
            }
        }

        self.totalcharacters += self.length as f64;
        self.wdp_count = 0.0;

        if self.params.multisequencefile {
            let text = format!("Scanning Sequence {}...", self.params.sequenceordinal);
            self.progress(&text);
        } else {
            self.progress("Scanning...");
        }

        self.newtupbo();

        self.cell_total = 0.0;
        let mut i = self.maxdistance;
        while i >= 1 {
            self.cell_total += self.cell_count[i as usize];
            if i <= SMALLDISTANCE {
                self.wasted_total += (self.criteria_count[i as usize]
                    + self.consensus_count[i as usize]) as f64
                    * (i * i * 2) as f64;
            } else {
                self.wasted_total += (self.criteria_count[i as usize]
                    + self.consensus_count[i as usize]) as f64
                    * i as f64
                    * (2 * self.d_range(i) + 1) as f64
                    * 2.0;
            }
            i -= 1;
        }

        if !self.params.htmloff {
            self.msg("Done.\n");
            self.msg("</PRE>");
            self.msg("</BODY>");
            self.msg("</HTML>");
        }
        self.fptxt = None;

        self.progress("Freeing Memory...");
        self.progress("Resolving output...");

        let maxperiod = self.maxperiod;
        self.trf_clean(&htmlstring, &txtstring, maxperiod);

        self.params.outputfilename = make_file_name(&htmlstring, 1);

        self.progress("Done.");
    }
}

/// Formats the `.dat` records of the current sequence.
fn format_dat_records(trf: &Trf) -> String {
    let mut out = String::new();
    for r in &trf.records {
        out.push_str(&format!(
            "{} {} {} {:.1} {} {} {} {} {} {} {} {} {:.2} {} ",
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
            r.entropy,
            r.pattern
        ));
        for c in r.first..=r.last {
            out.push(trf.sequence[c as usize] as char);
        }

        /* short flanks are appended in -ngs mode */
        if trf.params.ngs {
            let flankstart = std::cmp::max(1, r.first - 50);
            let flankend = std::cmp::min(trf.length, r.last + 50);
            out.push(' ');
            if r.first == 1 {
                out.push('.');
            } else {
                for c in flankstart..r.first {
                    out.push(trf.sequence[c as usize] as char);
                }
            }
            out.push(' ');
            if r.last == trf.length {
                out.push('.');
            } else {
                for c in (r.last + 1)..=flankend {
                    out.push(trf.sequence[c as usize] as char);
                }
            }
        }
        out.push('\n');
    }
    out
}

fn dat_program_header() -> String {
    format!(
        "Tandem Repeats Finder Program written by:\n\nGary Benson\nProgram in Bioinformatics\nBoston University\nVersion {}\n",
        VERSION
    )
}

/// Header of the multi-sequence `*.summary.html` table.
fn summary_header() -> String {
    let mut header = String::new();
    header.push_str("<HTML>");
    header.push_str("<HEAD>");
    header.push_str("<TITLE>Output Summary</TITLE>");
    header.push_str("</HEAD>");
    header.push_str("<BODY bgcolor=\"#FBF8BC\">");
    header.push_str("<PRE>");
    header.push_str("\nTandem Repeats Finder Program written by:<CENTER>");
    header.push_str("\nGary Benson");
    header.push_str("\nProgram in Bioinformatics");
    header.push_str("\nBoston University");
    header.push_str(&format!("\nVersion {}</CENTER>\n", VERSION));
    header.push_str("\nPlease cite:\nG. Benson,\n\"Tandem repeats finder: a program to analyze DNA sequences\"\nNucleic Acid Research(1999)\nVol. 27, No. 2, pp. 573-580.\n");
    header.push_str("\n\n<B>Multiple Sequence Summary</B>\n\n");
    header.push_str("Only sequences containing repeats are shown!\n\n");
    header.push_str("Click on sequence description to view repeat table.\n\n");
    header.push_str("<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");
    header.push_str(
        "<TR><TD WIDTH=80><CENTER>Sequence\nIndex</CENTER></TD><TD WIDTH=400><CENTER>Sequence\nDescription</CENTER></TD><TD WIDTH=80><CENTER>Number of\nRepeats</CENTER></TD></TR>\n",
    );
    header
}

/// Footer of the multi-sequence `*.summary.html` table.
fn summary_footer(foundsome: bool) -> String {
    let mut footer = String::new();
    footer.push_str("\n</TABLE>\n");
    if !foundsome {
        footer.push_str("\nNo Repeats Found!<BR>");
    }
    footer.push_str("\n</BODY></HTML>\n");
    footer
}

fn dat_sequence_header(params: &ParamSet, seq: &FastaSequence) -> String {
    format!(
        "\n\nSequence: {}\n\n\n\nParameters: {} {} {} {} {} {} {}\n\n\n",
        seq.name,
        params.match_weight,
        params.mismatch,
        params.indel,
        params.pm,
        params.pi,
        params.minscore,
        params.maxperiod
    )
}

/// Writes the data file. In `-ngs` mode the reference implementation writes to
/// stdout and then closes it, which is modelled by `params.stdout_closed`.
fn emit_dat(dest: Option<&str>, content: &str, params: &mut ParamSet) {
    match dest {
        Some(path) => {
            if let Err(e) = std::fs::write(path, content) {
                eprintln!("Error: unable to write {}: {}", path, e);
                std::process::exit(-1);
            }
        }
        None => {
            print!("{}", content);
            let _ = std::io::stdout().flush();
            /* fclose(stdout): messages printed afterwards are lost */
            params.stdout_closed = true;
        }
    }
}

/// Renders an I/O error the way `strerror()` does.
fn os_error_message(e: &std::io::Error) -> String {
    let text = e.to_string();
    match text.rfind(" (os error ") {
        Some(pos) => text[..pos].to_string(),
        None => text,
    }
}

/// `TRFControlRoutine()`: handles single and multiple sequence files.
pub fn trf_control_routine(params: &mut ParamSet) {
    let source = params.inputfilename.clone();
    let prefix = params.outputprefix.clone();

    let mut reader = if params.use_stdin {
        ByteReader::new(Box::new(BufReader::new(std::io::stdin())))
    } else {
        match File::open(&source) {
            Ok(f) => ByteReader::new(Box::new(BufReader::new(f))),
            Err(e) => {
                eprintln!(
                    "Error: Error opening {}: {}\n",
                    source,
                    os_error_message(&e)
                );
                params.endstatus = Some("Bad filename".into());
                params.running = false;
                return;
            }
        }
    };

    if !params.ngs {
        eprint!("\nLoading sequence...");
    }

    let loaded = match load_sequence(&mut reader) {
        Ok(v) => v,
        Err(e) => {
            if let LoadError::Message(m) = e {
                eprintln!("Error: {}", m);
            }
            eprintln!("Error: Could not load sequence. Empty file or bad format.");
            params.endstatus = Some("Bad format.".into());
            params.running = false;
            return;
        }
    };
    let (mut seq, mut more) = loaded;

    let paramstring = params.param_string();
    let hparameters = params.parameters_header();

    if !more {
        /* only one sequence in the file */
        params.multisequencefile = false;
        params.sequenceordinal = 1;

        let mut trf = Trf::new(params.clone(), &seq);
        trf.hsequence = format!("Sequence: {}\n", seq.name);
        trf.hparameters = hparameters.clone();
        trf.hlength = format!("Length:  {}", seq.length);
        trf.run();

        params.outputcount = trf.params.outputcount;
        params.outputfilename = trf.params.outputfilename.clone();

        if params.endstatus.is_some() {
            return;
        }

        if params.datafile {
            let mut content = String::new();
            if !params.ngs {
                content.push_str(&dat_program_header());
            }
            if params.ngs {
                if !trf.records.is_empty() {
                    content.push_str(&format!("@{}\n", seq.name));
                }
            } else {
                content.push_str(&dat_sequence_header(params, &seq));
            }
            content.push_str(&format_dat_records(&trf));
            let dest = if params.ngs {
                None
            } else {
                Some(format!("{}.{}.dat", prefix, paramstring))
            };
            emit_dat(dest.as_deref(), &content, params);
        }

        /* masked file is written after the .dat output so that the sequence is
        not "ruined" by Ns beforehand */
        let maskstring = format!("{}.{}.mask", params.outputprefix, paramstring);
        trf.make_masked_file(params.maskedfile, &maskstring);

        params.endstatus = None;
        params.running = false;
        return;
    }

    params.multisequencefile = true;
    params.sequenceordinal = 1;

    /*********************************************************
     * more than one sequence: produce summary-style output
     **********************************************************/
    let desth = format!("{}.{}.summary.html", prefix, paramstring);

    /* The reference implementation opens the summary, masked and data files
    before the first sequence is processed and writes to them incrementally.
    Doing the same keeps those files present (and partially filled) if the
    process terminates early. */
    let mut summary_fp: Option<File> = None;
    if !params.htmloff {
        match File::create(&desth) {
            Ok(mut f) => {
                let _ = f.write_all(summary_header().as_bytes());
                summary_fp = Some(f);
            }
            Err(_) => {
                eprintln!(
                    "Error: Unable to open summary file for writing in TRFControlRoutine routine!"
                );
                std::process::exit(-1);
            }
        }
    }

    let mut mask_fp: Option<File> = None;
    if params.maskedfile {
        match File::create(format!("{}.{}.mask", prefix, paramstring)) {
            Ok(f) => mask_fp = Some(f),
            Err(_) => {
                eprintln!(
                    "Error: Unable to open masked file for writing in TRFControlRoutine routine!"
                );
                std::process::exit(-1);
            }
        }
    }

    let mut dat_fp: Option<File> = None;
    if params.datafile && !params.ngs {
        match File::create(format!("{}.{}.dat", prefix, paramstring)) {
            Ok(f) => dat_fp = Some(f),
            Err(_) => {
                eprintln!(
                    "Error: Unable to open data file for writing in TRFControlRoutine routine!"
                );
                std::process::exit(-1);
            }
        }
    }

    let mut foundsome = false;
    let mut i = 1i32;

    loop {
        /* prefix used for the output files of this sequence */
        let input = format!("{}.s{}", prefix, i);
        params.inputfilename = input.clone();
        params.outputprefix = input.clone();

        let mut trf = Trf::new(params.clone(), &seq);
        trf.hsequence = format!("Sequence: {}\n", seq.name);
        trf.hparameters = hparameters.clone();
        trf.hlength = format!("Length:  {}", seq.length);
        trf.run();

        params.outputcount = trf.params.outputcount;
        params.outputfilename = trf.params.outputfilename.clone();

        if params.datafile {
            let mut content = String::new();
            if i == 1 && !params.ngs {
                content.push_str(&dat_program_header());
            }
            if params.ngs {
                if !trf.records.is_empty() {
                    content.push_str(&format!("@{}\n", seq.name));
                }
            } else {
                content.push_str(&dat_sequence_header(params, &seq));
            }
            content.push_str(&format_dat_records(&trf));
            match dat_fp.as_mut() {
                Some(f) => {
                    let _ = f.write_all(content.as_bytes());
                }
                None => {
                    print!("{}", content);
                    let _ = std::io::stdout().flush();
                }
            }
        }

        if let Some(f) = summary_fp.as_mut() {
            let outh = format!("{}.{}.1.html", input, paramstring);
            if params.outputcount > 0 {
                let row = format!(
                    "<TR><TD><CENTER>{}</CENTER></TD><TD><CENTER><A TARGET=\"{}\" HREF=\"{}\">{}</A></CENTER></TD><TD><CENTER>{}</CENTER></TD></TR>",
                    i, outh, outh, seq.name, params.outputcount
                );
                let _ = f.write_all(row.as_bytes());
                foundsome = true;
            } else {
                /* remove html files that contain no output */
                let _ = std::fs::remove_file(&outh);
                let line = format!("{}.{}.1.txt.html", input, paramstring);
                let _ = std::fs::remove_file(&line);
            }
        }

        let maskstring = format!("{}.s{}.{}.mask", prefix, i, paramstring);
        trf.make_masked_file(params.maskedfile, &maskstring);

        if params.maskedfile {
            match std::fs::read(&maskstring) {
                Ok(bytes) => {
                    if let Some(f) = mask_fp.as_mut() {
                        let _ = f.write_all(&bytes);
                    }
                }
                Err(_) => {
                    eprintln!("Error: Unable to open masked file for reading in TRFControlRoutine routine!");
                    std::process::exit(-1);
                }
            }
            let _ = std::fs::remove_file(&maskstring);
        }

        if more {
            if !params.ngs {
                eprint!("\nLoading sequence file...");
            }
            match load_sequence(&mut reader) {
                Ok((next_seq, next_more)) => {
                    seq = next_seq;
                    more = next_more;
                }
                Err(_) => {
                    eprintln!("Error: Could not load sequence. Empty file or bad format.");
                    params.endstatus = Some("Bad format.".into());
                    params.running = false;
                    return;
                }
            }
            params.sequenceordinal += 1;
            i += 1;
        } else {
            break;
        }
    }

    if let Some(mut f) = summary_fp {
        let _ = f.write_all(summary_footer(foundsome).as_bytes());
        let _ = f.flush();
    }
    if let Some(mut f) = dat_fp {
        let _ = f.flush();
    }
    if let Some(mut f) = mask_fp {
        let _ = f.flush();
    }
    if params.datafile && params.ngs {
        /* the reference implementation closes stdout here */
        params.stdout_closed = true;
    }

    params.outputfilename = desth;
    params.endstatus = None;
    params.running = false;
}
