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
use std::path::Path;

use crate::clean::make_file_name;
use crate::fasta::{load_sequence, ByteReader, FastaSequence, LoadError};
use crate::params::{ParamSet, VERSION};
use crate::state::{format_c_float, Trf, SMALL_DISTANCE};

/// The output files derived from one prefix: `<dir>/<prefix>.<ext>`.
///
/// `prefix` is the input file name without its directory part (single sequence
/// and the root of a multi-sequence run), or `<prefix>.s<ordinal>` for one
/// sequence of a multi-sequence run.  The search parameters are not part of
/// the names of this port; runs that should not overwrite each other write
/// into separate `--output` directories.
struct OutputPaths {
    /// Prefix written back to `ParamSet` for a sequence of a multi-sequence run.
    prefix: String,
    /// `<dir>/<prefix>.html` (page 1 is derived by [`Self::html_page`]).
    html: String,
    /// `<dir>/<prefix>.txt.html` (the alignment file).
    txt: String,
    /// `<dir>/<prefix>.dat`.
    dat: String,
    /// `<dir>/<prefix>.csv`.
    csv: String,
    /// `<dir>/<prefix>.mask`.
    mask: String,
    /// `<dir>/<prefix>.summary.html` (only read from the root paths).
    summary: String,
}

impl OutputPaths {
    fn new(dir: &Path, prefix: &str) -> OutputPaths {
        /* the default directory is `.`, which would leak a `./` prefix into
        the file names (and the HTML links) */
        let join = |rest: String| {
            if dir == Path::new(".") {
                rest
            } else {
                dir.join(rest).to_string_lossy().into_owned()
            }
        };
        OutputPaths {
            prefix: prefix.to_string(),
            html: join(format!("{prefix}.html")),
            txt: join(format!("{prefix}.txt.html")),
            dat: join(format!("{prefix}.dat")),
            csv: join(format!("{prefix}.csv")),
            mask: join(format!("{prefix}.mask")),
            summary: join(format!("{prefix}.summary.html")),
        }
    }

    /// Paths of sequence `ordinal` of a multi-sequence run.
    fn for_ordinal(dir: &Path, prefix: &str, ordinal: i32) -> OutputPaths {
        OutputPaths::new(dir, &format!("{prefix}.s{ordinal}"))
    }

    /// Table page name, using the same rule as `clean::make_file_name`.
    fn html_page(&self, page: i32) -> String {
        make_file_name(&self.html, page)
    }

    /// Alignment page name, using the same rule as `clean::make_file_name`.
    fn txt_page(&self, page: i32) -> String {
        make_file_name(&self.txt, page)
    }
}

impl Trf {
    /// `TRF()`: runs the algorithm on a single sequence and writes the output.
    pub fn run(&mut self) {
        let paths = OutputPaths::new(
            Path::new(&self.params.output_directory),
            &self.params.output_prefix,
        );

        self.prepare_sequence();
        self.open_alignment_output(&paths);
        self.initialize_tables();
        self.write_sequence_statistics();
        self.scan_sequence_with_progress();
        self.accumulate_statistics();
        self.close_alignment_output();
        self.resolve_output(&paths);

        self.progress("Done.");
    }

    /// Resets the per-run counters before the scan.
    fn prepare_sequence(&mut self) {
        self.init_best_period_list();
        self.print_flanking = self.params.flanking_sequence;
        self.progress("Allocating Memory...");

        self.dp.report_min = 0;
        self.dp.rows = 0.0;
        self.total_characters = 0.0;
    }

    /// Starts the alignment (txt) file and writes the HTML header.
    fn open_alignment_output(&mut self, paths: &OutputPaths) {
        if self.params.html_output_off {
            return;
        }
        match File::create(&paths.txt) {
            Ok(f) => self.alignment_file = Some(BufWriter::new(f)),
            Err(_) => {
                self.print_error("Unable to open alignment file for writing in TRF() routine!");
                std::process::exit(-1);
            }
        }
        let mut head: Vec<u8> = Vec::new();
        head.extend_from_slice(b"<HTML>");
        head.extend_from_slice(b"<HEAD>");
        head.extend_from_slice(format!("<TITLE>{}</TITLE>", paths.txt).as_bytes());
        head.extend_from_slice(b"</HEAD>");
        head.extend_from_slice(b"<BODY bgcolor=\"#FBF8BC\">");
        head.extend_from_slice(b"<PRE>");
        head.extend_from_slice(b"\nTandem Repeats Finder Program written by:");
        head.extend_from_slice(b"\n\n                 Gary Benson");
        head.extend_from_slice(b"\n      Program in Bioinformatics");
        head.extend_from_slice(b"\n          Boston University");
        head.extend_from_slice(format!("\n\nVersion {}", VERSION).as_bytes());
        head.extend_from_slice(b"\n\nSequence: ");
        head.extend_from_slice(&self.fasta.name);
        head.extend_from_slice(
            format!(
                "\n\nParameters: {} {} {} {} {} {} {}\n",
                self.params.match_weight,
                self.params.mismatch,
                self.params.indel,
                self.params.match_probability_percent,
                self.params.indel_probability_percent,
                self.params.minscore,
                self.params.maxperiod
            )
            .as_bytes(),
        );
        self.write_bytes_to_alignment_file(&head);
    }

    /// Initializes the criteria tables used by the scan.
    fn initialize_tables(&mut self) {
        self.progress("Initializing data structures...");

        self.progress("Computing TR Model Statistics...");
        self.init_criteria_tables();
    }

    /// Writes the sequence length and base composition to the alignment file.
    fn write_sequence_statistics(&mut self) {
        if self.params.html_output_off {
            return;
        }
        let text = format!("\n\nLength: {}", self.length);
        self.write_to_alignment_file(&text);
        let text = format!(
            "\nACGTcount: A:{}, C:{}, G:{}, T:{}\n\n",
            format_c_float(self.fasta.composition[0] as f64 / self.length as f64, 3, 2),
            format_c_float(self.fasta.composition[2] as f64 / self.length as f64, 3, 2),
            format_c_float(self.fasta.composition[6] as f64 / self.length as f64, 3, 2),
            format_c_float(self.fasta.composition[19] as f64 / self.length as f64, 3, 2)
        );
        self.write_to_alignment_file(&text);
        if self.fasta.length as i32 - self.fasta.nucleotides > 0 {
            let text = format!(
                "Warning! {} characters in sequence are not A, C, G, or T\n\n",
                self.fasta.length as i32 - self.fasta.nucleotides
            );
            self.write_to_alignment_file(&text);
        }
    }

    /// Updates the per-run character counters and scans the sequence.
    fn scan_sequence_with_progress(&mut self) {
        self.total_characters += self.length as f64;
        self.wraparound_dp_count = 0.0;

        if self.params.multi_sequence_file {
            let text = format!("Scanning Sequence {}...", self.params.sequence_ordinal);
            self.progress(&text);
        } else {
            self.progress("Scanning...");
        }

        self.scan_sequence();
    }

    /// Aggregates the cell and wasted-cell counters of the scan.
    fn accumulate_statistics(&mut self) {
        self.cell_total = 0.0;
        let mut distance = self.max_distance;
        while distance >= 1 {
            self.cell_total += self.cell_counts[distance as usize];
            if distance <= SMALL_DISTANCE {
                self.wasted_total += (self.criteria_counts[distance as usize]
                    + self.consensus_counts[distance as usize])
                    as f64
                    * (distance * distance * 2) as f64;
            } else {
                self.wasted_total += (self.criteria_counts[distance as usize]
                    + self.consensus_counts[distance as usize])
                    as f64
                    * distance as f64
                    * (2 * self.indel_random_walk_range(distance) + 1) as f64
                    * 2.0;
            }
            distance -= 1;
        }
    }

    /// Writes the closing HTML tags and drops the alignment writer.
    fn close_alignment_output(&mut self) {
        if !self.params.html_output_off {
            self.write_to_alignment_file("Done.\n");
            self.write_to_alignment_file("</PRE>");
            self.write_to_alignment_file("</BODY>");
            self.write_to_alignment_file("</HTML>");
        }
        self.alignment_file = None;
    }

    /// Removes redundant repeats and writes the table files.
    fn resolve_output(&mut self, paths: &OutputPaths) {
        self.progress("Freeing Memory...");
        self.progress("Resolving output...");

        self.trf_clean(&paths.html, &paths.txt);

        self.params.output_file_name = paths.html_page(1);
    }
}

/// Formats the `.dat` records of the current sequence.
fn format_dat_records(trf: &Trf) -> String {
    let mut out = String::new();
    for record in &trf.records {
        out.push_str(&format!(
            "{} {} {} {:.1} {} {} {} {} {} {} {} {} {:.2} {} ",
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
            record.entropy,
            record.pattern
        ));
        for position in record.first..=record.last {
            out.push(trf.sequence[position as usize] as char);
        }

        /* short flanks are appended in -ngs mode */
        if trf.params.ngs {
            let flank_start = std::cmp::max(1, record.first - 50);
            let flank_end = std::cmp::min(trf.length, record.last + 50);
            out.push(' ');
            if record.first == 1 {
                out.push('.');
            } else {
                for position in flank_start..record.first {
                    out.push(trf.sequence[position as usize] as char);
                }
            }
            out.push(' ');
            if record.last == trf.length {
                out.push('.');
            } else {
                for position in (record.last + 1)..=flank_end {
                    out.push(trf.sequence[position as usize] as char);
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

/// Header of the `.csv` output.  One column per `.dat` field of a repeat,
/// preceded by the index and the name of the sequence that holds it, so that
/// the table of a multi sequence file stays self contained.
const CSV_HEADER: &str = "sequence_index,sequence_name,start,end,period,\
copy_number,pattern_size,percent_matches,percent_indels,score,percent_a,\
percent_c,percent_g,percent_t,entropy,repeat_unit,repeat_sequence\n";

/// Appends one CSV field, quoting it when it contains a delimiter, a quote or
/// a line break (`RFC 4180`).  The sequence name is written as raw bytes, like
/// every other output of the program.
fn push_csv_field(out: &mut Vec<u8>, field: &[u8]) {
    if !field.iter().any(|b| matches!(b, b',' | b'"' | b'\n' | b'\r')) {
        out.extend_from_slice(field);
        return;
    }
    out.push(b'"');
    for &byte in field {
        if byte == b'"' {
            out.push(b'"');
        }
        out.push(byte);
    }
    out.push(b'"');
}

/// Formats the `.csv` records of the current sequence.  `sequence_index` is
/// the ordinal of the sequence in a multi sequence file (1 otherwise).
fn format_csv_records(trf: &Trf, sequence_index: i32, name: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for record in &trf.records {
        out.extend_from_slice(sequence_index.to_string().as_bytes());
        out.push(b',');
        push_csv_field(&mut out, name);
        out.push(b',');
        out.extend_from_slice(
            format!(
                "{},{},{},{:.1},{},{},{},{},{},{},{},{},{:.2},",
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
            )
            .as_bytes(),
        );
        push_csv_field(&mut out, record.pattern.as_bytes());
        out.push(b',');
        out.extend_from_slice(&trf.sequence[record.first as usize..=record.last as usize]);
        out.push(b'\n');
    }
    out
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

/// `header_sequence`: `"Sequence: %s\n"`, the header line shared by the HTML, data
/// and masked output.  Kept as bytes: the name may not be valid UTF-8.
fn sequence_header(seq: &FastaSequence) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.name.len() + 11);
    out.extend_from_slice(b"Sequence: ");
    out.extend_from_slice(&seq.name);
    out.push(b'\n');
    out
}

fn dat_sequence_header(params: &ParamSet, seq: &FastaSequence) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(b"\n\nSequence: ");
    out.extend_from_slice(&seq.name);
    out.extend_from_slice(
        format!(
            "\n\n\n\nParameters: {} {} {} {} {} {} {}\n\n\n",
            params.match_weight,
            params.mismatch,
            params.indel,
            params.match_probability_percent,
            params.indel_probability_percent,
            params.minscore,
            params.maxperiod
        )
        .as_bytes(),
    );
    out
}

/// Writes the data file. In `-ngs` mode the reference implementation writes to
/// stdout and then closes it, which is modelled by `params.stdout_closed`.
fn emit_dat(dest: Option<&str>, content: &[u8], params: &mut ParamSet) {
    match dest {
        Some(path) => {
            if let Err(e) = std::fs::write(path, content) {
                eprintln!("Error: unable to write {}: {}", path, e);
                std::process::exit(-1);
            }
        }
        None => {
            let mut stdout = std::io::stdout();
            let _ = stdout.write_all(content);
            let _ = stdout.flush();
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
    let source = params.input_file_name.clone();
    let prefix = params.output_prefix.clone();
    let output_dir = std::path::PathBuf::from(&params.output_directory);

    /* the report and data files of a non `-ngs` run live in the output
    directory; create it up front so a bad path fails before any work */
    if !params.ngs {
        if let Err(e) = std::fs::create_dir_all(&output_dir) {
            eprintln!(
                "Error: unable to create output directory {}: {}",
                output_dir.display(),
                os_error_message(&e)
            );
            params.end_status = Some(format!(
                "Cannot create output directory {}",
                output_dir.display()
            ));
            params.running = false;
            return;
        }
    }

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
                params.end_status = Some("Bad filename".into());
                params.running = false;
                return;
            }
        }
    };

    if params.verbose && !params.ngs {
        eprint!("\nLoading sequence...");
    }

    let loaded = match load_sequence(&mut reader) {
        Ok(v) => v,
        Err(e) => {
            if let LoadError::Message(m) = e {
                eprintln!("Error: {}", m);
            }
            eprintln!("Error: Could not load sequence. Empty file or bad format.");
            params.end_status = Some("Bad format.".into());
            params.running = false;
            return;
        }
    };
    let (mut seq, mut more) = loaded;

    let root_paths = OutputPaths::new(&output_dir, &prefix);
    let header_parameters = params.parameters_header();

    if !more {
        /* only one sequence in the file */
        params.multi_sequence_file = false;
        params.sequence_ordinal = 1;

        let mut trf = Trf::new(params.clone(), &seq);
        trf.header_sequence = sequence_header(&seq);
        trf.header_parameters = header_parameters.clone();
        trf.header_length = format!("Length:  {}", seq.length);
        trf.run();

        params.output_count = trf.params.output_count;
        params.output_file_name = trf.params.output_file_name.clone();

        if params.end_status.is_some() {
            return;
        }

        if params.data_file {
            let mut content: Vec<u8> = Vec::new();
            if !params.ngs {
                content.extend_from_slice(dat_program_header().as_bytes());
            }
            if params.ngs {
                if !trf.records.is_empty() {
                    content.push(b'@');
                    content.extend_from_slice(&seq.name);
                    content.push(b'\n');
                }
            } else {
                content.extend_from_slice(&dat_sequence_header(params, &seq));
            }
            content.extend_from_slice(format_dat_records(&trf).as_bytes());
            let dest = if params.ngs {
                None
            } else {
                Some(root_paths.dat.clone())
            };
            emit_dat(dest.as_deref(), &content, params);
        }

        if params.csv_file {
            let mut content = CSV_HEADER.as_bytes().to_vec();
            content.extend_from_slice(&format_csv_records(&trf, 1, &seq.name));
            if let Err(e) = std::fs::write(&root_paths.csv, &content) {
                eprintln!("Error: unable to write {}: {}", root_paths.csv, e);
                std::process::exit(-1);
            }
        }

        /* masked file is written after the .dat output so that the sequence is
        not "ruined" by Ns beforehand */
        let maskstring = root_paths.mask.clone();
        trf.make_masked_file(params.masked_file, &maskstring);

        params.end_status = None;
        params.running = false;
        return;
    }

    params.multi_sequence_file = true;
    params.sequence_ordinal = 1;

    /*********************************************************
     * more than one sequence: produce summary-style output
     **********************************************************/
    let desth = root_paths.summary.clone();

    /* The reference implementation opens the summary, masked and data files
    before the first sequence is processed and writes to them incrementally.
    Doing the same keeps those files present (and partially filled) if the
    process terminates early. */
    let mut summary_fp: Option<File> = None;
    if !params.html_output_off {
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
    if params.masked_file {
        match File::create(&root_paths.mask) {
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
    if params.data_file && !params.ngs {
        match File::create(&root_paths.dat) {
            Ok(f) => dat_fp = Some(f),
            Err(_) => {
                eprintln!(
                    "Error: Unable to open data file for writing in TRFControlRoutine routine!"
                );
                std::process::exit(-1);
            }
        }
    }

    let mut csv_fp: Option<BufWriter<File>> = None;
    if params.csv_file {
        match File::create(&root_paths.csv) {
            Ok(f) => {
                let mut f = BufWriter::new(f);
                let _ = f.write_all(CSV_HEADER.as_bytes());
                csv_fp = Some(f);
            }
            Err(_) => {
                eprintln!(
                    "Error: Unable to open csv file for writing in TRFControlRoutine routine!"
                );
                std::process::exit(-1);
            }
        }
    }

    let mut foundsome = false;
    let mut ordinal = 1i32;

    loop {
        /* prefix used for the output files of this sequence */
        let seq_paths = OutputPaths::for_ordinal(&output_dir, &prefix, ordinal);
        params.input_file_name = seq_paths.prefix.clone();
        params.output_prefix = seq_paths.prefix.clone();

        let mut trf = Trf::new(params.clone(), &seq);
        trf.header_sequence = sequence_header(&seq);
        trf.header_parameters = header_parameters.clone();
        trf.header_length = format!("Length:  {}", seq.length);
        trf.run();

        params.output_count = trf.params.output_count;
        params.output_file_name = trf.params.output_file_name.clone();

        if params.data_file {
            let mut content: Vec<u8> = Vec::new();
            if ordinal == 1 && !params.ngs {
                content.extend_from_slice(dat_program_header().as_bytes());
            }
            if params.ngs {
                if !trf.records.is_empty() {
                    content.push(b'@');
                    content.extend_from_slice(&seq.name);
                    content.push(b'\n');
                }
            } else {
                content.extend_from_slice(&dat_sequence_header(params, &seq));
            }
            content.extend_from_slice(format_dat_records(&trf).as_bytes());
            match dat_fp.as_mut() {
                Some(f) => {
                    let _ = f.write_all(&content);
                }
                None => {
                    let mut stdout = std::io::stdout();
                    let _ = stdout.write_all(&content);
                    let _ = stdout.flush();
                }
            }
        }

        if let Some(f) = csv_fp.as_mut() {
            let _ = f.write_all(&format_csv_records(&trf, ordinal, &seq.name));
        }

        if let Some(f) = summary_fp.as_mut() {
            let outh = seq_paths.html_page(1);
            if params.output_count > 0 {
                let mut row: Vec<u8> = Vec::new();
                row.extend_from_slice(
                    format!(
                        "<TR><TD><CENTER>{}</CENTER></TD><TD><CENTER><A TARGET=\"{}\" HREF=\"{}\">",
                        ordinal, outh, outh
                    )
                    .as_bytes(),
                );
                row.extend_from_slice(&seq.name);
                row.extend_from_slice(
                    format!(
                        "</A></CENTER></TD><TD><CENTER>{}</CENTER></TD></TR>",
                        params.output_count
                    )
                    .as_bytes(),
                );
                let _ = f.write_all(&row);
                foundsome = true;
            } else {
                /* remove html files that contain no output */
                let _ = std::fs::remove_file(&outh);
                let line = seq_paths.txt_page(1);
                let _ = std::fs::remove_file(&line);
            }
        }

        let maskstring = seq_paths.mask.clone();
        trf.make_masked_file(params.masked_file, &maskstring);

        if params.masked_file {
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
            if params.verbose && !params.ngs {
                eprint!("\nLoading sequence file...");
            }
            match load_sequence(&mut reader) {
                Ok((next_seq, next_more)) => {
                    seq = next_seq;
                    more = next_more;
                }
                Err(_) => {
                    eprintln!("Error: Could not load sequence. Empty file or bad format.");
                    params.end_status = Some("Bad format.".into());
                    params.running = false;
                    return;
                }
            }
            params.sequence_ordinal += 1;
            ordinal += 1;
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
    if let Some(mut f) = csv_fp {
        let _ = f.flush();
    }
    if let Some(mut f) = mask_fp {
        let _ = f.flush();
    }
    if params.data_file && params.ngs {
        /* the reference implementation closes stdout here */
        params.stdout_closed = true;
    }

    params.output_file_name = desth;
    params.end_status = None;
    params.running = false;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_paths_use_the_output_directory() {
        let paths = OutputPaths::new(Path::new("out"), "seq.fasta");
        assert_eq!(paths.prefix, "seq.fasta");
        assert_eq!(paths.html, "out/seq.fasta.html");
        assert_eq!(paths.txt, "out/seq.fasta.txt.html");
        assert_eq!(paths.dat, "out/seq.fasta.dat");
        assert_eq!(paths.csv, "out/seq.fasta.csv");
        assert_eq!(paths.mask, "out/seq.fasta.mask");
        assert_eq!(paths.summary, "out/seq.fasta.summary.html");
        assert_eq!(paths.html_page(1), "out/seq.fasta.1.html");
        assert_eq!(paths.txt_page(1), "out/seq.fasta.txt.1.html");
    }

    #[test]
    fn sequence_paths_carry_the_multi_sequence_prefix() {
        let paths = OutputPaths::for_ordinal(Path::new("."), "multi.fasta", 3);
        assert_eq!(paths.prefix, "multi.fasta.s3");
        assert_eq!(paths.html_page(1), "multi.fasta.s3.1.html");
        assert_eq!(paths.txt_page(1), "multi.fasta.s3.txt.1.html");
        assert_eq!(paths.mask, "multi.fasta.s3.mask");
    }

    #[test]
    fn csv_fields_are_quoted_only_when_needed() {
        let mut out: Vec<u8> = Vec::new();
        push_csv_field(&mut out, b"plain");
        out.push(b'|');
        push_csv_field(&mut out, b"a,b");
        out.push(b'|');
        push_csv_field(&mut out, b"say \"hi\"");
        out.push(b'|');
        push_csv_field(&mut out, b"two\nlines");
        assert_eq!(
            out,
            b"plain|\"a,b\"|\"say \"\"hi\"\"\"|\"two\nlines\"".to_vec()
        );
    }
}
