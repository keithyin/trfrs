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

//! Program parameters (`TRFPARAMSET` in the reference implementation).

/// Maximum length of a FASTA sequence name, as stored in the C struct.
pub const MAX_SEQ_NAME_LEN: usize = 200;

/// Version reported by the program.
pub const VERSION: &str = "4.10.0";

/// Global knobs controlling the algorithm and the output.
#[derive(Clone, Debug)]
pub struct ParamSet {
    pub match_weight: u32,
    pub mismatch: u32,
    pub indel: u32,
    pub minscore: u32,
    pub maxperiod: u32,
    pub pm: u32,
    pub pi: u32,

    pub datafile: bool,
    pub maskedfile: bool,
    pub flankingsequence: bool,
    pub flankinglength: u32,
    pub htmloff: bool,
    pub redundoff: bool,
    pub ngs: bool,
    pub use_stdin: bool,
    pub maxwraplength: u32,

    pub inputfilename: String,
    pub outputprefix: String,
    #[allow(dead_code)]
    pub outputdirectory: String,
    pub outputfilename: String,

    /// Flags if the input file has more than one sequence.
    pub multisequencefile: bool,
    /// Sequence index, starting at 1.
    pub sequenceordinal: i32,
    /// Number of repeats found.
    pub outputcount: i32,
    pub running: bool,
    pub endstatus: Option<String>,
    /// Set when the reference implementation would have closed stdout
    /// (`fclose(destdfp)` in `-ngs` mode), which silences later messages.
    pub stdout_closed: bool,
    /// Progress percentage, `-1` means "hide the indicator".
    pub percent: i32,
}

impl Default for ParamSet {
    fn default() -> Self {
        ParamSet {
            match_weight: 0,
            mismatch: 0,
            indel: 0,
            minscore: 0,
            maxperiod: 0,
            pm: 0,
            pi: 0,
            datafile: false,
            maskedfile: false,
            flankingsequence: false,
            flankinglength: 500,
            htmloff: false,
            redundoff: false,
            ngs: false,
            use_stdin: false,
            maxwraplength: 2_000_000,
            inputfilename: String::new(),
            outputprefix: String::new(),
            outputdirectory: String::new(),
            outputfilename: String::new(),
            multisequencefile: false,
            sequenceordinal: 1,
            outputcount: 0,
            running: false,
            endstatus: None,
            stdout_closed: false,
            percent: 0,
        }
    }
}

impl ParamSet {
    /// `"%d.%d.%d.%d.%d.%d.%d"` parameter string used in output file names.
    pub fn param_string(&self) -> String {
        format!(
            "{}.{}.{}.{}.{}.{}.{}",
            self.match_weight,
            self.mismatch,
            self.indel,
            self.pm,
            self.pi,
            self.minscore,
            self.maxperiod
        )
    }

    /// `"Parameters: %d %d %d %d %d %d %d\n"` header line.
    pub fn parameters_header(&self) -> String {
        format!(
            "Parameters: {} {} {} {} {} {} {}\n",
            self.match_weight,
            self.mismatch,
            self.indel,
            self.pm,
            self.pi,
            self.minscore,
            self.maxperiod
        )
    }
}
