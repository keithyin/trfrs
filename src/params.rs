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

/// Lower bound of the period range selected by the STR mode (`--str`).
pub const STR_MIN_PERIOD: u32 = 1;

/// Upper bound of the period range selected by the STR mode (`--str`).
pub const STR_MAX_PERIOD: u32 = 15;

/// Global knobs controlling the algorithm and the output.
#[derive(Clone, Debug)]
pub struct ParamSet {
    pub match_weight: u32,
    pub mismatch: u32,
    pub indel: u32,
    pub minscore: u32,
    /// Lower bound of the reported repeat unit length (period).
    pub min_period: u32,
    pub maxperiod: u32,
    pub match_probability_percent: u32,
    pub indel_probability_percent: u32,

    pub data_file: bool,
    pub masked_file: bool,
    pub flanking_sequence: bool,
    pub flankinglength: u32,
    pub html_output_off: bool,
    pub redundancy_removal_off: bool,
    pub ngs: bool,
    pub use_stdin: bool,
    pub max_wrap_length: u32,

    pub input_file_name: String,
    pub output_prefix: String,
    #[allow(dead_code)]
    pub output_directory: String,
    pub output_file_name: String,

    /// Flags if the input file has more than one sequence.
    pub multi_sequence_file: bool,
    /// Sequence index, starting at 1.
    pub sequence_ordinal: i32,
    /// Number of repeats found.
    pub output_count: i32,
    pub running: bool,
    pub end_status: Option<String>,
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
            min_period: STR_MIN_PERIOD,
            maxperiod: 0,
            match_probability_percent: 0,
            indel_probability_percent: 0,
            data_file: false,
            masked_file: false,
            flanking_sequence: false,
            flankinglength: 500,
            html_output_off: false,
            redundancy_removal_off: false,
            ngs: false,
            use_stdin: false,
            max_wrap_length: 2_000_000,
            input_file_name: String::new(),
            output_prefix: String::new(),
            output_directory: String::new(),
            output_file_name: String::new(),
            multi_sequence_file: false,
            sequence_ordinal: 1,
            output_count: 0,
            running: false,
            end_status: None,
            stdout_closed: false,
            percent: 0,
        }
    }
}

impl ParamSet {
    /// `"%d.%d.%d.%d.%d.%d.%d"` parameter string used in output file names.
    ///
    /// The lower period bound of the STR mode is not part of the string: the
    /// reference implementation has no such parameter and both header and file
    /// names have to stay byte identical to the C program.
    pub fn param_string(&self) -> String {
        format!(
            "{}.{}.{}.{}.{}.{}.{}",
            self.match_weight,
            self.mismatch,
            self.indel,
            self.match_probability_percent,
            self.indel_probability_percent,
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
            self.match_probability_percent,
            self.indel_probability_percent,
            self.minscore,
            self.maxperiod
        )
    }
}
