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

//! Binary entry point.
//!
//! The driver itself lives in the library ([`trfrs::run_cli`]) so that the
//! algorithm can be embedded into other programs without duplicating modules.

fn main() {
    /* The reference program does not ignore SIGPIPE: writing into a closed
    pipe terminates it.  Rust ignores SIGPIPE by default, which would turn the
    same situation into a panic, so restore the default disposition before any
    output is written. */
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    std::process::exit(trfrs::run_cli());
}
