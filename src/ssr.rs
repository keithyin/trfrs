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

//! SSR filtering: turns the reported STR records into clean `motif ×
//! repeats` tracts (`--ssr`).
//!
//! [`replace_with_ssr_tracts`] is a pure filter: it takes the records that
//! survived the regular cleanup and rewrites each one into the clean tract(s)
//! found inside its region.  The consensus of an STR record is reused as-is
//! (it is not re-derived from the sequence); only two things are redone on
//! the record's region:
//!
//! * normalization: the pattern is collapsed to its minimal period, but only
//!   when it is a whole-number repetition of it (`TCTC` -> `TC`); a primitive
//!   pattern such as `ATATA` (minimal period 2 does not divide 5) keeps its
//!   full length, because the region is not a repetition of the shorter unit;
//! * phase: every rotation of the normalized unit is anchored inside the
//!   region and extended copy by copy, so a tract is the leftmost-anchored
//!   window with the most complete copies.  `GAGAGAGAG` reported as `AG`
//!   yields `GA*4` (the trailing single base is not a complete copy and is
//!   dropped); between rotations with the same number of copies the leftmost
//!   one wins.
//!
//! One STR region can contain several runs of the motif, each of which
//! becomes its own tract.  A tract needs at least [`SSR_MIN_COPIES`] complete
//! copies; there is no per-tract score threshold because exact repeats have
//! no mismatches by construction (the source regions were still subject to
//! `--minscore` during the scan).
//!
//! Tracts that overlap substantially (at least half of the smaller span) are
//! resolved by the extraction priority (most complete copies, leftmost on
//! ties); a mere boundary touch is resolved by trimming the later tract by
//! whole motif periods, so neighbouring loci are not lost to TRF's boundary
//! drift.  See [`resolve_candidates`].

use crate::clean::RepeatRecord;

/// Minimum number of complete motif copies a clean SSR tract must contain.
pub(crate) const SSR_MIN_COPIES: u32 = 3;

/// A candidate clean tract (1-based inclusive coordinates).
struct Tract {
    first: i32,
    last: i32,
    period: i32,
    copies: u32,
    pattern: String,
}

/// Rewrites the STR `records` into clean SSR tracts.
///
/// `sequence` is the 1-based scan sequence (index 0 is the sentinel) the
/// record coordinates refer to; `min_period`/`max_period` are the reported
/// period range and `match_weight` the matching weight used for the score.
/// The returned records are sorted by position and carry fresh `count`s.
pub(crate) fn replace_with_ssr_tracts(
    sequence: &[u8],
    records: &[RepeatRecord],
    min_period: u32,
    max_period: u32,
    match_weight: u32,
) -> Vec<RepeatRecord> {
    let mut candidates = Vec::new();
    for record in records {
        extract_record_tracts(sequence, record, min_period, max_period, &mut candidates);
    }
    let accepted = resolve_candidates(candidates);
    accepted
        .iter()
        .enumerate()
        .map(|(index, tract)| build_record(tract, sequence, index as i32 + 1, match_weight))
        .collect()
}

/// Minimal period of `pattern`: the smallest `mp` with `pattern[i] ==
/// pattern[i + mp]` for every valid `i`.
fn minimal_period(pattern: &[u8]) -> usize {
    let len = pattern.len();
    for mp in 1..=len {
        if pattern[..len - mp] == pattern[mp..] {
            return mp;
        }
    }
    /* only an empty pattern gets here: with len >= 1 the mp == len step
    always matches (both slices are empty), so the loop returns len itself */
    len
}

/// True when `a` and `b` are rotations of the same cyclic word (`GA` and
/// `AG`); equal-length motifs of the same period are always rotations of each
/// other when they are exact powers of one unit.
fn same_motif(a: &str, b: &str) -> bool {
    /* the length guard is what keeps unequal periods apart (`GA` is a
    substring of `GAGA` doubled but not a rotation of it): only equal-length
    motifs can be cyclic shifts of one another */
    a.len() == b.len() && b.repeat(2).contains(a)
}

/// Collects the clean tracts of one STR record.
fn extract_record_tracts(
    sequence: &[u8],
    record: &RepeatRecord,
    min_period: u32,
    max_period: u32,
    out: &mut Vec<Tract>,
) {
    let pattern = record.pattern.as_bytes();
    /* an empty consensus would make `minimal_period` return 0 and the
    modulo below panic, and there is no unit to anchor anyway */
    if pattern.is_empty() {
        return;
    }

    let mp = minimal_period(pattern);
    /* collapse only whole-number repetitions (`TCTC` -> `TC`): when the
    minimal period does not divide the length (`ATATA`) the region is not a
    repetition of the shorter unit, so the full pattern stays the unit */
    let (unit, period) = if pattern.len() % mp == 0 {
        (&pattern[..mp], mp)
    } else {
        (pattern, pattern.len())
    };
    /* the scan already filtered the *reported* period, but normalization
    changes it (`TCTC` collapses to period 2, a primitive pattern keeps its
    full length), so the effective unit is re-checked against the range */
    if (period as u32) < min_period || (period as u32) > max_period {
        return;
    }

    let first = record.first as usize;
    let region = &sequence[first..=record.last as usize];

    /* anchor every rotation of the unit inside the region and extend copy by
    copy.  After a run the scan resumes at `e - period + 1`: anchors whose
    window lies inside the found run only yield contained sub-windows and are
    skipped, but the last `period - 1` positions must be retried -- their
    windows reach past the run and can start a longer, phase-shifted run. */
    for k in 0..period {
        let motif: Vec<u8> = unit[k..].iter().chain(unit[..k].iter()).copied().collect();
        let mut s = 0usize;
        while s + period <= region.len() {
            if region[s..s + period] != motif[..] {
                s += 1;
                continue;
            }
            let mut e = s + period;
            while e + period <= region.len() && region[e..e + period] == region[s..s + period] {
                e += period;
            }
            let copies = ((e - s) / period) as u32;
            /* below SSR_MIN_COPIES complete copies a run is not a usable
            tract (a few bases of any short unit match somewhere by chance);
            the window arithmetic already excludes a trailing partial copy,
            so `copies` counts only complete ones */
            if copies >= SSR_MIN_COPIES {
                out.push(Tract {
                    first: (first + s) as i32,
                    last: (first + s + copies as usize * period - 1) as i32,
                    period: period as i32,
                    copies,
                    pattern: String::from_utf8_lossy(&motif).into_owned(),
                });
            }
            s = e - period + 1;
        }
    }
}

impl Tract {
    /// span in bases of the 1-based inclusive interval
    fn span(&self) -> i64 {
        self.last as i64 - self.first as i64 + 1
    }
}

/// shared bases of two 1-based inclusive spans; 0 when disjoint (adjacent
/// intervals do not overlap)
fn overlap_bases(a: &Tract, b: &Tract) -> i64 {
    let first = a.first.max(b.first) as i64;
    let last = a.last.min(b.last) as i64;
    (last - first + 1).max(0)
}

/// true when the overlap covers at least half of the smaller tract, the
/// threshold of smicing's `str_polish::detector` (`IOI_MIN_OVERLAP_FRACTION`).
/// a boundary touch (trfrs routinely extends a repeat 1-2 bp into its
/// neighbour) is not a conflict and is trimmed by
/// [`trim_boundary_overlaps`] instead.
fn regions_substantially_overlap(a: &Tract, b: &Tract) -> bool {
    let overlap = overlap_bases(a, b);
    if overlap == 0 {
        return false;
    }
    let smaller_span = a.span().min(b.span());
    /* overlap >= 0.5 * smaller_span, in exact integer arithmetic */
    overlap * 2 >= smaller_span
}

/// Keep the higher-priority tract of every substantial conflict; priority is
/// unchanged from the historical greedy step (most complete copies first,
/// leftmost start on ties).  Boundary touches survive and are trimmed by
/// [`trim_boundary_overlaps`].
fn drop_substantially_overlapping(mut tracts: Vec<Tract>) -> Vec<Tract> {
    tracts.sort_by_key(|tract| (std::cmp::Reverse(tract.copies), tract.first));
    let mut kept: Vec<Tract> = Vec::with_capacity(tracts.len());
    for tract in tracts {
        if kept
            .iter()
            .any(|kept| regions_substantially_overlap(kept, &tract))
        {
            continue;
        }
        kept.push(tract);
    }
    kept
}

/// Make the surviving tracts disjoint by trimming the LATER one (by start)
/// instead of dropping it: `first` moves right by `units = ceil(overlap /
/// period)` whole periods and `copies` shrinks by `units`.  The shift is a
/// multiple of the period, so the phase-aligned `pattern` still rebuilds the
/// tract exactly; `last`, `period` and `pattern` are untouched.  A tract left
/// with fewer than [`SSR_MIN_COPIES`] copies is dropped rather than reported
/// as a truncated repeat.  Ported from smicing's `trim_boundary_overlaps`.
fn trim_boundary_overlaps(mut tracts: Vec<Tract>) -> Vec<Tract> {
    tracts.sort_by_key(|tract| (tract.first, tract.last));
    let mut kept: Vec<Tract> = Vec::with_capacity(tracts.len());
    for mut tract in tracts {
        let Some(prev) = kept.last() else {
            kept.push(tract);
            continue;
        };
        /* substantial conflicts were resolved above, so containment is gone:
        each kept tract ends after the previous one and the last kept tract
        carries the maximum end (trimming only moves `first` right) */
        let prev_last = prev.last;
        if tract.first <= prev_last {
            let overlap = prev_last as i64 - tract.first as i64 + 1;
            let period = tract.period as i64;
            /* overlap > 0 and period > 0, so ceiling division is exact */
            let units = ((overlap + period - 1) / period) as u64;
            let copies = tract.copies as u64;
            if units >= copies || copies - units < SSR_MIN_COPIES as u64 {
                continue;
            }
            tract.first = (tract.first as i64 + units as i64 * period) as i32;
            tract.copies = (copies - units) as u32;
            debug_assert_eq!(
                tract.first as i64 + tract.copies as i64 * period - 1,
                tract.last as i64
            );
        }
        kept.push(tract);
    }
    kept
}

/// Deduplicates and resolves the candidates of all records.
///
/// Identical spans (they can arise from overlapping source records) collapse;
/// a tract fully contained in a *same cyclic word* tract is a phase-shifted
/// nested artifact and is dropped, while containment by a different motif
/// (for example an `AT*4` window inside an exact `ATATATATC*3` region) is a
/// genuine finer repeat and is kept.  Remaining overlaps are resolved in two
/// steps, ported from smicing's `str_polish::detector`:
///
/// * a **substantial** conflict (the overlap covers at least half of the
///   smaller tract) keeps the higher-priority tract -- most copies first,
///   leftmost start on ties -- so in the different-motif case the finer tract
///   usually wins by copy number and its container is dropped;
/// * a mere **boundary touch** is not a conflict: the later tract is trimmed
///   by whole motif periods (see [`trim_boundary_overlaps`]), so trfrs
///   extending a repeat 1-2 bp into its neighbour no longer loses a genuine
///   locus.
fn resolve_candidates(mut candidates: Vec<Tract>) -> Vec<Tract> {
    candidates.sort_by_key(|tract| (tract.first, tract.last));
    candidates.dedup_by(|a, b| a.first == b.first && a.last == b.last);

    candidates.sort_by_key(|tract| (tract.first, std::cmp::Reverse(tract.last)));
    let mut outer: Vec<Tract> = Vec::with_capacity(candidates.len());
    for tract in candidates {
        let nested = outer.iter().any(|kept| {
            kept.first <= tract.first
                && tract.last <= kept.last
                && same_motif(&kept.pattern, &tract.pattern)
        });
        if !nested {
            outer.push(tract);
        }
    }

    let mut accepted = trim_boundary_overlaps(drop_substantially_overlapping(outer));

    accepted.sort_by_key(|tract| (tract.first, tract.last));
    accepted
}

/// Builds the reported record of a clean tract.
fn build_record(tract: &Tract, sequence: &[u8], count: i32, match_weight: u32) -> RepeatRecord {
    let bases = &sequence[tract.first as usize..=tract.last as usize];
    let total = bases.len() as f64;

    let mut acgt = [0i32; 4];
    for &base in bases {
        match base {
            b'A' => acgt[0] += 1,
            b'C' => acgt[1] += 1,
            b'G' => acgt[2] += 1,
            b'T' => acgt[3] += 1,
            /* the scan sequence may hold non-ACGT bytes (the scanner skips
            them): they count toward the denominator `total` but toward no
            base percentage, mirroring the alignment report */
            _ => {}
        }
    }

    /* base composition and entropy, computed like the alignment report */
    let mut diversity = [
        acgt[0] as f64 / total,
        acgt[1] as f64 / total,
        acgt[2] as f64 / total,
        acgt[3] as f64 / total,
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
    /* every accumulated term `frequency * log2(frequency)` is <= 0, so the
    running sum is the negated entropy and one sign flip yields the reported
    value; an exact 0.0 stays 0.0 instead of becoming -0.0 */
    if entropy < 0.0 {
        entropy = -entropy;
    }

    /* the score a perfect repeat reaches: every complete copy beyond the
    first contributes `period` matched bases at the matching weight; the
    product is built in i64 and clamped because a very long tract with a
    large `--matchweight` would overflow the reported i32 field */
    let score = ((tract.copies as i64 - 1) * tract.period as i64 * match_weight as i64)
        .min(i32::MAX as i64) as i32;

    let copies = tract.copies as f32;
    RepeatRecord {
        count,
        reference: format!(
            "{}--{},{},{:3.1},{},{}",
            tract.first, tract.last, tract.period, copies, tract.period, count
        ),
        first: tract.first,
        last: tract.last,
        period: tract.period,
        copies,
        size: tract.period,
        matches: 100,
        indels: 0,
        score,
        a_count: (100.0 * acgt[0] as f64 / total) as i32,
        c_count: (100.0 * acgt[1] as f64 / total) as i32,
        g_count: (100.0 * acgt[2] as f64 / total) as i32,
        t_count: (100.0 * acgt[3] as f64 / total) as i32,
        entropy: entropy as f32,
        pattern: tract.pattern.clone(),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// A 1-based scan sequence (index 0 is the sentinel).
    fn seq(bases: &str) -> Vec<u8> {
        let mut sequence = vec![0u8];
        sequence.extend_from_slice(bases.as_bytes());
        sequence
    }

    /// A synthetic STR record over `first..=last`.
    fn str_record(first: i32, last: i32, pattern: &str) -> RepeatRecord {
        RepeatRecord {
            count: 1,
            reference: String::new(),
            first,
            last,
            period: pattern.len() as i32,
            copies: 2.0,
            size: pattern.len() as i32,
            matches: 90,
            indels: 5,
            score: 50,
            a_count: 0,
            c_count: 0,
            g_count: 0,
            t_count: 0,
            entropy: 0.0,
            pattern: pattern.to_string(),
        }
    }

    /// Runs the filter and returns `(first, last, period, pattern, copies)`.
    fn tracts(sequence: &[u8], records: &[RepeatRecord]) -> Vec<(i32, i32, i32, String, u32)> {
        replace_with_ssr_tracts(sequence, records, 1, 15, 2)
            .into_iter()
            .map(|record| {
                (
                    record.first,
                    record.last,
                    record.period,
                    record.pattern,
                    record.copies as u32,
                )
            })
            .collect()
    }

    #[test]
    fn minimal_period_collapses_repetitions() {
        assert_eq!(minimal_period(b"TC"), 2);
        assert_eq!(minimal_period(b"TCTC"), 2);
        assert_eq!(minimal_period(b"GAGAGAGAGA"), 2);
        assert_eq!(minimal_period(b"ATATA"), 2);
        assert_eq!(minimal_period(b"TATATATATATAG"), 13);
        assert_eq!(minimal_period(b"A"), 1);
        assert_eq!(minimal_period(b"AAAAAAAA"), 1);
    }

    #[test]
    fn prefers_the_rotation_with_more_complete_copies() {
        /* the whole region is AG*7; GA*6 would leave one base out */
        let sequence = seq("AGAGAGAGAGAGAG");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 14, "GA")]),
            vec![(1, 14, 2, "AG".to_string(), 7)]
        );
    }

    #[test]
    fn leftmost_rotation_wins_a_tie() {
        /* GAGAGAGAG: GA*4 at offset 0 and AG*4 at offset 1 are equal, the
        leftmost one is reported and the trailing base is dropped */
        let sequence = seq("GAGAGAGAG");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 9, "AG")]),
            vec![(1, 8, 2, "GA".to_string(), 4)]
        );
    }

    #[test]
    fn reanchors_the_phase_after_a_broken_prefix() {
        let sequence = seq("CGAGAGAGAG");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 10, "AG")]),
            vec![(2, 9, 2, "GA".to_string(), 4)]
        );
    }

    #[test]
    fn a_short_anchor_run_does_not_hide_a_longer_one() {
        /* "GCG" at offset 0 is a one-copy run; skipping it entirely used to
        skip the offset 2 anchor, where the real GCG*5 run starts (regression
        test for the run-skip rule) */
        let sequence = seq("GCGCGGCGGCGGCGGCGG");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 18, "GGC")]),
            vec![(3, 17, 3, "GCG".to_string(), 5)]
        );
    }

    #[test]
    fn normalizes_the_pattern_to_its_minimal_period() {
        let sequence = seq("TCTCTCTCTC");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 10, "TCTC")]),
            vec![(1, 10, 2, "TC".to_string(), 5)]
        );
    }

    #[test]
    fn keeps_primitive_patterns_at_full_length() {
        /* ATATA has minimal period 2, but 2 does not divide 5: the region is
        not a repetition of AT, so the motif stays ATATA */
        let sequence = seq("ATATAATATAATATA");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 15, "ATATA")]),
            vec![(1, 15, 5, "ATATA".to_string(), 3)]
        );
    }

    #[test]
    fn reports_one_tract_per_run() {
        /* (GA)4 + junk + (GA)3 inside one STR region */
        let sequence = seq("GAGAGAGACTTGAGAGA");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 17, "GA")]),
            vec![
                (1, 8, 2, "GA".to_string(), 4),
                (12, 17, 2, "GA".to_string(), 3)
            ]
        );
    }

    #[test]
    fn drops_tracts_below_the_minimum_copy_number() {
        /* exactly SSR_MIN_COPIES copies are still reported */
        let sequence = seq("GAGAGA");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 6, "GA")]),
            vec![(1, 6, 2, "GA".to_string(), 3)]
        );
        /* (GA)2 has only two complete copies and is dropped */
        let sequence = seq("GAGA");
        assert!(tracts(&sequence, &[str_record(1, 4, "GA")]).is_empty());
    }

    #[test]
    fn skips_patterns_outside_the_period_range() {
        let sequence = seq("ACGTACGTACGTACGTAC");
        /* minimal period 4 does not divide 18, so the pattern stays an
        18-mer, which exceeds --maxperiod 15 */
        assert!(tracts(&sequence, &[str_record(1, 18, "ACGTACGTACGTACGTAC")]).is_empty());
    }

    #[test]
    fn overlapping_source_records_yield_one_tract() {
        /* the same locus reported twice: the identical span collapses and the
        remaining overlaps lose against the leftmost, longest tract */
        let sequence = seq("GAGAGAGAGAGA");
        let records = [str_record(1, 10, "GA"), str_record(3, 12, "GA")];
        assert_eq!(
            tracts(&sequence, &records),
            vec![(1, 10, 2, "GA".to_string(), 5)]
        );
    }

    #[test]
    fn a_nested_tract_loses_against_its_container() {
        /* (ATG)3 sits inside (GAT)4: the container covers it completely */
        let sequence = seq("GATGATGATGAT");
        let records = [str_record(1, 12, "GAT"), str_record(2, 10, "ATG")];
        assert_eq!(
            tracts(&sequence, &records),
            vec![(1, 12, 3, "GAT".to_string(), 4)]
        );
    }

    #[test]
    fn a_finer_period_tract_survives_a_different_motif_container() {
        /* every AT*4 window lies inside the exact ATATATATC*3 region; the
        coarser region must not swallow the finer, different-motif tracts */
        let sequence = seq("ATATATATCATATATATCATATATATCA");
        let records = [str_record(1, 27, "AT"), str_record(1, 27, "ATATATATC")];
        assert_eq!(
            tracts(&sequence, &records),
            vec![
                (1, 8, 2, "AT".to_string(), 4),
                (10, 17, 2, "AT".to_string(), 4),
                (19, 26, 2, "AT".to_string(), 4)
            ]
        );
    }

    #[test]
    fn a_minor_boundary_overlap_trims_the_later_tract() {
        /* (AT)^7; the [7,14] record shares two bases with [1,8] (25% of the
        smaller span): not redundant, the later tract loses one whole period */
        let sequence = seq("ATATATATATATAT");
        assert_eq!(
            tracts(
                &sequence,
                &[str_record(1, 8, "AT"), str_record(7, 14, "AT")]
            ),
            vec![
                (1, 8, 2, "AT".to_string(), 4),
                /* first moved 7 -> 9, copies 4 -> 3, last and pattern intact */
                (9, 14, 2, "AT".to_string(), 3),
            ]
        );
    }

    #[test]
    fn a_trim_that_leaves_fewer_than_three_copies_drops_the_tract() {
        /* (AT)^6: the [7,12] tract only has three copies and loses one to the
        boundary overlap, so it is dropped instead of reported truncated */
        let sequence = seq("ATATATATATAT");
        assert_eq!(
            tracts(
                &sequence,
                &[str_record(1, 8, "AT"), str_record(7, 12, "AT")]
            ),
            vec![(1, 8, 2, "AT".to_string(), 4)]
        );
    }

    #[test]
    fn an_exact_half_overlap_is_substantial() {
        /* 4 of the 8 bases overlap: exactly 50%, which is still a substantial
        conflict, so the later tract is dropped rather than trimmed */
        let sequence = seq("ATATATATATAT");
        assert_eq!(
            tracts(
                &sequence,
                &[str_record(1, 8, "AT"), str_record(5, 12, "AT")]
            ),
            vec![(1, 8, 2, "AT".to_string(), 4)]
        );
    }

    #[test]
    fn a_substantial_conflict_uses_priority_not_position() {
        /* exactly 50% overlap again: the 12-copy tract wins although it starts
        later, so priority (copies desc) still decides substantial conflicts */
        let sequence = seq("AAAAAAAAAAAAAAAA");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 8, "A"), str_record(5, 16, "A")]),
            vec![(5, 16, 1, "A".to_string(), 12)]
        );
    }

    #[test]
    fn a_minor_overlap_trims_the_later_tract_even_when_it_has_more_copies() {
        /* 3/8 = 37.5% overlap: position decides (the later tract is trimmed),
        not priority, so the smaller tract keeps its full span */
        let sequence = seq("AAAAAAAAAAAAAAAA");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 8, "A"), str_record(6, 16, "A")]),
            vec![
                (1, 8, 1, "A".to_string(), 8),
                (9, 16, 1, "A".to_string(), 8),
            ]
        );
    }

    #[test]
    fn successive_boundary_trims_use_the_unchanged_end() {
        /* (AT)^10; the third trim is computed from the second tract's
        unchanged `last`, not from its shifted `first` */
        let sequence = seq("ATATATATATATATATATAT");
        assert_eq!(
            tracts(
                &sequence,
                &[
                    str_record(1, 8, "AT"),
                    str_record(7, 14, "AT"),
                    str_record(13, 20, "AT"),
                ]
            ),
            vec![
                (1, 8, 2, "AT".to_string(), 4),
                (9, 14, 2, "AT".to_string(), 3),
                (15, 20, 2, "AT".to_string(), 3),
            ]
        );
    }

    #[test]
    fn touching_tracts_are_not_trimmed() {
        /* `prev_last + 1 == first` is adjacency, not an overlap */
        let sequence = seq("ATATATATATATAT");
        assert_eq!(
            tracts(
                &sequence,
                &[str_record(1, 8, "AT"), str_record(9, 14, "AT")]
            ),
            vec![
                (1, 8, 2, "AT".to_string(), 4),
                (9, 14, 2, "AT".to_string(), 3),
            ]
        );
    }

    #[test]
    fn builds_the_record_fields_of_a_clean_tract() {
        let sequence = seq("GAGAGAGAG");
        let records = replace_with_ssr_tracts(&sequence, &[str_record(1, 9, "AG")], 1, 15, 2);
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(record.count, 1);
        assert_eq!(record.first, 1);
        assert_eq!(record.last, 8);
        assert_eq!(record.period, 2);
        assert_eq!(record.size, 2);
        assert_eq!(record.copies, 4.0);
        assert_eq!(record.matches, 100);
        assert_eq!(record.indels, 0);
        assert_eq!(record.score, 12); /* (4 - 1) * 2 * 2 */
        /* GA*4: half G, half A */
        assert_eq!(record.a_count, 50);
        assert_eq!(record.g_count, 50);
        assert_eq!(record.c_count, 0);
        assert_eq!(record.t_count, 0);
        assert_eq!(record.entropy, 1.0);
        assert_eq!(record.pattern, "GA");
    }

    #[test]
    fn mono_nucleotide_tracts_are_found() {
        let sequence = seq("AAAAAAAAAAAA");
        assert_eq!(
            tracts(&sequence, &[str_record(1, 12, "A")]),
            vec![(1, 12, 1, "A".to_string(), 12)]
        );
    }
}
