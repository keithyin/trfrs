#!/usr/bin/env python3
"""SSR reference implementation, consistency checker and HTML report for `trfrs --ssr`.

The reference algorithm (ground truth) is an exhaustive exact scan that is
independent of TRF's heuristics: for every period p in 1..MAX_PERIOD it finds
the maximal p-periodic stretches of the sequence, normalizes the stretch unit
to its minimal period (skipping stretches already covered by a finer period)
and keeps the stretches with at least MIN_COPIES complete copies.  For every
stretch it records the *span* (the whole exact-periodic window, including a
possible partial copy) and the *tract* (the leftmost-anchored whole copies,
the canonical form the `--ssr` filter emits).

The checker verifies three levels:

1. self-consistency: every `--ssr` record is, checked directly on the
   sequence, an exact `motif x copies` region with copies >= MIN_COPIES and a
   period inside the reported range (catches filter bugs without the
   reference);
2. soundness: every `--ssr` tract lies inside a reference stretch of the same
   period (a violation would be a fabricated tract);
3. completeness: every reference stretch is accounted for:
   - when an STR region with a rotation-equivalent motif fully covers the
     stretch span and no larger tract interferes, the `--ssr` output must
     contain exactly the reference tract (a violation is a filter bug);
   - when a larger (or equal-copies leftmost) tract overlaps the span, the
     greedy overlap resolution may consume it (expected);
   - when the covering regions leave room for MIN_COPIES copies of the
     stretch but nothing is reported, that is a filter bug as well;
   - otherwise (unrelated region motifs / partial overlaps / no region at
     all) the stretch is not reachable for the filter -- counted as expected.

The exit status is non-zero when level 1/2 or the "filter bug" class of
level 3 is non-empty.

With `--html PATH` the same verification is rendered as a single-file HTML
report (Chinese, inline CSS, standard library only): a summary header with
the three-level verdict, and one section per sequence showing every `--ssr`
record with its surrounding bases (the SSR region highlighted in red, shaded
per motif copy, the motif spelled out), the reference stretches of the same
window for comparison, and the full reference result table.

Usage:
    python3 scripts/ssr_check.py <input.fa> <ssr.csv> [--str-regions <str.csv>]
                                 [--flank N] [--html report.html]

`ssr.csv` is the `--csv` output of `trfrs --ssr` on the same FASTA input,
`str.csv` optionally the `--csv` output of `trfrs --str`.  Coordinates in
both files are 1-based inclusive, `sequence_index` is the 1-based ordinal of
the FASTA record.
"""

import argparse
import csv
import html
import sys
from datetime import datetime

MAX_PERIOD = 15
MIN_COPIES = 3
DEFAULT_FLANK = 30

# Number of differences printed per level before the list is truncated.
MAX_LISTED = 20

STATUS_CN = {
    "accounted": "已报告",
    "resolved": "被更大 tract 消解",
    "mismatch": "区域 motif 不相关 / 部分重叠",
    "outside": "STR 区域外",
    "bug": "未报告（filter bug）",
}


def read_fasta(path):
    """Returns `(header, sequence)` pairs; sequences are uppercase
    letter-only strings, i.e. the same normalization trfrs applies to its
    input (only ASCII letters are kept)."""
    records = []
    header = None
    chunks = []
    with open(path, "rb") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if header is not None:
                    records.append((header, b"".join(chunks)))
                header = raw[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(bytes(c for c in line if 65 <= c <= 90 or 97 <= c <= 122))
    if header is not None:
        records.append((header, b"".join(chunks)))
    return [
        (head.decode("ascii", "replace"), chunk.decode("ascii", "replace").upper())
        for head, chunk in records
    ]


def read_csv_tracts(path):
    """Reads a trfrs `--csv` file into per-sequence tract lists.

    Returns a dict mapping the 1-based sequence index to a list of dicts with
    0-based half-open `start`/`end`, `period`, `copies`, `score`,
    `matches` and the normalized (`motif`, `motif_period`) pair."""
    per_sequence = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            index = int(row["sequence_index"])
            motif, motif_period = normalize_motif(row["repeat_unit"])
            per_sequence.setdefault(index, []).append(
                {
                    "start": int(row["start"]) - 1,
                    "end": int(row["end"]),
                    "period": int(row["period"]),
                    "copies": float(row["copy_number"]),
                    "score": int(row["score"]),
                    "matches": int(row["percent_matches"]),
                    "motif": motif,
                    "motif_period": motif_period,
                }
            )
    return per_sequence


def minimal_period(unit):
    """Smallest `mp` with `unit[i] == unit[i + mp]` for every valid `i`."""
    n = len(unit)
    for mp in range(1, n + 1):
        if unit[: n - mp] == unit[mp:]:
            return mp
    return n


def normalize_motif(pattern):
    """Collapses a motif to its minimal period, but only when the motif is a
    whole-number repetition of it (`TCTC` -> `TC`); a primitive motif such as
    `ATATA` keeps its full length.  Mirrors the Rust filter."""
    mp = minimal_period(pattern)
    if len(pattern) % mp == 0:
        return pattern[:mp], mp
    return pattern, len(pattern)


def same_motif(a, b):
    """True when two normalized motifs are rotations of each other."""
    return len(a) == len(b) and a in b + b


def overlap(a_start, a_end, b_start, b_end):
    """Length of the intersection of two half-open intervals."""
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def reference_stretches(seq):
    """All clean SSR stretches of one sequence (ground truth).

    Returns a list of dicts with 0-based half-open `span_start`/`span_end`
    (the maximal exact-periodic window), `tract_start`/`tract_end` (its
    leftmost-anchored whole copies), `period`, `unit` and `copies`."""
    n = len(seq)
    stretches = {}
    for period in range(1, MAX_PERIOD + 1):
        if MIN_COPIES * period > n:
            break
        match = [False] * n
        for i in range(period, n):
            match[i] = seq[i] == seq[i - period]
        run_start = None
        for i in range(n + 1):
            in_run = i < n and match[i]
            if in_run and run_start is None:
                run_start = i
            elif not in_run and run_start is not None:
                run_end = i  # exclusive; the match run covers [run_start, run_end)
                span_start = run_start - period
                span_end = run_end
                unit = seq[span_start : span_start + period]
                mp = minimal_period(unit)
                if mp < period and period % mp == 0:
                    run_start = None
                    continue  # the finer-period scan reports this stretch
                span_len = span_end - span_start
                copies = span_len // period
                if copies >= MIN_COPIES:
                    tract_end = span_start + copies * period
                    stretches[(span_start, span_end)] = {
                        "span_start": span_start,
                        "span_end": span_end,
                        "period": period,
                        "unit": unit,
                        "copies": copies,
                        "tract_start": span_start,
                        "tract_end": tract_end,
                    }
                run_start = None
    return sorted(stretches.values(), key=lambda s: (s["span_start"], s["span_end"]))


def annotate_tracts(ssr_tracts, str_regions, reference):
    """Annotates every `--ssr` tract with

    * `relation`: `exact` (same coordinates and motif as a reference tract),
      `in_span` (inside a same-period reference span -- region clipping or
      phase difference) or `none`;
    * `relation_stretch`: the reference stretch behind the relation;
    * `source`: the smallest `--str` region containing the tract, if any."""
    for index in sorted(ssr_tracts):
        for tract in ssr_tracts[index]:
            exact = None
            in_span = None
            for stretch in reference.get(index, []):
                if stretch["period"] != tract["motif_period"]:
                    continue
                if (
                    tract["start"] == stretch["tract_start"]
                    and tract["end"] == stretch["tract_end"]
                    and same_motif(tract["motif"], stretch["unit"])
                ):
                    exact = stretch
                elif (
                    stretch["span_start"] <= tract["start"]
                    and tract["end"] <= stretch["span_end"]
                    and (
                        in_span is None
                        or stretch["span_end"] - stretch["span_start"]
                        < in_span["span_end"] - in_span["span_start"]
                    )
                ):
                    in_span = stretch
            if exact is not None:
                tract["relation"], tract["relation_stretch"] = "exact", exact
            elif in_span is not None:
                tract["relation"], tract["relation_stretch"] = "in_span", in_span
            else:
                tract["relation"], tract["relation_stretch"] = "none", None

            source = None
            if str_regions is not None:
                for region in str_regions.get(index, []):
                    if (
                        region["start"] <= tract["start"]
                        and tract["end"] <= region["end"]
                        and (
                            source is None
                            or region["end"] - region["start"]
                            < source["end"] - source["start"]
                        )
                    ):
                        source = region
            tract["source"] = source


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fasta", help="FASTA input given to trfrs")
    parser.add_argument("ssr_csv", help="`trfrs --ssr --csv` output for the same input")
    parser.add_argument(
        "--str-regions",
        help="`trfrs --str --csv` output for the same input (enables the "
        "level 3 classification of missed stretches)",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=DEFAULT_FLANK,
        help=f"surrounding bases shown around each SSR region in the HTML "
        f"report (default {DEFAULT_FLANK})",
    )
    parser.add_argument(
        "--html",
        metavar="PATH",
        help="write the verification report as a single-file HTML page",
    )
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    names = [name for name, _ in records]
    sequences = [seq for _, seq in records]
    ssr_tracts = read_csv_tracts(args.ssr_csv)
    str_regions = read_csv_tracts(args.str_regions) if args.str_regions else None

    errors = 0

    # Reference stretches per sequence.
    reference = {}
    for index, seq in enumerate(sequences, start=1):
        reference[index] = reference_stretches(seq)
    annotate_tracts(ssr_tracts, str_regions, reference)

    # Level 1: self-consistency of the --ssr records.
    level1 = []
    for index in sorted(ssr_tracts):
        seq = sequences[index - 1]
        for tract in ssr_tracts[index]:
            region = seq[tract["start"] : tract["end"]]
            length = tract["end"] - tract["start"]
            copies = length // tract["motif_period"]
            exact = (
                length % tract["motif_period"] == 0
                and copies >= MIN_COPIES
                and 1 <= tract["motif_period"] <= MAX_PERIOD
                and region == tract["motif"] * copies
            )
            tract["exact_ok"] = exact
            if not exact:
                level1.append((index, tract, region))
    errors += len(level1)

    # Level 2: every --ssr tract lies inside the span of a reference stretch
    # of the same period.
    level2 = []
    for index in sorted(ssr_tracts):
        for tract in ssr_tracts[index]:
            inside = any(
                stretch["period"] == tract["motif_period"]
                and stretch["span_start"] <= tract["start"]
                and tract["end"] <= stretch["span_end"]
                for stretch in reference[index]
            )
            if not inside:
                level2.append((index, tract))
    errors += len(level2)

    # Level 3: every reference stretch is accounted for.
    accounted, bug, resolved, mismatch, outside = [], [], [], [], []
    for index in sorted(reference):
        seq_tracts = ssr_tracts.get(index, [])
        regions = str_regions.get(index, []) if str_regions else None
        for stretch in reference[index]:
            span_start, span_end = stretch["span_start"], stretch["span_end"]
            period, unit = stretch["period"], stretch["unit"]

            interference = any(
                overlap(tract["start"], tract["end"], span_start, span_end) > 0
                and (
                    tract["copies"] > stretch["copies"]
                    or (
                        tract["copies"] == stretch["copies"]
                        and tract["start"] < stretch["tract_start"]
                    )
                )
                for tract in seq_tracts
            )
            fully_covered = regions is not None and any(
                region["start"] <= span_start
                and span_end <= region["end"]
                and same_motif(unit, region["motif"])
                for region in regions
            )

            if fully_covered and not interference:
                exact = any(
                    tract["start"] == stretch["tract_start"]
                    and tract["end"] == stretch["tract_end"]
                    and same_motif(tract["motif"], unit)
                    for tract in seq_tracts
                )
                if exact:
                    stretch["status"] = "accounted"
                    accounted.append((index, stretch))
                else:
                    stretch["status"] = "bug"
                    stretch["status_note"] = "canonical tract missing"
                    bug.append((index, stretch, stretch["status_note"]))
            elif fully_covered and interference:
                stretch["status"] = "resolved"
                resolved.append((index, stretch))
            else:
                same_motif_overlap = any(
                    tract["motif_period"] == period
                    and same_motif(tract["motif"], unit)
                    and overlap(tract["start"], tract["end"], span_start, span_end) > 0
                    for tract in seq_tracts
                )
                if same_motif_overlap:
                    stretch["status"] = "accounted"
                    accounted.append((index, stretch))
                elif interference:
                    stretch["status"] = "resolved"
                    resolved.append((index, stretch))
                elif regions is not None:
                    reachable = any(
                        region["start"] <= span_start
                        and span_end <= region["end"]
                        and same_motif(unit, region["motif"])
                        for region in regions
                    )
                    room = any(
                        same_motif(unit, region["motif"])
                        and overlap(region["start"], region["end"], span_start, span_end)
                        >= MIN_COPIES * period
                        for region in regions
                    )
                    if reachable or room:
                        stretch["status"] = "bug"
                        stretch["status_note"] = (
                            "matching region leaves room but nothing reported"
                        )
                        bug.append((index, stretch, stretch["status_note"]))
                    elif any(
                        overlap(region["start"], region["end"], span_start, span_end) > 0
                        for region in regions
                    ):
                        stretch["status"] = "mismatch"
                        mismatch.append((index, stretch))
                    else:
                        stretch["status"] = "outside"
                        outside.append((index, stretch))
                else:
                    stretch["status"] = "outside"
                    outside.append((index, stretch))
    errors += len(bug)

    total_ssr = sum(len(v) for v in ssr_tracts.values())
    print(f"sequences: {len(sequences)}")
    print(f"--ssr records: {total_ssr}")
    print(f"reference stretches: {sum(len(v) for v in reference.values())}")
    print(f"accounted reference stretches: {len(accounted)}")
    print()
    print(f"level 1 violations (record is not an exact motif*copies region): {len(level1)}")
    for index, tract, region in level1[:MAX_LISTED]:
        print(
            f"  seq {index} [{tract['start']}, {tract['end']}) period={tract['motif_period']} "
            f"motif={tract['motif']!r} copies={tract['copies']} region={region[:40]!r}"
        )
    print(f"level 2 violations (tract outside any same-period reference span): {len(level2)}")
    for index, tract in level2[:MAX_LISTED]:
        print(
            f"  seq {index} [{tract['start']}, {tract['end']}) period={tract['motif_period']} "
            f"motif={tract['motif']!r}"
        )
    print(f"level 3 filter bugs (stretch reachable for the filter but not reported): {len(bug)}")
    for index, stretch, why in bug[:MAX_LISTED]:
        print(
            f"  seq {index} [{stretch['span_start']}, {stretch['span_end']}) "
            f"tract [{stretch['tract_start']}, {stretch['tract_end']}) period={stretch['period']} "
            f"unit={stretch['unit']!r}: {why}"
        )
    print(f"level 3 resolved by a larger tract (greedy overlap resolution): {len(resolved)}")
    for index, stretch in resolved[:MAX_LISTED]:
        print(
            f"  seq {index} [{stretch['span_start']}, {stretch['span_end']}) "
            f"period={stretch['period']} unit={stretch['unit']!r}"
        )
    print(f"level 3 motif mismatch / partial region overlap (expected): {len(mismatch)}")
    for index, stretch in mismatch[:MAX_LISTED]:
        print(
            f"  seq {index} [{stretch['span_start']}, {stretch['span_end']}) "
            f"period={stretch['period']} unit={stretch['unit']!r}"
        )
    print(f"level 3 outside STR regions (TRF detection stage, expected): {len(outside)}")
    for index, stretch in outside[:MAX_LISTED]:
        print(
            f"  seq {index} [{stretch['span_start']}, {stretch['span_end']}) "
            f"period={stretch['period']} unit={stretch['unit']!r}"
        )

    if args.html:
        render_html(
            args.html,
            args,
            names,
            sequences,
            ssr_tracts,
            str_regions,
            reference,
            len(level1),
            len(level2),
            bug,
            len(accounted),
            len(resolved),
            len(mismatch),
            len(outside),
        )
        print(f"\nHTML report written to {args.html}")

    if errors:
        print(f"\nFAILED: {errors} violation(s)", file=sys.stderr)
        return 1
    print("\nOK: --ssr output is consistent with the reference")
    return 0


# --------------------------------------------------------------------------
# HTML report
# --------------------------------------------------------------------------

_CSS = """
body { font-family: "PingFang SC", "Microsoft YaHei", system-ui, sans-serif;
       margin: 2rem auto; max-width: 75rem; color: #1c1c1c; background: #fff; }
h1 { font-size: 1.5rem; margin: 0 0 0.5rem; }
h2 { font-size: 1.15rem; margin: 2rem 0 0.6rem; padding: 0.4rem 0.7rem;
     background: #eef3f8; border-left: 4px solid #3b6ea5; }
h3 { font-size: 1rem; margin: 1.2rem 0 0.5rem; color: #3b6ea5; }
table { border-collapse: collapse; font-size: 0.82rem; margin: 0.4rem 0; }
th, td { border: 1px solid #d9d9d9; padding: 0.15rem 0.55rem; text-align: left; }
th { background: #f2f2f2; }
td.mono, .mono { font-family: ui-monospace, Menlo, Consolas, monospace; }
.meta-table td { padding: 0.2rem 0.9rem 0.2rem 0.2rem; }
.verdict { display: inline-block; padding: 0.3rem 0.9rem; border-radius: 4px;
           font-weight: 600; margin: 0.6rem 0; }
.verdict.ok { background: #e6f6e6; color: #1a7a1a; border: 1px solid #9fd39f; }
.verdict.failed { background: #fde8e8; color: #a11; border: 1px solid #eaa; }
.legend { font-size: 0.85rem; color: #444; background: #f7f7f7; padding: 0.5rem 0.8rem;
          border-radius: 4px; margin: 0.8rem 0; }
.tract { margin: 0.7rem 0 0.9rem; padding: 0.4rem 0.6rem; border: 1px solid #eee;
         border-radius: 4px; }
.tract-head { font-size: 0.85rem; margin-bottom: 0.25rem; }
.tract-head .meta { color: #555; margin-right: 0.9rem; }
.motif-desc b { color: #a11; font-family: ui-monospace, monospace; font-size: 0.95rem; }
.badge { display: inline-block; padding: 0.05rem 0.5rem; border-radius: 10px;
         font-size: 0.78rem; font-weight: 600; margin-right: 0.6rem; }
.badge.ok { background: #e6f6e6; color: #1a7a1a; }
.badge.warn { background: #fff4dd; color: #8a6d00; }
.badge.err { background: #fde8e8; color: #a11; }
.badge.info { background: #e8f0fb; color: #2b5aa0; }
.seq-line { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 0.82rem;
            white-space: pre; overflow-x: auto; background: #fafafa;
            padding: 0.25rem 0.4rem; border-radius: 4px; margin: 0.15rem 0;
            line-height: 1.5; }
.line-label { font-size: 0.75rem; color: #666; margin-top: 0.35rem; }
.flank { color: #777; }
.ssr .c0 { background: #e53935; color: #fff; }
.ssr .c1 { background: #f6999a; color: #501010; }
.ref-hit { background: #a8d4ff; color: #0b2a4a; }
.ref-left { background: #dcebff; color: #33567a; }
.ref-other { background: #e4e4e4; color: #444; }
.summary-line { color: #666; font-size: 0.9rem; }
details { margin: 0.3rem 0; }
summary { cursor: pointer; color: #3b6ea5; font-size: 0.88rem; }
.note { color: #666; font-size: 0.88rem; }
tr.bug-row td { background: #fde8e8; }
footer { margin-top: 2.5rem; color: #888; font-size: 0.8rem; }
"""


def _render_context_trfrs(seq, start, end, period, flank):
    """trfrs view: flanks plus the SSR region in red, shaded per motif copy.

    Returns the HTML fragment and the window bounds `[w0, w1)`."""
    w0 = max(0, start - flank)
    w1 = min(len(seq), end + flank)
    parts = []
    if w0 < start:
        parts.append(f'<span class="flank">{html.escape(seq[w0:start])}</span>')
    copies = (end - start) // period
    for i in range(copies):
        chunk = html.escape(seq[start + i * period : start + (i + 1) * period])
        shade = "c0" if i % 2 == 0 else "c1"
        parts.append(f'<span class="{shade}" title="第 {i + 1} 个完整拷贝">{chunk}</span>')
    leftover = end - (start + copies * period)
    if leftover:
        parts.append(f'<span class="flank">{html.escape(seq[end - leftover:end])}</span>')
    if end < w1:
        parts.append(f'<span class="flank">{html.escape(seq[end:w1])}</span>')
    return "".join(parts), w0, w1


def _render_context_reference(seq, w0, w1, stretches, tract):
    """Reference view: per-base painting of the reference stretches inside the
    same window.  For stretches overlapping the tract the canonical tract
    (whole copies) is painted solid blue and the span leftover (the partial
    copy the reference keeps but `--ssr` does not report) light blue; the
    remaining stretches fill unpainted bases (grey)."""
    n = w1 - w0
    classes = [None] * n
    titles = [None] * n
    ordered = sorted(
        stretches,
        key=lambda st: (
            0 if st["span_start"] < tract["end"] and tract["start"] < st["span_end"] else 1,
            st["span_start"],
            -st["span_end"],
        ),
    )
    for st in ordered:
        s0 = max(st["span_start"], w0)
        s1 = min(st["span_end"], w1)
        if s0 >= s1:
            continue
        title = (
            f"参考 stretch：1-based span [{st['span_start'] + 1}..{st['span_end']}]，"
            f"tract [{st['tract_start'] + 1}..{st['tract_end']}]，"
            f"period {st['period']}，motif {st['unit']} × {st['copies']}"
            f"（{STATUS_CN.get(st.get('status'), '?')}）"
        )
        if not (s0 < tract["end"] and tract["start"] < s1):
            ranges = [(s0, s1, "ref-other", title)]
        else:
            # 与该 tract 相交：完整拷贝部分实心蓝，span 残余（不完整拷贝）浅蓝
            ranges = []
            ta = max(st["tract_start"], s0)
            tb = min(st["tract_end"], s1)
            if ta < tb:
                ranges.append((ta, tb, "ref-hit", title))
            if st["span_start"] < ta:
                ranges.append((max(st["span_start"], s0), ta, "ref-left",
                               title + "；浅蓝 = span 残余（不完整拷贝，--ssr 不报告）"))
            if tb < st["span_end"]:
                ranges.append((tb, min(st["span_end"], s1), "ref-left",
                               title + "；浅蓝 = span 残余（不完整拷贝，--ssr 不报告）"))
        for a, b, cls, tip in ranges:
            for i in range(a, b):
                j = i - w0
                if classes[j] is None:
                    classes[j] = cls
                    titles[j] = tip
    parts = []
    i = 0
    while i < n:
        j = i
        while j < n and classes[j] == classes[i] and titles[j] == titles[i]:
            j += 1
        text = html.escape(seq[w0 + i : w0 + j])
        if classes[i] is None:
            parts.append(text)
        else:
            parts.append(
                f'<span class="{classes[i]}" title="{html.escape(titles[i], quote=True)}">{text}</span>'
            )
        i = j
    return "".join(parts)


def _tract_badge(tract):
    if not tract.get("exact_ok", True):
        return '<span class="badge err">✗ level 1 违例：不是精确的 motif × copies 区域</span>'
    if tract["relation"] == "exact":
        return '<span class="badge ok">✓ 与参考完全一致</span>'
    if tract["relation"] == "in_span":
        return '<span class="badge warn">≈ 位于参考 span 内（区域裁剪 / 相位差异）</span>'
    return '<span class="badge err">✗ 无同周期参考 span</span>'


def _render_tract_block(seq, tract, stretches, flank):
    start, end = tract["start"], tract["end"]
    period = tract["motif_period"]
    motif = tract["motif"]
    copies = (end - start) // period
    bp = end - start

    trfrs_line, w0, w1 = _render_context_trfrs(seq, start, end, period, flank)
    ref_line = _render_context_reference(seq, w0, w1, stretches, tract)

    source_html = ""
    source = tract.get("source")
    if source is not None:
        source_html = (
            f'<span class="meta">来源 STR 区域 1-based [{source["start"] + 1}..{source["end"]}]'
            f'（STR motif {html.escape(source["motif"])} × {source["copies"]:g}）</span>'
        )

    badge = _tract_badge(tract)
    head = (
        f'<div class="tract-head">{badge}'
        f'<span class="meta">1-based [{start + 1}..{end}]</span>'
        f'<span class="motif-desc">motif：<b>{html.escape(motif)}</b>'
        f'（period {period}）× {copies} 次完整重复 = {bp} bp</span>'
        f'<span class="meta">score {tract["score"]} · matches {tract["matches"]}%</span>'
        f"{source_html}</div>"
    )
    return (
        f'<div class="tract">{head}'
        f'<div class="line-label">trfrs --ssr（红色 = 检测到的 SSR 区域，深浅交替 = motif 的逐拷贝平铺）：</div>'
        f'<div class="seq-line ssr">{trfrs_line}</div>'
        f'<div class="line-label">Python 参考实现（实心蓝 = 参考 tract 的完整拷贝，浅蓝 = 参考 span 残余碱基'
        f'——不完整拷贝，--ssr 不报告，灰色 = 其他参考 stretch，悬浮查看详情）：</div>'
        f'<div class="seq-line ref">{ref_line}</div>'
        f"</div>"
    )


def _render_stretch_table(stretches, categories):
    rows = []
    for st in sorted(
        stretches,
        key=lambda s: (s.get("status") != "bug", s["span_start"], s["span_end"]),
    ):
        if st.get("status") not in categories:
            continue
        cls = "bug-row" if st.get("status") == "bug" else ""
        rows.append(
            f'<tr class="{cls}">'
            f'<td class="mono">[{st["span_start"] + 1}..{st["span_end"]}]</td>'
            f"<td>{st['period']}</td>"
            f"<td class=\"mono\">{html.escape(st['unit'])}</td>"
            f"<td>{st['copies']}</td>"
            f'<td class="mono">[{st["tract_start"] + 1}..{st["tract_end"]}]</td>'
            f"<td>{STATUS_CN.get(st.get('status'), '?')}"
            + (
                f'：{html.escape(st.get("status_note", ""))}'
                if st.get("status") == "bug" and st.get("status_note")
                else ""
            )
            + "</td></tr>"
        )
    if not rows:
        return None
    return (
        "<table><tr><th>参考 span（1-based）</th><th>period</th><th>motif</th>"
        "<th>完整拷贝</th><th>tract（1-based）</th><th>状态</th></tr>"
        + "".join(rows)
        + "</table>"
    )


def _render_sequence_section(seq_no, name, seq, tracts, stretches, flank, has_regions):
    counts = f"--ssr 记录 {len(tracts)} · 参考 stretch {len(stretches)}"
    head = f'<h2 id="seq-{seq_no}">序列 {seq_no}：{html.escape(name)} <span class="summary-line">（{counts}）</span></h2>'

    if tracts:
        blocks = "".join(_render_tract_block(seq, tract, stretches, flank) for tract in tracts)
        tracts_html = f"<h3>trfrs --ssr 检测结果（每条记录：区域标红 + motif 说明 + 前后侧翼）</h3>{blocks}"
    else:
        tracts_html = '<p class="note">该序列没有 --ssr 记录。</p>'

    main_table = _render_stretch_table(stretches, {"accounted", "resolved", "bug"})
    other_table = _render_stretch_table(stretches, {"mismatch", "outside"})
    bug_count = sum(1 for st in stretches if st.get("status") == "bug")

    ref_html = "<h3>Python 参考实现（穷举精确扫描）</h3>"
    if main_table:
        ref_html += main_table
    else:
        ref_html += '<p class="note">没有与 --ssr 结果相关的参考 stretch。</p>'
    if other_table:
        other_count = sum(1 for st in stretches if st.get("status") in ("mismatch", "outside"))
        ref_html += (
            f"<details><summary>其余参考 stretch：{other_count} 条"
            f"（STR 区域外 / 区域 motif 不相关，对过滤器不可达，属预期）</summary>{other_table}</details>"
        )
    if bug_count:
        ref_html = (
            f'<p class="note" style="color:#a11">⚠ 该序列存在 {bug_count} 条 filter bug（见上表红色行）。</p>'
            + ref_html
        )
    if not has_regions:
        ref_html += (
            '<p class="note">未提供 --str-regions，无法区分“区域外”与“区域 motif 不相关”，'
            "未报告的 stretch 均计入 STR 区域外。</p>"
        )

    return head + tracts_html + ref_html


def render_html(
    path,
    args,
    names,
    sequences,
    ssr_tracts,
    str_regions,
    reference,
    level1_count,
    level2_count,
    bug,
    accounted_count,
    resolved_count,
    mismatch_count,
    outside_count,
):
    total_ssr = sum(len(v) for v in ssr_tracts.values())
    total_ref = sum(len(v) for v in reference.values())
    errors = level1_count + level2_count + len(bug)
    verdict = (
        '<span class="verdict ok">OK：--ssr 输出与参考实现一致</span>'
        if errors == 0
        else f'<span class="verdict failed">FAILED：{errors} 条违例（level 1 {level1_count} · '
        f"level 2 {level2_count} · filter bug {len(bug)}）</span>"
    )

    parts = [
        "<!DOCTYPE html>\n<html lang=\"zh-CN\">\n<head>\n<meta charset=\"utf-8\">",
        "<title>trfrs --ssr 一致性验证报告</title>",
        f"<style>{_CSS}</style>\n</head>\n<body>",
        "<h1>trfrs --ssr 一致性验证报告</h1>",
        f'<p class="summary-line">生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}</p>',
        verdict,
        "<table class=\"meta-table\">",
        f"<tr><td>输入 FASTA</td><td class=\"mono\">{html.escape(args.fasta)}</td></tr>",
        f"<tr><td>--ssr CSV</td><td class=\"mono\">{html.escape(args.ssr_csv)}</td></tr>",
        "<tr><td>--str CSV</td><td class=\"mono\">"
        + (html.escape(args.str_regions) if args.str_regions else "未提供")
        + "</td></tr>",
        f"<tr><td>参数</td><td>period 1..{MAX_PERIOD} · 最少完整拷贝 {MIN_COPIES} · "
        f"侧翼 ±{args.flank} bp</td></tr>",
        f"<tr><td>规模</td><td>序列 {len(sequences)} 条 · --ssr 记录 {total_ssr} 条 · "
        f"参考 stretch {total_ref} 个</td></tr>",
        f"<tr><td>三层比对</td><td>level 1 违例 {level1_count} · level 2 违例 {level2_count} · "
        f"filter bug {len(bug)} · <a href=\"#legend\">参考 stretch 分类</a>：已报告 "
        f"{accounted_count} · 被更大 tract 消解 {resolved_count} · 区域 motif 不相关 "
        f"{mismatch_count} · STR 区域外 {outside_count}</td></tr>",
        "</table>",
        '<div class="legend" id="legend"><b>图例：</b>'
        "trfrs 行中<b>红色</b>为检测到的 SSR 区域 = motif × 完整拷贝数（深浅交替 = motif 的逐拷贝平铺），"
        "灰色为前后侧翼；参考行中<b>实心蓝</b>为参考 tract 的完整拷贝部分（与红色对应），"
        "<b>浅蓝</b>为参考 span 的残余碱基——极大精确周期窗口末尾的不完整拷贝，--ssr 刻意不在"
        "完整拷贝处截断故不报告（悬浮显示坐标、motif 与状态），灰色为窗口内其他参考 stretch。徽章："
        "<span class=\"badge ok\">✓ 与参考完全一致</span>"
        "<span class=\"badge warn\">≈ 位于参考 span 内</span>"
        "<span class=\"badge err\">✗ 违例</span>。"
        "参考 stretch 状态：已报告 / 被更大 tract 消解（贪心重叠消解，属预期）/ "
        "区域 motif 不相关或部分重叠（过滤器不可达，属预期）/ STR 区域外（TRF 检测阶段未发现，属预期）。"
        "</div>",
    ]

    for seq_no, (name, seq) in enumerate(zip(names, sequences), start=1):
        parts.append(
            _render_sequence_section(
                seq_no,
                name,
                seq,
                ssr_tracts.get(seq_no, []),
                reference.get(seq_no, []),
                args.flank,
                str_regions is not None,
            )
        )

    parts.append(
        "<footer>trfrs --ssr 一致性验证报告 · 由 scripts/ssr_check.py 生成 · "
        f"{datetime.now():%Y-%m-%d %H:%M:%S}</footer>"
    )
    parts.append("</body>\n</html>\n")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(parts))


if __name__ == "__main__":
    sys.exit(main())
