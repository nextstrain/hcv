#! /usr/bin/env python
"""
Identify regions of a reference-based alignment (e.g. the output of
`nextclade run --output-fasta`) where more than a given fraction of the
sequences have a gap.

This is used during set-up of the dataset to determine the regions that are
masked before tree building (`mask_ranges` in config.yml) and ignored during
placement (`placementMaskRanges` in pathogen.json). It is not part of the
workflow itself.

Ranges are reported in 0-based, end-exclusive reference coordinates, which is
the convention of both `placementMaskRanges` and BED files.

With --config and/or --pathogen-json, the ranges are written into the
respective files, replacing `mask_ranges` / `placementMaskRanges`. Otherwise
they are printed to STDOUT.
"""

import argparse
import json
import re
import sys

import numpy as np
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True, help="Aligned FASTA (all sequences of reference length)")
    parser.add_argument("--max-gap-fraction", type=float, default=0.3,
        help="Report alignment columns in which more than this fraction of sequences have a gap (default: %(default)s)")
    parser.add_argument("--count-n", action="store_true",
        help="Count N as a gap (default: only '-' counts as gap)")
    parser.add_argument("--merge-distance", type=int, default=0,
        help="Merge ranges separated by at most this many non-gappy columns (default: %(default)s)")
    parser.add_argument("--min-length", type=int, default=1,
        help="Only report ranges of at least this length after merging (default: %(default)s)")
    parser.add_argument("--keep-ranges", type=str, nargs="*", default=[],
        help="Additional ranges to always include, as BEGIN-END (0-based, end-exclusive), e.g. 0-336 for the 5'UTR")
    parser.add_argument("--config", help="config.yml in which to replace `mask_ranges`")
    parser.add_argument("--pathogen-json", help="pathogen.json in which to replace `placementMaskRanges`")
    return parser.parse_args()


def gap_fractions(alignment, count_n):
    gap_chars = b"-N" if count_n else b"-"
    counts = None
    n_seqs = 0
    for record in SeqIO.parse(alignment, "fasta"):
        seq = np.frombuffer(str(record.seq).upper().encode(), dtype=np.uint8)
        if counts is None:
            counts = np.zeros(len(seq), dtype=int)
        elif len(seq) != len(counts):
            sys.exit(f"ERROR: sequence {record.id} has length {len(seq)}, expected {len(counts)}. Is this an alignment?")
        counts += np.isin(seq, np.frombuffer(gap_chars, dtype=np.uint8))
        n_seqs += 1
    return counts / n_seqs, n_seqs


def to_ranges(mask):
    """Convert a boolean array into a list of [begin, end) ranges of True values."""
    padded = np.concatenate([[False], mask, [False]]).astype(int)
    diff = np.diff(padded)
    return list(zip(np.flatnonzero(diff == 1), np.flatnonzero(diff == -1)))


def merge_ranges(ranges, distance):
    merged = []
    for begin, end in sorted(ranges):
        if merged and begin - merged[-1][1] <= distance:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([begin, end])
    return merged


def write_config(path, ranges):
    """Replace the `mask_ranges` entry of the YAML config, keeping comments and formatting elsewhere."""
    with open(path) as fh:
        text = fh.read()
    block = "mask_ranges:\n" + "".join(f"  - {{begin: {b}, end: {e}}}\n" for b, e in ranges) if ranges else "mask_ranges: []\n"
    pattern = re.compile(r"^mask_ranges:.*?(?=^\S|\Z)", re.MULTILINE | re.DOTALL)
    if not pattern.search(text):
        sys.exit(f"ERROR: no `mask_ranges` entry found in {path}")
    text = pattern.sub(lambda _: block + "\n", text, count=1)
    with open(path, "w") as fh:
        fh.write(text)


def write_pathogen_json(path, ranges):
    with open(path) as fh:
        pathogen = json.load(fh)
    pathogen["placementMaskRanges"] = [{"begin": b, "end": e} for b, e in ranges]
    with open(path, "w") as fh:
        json.dump(pathogen, fh, indent=2)
        fh.write("\n")


def main():
    args = parse_args()
    fractions, n_seqs = gap_fractions(args.alignment, args.count_n)

    ranges = to_ranges(fractions > args.max_gap_fraction)
    ranges += [tuple(int(x) for x in r.split("-")) for r in args.keep_ranges]
    ranges = [(int(b), int(e)) for b, e in merge_ranges(ranges, args.merge_distance) if e - b >= args.min_length]

    total = sum(e - b for b, e in ranges)
    print(f"{n_seqs} sequences, alignment length {len(fractions)}", file=sys.stderr)
    print(f"{len(ranges)} ranges with gap fraction > {args.max_gap_fraction} covering {total} sites:", file=sys.stderr)
    for b, e in ranges:
        print(f"  {b}-{e}\tlength {e - b}\tmean gap fraction {fractions[b:e].mean():.2f}", file=sys.stderr)

    if args.config:
        write_config(args.config, ranges)
        print(f"Wrote mask_ranges to {args.config}", file=sys.stderr)
    if args.pathogen_json:
        write_pathogen_json(args.pathogen_json, ranges)
        print(f"Wrote placementMaskRanges to {args.pathogen_json}", file=sys.stderr)
    if not (args.config or args.pathogen_json):
        json.dump([{"begin": b, "end": e} for b, e in ranges], sys.stdout, indent=2)
        print()


if __name__ == "__main__":
    main()
