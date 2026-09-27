#! /usr/bin/env python
"""
Generate an Auspice colors TSV for HCV genotypes and subtypes.

Each genotype is assigned a distinct base hue. Subtypes of a genotype get hues
within a bracket around the genotype's hue, alternating in lightness so that
neighbouring subtypes remain distinguishable. Colors are constructed in OKLCH
(perceptually uniform lightness/chroma/hue) and converted to sRGB, reducing
chroma where needed to stay within the sRGB gamut.

Subtypes are read from the given metadata column(s), so the scale covers all
subtypes present in the data, e.g. the ingest output. Colors for a given set
of subtypes are deterministic, but adding a subtype to a genotype shifts the
hues of the other subtypes of that genotype.

The same colors are written for each of the given genotype and subtype traits
(e.g. both `genotype` and `ncbi_genotype`).
"""

import argparse
import csv
import re
import sys

import numpy as np

# Base color (OKLCH hue in degrees, lightness, chroma) for each genotype.
# Hues are spread around the colour wheel and lightness varies between
# genotypes, which improves separation for colour vision deficiencies. These
# values were optimised to maximise the worst-pair color distance between all
# eight genotypes under protan/deutan simulation (Delta E 6.8, OKLab x100),
# while keeping the normal-vision distance between all pairs >= 15.
GENOTYPE_COLORS = {
    "1": (242, 0.51, 0.13),  # blue
    "2": (158, 0.50, 0.14),  # green
    "3": (26, 0.72, 0.17),   # salmon red
    "4": (312, 0.69, 0.16),  # lilac
    "5": (46, 0.53, 0.15),   # dark orange
    "6": (183, 0.67, 0.14),  # teal
    "7": (96, 0.65, 0.14),   # olive yellow
    "8": (331, 0.46, 0.15),  # plum
}
RECOMBINANT_COLOR = "#8c8c8c"

# Subtypes: hues are spread evenly within a bracket of at most +/- HUE_BRACKET
# around the genotype hue (narrower for genotypes with close neighbours),
# and neighbouring subtypes are offset in lightness from the genotype color,
# cycling through SUBTYPE_LIGHTNESS_OFFSETS (ordered such that consecutive
# subtypes differ by >= 0.16 in lightness).
HUE_BRACKET = 16
SUBTYPE_LIGHTNESS_OFFSETS = [-0.16, 0.0, 0.16, -0.08, 0.08]
SUBTYPE_LIGHTNESS_RANGE = (0.38, 0.85)
# Minimal OKLab distance (x100) of a subtype color to the colors of other genotypes
MIN_DELTA_E_OTHER_GENOTYPE = 10
# Minimal OKLab distance (x100) between consecutive subtypes of a genotype
MIN_DELTA_E_SIBLING = 10


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metadata", required=True, help="Metadata TSV with subtype annotations")
    parser.add_argument("--subtype-columns", nargs="+", default=["ncbi_subtype"],
        help="Columns of the metadata to read subtypes from (default: %(default)s)")
    parser.add_argument("--extra-subtypes", nargs="*", default=["7a", "8a"],
        help="Subtypes to include in addition to those found in the metadata (default: %(default)s)")
    parser.add_argument("--genotype-traits", nargs="+", default=["clade_membership", "genotype", "ncbi_genotype"],
        help="Traits (coloring keys) to write genotype colors for (default: %(default)s)")
    parser.add_argument("--subtype-traits", nargs="+", default=["subtype", "ncbi_subtype"],
        help="Traits (coloring keys) to write subtype colors for (default: %(default)s)")
    parser.add_argument("--output", required=True, help="Output colors TSV")
    return parser.parse_args()


def oklch_to_srgb(L, C, h):
    """Convert OKLCH to sRGB in [0, 1] (may be out of gamut)."""
    a = C * np.cos(np.radians(h))
    b = C * np.sin(np.radians(h))
    l_ = L + 0.3963377774 * a + 0.2158037573 * b
    m_ = L - 0.1055613458 * a - 0.0638541728 * b
    s_ = L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    rgb_linear = np.array([
        +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s,
    ])
    return np.where(
        rgb_linear <= 0.0031308,
        12.92 * rgb_linear,
        1.055 * np.sign(rgb_linear) * np.abs(rgb_linear) ** (1 / 2.4) - 0.055,
    )


def oklch_to_hex(L, C, h):
    """Convert OKLCH to hex, reducing chroma until the color is within the sRGB gamut."""
    while C > 0:
        rgb = oklch_to_srgb(L, C, h)
        if np.all(rgb >= -1e-6) and np.all(rgb <= 1 + 1e-6):
            break
        C -= 0.005
    rgb = np.clip(oklch_to_srgb(L, max(C, 0), h), 0, 1)
    return "#" + "".join(f"{round(x * 255):02x}" for x in rgb)


def subtype_sort_key(subtype):
    # e.g. "1a" < "1b" < ... < "1z" < "6xa"
    match = re.match(r"^(\d+)([a-z]*)$", subtype)
    return (int(match.group(1)), len(match.group(2)), match.group(2))


def read_subtypes(metadata, columns):
    subtypes = set()
    with open(metadata) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for column in columns:
                value = (row.get(column) or "").strip().lower()
                # skip recombinants (e.g. "2k/1b") and values not of the form <genotype><letters>
                if re.match(r"^[1-8][a-z]+$", value):
                    subtypes.add(value)
    return sorted(subtypes, key=subtype_sort_key)


def hue_bracket(genotype):
    """
    Half-width of the hue bracket of a genotype: at most HUE_BRACKET, and
    narrower for genotypes with a close neighbouring hue, such that the
    brackets of different genotypes do not overlap.
    """
    hue = GENOTYPE_COLORS[genotype][0]
    nearest = min(
        min(abs(hue - other[0]), 360 - abs(hue - other[0]))
        for g, other in GENOTYPE_COLORS.items() if g != genotype
    )
    return min(HUE_BRACKET, 0.45 * nearest)


def hex_to_oklab(color):
    rgb = np.array([int(color[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    l = 0.4122214708 * lin[0] + 0.5363325363 * lin[1] + 0.0514459929 * lin[2]
    m = 0.2119034982 * lin[0] + 0.6806995451 * lin[1] + 0.1073969566 * lin[2]
    s = 0.0883024619 * lin[0] + 0.2817188376 * lin[1] + 0.6299787005 * lin[2]
    l_, m_, s_ = np.cbrt([l, m, s])
    return np.array([
        0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
        1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
        0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_,
    ])


def delta_e(color, others):
    """Smallest OKLab distance (x100) between a color and a list of colors."""
    if not others:
        return np.inf
    lab = hex_to_oklab(color)
    return min(100 * np.linalg.norm(lab - hex_to_oklab(other)) for other in others)


def subtype_colors(subtypes, genotype_colors):
    """
    Assign subtype colors. For each subtype, the lightness levels are tried
    starting from its position in the SUBTYPE_LIGHTNESS_OFFSETS cycle; the first
    level that keeps the color at least MIN_DELTA_E_OTHER_GENOTYPE away from the
    colors of other genotypes (and their subtypes), and MIN_DELTA_E_SIBLING
    away from the previous subtype of the same genotype, is used.
    """
    colors = {}
    colors_by_genotype = {g: [c] for g, c in genotype_colors.items()}
    for genotype, (hue, genotype_lightness, chroma) in GENOTYPE_COLORS.items():
        members = [s for s in subtypes if s[0] == genotype]
        n = len(members)
        if n == 1:
            colors[members[0]] = genotype_colors[genotype]
            continue
        bracket = hue_bracket(genotype)
        for i, subtype in enumerate(members):
            offset = -bracket + 2 * bracket * i / (n - 1)
            other_colors = [c for g, cs in colors_by_genotype.items() if g != genotype for c in cs]
            candidates = []
            for k in range(len(SUBTYPE_LIGHTNESS_OFFSETS)):
                lightness = genotype_lightness + SUBTYPE_LIGHTNESS_OFFSETS[(i + k) % len(SUBTYPE_LIGHTNESS_OFFSETS)]
                lightness = min(max(lightness, SUBTYPE_LIGHTNESS_RANGE[0]), SUBTYPE_LIGHTNESS_RANGE[1])
                candidates.append(oklch_to_hex(lightness, chroma, (hue + offset) % 360))
            previous = colors_by_genotype[genotype][-1:] if i > 0 else []
            color = next(
                (c for c in candidates
                 if delta_e(c, other_colors) >= MIN_DELTA_E_OTHER_GENOTYPE
                 and delta_e(c, previous) >= MIN_DELTA_E_SIBLING),
                max(candidates, key=lambda c: min(delta_e(c, other_colors), delta_e(c, previous))),
            )
            colors[subtype] = color
            colors_by_genotype[genotype].append(color)
    return colors


def main():
    args = parse_args()
    subtypes = sorted(set(read_subtypes(args.metadata, args.subtype_columns)) | set(args.extra_subtypes), key=subtype_sort_key)

    genotype_colors = {g: oklch_to_hex(L, C, h) for g, (h, L, C) in GENOTYPE_COLORS.items()}
    genotype_colors["recombinant"] = RECOMBINANT_COLOR
    sub_colors = subtype_colors(subtypes, genotype_colors)

    with open(args.output, "w") as fh:
        for trait in args.genotype_traits:
            for value, color in genotype_colors.items():
                fh.write(f"{trait}\t{value}\t{color}\n")
        for trait in args.subtype_traits:
            for value, color in sub_colors.items():
                fh.write(f"{trait}\t{value}\t{color}\n")

    print(f"Wrote colors for {len(genotype_colors)} genotypes and {len(sub_colors)} subtypes to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
