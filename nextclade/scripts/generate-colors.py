#! /usr/bin/env python
"""
Generate an Auspice colors TSV for HCV genotypes and subtypes.

Each genotype is assigned a distinct base hue. Subtypes of a genotype are
picked from a region of OKLCH space (perceptually uniform lightness/chroma/hue)
around the genotype color: hues between the genotype's hue and part of the way
to its neighbouring genotype hues (more for genotypes with more subtypes), and
a range of lightness and chroma. Colors are
picked greedily such that each new subtype color is as far as possible from the
colors already picked for the genotype (farthest point sampling), starting with
the genotype color itself, while keeping a minimal distance to the colors of
other genotypes. The first subtypes (e.g. 1a, 1b, ...) thus get the most
distinct colors. Genotypes pick their colors in turns, so adding a subtype
can shift the colors of subtypes of other genotypes that sort after it.

Subtypes are read from the given metadata column(s), so the scale covers all
subtypes present in the data.

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

# Subtype colors are picked within +/- HUE_BRACKET (at most) around the genotype
# hue. The hue gap to a neighbouring genotype is split in proportion to the
# number of subtypes of the two genotypes (bounded by HUE_SPLIT_BOUNDS),
# leaving a margin of HUE_MARGIN.
HUE_BRACKET = 40
HUE_SPLIT_BOUNDS = (0.3, 0.7)
HUE_MARGIN = 3
SUBTYPE_LIGHTNESS_RANGE = (0.35, 0.85)
SUBTYPE_CHROMA_LEVELS = (0.07, 0.10, 0.13, 0.16)
# Minimal OKLab distance (x100) of a subtype color to the colors of other genotypes
MIN_DELTA_E_OTHER_GENOTYPE = 6

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


def hue_range(genotype, counts):
    """
    Hue range (degrees, relative to the genotype hue) for the subtypes of a
    genotype: up to HUE_BRACKET on either side, and up to a fraction of the hue
    gap to the nearest genotype on that side, proportional to the number of
    subtypes (`counts`) of the two genotypes, minus HUE_MARGIN.
    """
    hue = GENOTYPE_COLORS[genotype][0]
    limits = []
    for direction in (-1, 1):
        neighbour, gap = min(
            ((g, (direction * (other[0] - hue)) % 360) for g, other in GENOTYPE_COLORS.items() if g != genotype),
            key=lambda x: x[1],
        )
        total = counts.get(genotype, 0) + counts.get(neighbour, 0)
        fraction = counts.get(genotype, 0) / total if total else 0.5
        fraction = min(max(fraction, HUE_SPLIT_BOUNDS[0]), HUE_SPLIT_BOUNDS[1])
        limits.append(min(HUE_BRACKET, fraction * gap - HUE_MARGIN))
    return -limits[0], limits[1]


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


def candidate_colors(genotype, counts):
    """Grid of in-gamut candidate colors (hex, OKLab) for the subtypes of a genotype."""
    hue = GENOTYPE_COLORS[genotype][0]
    lower, upper = hue_range(genotype, counts)
    candidates = {}
    for dh in np.arange(lower, upper + 1e-9, 2):
        for L in np.arange(SUBTYPE_LIGHTNESS_RANGE[0], SUBTYPE_LIGHTNESS_RANGE[1] + 1e-9, 0.025):
            for C in SUBTYPE_CHROMA_LEVELS:
                color = oklch_to_hex(L, C, (hue + dh) % 360)
                candidates[color] = hex_to_oklab(color)
    return candidates


def subtype_colors(subtypes, genotype_colors):
    """
    Assign subtype colors by farthest point sampling: the first subtype of a
    genotype gets the genotype color, each further subtype gets the candidate
    color that is farthest from the colors already assigned to the genotype,
    excluding candidates closer than MIN_DELTA_E_OTHER_GENOTYPE to the colors
    assigned to other genotypes. Genotypes pick in turns (round robin), such
    that genotypes with many subtypes aren't crowded out by those processed
    before them.
    """
    members = {
        g: [s for s in subtypes if re.match(r"^\d+", s).group() == g]
        for g in GENOTYPE_COLORS
    }
    counts = {g: len(m) for g, m in members.items()}
    colors = {m[0]: genotype_colors[g] for g, m in members.items() if m}
    assigned = [(g, hex_to_oklab(c)) for g, c in genotype_colors.items()]
    state = {}
    for g in GENOTYPE_COLORS:
        hexes, labs = zip(*candidate_colors(g, counts).items())
        labs = np.array(labs)
        # distance of each candidate to the nearest color assigned to this genotype
        nearest = 100 * np.linalg.norm(labs - hex_to_oklab(genotype_colors[g]), axis=1)
        state[g] = (list(hexes), labs, nearest)

    for k in range(1, max(len(m) for m in members.values())):
        for g, m in members.items():
            if k >= len(m):
                continue
            hexes, labs, nearest = state[g]
            others = np.array([lab for og, lab in assigned if og != g])
            to_others = 100 * np.min(np.linalg.norm(labs[:, None, :] - others[None, :, :], axis=2), axis=1)
            score = np.where(to_others >= MIN_DELTA_E_OTHER_GENOTYPE, nearest, -np.inf)
            i = int(np.argmax(score))
            colors[m[k]] = hexes[i]
            assigned.append((g, labs[i]))
            state[g] = (hexes, labs, np.minimum(nearest, 100 * np.linalg.norm(labs - labs[i], axis=1)))
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
