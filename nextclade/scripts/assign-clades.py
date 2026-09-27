#! /usr/bin/env python
"""
Assign clade-like annotations (e.g. genotype or subtype) to all nodes of a
tree, based on the annotations of the tips in the metadata, and write them as
an augur node data JSON.

For each value of the metadata field (e.g. subtype "1a"), the most recent
common ancestor (MRCA) of all tips with that value is identified. The branch
leading to the MRCA is labeled with the value and the MRCA and all nodes
downstream of it are assigned the value. Nested groups (e.g. a subtype within
a genotype, or a group whose MRCA lies within another group) override the
assignment of the enclosing group, since assignments are made in preorder.

Tips that have an annotation in the metadata keep it, even when it differs
from the clade they fall in (reported as a conflict). Tips without annotation
are assigned the value of the clade they fall in.

Output structure, for attribute `subtype`:

    {
      "nodes":    {"<node name>": {"subtype": "1a"}, ...},
      "branches": {"<MRCA name>": {"labels": {"subtype": "1a"}}, ...}
    }

For the standard clade attribute `clade_membership`, the branch label key is
`clade`, as expected by augur/auspice/nextclade.
"""

import argparse
import json
import sys
from collections import defaultdict

import pandas as pd
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", required=True, help="Newick tree file")
    parser.add_argument("--metadata", required=True, help="Metadata TSV with the annotations of the tips")
    parser.add_argument("--metadata-id-column", default="accession", help="Column with the tip names (default: %(default)s)")
    parser.add_argument("--field", required=True, help="Metadata column with the annotation, e.g. genotype or subtype")
    parser.add_argument("--attribute",
        help="Name of the node attribute to write (default: same as --field). "
             "Use `clade_membership` for the standard clade attribute.")
    parser.add_argument("--label-key",
        help="Key of the branch label (default: `clade` for attribute `clade_membership`, otherwise the attribute name)")
    parser.add_argument("--min-count", type=int, default=1,
        help="Ignore values with fewer than this many tips (default: %(default)s)")
    parser.add_argument("--ignore-values", nargs="*", default=["recombinant"],
        help="Values of the field that do not define clades (default: %(default)s)")
    parser.add_argument("--unassigned-value", default=None,
        help="Value for nodes that are not within any clade (default: no attribute is written for such nodes)")
    parser.add_argument("--output", required=True, help="Output node data JSON")
    return parser.parse_args()


def main():
    args = parse_args()
    attribute = args.attribute or args.field
    label_key = args.label_key or ("clade" if attribute == "clade_membership" else attribute)

    tree = Phylo.read(args.tree, "newick")
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str, keep_default_na=False, usecols=[args.metadata_id_column, args.field])
    tip_annotations = {
        name: value.strip()
        for name, value in zip(metadata[args.metadata_id_column], metadata[args.field])
        if value.strip() and value.strip() not in args.ignore_values
    }

    tips_by_value = defaultdict(list)
    for tip in tree.get_terminals():
        if tip.name in tip_annotations:
            tips_by_value[tip_annotations[tip.name]].append(tip)

    # Identify the MRCA of each group
    mrca_by_value = {}
    for value, tips in sorted(tips_by_value.items()):
        if len(tips) < args.min_count:
            print(f"Ignoring {args.field} {value!r}: only {len(tips)} tip(s)", file=sys.stderr)
            continue
        mrca = tree.common_ancestor(tips)
        clade_tips = [t.name for t in mrca.get_terminals()]
        intruders = [t for t in clade_tips if t in tip_annotations and tip_annotations[t] != value]
        unannotated = [t for t in clade_tips if t not in tip_annotations]
        if intruders:
            print(f"{args.field} {value!r} is not monophyletic: its MRCA also contains "
                  f"{len(intruders)} tip(s) with other annotations: "
                  + ", ".join(f"{t}({tip_annotations[t]})" for t in intruders[:10])
                  + (" ..." if len(intruders) > 10 else ""), file=sys.stderr)
        if unannotated:
            print(f"{args.field} {value!r} is assigned to {len(unannotated)} unannotated tip(s): "
                  + ", ".join(unannotated[:10]) + (" ..." if len(unannotated) > 10 else ""), file=sys.stderr)
        mrca_by_value[value] = mrca

    # Resolve groups that share the same MRCA: keep the larger group
    value_by_mrca = {}
    for value, mrca in mrca_by_value.items():
        other = value_by_mrca.get(mrca)
        if other is not None:
            keep, drop = sorted([value, other], key=lambda v: -len(tips_by_value[v]))
            print(f"{args.field} {value!r} and {other!r} share the same MRCA; labeling it as {keep!r}", file=sys.stderr)
            value_by_mrca[mrca] = keep
        else:
            value_by_mrca[mrca] = value

    # Assign values in preorder, such that nested groups override enclosing ones
    assignment = {}
    for node in tree.find_clades(order="preorder"):
        if node in value_by_mrca:
            for descendant in node.find_clades():
                assignment[descendant] = value_by_mrca[node]

    node_data = {"nodes": {}, "branches": {}}
    n_conflicts = 0
    for node in tree.find_clades(order="preorder"):
        value = assignment.get(node, args.unassigned_value)
        if node.is_terminal() and node.name in tip_annotations:
            if value != tip_annotations[node.name]:
                print(f"Tip {node.name} has {args.field} {tip_annotations[node.name]!r} "
                      f"but falls in clade {value!r}; keeping {tip_annotations[node.name]!r}", file=sys.stderr)
                n_conflicts += 1
            value = tip_annotations[node.name]
        if value is not None:
            node_data["nodes"][node.name] = {attribute: value}
        if node in value_by_mrca:
            node_data["branches"][node.name] = {"labels": {label_key: value_by_mrca[node]}}

    print(f"Labeled {len(node_data['branches'])} clades by {args.field}; "
          f"{n_conflicts} tip(s) conflict with their clade", file=sys.stderr)

    with open(args.output, "w") as fh:
        json.dump(node_data, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
