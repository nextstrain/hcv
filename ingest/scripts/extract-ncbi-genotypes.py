#! /usr/bin/env python
"""
Parses the NCBI Datasets virus data report (JSON lines, provided on STDIN) to
extract the HCV genotype and subtype from the NCBI taxonomy assigned to each
record. Outputs TSV to STDOUT.

NCBI assigns some HCV records to taxa below the species level, e.g.
"hepatitis C virus genotype 1a", "Hepatitis C virus subtype 2k/1b",
"Hepatitis C virus genotype 3", or named isolates such as
"Hepatitis C virus (isolate H77)" that are nested within a subtype taxon.
The full taxonomic lineage of each record is included in the data report,
so we walk it from the most specific taxon upwards and take the first
subtype and genotype designations we find.

Intra- and inter-genotypic recombinants (e.g. subtype "2k/1b") and
lab-generated chimeric viruses (taxa named "Recombinant ...") are given
the genotype "recombinant".

Records that are only assigned to the species (the majority) get empty
genotype and subtype values.

Intended to be merged with other TSV-formatted data in subsequent
`augur merge` step in workflow.
"""

import csv
import json
import re
import sys

SUBTYPE_PATTERN = re.compile(
    r"(?:genotype|subtype)\s+([1-8][a-z](?:/[1-8][a-z])*)\b", re.IGNORECASE
)
GENOTYPE_PATTERN = re.compile(r"genotype\s+([1-8])\b", re.IGNORECASE)


def main():
    tsv_writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    tsv_writer.writerow(["accession", "ncbi_genotype", "ncbi_subtype"])

    for line in sys.stdin:
        if not line.strip():
            continue
        report = json.loads(line)
        accession = report["accession"].split(".")[0]
        genotype, subtype = parse_lineage(report.get("virus", {}).get("lineage", []))
        tsv_writer.writerow([accession, genotype, subtype])


def parse_lineage(lineage):
    """
    Given the list of taxa (dicts with a "name" key, ordered from root to
    leaf), return a (genotype, subtype) tuple of strings.
    """
    genotype = ""
    subtype = ""

    for taxon in reversed(lineage):
        name = taxon.get("name", "")

        if name.lower().startswith("recombinant"):
            genotype = "recombinant"

        if not subtype:
            match = SUBTYPE_PATTERN.search(name)
            if match:
                subtype = match.group(1).lower()

        if not genotype:
            match = GENOTYPE_PATTERN.search(name)
            if match:
                genotype = match.group(1)

    if "/" in subtype:
        genotype = "recombinant"
    elif subtype and not genotype:
        genotype = subtype[0]

    return genotype, subtype


if __name__ == "__main__":
    main()
