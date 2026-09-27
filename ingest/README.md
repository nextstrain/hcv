# Ingest workflow for hepatitis C virus

This workflow ingests public data from NCBI and outputs curated
metadata and sequences that can be used as input for the HCV
`nextclade` and `phylogenetic` workflows.

If you have another data source or private data that needs to be
formatted for the phylogenetic workflow, then you can use a similar
workflow to curate your own data.

## Usage

From within the `ingest` directory:

```
nextstrain build .
```

or, with a local Snakemake installation, `snakemake --cores 4`.

This fetches all records for the HCV species *Orthohepacivirus hominis*
([taxid 3052230][]) via [NCBI Datasets][] and produces

- `results/sequences.fasta`
- `results/metadata.tsv`

To run on the small bundled example dataset instead of fetching from NCBI:

```
nextstrain build . --configfile build-configs/ci/config.yaml
```

## Genotypes

NCBI assigns some records to taxa below the species level (genotypes,
subtypes, or named isolates nested within a subtype). The script
[scripts/extract-ncbi-genotypes.py][] walks the taxonomic lineage in the
NCBI Datasets data report and outputs the `ncbi_genotype` and
`ncbi_subtype` columns. Recombinants (e.g. subtype `2k/1b`) and lab
chimeras are given the genotype `recombinant`. Note that only a minority
of records (~20%) are assigned to such sub-taxa; the remaining records
have empty genotype/subtype values.

## Config

The config directory contains all of the default configurations for
the ingest workflow.

[defaults/config.yaml][] contains all of the default configuration
parameters used for the ingest workflow. Use Snakemake's
`--configfile`/`--config` options to override these default values.

## Snakefile and rules

The rules directory contains separate Snakefiles (`*.smk`) as modules
of the core ingest workflow. The modules of the workflow are in
separate files to keep the main ingest [Snakefile][] succinct and
organized. Modules are all [included][] in the main Snakefile in the
order that they are expected to run.

## Vendored

This repository uses [`git subrepo`][] to manage copies of shared
scripts in [shared/vendored][], from [nextstrain/shared][].

See [shared/vendored/README.md][] for instructions on how to update the
vendored scripts.

[defaults/config.yaml]: ./defaults/config.yaml
[`git subrepo`]: https://github.com/ingydotnet/git-subrepo
[included]: https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes
[NCBI Datasets]: https://www.ncbi.nlm.nih.gov/datasets/
[nextstrain/shared]: https://github.com/nextstrain/shared
[scripts/extract-ncbi-genotypes.py]: ./scripts/extract-ncbi-genotypes.py
[Snakefile]: ./Snakefile
[shared/vendored]: ../shared/vendored
[shared/vendored/README.md]: ../shared/vendored/README.md#vendoring
[taxid 3052230]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3052230
