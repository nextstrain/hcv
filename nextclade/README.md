# Nextclade dataset workflow for hepatitis C virus

This workflow builds a reference tree spanning all HCV genotypes and
subtypes, to be used for a Nextclade dataset aligned to the H77 reference
(NC_038882.1).

## Inputs

The curated output of the [ingest](../ingest) workflow
(`../ingest/results/{metadata.tsv,sequences.fasta}`), see
[defaults/config.yml][].

Genotype and subtype are taken from the NCBI annotations (`ncbi_genotype`,
`ncbi_subtype`) and augmented/overridden with
[defaults/genotype_annotations.tsv][]. Accessions in [defaults/include.txt][]
(reference and minimizer references) are always included.

## Workflow

1. Subsample near full length genomes (`subsample.min_length`) evenly across
   genotypes and subtypes (`augur filter --group-by genotype subtype --sequences-per-group`).
2. Align the subsample with Nextclade and translate the CDSs.
3. Mask the gappy regions listed in `mask_ranges`.
4. Build the tree (`augur tree`, midpoint rooted by `augur refine`), reconstruct
   ancestral sequences and mutations (`augur ancestral`, using the Nextclade
   translations) and export to `auspice/tree.json`.
5. Annotate genotypes and subtypes on the tree (`scripts/assign-clades.py`,
   configured in `assign_clades`): the MRCA of all tips of a genotype/subtype
   gets a branch label, and the value is assigned to all downstream nodes.
   Genotypes are written as `clade_membership` (branch label `clade`), subtypes
   as `subtype` (branch label `subtype`).
6. Assemble the Nextclade dataset in `dataset/` from [defaults/pathogen.json][],
   the reference, annotation, dataset README/CHANGELOG, the tree, example
   sequences (near full length genomes not in the tree, config `examples`) and
   the minimizer references (config `minimizer_references`, written to
   `additional_references.fasta` and listed in `minimizerIndex` of the
   pathogen.json), and test it on the examples (`test_output/`).

Run with

```
nextstrain build .
```

or, with a local Snakemake installation, `snakemake --cores 8`.

## Set-up: finding regions to mask

HCV is very diverse and the alignment against a single reference contains
regions with many gaps. These regions are masked before tree building and
excluded from placement. They are identified once during set-up (this is not
part of the workflow):

```
snakemake --cores 8 results/aligned.fasta
./scripts/find-gappy-regions.py \
    --alignment results/aligned.fasta \
    --max-gap-fraction 0.3 \
    --keep-ranges 0-336 \
    --config defaults/config.yml \
    --pathogen-json defaults/pathogen.json
```

This writes the ranges (0-based, end-exclusive) to `mask_ranges` in
[defaults/config.yml][] and to `placementMaskRanges` in
[defaults/pathogen.json][]. `--keep-ranges 0-336` keeps the 5'UTR in the
placement mask. Run without `--config`/`--pathogen-json` to only print the
ranges.

## Colors

The `colors` rule generates `results/colors.tsv` with
[scripts/generate-colors.py][] from the subtypes in the annotated metadata
(`subtype` from the reference alignment/manual annotations and
`ncbi_subtype`), so every subtype in the tree gets a color.

Each genotype has a distinct base color. Its subtypes are picked from a region
of hue, lightness and chroma around the genotype color by farthest point
sampling, so that subtypes within a genotype are as distinguishable as
possible (the first subtype gets the genotype color, and 1a, 1b, ... get the
most distinct colors). With more than 30 subtypes of genotype 6 in a small
region of the sRGB gamut, some of its subtypes will still look similar.

## Notes

- `minSeedCover` is lowered to 0.1 in [defaults/pathogen.json][], as the seed
  matches of genotypes 2–8 cover only 13–21% of the H77 reference (the
  default is 0.33).
- [defaults/annotation.gff][] has `gene` rows as parents of each CDS, as
  `augur` reads gene names from `gene` rows only.

[defaults/annotation.gff]: ./defaults/annotation.gff
[defaults/config.yml]: ./defaults/config.yml
[scripts/generate-colors.py]: ./scripts/generate-colors.py
[defaults/genotype_annotations.tsv]: ./defaults/genotype_annotations.tsv
[defaults/include.txt]: ./defaults/include.txt
[defaults/pathogen.json]: ./defaults/pathogen.json
