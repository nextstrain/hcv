# HCV Nextclade dataset: run on the full GenBank ingest

Dataset `dataset/` (built 2026-11-13) run with Nextclade 3.16.0 on all 279,709
sequences in `../ingest/results/sequences.fasta`. Calls are compared to NCBI
genotype/subtype annotations and to the ICTV reference alignment.

Files: `nextclade.tsv` (Nextclade output), `comparison.tsv` (merged calls,
annotations and genomic region), `compare.py` (comparison script).

## Coverage

| Sequence length | n | Aligned | Genotype | Subtype |
|---|---:|---:|---:|---:|
| ≤300 nt | 68,654 | 68% | 68% | 63% |
| 300–500 | 87,351 | 96% | 96% | 93% |
| 500–1,000 | 58,997 | 99% | 99% | 98% |
| 1,000–8,000 | 57,450 | 100% | 100% | >99% |
| 8,000–10,000 | 7,183 | 99% | 99% | 95% |

- 25,953 sequences failed to align, mostly short fragments (median 81 nt).
- 74 near-complete genomes failed, 71 of them patent sequences (`FV…`). 23 of
  the 74 align with `--retry-reverse-complement`.
- Some full-length genomes get a genotype but no subtype because they are
  placed above the base of a subtype's clade.

## Agreement with existing annotations

| Comparison | n | Agreement |
|---|---:|---:|
| Genotype vs NCBI | 49,672 | 99.3% (332 disagree) |
| Subtype vs NCBI | 47,585 | 98.3% (791 disagree) |
| Genotype/subtype vs reference alignment | 238 / 184 | 100%* |

\* Not independent: the reference sequences are always included in the tree.

### Main sources of disagreement

- **5′UTR-only sequences:** 29% genotype disagreement (165/560, mostly NCBI
  genotype 1 called 6) and 48% subtype disagreement (115/240, e.g. 1b→6k,
  1a→1o), all passing QC.
- **NCBI "recombinant" called 2a (57):** JFH1-based lab chimeras (Gottwein,
  Scheel, Li et al.), so 2a is expected.
- **NCBI 2b called 1a (68, NS5B):** all from Kageyama et al. (Philippines);
  consistent with the known 2b/1a recombinant, which is 1a in NS5B.
- **1a↔1b (501, mostly E2 and NS5A):** concentrated in Torres-Puente et al.
  (Spain, 298 of the 501), which suggests annotation or mixed-infection
  issues in that study.
- **Other genomic regions:** 0–3% subtype disagreement. Full genomes: 0.3%
  (5/1,661).

## Recommendations

- 5′UTR-only sequences are still placed (masking the UTR would leave them
  unplaced); the dataset README notes that their assignments are unreliable.
- `retryReverseComplement` is now enabled in `alignmentParams` of the
  pathogen.json (rescues 23 of the 74 failed near-complete genomes).
- `mixedSitesThreshold` raised from 20 to 50 (see QC parameter options).

## Quality control

Of the 253,756 aligned sequences, 235,431 (92.8%) are `good`, 6,676 (2.6%)
`mediocre` and 11,649 (4.6%) `bad`. Sequences shorter than 1 kb are almost
always `good`; 15–19% of sequences 1–8 kb long are `bad`. Scripts:
`qc_analysis.py` (writes `qc.tsv`).

| QC rule | Mediocre | Bad | Sequences non-good from this rule alone |
|---|---:|---:|---:|
| Mixed sites (threshold 20) | 12,939 | 9,394 | 14,766 |
| Frameshifts | 3,170 | 1,937 | 1,136 |
| Stop codons | 3,984 | 156 | 470 |
| Frameshifts + stop codons | | | 1,201 |
| Private mutations (typical 500, cutoff 1000) | 616 | 96 | 291 |
| Missing data | 30 | 59 | 54 |

### Mixed sites

The main reason sequences are not `good`. Most are population (Sanger)
sequences from clinical/drug resistance studies (Berger, Bartels, Hezode et
al.) with a median of ~2% ambiguous sites among flagged sequences, which
reflects within-host diversity. Ambiguities don't make genotype and subtype
calls worse: none of the 120 annotated sequences with ≥20 mixed sites
disagree with the GenBank subtype. Because the threshold is absolute, longer
sequences are flagged more often.

### Frameshifts (5,107 sequences)

- 992 sequences of a single submission (SHAW et al., genotype 1, ~4.9 kb)
  share an identical set of ~15 frameshifts and a stop at E2:168. The raw
  sequences have no intact reading frame, so these are real defects of the
  submitted sequences.
- 563 sequences of one study (Li et al., mostly genotype 3) share an
  insertion at position 5048 (`NS3:544-631`) that is absent from other
  sequences, so it's study-specific.
- Most remaining frameshifts in sequences with an otherwise intact reading
  frame lie within 15 nt of the alignment ends (82%). These are errors at
  the ends of Sanger reads or misalignments of the terminal nucleotides.

### Stop codons (4,140 sequences)

- 1,589 stop codons are in F (all genotypes). F is not conserved, so these
  don't indicate poor quality. Ignoring all stops in F would change the
  status of only 223 sequences, because a single stop codon doesn't make a
  sequence non-`good`.
- Outside F, stops in sequences with an intact reading frame are mostly at
  the alignment ends (85%), e.g. a few mismatched terminal nucleotides
  instead of an indel. Some are real, e.g. stop codons appended to
  expression constructs.

### QC parameter options

| Change | Good | Mediocre | Bad |
|---|---:|---:|---:|
| Current | 235,431 | 6,676 | 11,649 |
| Ignore all stops in F | 235,654 | 6,568 | 11,534 |
| + `mixedSitesThreshold` 30 | 239,801 | 5,348 | 8,607 |
| + `mixedSitesThreshold` 50 | 243,653 | 4,862 | 5,241 |
| + `mixedSitesThreshold` 100 | 248,069 | 3,077 | 2,610 |

Frameshift, stop codon and private mutation flags mostly point to real
problems in the submitted sequences, so their parameters look appropriate.
Private mutations in full genomes: median 452, 99th percentile 1,564
(higher for genotypes 5 and 7, which have few relatives in the tree).
