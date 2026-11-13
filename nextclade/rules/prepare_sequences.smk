"""
This part of the workflow selects near full length genomes evenly across
genotypes and subtypes, aligns them to the reference and masks gappy
regions of the alignment.

REQUIRED INPUTS:

    metadata    = ../ingest/results/metadata.tsv
    sequences   = ../ingest/results/sequences.fasta

OUTPUTS:

    alignment           = results/aligned.fasta
    translations        = results/translations/{cds}.fasta
    masked alignment    = results/masked.fasta

"""

rule parse_reference_alignment:
    """
    Extract genotype/subtype assignments from the ICTV reference alignment
    (headers like `1a_AF009606` or `1_AJ851228` for unassigned subtypes),
    merge them with the manual annotations and add the reference accessions
    to the include list.
    """
    input:
        alignment=config["files"]["reference_alignment"],
        annotations=config["files"]["genotype_annotations"],
        include=config["files"]["include"],
    output:
        annotations="results/genotype_annotations.tsv",
        include="results/include.txt",
    run:
        import re
        import pandas as pd
        from Bio import SeqIO

        records = []
        for record in SeqIO.parse(input.alignment, "fasta"):
            match = re.fullmatch(r"(\d+)([a-z]*)\??_(\S+)", record.id)
            if match is None:
                print(f"WARNING: skipping unparsable reference header {record.id!r}")
                continue
            genotype, subtype_letters, accession = match.groups()
            records.append(
                {
                    "accession": accession,
                    "genotype": genotype,
                    # genotype-only labels (e.g. `1_AJ851228`) have no subtype
                    "subtype": genotype + subtype_letters if subtype_letters else "",
                }
            )
        reference_df = pd.DataFrame(records, columns=["accession", "genotype", "subtype"])
        reference_df = reference_df.replace("", None)

        # Non-empty manual annotations take precedence over the reference alignment
        manual_df = pd.read_csv(input.annotations, sep="\t", dtype=str).set_index("accession")
        annotations_df = manual_df.combine_first(reference_df.set_index("accession"))
        annotations_df.to_csv(output.annotations, sep="\t", columns=["genotype", "subtype"])

        with open(input.include) as fh:
            include = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        with open(output.include, "w") as fh:
            for accession in dict.fromkeys(include + list(reference_df["accession"])):
                fh.write(f"{accession}\n")


rule annotate_genotypes:
    """
    Add `genotype` and `subtype` columns based on the NCBI annotations,
    augmented/overridden by manual annotations.
    """
    input:
        metadata=config["files"]["metadata"],
        annotations="results/genotype_annotations.tsv",
    output:
        metadata="data/metadata.tsv",
    params:
        strain_id=config["strain_id_field"],
    log:
        "logs/annotate_genotypes.txt",
    benchmark:
        "benchmarks/annotate_genotypes.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk mutate -t -f ncbi_genotype -n genotype {input.metadata:q} \
            | csvtk mutate -t -f ncbi_subtype -n subtype \
            > data/metadata_with_genotypes.tsv

        augur merge \
            --metadata metadata=data/metadata_with_genotypes.tsv annotations={input.annotations:q} \
            --metadata-id-columns {params.strain_id:q} \
            --no-source-columns \
            --output-metadata {output.metadata:q}

        rm data/metadata_with_genotypes.tsv
        """


rule subsample:
    """Select near full length genomes evenly across genotypes and subtypes"""
    input:
        metadata="data/metadata.tsv",
        sequences=config["files"]["sequences"],
        include="results/include.txt",
        exclude=config["files"]["exclude"],
    output:
        sequences="results/subsampled.fasta",
        metadata="results/subsampled_metadata.tsv",
    params:
        strain_id=config["strain_id_field"],
        min_length=config["subsample"]["min_length"],
        query=config["subsample"]["query"],
        group_by=config["subsample"]["group_by"],
        sequences_per_group=config["subsample"]["sequences_per_group"],
        seed=config["subsample"]["seed"],
    log:
        "logs/subsample.txt",
    benchmark:
        "benchmarks/subsample.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur filter \
            --sequences {input.sequences:q} \
            --metadata {input.metadata:q} \
            --metadata-id-columns {params.strain_id:q} \
            --query {params.query:q} \
            --min-length {params.min_length:q} \
            --group-by {params.group_by:q} \
            --sequences-per-group {params.sequences_per_group:q} \
            --subsample-seed {params.seed:q} \
            --include {input.include:q} \
            --exclude {input.exclude:q} \
            --output-sequences {output.sequences:q} \
            --output-metadata {output.metadata:q}
        """


rule align:
    """Align the subsampled genomes to the reference and translate CDSs"""
    input:
        sequences="results/subsampled.fasta",
        reference=config["files"]["reference"],
        annotation=config["files"]["annotation"],
        pathogen_json=config["files"]["pathogen_json"],
    output:
        alignment="results/aligned.fasta",
        nextclade="results/nextclade.tsv",
        translations=expand(
            "results/translations/{cds}.fasta", cds=config["ancestral"]["genes"]
        ),
    params:
        translations=lambda wildcards: "results/translations/{cds}.fasta",
    log:
        "logs/align.txt",
    benchmark:
        "benchmarks/align.txt"
    threads: workflow.cores
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade run \
            --jobs {threads} \
            --input-ref {input.reference:q} \
            --input-annotation {input.annotation:q} \
            --input-pathogen-json {input.pathogen_json:q} \
            --output-fasta {output.alignment:q} \
            --output-translations {params.translations:q} \
            --output-tsv {output.nextclade:q} \
            {input.sequences:q}
        """


rule mask_bed:
    """
    Write the gappy regions to mask (`placementMaskRanges` in the
    pathogen.json, 0-based, end-exclusive) as a BED file
    """
    input:
        pathogen_json=config["files"]["pathogen_json"],
        reference=config["files"]["reference"],
    output:
        bed="results/mask.bed",
    run:
        import json
        from Bio import SeqIO

        seq_id = SeqIO.read(input.reference, "fasta").id
        with open(input.pathogen_json) as fh:
            ranges = json.load(fh).get("placementMaskRanges", [])
        with open(output.bed, "w") as fh:
            for r in ranges:
                fh.write(f"{seq_id}\t{r['begin']}\t{r['end']}\n")


rule mask:
    """Mask gappy regions of the alignment before tree building"""
    input:
        alignment="results/aligned.fasta",
        bed="results/mask.bed",
    output:
        alignment="results/masked.fasta",
    log:
        "logs/mask.txt",
    benchmark:
        "benchmarks/mask.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        if [ -s {input.bed:q} ]; then
            augur mask \
                --sequences {input.alignment:q} \
                --mask {input.bed:q} \
                --output {output.alignment:q}
        else
            echo "No mask ranges configured; copying alignment unmasked."
            cp {input.alignment:q} {output.alignment:q}
        fi
        """
