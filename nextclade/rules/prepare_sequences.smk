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


rule annotate_genotypes:
    """
    Add `genotype` and `subtype` columns based on the NCBI annotations,
    augmented/overridden by manual annotations.
    """
    input:
        metadata=config["files"]["metadata"],
        annotations=config["files"]["genotype_annotations"],
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
        include=config["files"]["include"],
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
    """Write the gappy regions to mask from the config as a BED file"""
    output:
        bed="results/mask.bed",
    params:
        ranges=config["mask_ranges"],
        seq_id="NC_038882.1",
    run:
        with open(output.bed, "w") as fh:
            for r in params.ranges:
                fh.write(f"{params.seq_id}\t{r['begin']}\t{r['end']}\n")


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
