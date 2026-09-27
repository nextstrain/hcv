"""
This part of the workflow reconstructs ancestral nucleotide and amino acid
sequences on the tree.

The unmasked alignment is used here, such that the tree contains all
mutations relative to the reference.
"""


rule ancestral:
    input:
        tree="results/tree.nwk",
        alignment="results/aligned.fasta",
        translations=expand(
            "results/translations/{cds}.fasta", cds=config["ancestral"]["genes"]
        ),
        reference=config["files"]["reference"],
        annotation=config["files"]["annotation"],
    output:
        node_data="results/muts.json",
    params:
        inference=config["ancestral"]["inference"],
        genes=config["ancestral"]["genes"],
        translations="results/translations/%GENE.fasta",
    log:
        "logs/ancestral.txt",
    benchmark:
        "benchmarks/ancestral.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur ancestral \
            --tree {input.tree:q} \
            --alignment {input.alignment:q} \
            --root-sequence {input.reference:q} \
            --annotation {input.annotation:q} \
            --genes {params.genes:q} \
            --translations {params.translations:q} \
            --inference {params.inference:q} \
            --output-node-data {output.node_data:q}
        """


rule assign_clades:
    """
    Label the MRCA of each genotype/subtype in the tree and assign the label
    to all downstream nodes.
    """
    input:
        tree="results/tree.nwk",
        metadata="results/subsampled_metadata.tsv",
    output:
        node_data="results/clades_{field}.json",
    params:
        strain_id=config["strain_id_field"],
        attribute=lambda w: config["assign_clades"][w.field]["attribute"],
        min_count=lambda w: config["assign_clades"][w.field].get("min_count", 1),
    log:
        "logs/assign_clades_{field}.txt",
    benchmark:
        "benchmarks/assign_clades_{field}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        ./scripts/assign-clades.py \
            --tree {input.tree:q} \
            --metadata {input.metadata:q} \
            --metadata-id-column {params.strain_id:q} \
            --field {wildcards.field:q} \
            --attribute {params.attribute:q} \
            --min-count {params.min_count:q} \
            --output {output.node_data:q}
        """
