"""
This part of the workflow collects the phylogenetic tree and
annotations to export an Auspice JSON.
"""


rule export:
    input:
        tree="results/tree.nwk",
        metadata="results/subsampled_metadata.tsv",
        branch_lengths="results/branch_lengths.json",
        muts="results/muts.json",
        clades=expand("results/clades_{field}.json", field=config["assign_clades"]),
        auspice_config=config["files"]["auspice_config"],
        colors=config["files"]["colors"],
    output:
        auspice_json="auspice/tree.json",
    params:
        strain_id=config["strain_id_field"],
        metadata_columns=config["export"]["metadata_columns"],
    log:
        "logs/export.txt",
    benchmark:
        "benchmarks/export.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur export v2 \
            --tree {input.tree:q} \
            --metadata {input.metadata:q} \
            --metadata-id-columns {params.strain_id:q} \
            --node-data {input.branch_lengths:q} {input.muts:q} {input.clades:q} \
            --metadata-columns {params.metadata_columns:q} \
            --auspice-config {input.auspice_config:q} \
            --colors {input.colors:q} \
            --include-root-sequence-inline \
            --output {output.auspice_json:q}
        """
