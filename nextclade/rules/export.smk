"""
This part of the workflow collects the phylogenetic tree and
annotations to export an Auspice JSON.
"""


rule colors:
    """
    Generate colors for genotypes and subtypes, covering all subtypes in the
    (annotated) metadata, such that subtypes within a genotype are distinguishable.
    """
    input:
        metadata="data/metadata.tsv",
    output:
        colors="results/colors.tsv",
    log:
        "logs/colors.txt",
    benchmark:
        "benchmarks/colors.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        ./scripts/generate-colors.py \
            --metadata {input.metadata:q} \
            --subtype-columns subtype ncbi_subtype \
            --output {output.colors:q}
        """


rule string_traits:
    """
    Write categorical metadata columns as node data with string values.
    augur export infers numeric types from metadata (e.g. genotype `1` -> 1),
    but Nextclade expects clade-like attributes to be strings. Node data
    overrides metadata values of the same name in augur export.
    """
    input:
        metadata="results/subsampled_metadata.tsv",
    output:
        node_data="results/string_traits.json",
    params:
        strain_id=config["strain_id_field"],
        columns=config["export"]["string_columns"],
    run:
        import json
        import pandas as pd

        metadata = pd.read_csv(
            input.metadata, sep="\t", dtype=str, keep_default_na=False
        ).set_index(params.strain_id)
        nodes = {
            strain: {col: row[col] for col in params.columns if row[col]}
            for strain, row in metadata.iterrows()
        }
        with open(output.node_data, "w") as fh:
            json.dump({"nodes": nodes}, fh, indent=1)


rule export:
    input:
        tree="results/tree.nwk",
        metadata="results/subsampled_metadata.tsv",
        branch_lengths="results/branch_lengths.json",
        muts="results/muts.json",
        clades=expand("results/clades_{field}.json", field=config["assign_clades"]),
        string_traits="results/string_traits.json",
        auspice_config=config["files"]["auspice_config"],
        colors="results/colors.tsv",
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
            --node-data {input.branch_lengths:q} {input.muts:q} {input.string_traits:q} {input.clades:q} \
            --metadata-columns {params.metadata_columns:q} \
            --auspice-config {input.auspice_config:q} \
            --colors {input.colors:q} \
            --include-root-sequence-inline \
            --output {output.auspice_json:q}
        """
