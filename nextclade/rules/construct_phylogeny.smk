"""
This part of the workflow constructs the phylogenetic tree from the
masked alignment.
"""


rule tree:
    input:
        alignment="results/masked.fasta",
    output:
        tree="results/tree_raw.nwk",
    params:
        tree_builder_args=config["tree"]["tree_builder_args"],
    log:
        "logs/tree.txt",
    benchmark:
        "benchmarks/tree.txt"
    threads: workflow.cores
    shell:
        r"""
        exec &> >(tee {log:q})

        augur tree \
            --alignment {input.alignment:q} \
            --tree-builder-args={params.tree_builder_args:q} \
            --nthreads {threads} \
            --output {output.tree:q}
        """


rule refine:
    input:
        alignment="results/masked.fasta",
        tree="results/tree_raw.nwk",
    output:
        node_data="results/branch_lengths.json",
        tree="results/tree.nwk",
    log:
        "logs/refine.txt",
    benchmark:
        "benchmarks/refine.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur refine \
            --tree {input.tree:q} \
            --alignment {input.alignment:q} \
            --output-tree {output.tree:q} \
            --output-node-data {output.node_data:q} \
            --root mid_point
        """
