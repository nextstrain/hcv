"""
This part of the workflow assembles the Nextclade dataset and tests it on
the example sequences.

OUTPUTS:

    dataset     = dataset/
    test output = test_output/

"""


rule minimizer_references:
    """
    Extract the minimizer references (config `minimizer_references`) used to
    build the dataset's entry in the minimizer index for `nextclade sort`.
    """
    input:
        metadata="data/metadata.tsv",
        sequences=config["files"]["sequences"],
    output:
        include="results/minimizer_references.txt",
        sequences="results/additional_references.fasta",
    params:
        strain_id=config["strain_id_field"],
        references=config["minimizer_references"],
    log:
        "logs/minimizer_references.txt",
    benchmark:
        "benchmarks/minimizer_references.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        printf '%s\n' {params.references:q} > {output.include:q}

        augur filter \
            --sequences {input.sequences:q} \
            --metadata {input.metadata:q} \
            --metadata-id-columns {params.strain_id:q} \
            --exclude-all \
            --include {output.include:q} \
            --output-sequences {output.sequences:q}

        n_expected=$(wc -l < {output.include:q})
        n_found=$(grep -c '>' {output.sequences:q})
        if [ "$n_found" -ne "$n_expected" ]; then
            echo "ERROR: found $n_found of $n_expected minimizer references" >&2
            exit 1
        fi
        """


rule example_sequences:
    """
    Select example sequences for the dataset: near full length genomes
    evenly across genotypes and subtypes that are not part of the tree.
    """
    input:
        metadata="data/metadata.tsv",
        sequences=config["files"]["sequences"],
        tree_metadata="results/subsampled_metadata.tsv",
    output:
        exclude="results/examples_exclude.txt",
        sequences="results/example_sequences.fasta",
    params:
        strain_id=config["strain_id_field"],
        min_length=config["examples"]["min_length"],
        query=config["examples"]["query"],
        group_by=config["examples"]["group_by"],
        sequences_per_group=config["examples"]["sequences_per_group"],
        seed=config["examples"]["seed"],
    log:
        "logs/example_sequences.txt",
    benchmark:
        "benchmarks/example_sequences.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.strain_id:q} {input.tree_metadata:q} \
          | csvtk del-header -t \
          > {output.exclude:q}

        augur filter \
            --sequences {input.sequences:q} \
            --metadata {input.metadata:q} \
            --metadata-id-columns {params.strain_id:q} \
            --query {params.query:q} \
            --min-length {params.min_length:q} \
            --exclude {output.exclude:q} \
            --group-by {params.group_by:q} \
            --sequences-per-group {params.sequences_per_group:q} \
            --subsample-seed {params.seed:q} \
            --output-sequences {output.sequences:q}
        """


rule assemble_dataset:
    input:
        pathogen_json=config["files"]["pathogen_json"],
        reference=config["files"]["reference"],
        annotation=config["files"]["annotation"],
        tree="auspice/tree.json",
        readme=config["files"]["dataset_readme"],
        changelog=config["files"]["dataset_changelog"],
        examples="results/example_sequences.fasta",
        minimizer_references="results/additional_references.fasta",
    output:
        pathogen_json="dataset/pathogen.json",
        reference="dataset/reference.fasta",
        annotation="dataset/annotation.gff",
        tree="dataset/tree.json",
        readme="dataset/README.md",
        changelog="dataset/CHANGELOG.md",
        examples="dataset/sequences.fasta",
        minimizer_references="dataset/additional_references.fasta",
    log:
        "logs/assemble_dataset.txt",
    benchmark:
        "benchmarks/assemble_dataset.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        cp -v {input.pathogen_json:q} {output.pathogen_json:q}
        cp -v {input.reference:q} {output.reference:q}
        cp -v {input.annotation:q} {output.annotation:q}
        cp -v {input.tree:q} {output.tree:q}
        cp -v {input.readme:q} {output.readme:q}
        cp -v {input.changelog:q} {output.changelog:q}
        cp -v {input.examples:q} {output.examples:q}
        cp -v {input.minimizer_references:q} {output.minimizer_references:q}
        """


rule test_dataset:
    """Run the assembled dataset on the example sequences"""
    input:
        dataset=rules.assemble_dataset.output,
        sequences="dataset/sequences.fasta",
    output:
        outdir=directory("test_output"),
    params:
        dataset_dir="dataset",
    log:
        "logs/test_dataset.txt",
    benchmark:
        "benchmarks/test_dataset.txt"
    threads: workflow.cores
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade run \
            --jobs {threads} \
            --input-dataset {params.dataset_dir:q} \
            --output-all {output.outdir:q} \
            {input.sequences:q}
        """
