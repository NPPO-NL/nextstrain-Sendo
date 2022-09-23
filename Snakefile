rule all:
    input:
        auspice = "auspice/Sendo_20220801.json"

input_fasta = "data/Sendo_20220801.fa",
input_metadata = "data/Sendo_meta_20220801.tsv",
reference = "config/mtDNA_MB42.gb",
colors = "config/colors.tsv",
lat_longs = "config/lat_longs.tsv",
description = "config/description.md",
auspice_config = "config/auspice_config.json"


rule align:
    message:
        """
        Aligning sequences to {input}
          - filling gaps with N
        """
    input:
        sequences = input_fasta
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            --method mafft \
            --nthreads auto
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --method raxml \
            --output {output.tree}\
            --substitution-model GTR \
            --nthreads auto
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = input_metadata
    output:
        tree = "results/refined_tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        clock = "0.00000000151",
        date_inference = "marginal",
        root = "SeTr-135-001"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --clock-rate {params.clock} \
            --date-confidence \
            --date-inference {params.date_inference}\
            --root {params.root}
        """

rule traits:
    message: "Infer ancestral traits from an existing phylogenetic tree and the metadata annotating each tip of the tree"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        traits = "results/traits.json"
    params:
        columns = "continent country pathotype"
    shell:
        """
        augur traits \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --output {output.traits} \
        --columns {params.columns}\
        --confidence
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = "results/nt_muts.json",
        sequences = "results/nucleotide-mutations.fasta"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --output-sequences {output.sequences}\
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule export:
    message: "Exporting data files fo auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.traits,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = colors,
        lat_longs = lat_longs,
        description = description,
        auspice_config = auspice_config
    output:
        auspice = rules.all.input.auspice
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --lat-longs {input.lat_longs} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --output {output.auspice}
        """
