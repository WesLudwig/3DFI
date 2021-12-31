"""
TODO: 
    Write a descriptive filestring
    Get it working + debug
"""

configfile: "/opt/3DFI/config.yaml"
workdir: config["workdir"]

import glob
import pandas as pd

# Checkpoint code to find dynamic files for downstream input.
class Checkpoint_FindFiles:
    def __init__(self, checkpoint_code):
        self.checkpoint_code = checkpoint_code
    def __call__(self, wildcards):
        global checkpoints

        # Checkpoint.get call to trigger exception until checkpoint is run.
        if self.checkpoint_code == 1:
            checkpoints.split_fastas_1.get(**wildcards)
            names = glob_wildcards("resources/1_split_fastas/{name}.fasta")[0]
            return expand("resources/1_split_fastas/{name}.fasta", name=names)

        elif self.checkpoint_code == 3:
            checkpoints.parse_af_results_3.get(**wildcards)
            names = glob_wildcards("resources/3_parsed_af_results/{name}.pdb")[0]
            return expand("resources/3_parsed_af_results/{name}.pdb", name=names)

        elif self.checkpoint_code == 4:
            checkpoints.run_gesamt_4.get(**wildcards)
            names = glob_wildcards("resources/4_gesamt_results/{name}.gesamt")[0]
            return expand("resources/4_gesamt_results/{name}.gesamt", name=names)

        elif self.checkpoint_code == 6:
            checkpoints.split_fastas_1.get(**wildcards)
            names = glob_wildcards("results/5_aligned_structures/{name}/")[0]
            # This scheme does not account for the .csx files that are created upon successful execution.
            return expand("results/5_aligned_structures/{name}/{name}.pdb", name=names)

rule all:
    input:
        "results/checkpoints/checkpoint_6.touch",
        Checkpoint_FindFiles(6)

""" 
This should be performed manually, the pipeline should start with a .csv of a subset of .fastas    
rule parse_jcast_results_0:
    input:
        jcast_dir = config["jcast_dir"],
        swissprot_tsv = "/mnt/databases/Human_SwissProt_Proteome_with_Subcellular_Locations.tsv"
    output:
        isoform_csv = "resources/0_parsed_jcast_results/"+config["run_name"]+"_isoform_dataframe.csv",
        plot_file = "Results/"+config["run_name"]+"_distribution_plot.pdf"
    script:
        config["3DFI_dir"]+"Misc_tools/parse_jcast_results.py"
"""

checkpoint split_fastas_1:
    input:
        isoform_csv = config["isoform_csv"],
        out_dir = "resources/1_split_fastas/",
    output:
        touch("results/checkpoints/checkpoint_1.touch")
    script:
        config["3DFI_dir"]+"Misc_tools/1_split_fastas.py"

checkpoint run_alphafold_2:
    input:
        "results/checkpoints/checkpoint_1.touch",
        fastas = Checkpoint_FindFiles(1),
        out_dir = "resources/2_run_alphafold/",
    output:
        touch("results/checkpoints/checkpoint_2.touch")
    run:
        shell(
            config['3DFI_dir'] + "Prediction/AlphaFold2/2_alphafold.pl"
                + " -f " + " ".join(input.fastas)
                + " -o " + input.out_dir + "/" # Snakemake removes trailing slashes.
        )

checkpoint parse_af_results_3:
    input:
        "results/checkpoints/checkpoint_2.touch",
        af_dir = "resources/2_run_alphafold/",
        out_dir = "resources/3_parsed_af_results/"
    output:
        touch("results/checkpoints/checkpoint_3.touch")
    run:
        shell(
            config['3DFI_dir'] + "Prediction/AlphaFold2/3_parse_af_results.pl"
                + " -a " + input.af_dir + "/"
                + " -o " + input.out_dir + "/" # Snakemake removes trailing slashes.
        )

checkpoint run_gesamt_4:
    input:
        "results/checkpoints/checkpoint_3.touch",
        pdbs = Checkpoint_FindFiles(3),
        out_dir = "resources/4_gesamt_results/"
    output:
        touch("results/checkpoints/checkpoint_4.touch")
    run:
        shell(
            config['3DFI_dir'] + "Homology_search/4_run_gesamt.pl"
                + " -c 28"
                + " -a " + config['gesamt_archive']
                + " -i " + " ".join(input.pdbs)
                + " -o " + input.out_dir + "/"
                + " --query"
        )

rule parse_gesamt_results_5:
    input:
        "results/checkpoints/checkpoint_4.touch",
        gesamt_files = Checkpoint_FindFiles(4)
    output:
        tsv = "resources/5_parsed_gesamt_results/"+config["run_name"]+"_parsed_gesamt_matches.tsv"
    script:
        config["3DFI_dir"]+"Homology_search/descriptive_GESAMT_matches.py"

checkpoint align_structures_6:
    priority: 1
    input:
        gesamt_tsv = rules.parse_gesamt_results_5.output.tsv,
        pdb_dir = "resources/3_parsed_af_results/",
        rcsb_dir = "/mnt/Databases/RCSB_PDB/",
        out_dir = "results/5_aligned_structures/",
    output:
        touch("results/checkpoints/checkpoint_6.touch")
    run:
        shell(
            config['3DFI_dir'] + "Visualization/prepare_visualizations.pl"
            + " -g " + input.gesamt_tsv
            + " -p " + input.pdb_dir + "/"
            + " -r " + input.rcsb_dir + "/"
            + " -o " +input.out_dir + "/" # Snakemake removes trailing slashes.
        )