configfile: "config.yaml"

import pandas as pd

rule all:
	input:
	
rule parse_jcast_results_0:
	input:
		config["jcast_dir"]+"psq_T1.fasta",
		config["jcast_dir"]+"psq_T2.fasta",
		config["jcast_dir"]+"psq_T3.fasta",
	output:
		isoform_csv = "resources/0_parsed_jcast_results/"+config["run_name"]+"_isoform_dataframe.csv"
	script:
		"scripts/parse_jcast_results.py"

checkpoint check_isoform_csv_1:
	input: 
		rules.parse_jcast_results_0.output.isoform_csv 
	output:
		touch(".check_isoform_csv.touch")
		

# Checkpoint code to find all the files and specify all the outputs
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_isoform_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.make_files.get(**w)

        # use glob_wildcards to find the (as-yet-unknown) new files.
        names = glob_wildcards('output-{rs}.txt')[0] # TODO: edit to match use.

        pattern = expand(self.pattern, name=names, **w)
        return pattern

checkpoint split_fastas_1:
	input:
		rules.parse_jcast_results_0.output.isoform_csv
	output:
		touch("split_fastas_1.touch")
	script:
		"split_fastas.py"

rule run_alphafold_2:
	input:

	output:

	shell:

rule parse_af_results_3:
	input:

	output:

	shell:

rule run_gesamt_4:
	input:

	output:

	shell:

rule parse_gesamt_results_5:
	input:

	output:

	shell:

rule align_structures_6:
	input:

	output:

	shell: