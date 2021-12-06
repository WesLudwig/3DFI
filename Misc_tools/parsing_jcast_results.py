"""Parses JCAST .fasta results files into a single DataFrame of splice isoform sequences and optionally annotates the subcellular location of the isoform's canonical sequence.

Splice isoforms resulting from JCAST are parsed from the tier 1-3 .fasta files and combined into a single DataFrame. Their canonical forms are queried against the SwissProt 
database (all non-orphan JCAST results should have SP Canonical forms), and the subcellular locations of the canonical forms are added to the isoform's DataFrame row.

Usage:
	> python parsing_jcast_results.py \
	  --input_directory /path/to/dir/containing/jcast.fastas \
	  --output_file /path/to/combined/dataframe.csv \
	  --swissprot_tsv /path/to/swissprot_with_subcellular_locations.tsv (/mnt/databases/Human_SwissProt_Proteome_with_Subcellular_Locations.tsv on Lam Lab GPU workstation) \
	  --plot_file /path/to/isoform_length_distributions_plot.png \
"""
import os
import glob
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def main(input_directory, output_file, swissprot_tsv, plot_file):
	"""Parses JCAST .fasta results files into a single DataFrame of splice isoform sequences and optionally annotates the subcellular location of the isoform's canonical sequence.

	Splice isoforms resulting from JCAST are parsed from the tier 1-3 .fasta files and combined into a single DataFrame. Their canonical forms are queried against the SwissProt 
	database (all non-orphan JCAST results should have SP Canonical forms), and the subcellular locations of the canonical forms are added to the isoform's DataFrame row.

	Args:
		input_directory (string): /path/to/dir/ containing JCAST .fasta files.
		output_file (string): /path/to/file.csv for combined isoform dataframe.
		swissprot_tsv (string): /path/to/file.tsv containing a local copy of the SwissProt database with subcellular locations.
		plot_file (string): /path/to/file.(png, jpeg, tiff, etc.) for plot of isoform length distribution.

	Returns:
		None

	"""
	# Globs tier 1 through 3 .fasta isoform file names from JCAST into a list.
	isoform_fastas = sorted(list(glob.glob(input_directory+"psq_T[1-3].fasta")))
	
	columns = [
	    "database", 
	    "accession_number", 
	    "prot_name", 
	    "ensembl_gene_id", 
	    "rmats_jxn_type", 
	    "input_file_row", 
	    "chromosome", 
	    "anchor_exon_start_end", 
	    "alt_exon_start_end", 
	    "translated_strand_and_phase", 
	    "msjc", 
	    "tier",
	    "prot_sequence",
	]

	main_df = pd.DataFrame(columns=columns)

	i = 0
	for fasta in isoform_fastas:
	    for rec in SeqIO.parse(fasta, "fasta"):
	        main_df.loc[i] = rec.id.split("|") + [str(rec.seq)]
	        i+=1

	# Adds locations of canonical proteoform to splice isoform rows in main_df.
	if swissprot_tsv is not None:
		sp = pd.read_csv(sp_tsv, sep = "\t")
		locations = []
		for i in range(len(main_df)):
		    locs = list(sp.loc[sp["Entry name"]==main_df.iloc[i]["prot_name"]]["Subcellular location [CC]"].values)
		    if len(locs) > 0:
		        if len(locs) > 1:
		            for j in range(len(locs)):
		                if locs[j] != "nan" and locs[j] != "NaN" and not np.isnan(locs[j]):
		                    locations.append(locs[j])
		        else:
		            if isinstance(locs[0], float):
		                locations.append("Not Found")
		            else:
		                locations.append(locs[0])
		    else:
		        locations.append("Not Found")

		# Add locations to output DataFrame.
		main_df["location"]=locations
	
	# Saves the main DataFrame as a .csv.
	if input_dir is not None:
		main_df.to_csv(output_file)

	# Saves a distribution plot of lengths of splice isoforms.
	if plot_file is not None:
		sns.displot([len(seq) for seq in main_df["prot_sequence"].values]).save
		plt.savefig(plot_file)

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="Parses JCAST output .fasta files into a single Pandas Dataframe, optionally annotates JCAST isoforms with subcellular locations.")
	parser.add_argument("input_directory", help="Directory containing the output .fasta files from JCAST.")
	parser.add_argument("-o", "--output_file", help="Path to .csv file to be made from a Pandas Dataframe of all JCAST tier 1-3 isoforms. Must include .csv extension.")
	parser.add_argument("-s", "--swissprot_tsv", help="Path to .tsv file containing the SwissProt database with subcellular locations, /mnt/databases/Human_SwissProt_Proteome_with_Subcellular_Locations.tsv on Lam Lab GPU workstation.")
	parser.add_argument("-p", "--plot_file", help="Path to image file to be made from distribution plot of isoform lengths. Must include extension: .png, .jpeg, etc.")
	args = parser.parse_arguments()

	main(
		input_directory = args.input_directory,
		output_file = args.output_file,
		swissprot_tsv = args.swissprot_tsv,
		plot_file = args.plot_file,
	)