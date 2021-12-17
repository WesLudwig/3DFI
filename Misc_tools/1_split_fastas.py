"""Saves all rows of a Pandas DataFrame containing splice isoforms identified by JCAST as single-sequence .fasta files.

Usage: 
	> python split_fastas.py /path/to/subset_of_all_isoforms.csv path/to/output/dir/
"""
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(input_file, out_dir):
	"""Reads the input.csv into a Pandas DataFrame and saves each row to a single .fasta file.

	Args:
		input_file (string): Path/to/file.csv Pandas DataFrame containing splice isoform sequences to be saved as single-sequence .fasta files.
		out_dir (string): Path/to/directory/ where single-sequence .fasta files will be saved."

	Returns:
		None

	"""
	df = pd.read_csv(input_file)
	# Creates a single .fasta for each isoform.
	for i in range(len(df)):
	    subdf = df.iloc[i]
	    with open(out_dir+"_".join([str(val) for val in list(df.loc[i][["accession_number", "prot_name"]].values)+[i]])+".fasta", "w") as out:
	        SeqIO.write(
	            SeqRecord(
	                Seq(subdf["prot_sequence"]),
	                id="|".join([str(val) for val in subdf.values[1:-1]])
	            ),
	            out,
	            "fasta"
	        )

if __name__ == "__main__":

	# If the snakemake global object is present, save expected arguments from snakemake to be passed to main().
	if "snakemake" in globals():
		isoform_csv = snakemake.input.isoform_csv
		out_dir = snakemake.input.out_dir

		main(
			input_file = isoform_csv,
			out_dir = out_dir,
		)

	# CLI context, set expected arguments with argparse module.
	else:
		parser = argparse.Parser(description="Takes a Pandas DataFrame with many protein sequences from parse_jcast_results.py and saves all rows to separate .fasta files.")
		parser.add_argument("input_file", help="Path/to/file.csv Pandas DataFrame containing splice isoform sequences to be saved as single-sequence .fasta files.")
		parser.add_argument("out_dir", help="Path/to/directory/ where single-sequence .fasta files will be saved.")
		args = parser.parse_args()

		main(
			input_file = args.input_file,
			out_dir = args.out_dir,
		)