import sys
import glob
import pypdb
import argparse
import pandas as pd

def main(input_directory=None, input_files=None, dataframe_output=None, tdfi_output, n_matches=5):

    # Columns expected for .gesamt files.
    columns = [
     'Hit_No.',
     'PDB_code',
     'Chain_Id',
     'Q-score',
     'r.m.s.d',
     'Seq._Id.',
     'Nalign',
     'nRes',
     'File_name'
    ]

    # Handles required inputs.
    if (input_directory is None) and (input_files is None):
        sys.exit("Missing required argument: input files or directory.")
    
    # Snakemake context.
    elif input_files is not None:
        gesamt_results = input_files
    
    # CLI context.
    elif input_directory is not None:
        gesamt_results = glob.glob(input_directory+"*.gesamt")

    main_df = pd.DataFrame(columns=columns+["parent_file", "function"])
    for fn in gesamt_results:
        
        # Opens .gesamt file and removes column header lines.
        with open(fn) as file:
            my_name = fn.split("/")[-1].split(".")[0]
            l1 = file.readline().strip().split()[1:]
            l2 = file.readline().strip().split()[1:]
            
            # Creates the dict that will be converted to a DataFrame for this file.
            df_source = dict.fromkeys(columns+["parent_file", "function"])
            for column in df_source:
                df_source[column] = []
            
            # Iterates over the first n_matches of the file, parses lines into dict.
            for i in range(n_matches):
                line = file.readline().strip().split()
                for j in range(len(line)):
                    df_source[columns[j]].append(line[j])
                df_source["parent_file"].append(my_name)

                # Searches PDB for match code and adds matching protein's function to dict.
                struct_info = pypdb.describe_pdb(df_source["PDB_code"][-1])["struct"]
                if len(struct_info.keys()) == 1:
                    df_source["function"].append(struct_info[list(struct_info.keys())[0]])
                else:
                    if "pdbx_descriptor" in struct_info.keys():
                        df_source["function"].append(struct_info["pdbx_descriptor"])
                    elif "title" in struct_info.keys():
                        df_source["function"].append(struct_info["title"])
                    else:
                        df_source["function"].append(struct_info[list(struct_info.keys())[-1]])
                        
            # Converts dict to DataFrame and concatenates to the master DF.
            df_in = pd.DataFrame(df_source)
            main_df = pd.concat([main_df, df_in])
    
    if dataframe_output is not None:
        main_df.to_csv(output_file)
    
    if tdfi_output is not None:
        with open(tdfi_output, "w+") as file:
            last_parent = ""
            for i in range(len(main_df)):
                my_line = main_df.iloc[i]
                my_parent = my_line["parent_file"]
                if my_parent != last_parent:
                    file.write("### " + my_parent + "; Query mode = normal\n")
                    file.write("\t".join([my_parent]+list(my_line[columns+["function"]].values)) + "\n")
                    last_parent = my_parent
                else:
                    file.write("\t".join([my_parent]+list(my_line[columns+["function"]].values)) + "\n")

if __name__ == "__main__":

    # If the snakemake global object is present, save expected arguments from snakemake to be passed to main().
    if "snakemake" in globals():
        input_files = snakemake.input.gesamt_files
        output_tsv = snakemake.output.to_csv

        main(
            input_files = input_files,
            tdfi_output = output_tsv,
        )

    # CLI context, set expected arguments with argparse module.
    else:
        parser = argparse.ArgumentParser(description="Combines the top n_matches for each result of run_GESAMT.pl with matching proteins' functions sourced from the RCSB PDB.")
        
        parser.add_argument("-i", "--input_directory", required=True, help="/path/to/directory/ containing .gesamt files to combine and annotate.")
        parser.add_argument("--dataframe_output", help="/path/to/DataFrame.csv, the combined DataFrame saved as a .csv file.")
        parser.add_argument("--tdfi_output", help="/path/to/descriptive_gesamt_matches.tsv, formatted to work with downstream structural alignment in 3DFI.")
        parser.add_argument("--n_matches", default=5, type=int, help="Number of top matches to include for each .gesamt file. Default: 5")
        
        args = parser.parse_args()
        
        main(
        input_directory=args.input_directory,
        dataframe_output=args.dataframe_output,
        tdfi_output=args.tdfi_output,
        n_matches=args.n_matches,
        )