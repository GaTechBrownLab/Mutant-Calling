#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os, glob
import sys
import math
import re

if __name__ == "__main__":
    length = sys.argv[1]
    blast_output = sys.argv[2]
    gene = sys.argv[3]
    gene_proportion = sys.argv[4]
    gene_difference = sys.argv[5]
    tolerance = sys.argv[6]
    split_tolerance = sys.argv[7]

complete=f"{gene}/Blast_gene_status/Complete"
incomplete=f"{gene}/Blast_gene_status/Incomplete"
missing=f"{gene}/Blast_gene_status/Missing"
split=f"{gene}/Blast_gene_status/Split"

# Print gene length and what gene proportion equates to
print(length, gene_proportion)

# Read in blast result
table = pd.read_csv(blast_output, sep = "\t", names = ['Gene', 'Prot_acc', 'Prot_acc_stripped', 'PID', 'nident', 'Gene_Length', 'Aln_length', 'evalue', 'slen', 'qstart', 'qend', 'Start', 'End', 'Strand'])

filename = os.path.basename(blast_output)
filename_no_ext = os.path.splitext(filename)[0]
filename_no_ext_clean = re.sub(f"_{gene}$", "", filename_no_ext)

def saving_funct(folder, type_st):
    """Function to save blast output to defined folder"""
    table.to_csv(f"{folder}/{filename}", index = False, header = False, sep = "\t")

    with open(f"{gene}/{type_st}_output.txt", "a") as f:
        f.write(filename_no_ext_clean + "\n")

# If blast output is empty, assign to missing
if len(table) == 0:
    saving_funct(missing, "Missing")

# If blast output yeilds one result, assign to complete
elif len(table) == 1:
    saving_funct(complete, "Complete")
else:
    # If gene length in the first line matches alignment length in the first line, save to complete
    if table.iloc[0]['Gene_Length'] == table.iloc[0]['Aln_length']:
        saving_funct(complete, "Complete")
    else:
        # Calculate length threshold
        length_threshold = math.floor(float(length) - (float(length) * float(gene_proportion)))

        # Calculate split threshold
        split_threshold = table.iloc[0]['slen'] - float(split_tolerance)

        # Establishes contig breaks (if a contig starts or ends in the middle of the gene of interest)
        if table.iloc[0]['Start'] == 1:
            saving_funct(split, "Split")
        elif table.iloc[0]['End'] == 1:
            saving_funct(split, "Split")
        elif table.iloc[0]['Start'] >= split_threshold:
            saving_funct(split, "Split")
        elif table.iloc[0]['End'] >= split_threshold:
            saving_funct(split, "Split")
        else:
            # If alignment length of the first line is greater than or equal to the length threshold, assign to complete
            if table.iloc[0]['Aln_length'] >= length_threshold:
                saving_funct(complete, "Complete")
            else:
                # Calculate gene difference threshold
                gene_diff_threshold = math.floor(float(length) * float(gene_difference))

                # If alignment length of the first line is less than the length threshold, assign to missing
                if table.iloc[0]['Aln_length'] < gene_diff_threshold:
                    saving_funct(missing, "Missing")
                else:
                    if table.iloc[0]['qstart'] == 1:
                        # If the query start in the first line equals one
                        # Calculate the differences between query start and end of the first two lines
                        diff_qstart = abs(table.iloc[1]['qstart'] - table.iloc[0]['qend'])
                        diff_qend = abs(table.iloc[1]['qend'] - float(length))

                        # If the differences are less than or equal to the tolerances and the contig accessions match, assign as incomplete (biologically split)
                        # If the contig accessions do not match, assign to contig split
                        if (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] == table.iloc[1]['Prot_acc']:
                            saving_funct(incomplete, "Incomplete")
                        elif (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] != table.iloc[1]['Prot_acc']:
                            saving_funct(split,  "Split")
                        else:
                            saving_funct(complete, "Complete")
                    elif table.iloc[1]['qstart'] == 1:
                        # If the query start in the second line equals one
                        # Calculate the differences between query start and end of the first two lines
                        diff_qstart = abs(table.iloc[0]['qstart'] - table.iloc[1]['qend'])
                        diff_qend = abs(table.iloc[0]['qend'] - float(length))
                        
                        # If the differences are less than or equal to the tolerances and the contig accessions match, assign as incomplete (biologically split)
                        # If the contig accessions do not match, assign to contig split
                        if (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] == table.iloc[1]['Prot_acc']:
                            saving_funct(incomplete, "Incomplete")
                        elif (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] != table.iloc[1]['Prot_acc']:
                            saving_funct(split, "Split")
                        else:
                            saving_funct(complete, "Complete")
                    else:
                        saving_funct(complete, "Complete")

