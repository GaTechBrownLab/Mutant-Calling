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

complete=f"{gene}/Blast_gene_status/Complete"
incomplete=f"{gene}/Blast_gene_status/Incomplete"
missing=f"{gene}/Blast_gene_status/Missing"

table = pd.read_csv(blast_output, sep = "\t", names = ['Gene', 'Prot_acc', 'Prot_acc_stripped', 'PID', 'nident', 'Gene_Length', 'Aln_length', 'evalue', 'slen', 'qstart', 'qend', 'Start', 'End', 'Strand'])

filename = os.path.basename(blast_output)
filename_no_ext = os.path.splitext(filename)[0]
filename_no_ext_clean = re.sub(f"_{gene}$", "", filename_no_ext)

def saving_funct(folder, type_st):
    table.to_csv(f"{folder}/{filename}", index = False, header = False, sep = "\t")

    with open(f"{gene}/{type_st}_output.txt", "a") as f:
        f.write(filename_no_ext_clean + "\n")

if len(table) == 0:
    saving_funct(missing, "Missing")
elif len(table) == 1:
    saving_funct(complete, "Complete")
else:
    if table.iloc[0]['Gene_Length'] == table.iloc[0]['Aln_length']:
        saving_funct(complete, "Complete")
    else:
        length_threshold = math.floor(float(length) - float(length) * float(gene_proportion))

        if table.iloc[0]['Aln_length'] >= length_threshold:
            saving_funct(complete, "Complete")
        else:
            gene_diff_threshold = math.floor(float(length) * float(gene_difference))

            if table.iloc[0]['Aln_length'] <= gene_diff_threshold:
                saving_funct(missing, "Missing")
            else:
                if table.iloc[0]['qstart'] == 1:
                    diff_qstart = abs(table.iloc[1]['qstart'] - table.iloc[0]['qend'])
                    diff_qend = abs(table.iloc[1]['qend'] - float(length))

                    if (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)):
                       saving_funct(incomplete, "Incomplete")
                    else:
                        saving_funct(complete, "Complete")
                elif table.iloc[1]['qstart'] == 1:
                    diff_qstart = abs(table.iloc[0]['qstart'] - table.iloc[1]['qend'])
                    diff_qend = abs(table.iloc[0]['qend'] - float(length))

                    if (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)):
                        saving_funct(incomplete, "Incomplete")
                    else:
                        saving_funct(complete, "Complete")
                else:
                    saving_funct(complete, "Complete")

