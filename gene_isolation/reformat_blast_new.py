#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

if __name__ == "__main__":
    input = sys.argv[1]
    data = sys.argv[2]
    outdir = sys.argv[3]

# Read in amino acid blast result table
table = pd.read_csv(input, sep = "\t", header = None)

# Rename columns
table.rename(columns = {0:'Prot_acc', 1: 'Locus Tag', 3:'PID', 6:'Aln_length', 8:'Gene_Length'}, inplace = True)
table_drop = table[["Prot_acc", "Locus Tag", "PID", "Aln_length", "Gene_Length"]]

# Calculate query coverage
table_drop["Query_Cov"] = (table_drop["Aln_length"]/table_drop["Gene_Length"])*table_drop["PID"]

table_drop = table_drop.sort_values('Query_Cov', ascending=False).drop_duplicates('Prot_acc').sort_index()

pao1_acc = pd.read_csv(f"{data}/Pseudomonas_aeruginosa_PAO1_107.csv", sep = ",")
pao1_acc_drop = pao1_acc[["Locus Tag", "Gene Name","Gene synonyns", "Product Name"]]

pao1_acc_merge = table_drop.merge(pao1_acc_drop, how = "left", on = "Locus Tag")

pao1_acc_merge.to_csv(f"{outdir}/PAO1_blast_result_new.csv", index = None)

pao1_acc_merge_sub = pao1_acc_merge[["Locus Tag"]]

pao1_acc_merge_sub.to_csv(f"{outdir}/PA_locus_tags.txt", sep = "\t", index = None, header = None)