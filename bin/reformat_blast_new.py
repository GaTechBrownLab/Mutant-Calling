#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

if __name__ == "__main__":
    input_table = sys.argv[1]
    base = sys.argv[2]
    WT_cutoff = sys.argv[3]
    mut_cutoff = sys.argv[4]

# Read in amino acid blast result table
table = pd.read_csv(input_table, sep = "\t", header = None)

# Rename columns
table.rename(columns = {0:'Prot_acc', 3:'PID', 6:'Aln_length', 8:'Gene_Length'}, inplace = True)
table_drop = table[["Prot_acc", "PID", "Aln_length", "Gene_Length"]]

# Calculate query coverage
table_drop["Query_Cov"] = (table_drop["Aln_length"]/table_drop["Gene_Length"])*table_drop["PID"]

table_drop = table_drop.sort_values('Query_Cov', ascending=False).drop_duplicates('Prot_acc').sort_index()

# Label wild-type as those that match to query coverage at 100%
WT = table_drop[(table_drop['Query_Cov'] == float(WT_cutoff))]

# Label mutant as those that match to query coverage between 40%-100%
mut = table_drop[(table_drop['Query_Cov'] > float(mut_cutoff)) & (table_drop['Query_Cov'] < float(WT_cutoff))]

mut.to_csv(f"{base}_muts_and_sig_diff.csv", index = False)
WT.to_csv(f"{base}_WT.csv", index = False)

# Get the accessions of just the mutants
mut1 = mut.copy()

mut1 = mut1.drop(mut1.columns[1:9], axis = 1)

mut1.to_csv(f"{base}_muts_accessions.txt", index = False, header = False)

# Assign deletion to those with less than 40% query coverage
deletions = table_drop.copy()

deletions = deletions[(deletions['Query_Cov'] <= float(mut_cutoff))]

deletions1 = deletions.copy()

deletions1 = deletions1.drop(deletions1.columns[1:9], axis = 1)

deletions1.to_csv(f"{base}_potential_deletions_accessions.txt", index = False, header = False)

# Create statuses for graphing table
WT['Status'] = "WT"
WT['Cheat'] = "WT"
mut['Status'] = "Mut"
mut['Cheat'] = "Mut"
deletions['Status'] = "Potential Deletion"
deletions['Cheat'] = "Potential Deletion"

las = pd.concat([WT, mut])

all_las = pd.concat([las, deletions])

all_las["qs_protein"] = f"{base}"

all_las.to_csv(f"{base}_muts_graphing.csv", index = False)

