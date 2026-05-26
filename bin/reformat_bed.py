#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

if __name__ == "__main__":
    input_file = sys.argv[1]
    base = sys.argv[2]

# Read in blast results
table = pd.read_csv(input_file, sep = "\t", names = ['Genome', 'Gene', 'Prot_acc', 'Prot_acc_stripped', 'PID', 'nident', 'Gene_Length', 'Aln_length', 'evalue', 'slen', 'qstart', 'qend', 'Start', 'End', 'Strand'])

# Reformat protein accession
split_genome = table['Genome'].str.split("_PA", expand = True)

table['Genome'] = split_genome[0]

table['Accession'] = table['Prot_acc'].str.extract(r'\|(.*?)\|')[0].fillna(table['Prot_acc'])

# Calculate query coverage
table["Query_Cov"] = (table["Aln_length"]/table["Gene_Length"])*table["PID"]

# Drop any rows with lower query coverage for each Accession
table = table.sort_values('Query_Cov', ascending=False).drop_duplicates('Genome').sort_index()

# Reformat into bed file
table_drop = table[['Accession', 'Start', 'End', 'Strand']]

table_drop.insert(3, 'col1', '.')
table_drop.insert(4, 'col2', '.')

table_drop['Strand'] = table_drop['Strand'].str.replace('minus', '-').str.replace('plus', '+')

def swap_if_negative(row):
    """Switch start and end if blast result is on the negative strand"""
    if row['Strand'] == '-':
        row['Start'], row['End'] = row['End'], row['Start']
    return row

# Apply the function to each row
table_drop = table_drop.apply(swap_if_negative, axis=1)

table_drop['Start'] = table_drop['Start'] - 1

table_drop.to_csv(f"{base}_bed.bed", sep="\t", index = False, header = False)

