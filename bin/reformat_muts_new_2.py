#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os, glob
import sys

if __name__ == "__main__":
    final_mut_table = sys.argv[1]
    gene = sys.argv[2]
    incomplete = sys.argv[3]
    missing = sys.argv[4]
    split = sys.argv[5]
    linking = sys.argv[6]
    muts_graphing = sys.argv[7]
    input_muts = sys.argv[8]

# Read in LasR final mut table
table = pd.read_table(final_mut_table, sep = ",")

# Fill in mutations for clustered strains
table['Mutation'] = table.groupby('Cluster')['Mutation'].transform(lambda x: x.ffill().bfill() if x.isna().any() else x)

##This will not apply to every run, but in this case - more than 5 seqential mutations will categorize the protein as frameshift
def check_sequential_numbers(mutation_str):
    if pd.isna(mutation_str):
        return False
    numbers = [int(m.group()) for m in re.finditer(r'\d+', mutation_str)]
    numbers.sort()
    longest_seq_len = 1
    current_seq_len = 1
    for i in range(1, len(numbers)):
        if numbers[i] == numbers[i - 1] + 1:
            current_seq_len += 1
            longest_seq_len = max(longest_seq_len, current_seq_len)
        else:
            current_seq_len = 1
    
    return longest_seq_len >= 5

frameshift_rows = table['Mutation'].apply(check_sequential_numbers)

frameshift = table[frameshift_rows].copy()
table = table[~frameshift_rows].copy()

frameshift['Mut_Type'] = 'Frameshift'

# Identify premature stop mutations
def filter_stop(mutation_str):
    if pd.isna(mutation_str):
        return False
    pattern = re.compile(r'\d+X')
    for m in mutation_str.split(', '):
        if pattern.search(m):
            return True
    return False

# Apply the filtering function to create the new dataframe
stop_rows = table['Mutation'].apply(filter_stop)

stop = table[stop_rows].copy()
table = table[~stop_rows].copy()

stop['Mut_Type'] = 'Premature stop'

# Specific mutations
with open(input_muts, 'r') as f:
    mutation_patterns = [line.strip() for line in f if line.strip()]

# Identify SNPs
def filter_mutations(mutation_str, patterns):
    if pd.isna(mutation_str):
        return False
    pattern = re.compile(r'\b(' + '|'.join([re.escape(p) + r'(?:[A-Za-z])?' for p in patterns]) + r')\b')
    for m in mutation_str.split(', '):
        if pattern.search(m):
            return True
    return False

# Apply function
snp_rows = table['Mutation'].apply(lambda x: filter_mutations(x, mutation_patterns))

# Filter and update DataFrames
snp = table[snp_rows].copy()
table = table[~snp_rows].copy()

snp['Mut_Type'] = 'SNP'

# Concatenate all loss of function dataframes
all = pd.concat([stop, frameshift, snp])

# Label their functions as non functional
all['Mut_Status'] = 'No function'

# Subset just ID, mutation, mut_status, and mut_type columnns
new_all = all[['ID', 'Mutation', 'Mut_Status', 'Mut_Type']]

# Rename ID column
new_all = new_all.rename(columns={'ID': 'Prot_acc'})

# Read in protein to genome mappiong file
mapping = pd.read_table(linking, sep = "\t", names=["Genome", "Prot_acc"])

# Read in the LasR graphing table
graphing_table = pd.read_table(muts_graphing, sep = ",")

# Reformat protein accession for graphing table
graphing_table['Prot_acc'] = graphing_table['Prot_acc'].str.replace(r'_1$', '', regex=True)

graphing_prot = graphing_table[['Prot_acc']]

graphing_genomes_merge = mapping.merge(graphing_prot, on="Prot_acc", how="inner")

graphing_genomes = graphing_genomes_merge[['Genome']]

# Check sizes
missing_size = os.stat(missing).st_size
incomplete_size = os.stat(incomplete).st_size
split_size = os.stat(split).st_size

# Check missing
if missing_size > 0:
    # Read in a list of missing (genomes with fully deleted gene)
    missing = pd.read_table(missing, header = None)

    # Rename column
    missing = missing.rename(columns={0: 'Genome'})

    # Remove extension
    missing['Genome'] = missing['Genome'].str.replace(f'_{gene}.txt', '')

    # Set full deletion status to all full deletions
    missing['Mutation'] = 'Full Deletion'
    missing['Mut_Status'] = 'No function'
    missing['Mut_Type'] = 'Full Deletion'

    missing_genomes = missing[['Genome']]

    all_ran = pd.concat([missing_genomes, graphing_genomes])
else:
    all_ran = graphing_genomes.copy()

# Check incomplete
if incomplete_size > 0:
    # Read in the list of incomplete LasR (these have IS insertions)
    incomplete = pd.read_table(incomplete, header = None)

    # Rename column
    incomplete = incomplete.rename(columns={0: 'Genome'})

    # Remove extension
    incomplete['Genome'] = incomplete['Genome'].str.replace(f'_{gene}.txt', '')

    # Set IS insertion status to all incomplete LasR genomes
    incomplete['Mutation'] = 'Insertion or disruption'
    incomplete['Mut_Status'] = 'No function'
    incomplete['Mut_Type'] = 'Insertion or disruption'

    incomplete_genomes = incomplete[['Genome']]

    all_ran = pd.concat([all_ran, incomplete_genomes])
else:
    all_ran = all_ran.copy()

# Check split
if split_size > 0:
    # Read in a list of genes split on contig
    split = pd.read_table(split, header = None)

    # Rename column
    split = split.rename(columns={0: 'Genome'})

    # Remove extension
    split['Genome'] = split['Genome'].str.replace(f'_{gene}.txt', '')

    # Set status to NA
    split['Mutation'] = 'Contig split'
    split['Mut_Status'] = 'NA'
    split['Mut_Type'] = 'NA'

    split_genomes = split[['Genome']]

    all_ran = pd.concat([all_ran, split_genomes])
else:
    all_ran = all_ran.copy()

# Set status to ran
all_ran['Run_status'] = "Ran"

# Merge with mapping file
all_ran_merge = all_ran.merge(mapping, on = "Genome", how = "outer")

# Any that are missing are a full deletion
all_ran_merge['Run_status'] = all_ran_merge['Run_status'].fillna('Full Deletion n to p')

# Subset all full deletions
full_del_rows = all_ran_merge['Run_status'].isin(['Full Deletion n to p'])

# Set full deletion status to all full deletions (between nucleotide to protein) rows
full_del = all_ran_merge[full_del_rows].copy()

full_del = full_del[['Genome']]

full_del['Mutation'] = 'Full Deletion n to p'
full_del['Mut_Status'] = 'No function'
full_del['Mut_Type'] = 'Full Deletion n to p'

# Merge the graphing table with the new complete file
graphing_merge = mapping.merge(graphing_table, how = 'inner', on = 'Prot_acc')

# Only keep columns genome, status, and cheat
graphing_merge = graphing_merge[['Genome', 'Status', 'Cheat']]

# Subset all WT from the complete and graphing file merge
wt_rows = graphing_merge['Status'].isin(['WT'])

# Set WT status to WT rows
wt = graphing_merge[wt_rows].copy()

wt = wt[['Genome']]

wt['Mutation'] = 'WT'
wt['Mut_Status'] = 'Functional'
wt['Mut_Type'] = 'WT'

# Subset all potential deletions (matching less than 40% to LasR)
pot_del_rows = graphing_merge['Status'].isin(['Potential Deletion'])

# Set potential deletion status to all potential deletion rows
pot_del = graphing_merge[pot_del_rows].copy()

pot_del = pot_del[['Genome']]

pot_del['Mutation'] = 'Potential Deletion'
pot_del['Mut_Status'] = 'No function'
pot_del['Mut_Type'] = 'Potential Deletion'

# Reformat protein accession of SNP/frameshift/stop codon genomes
new_all['Prot_acc'] = new_all['Prot_acc'].str.replace(r'_1$', '', regex=True)

# Merge the SNP/frameshift/stop codon genomes with the mapping file to get genome information
new_all_merge = mapping.merge(new_all, how='inner', on='Prot_acc')

# Only keep genome, mutation, mut_status, and mut_type columns
new_all_merge = new_all_merge[['Genome', 'Mutation', 'Mut_Status', 'Mut_Type']]

# Concatenate the SNP/frameshift/stop codon genomes with WT, potential deletion, full deletion, missing, incomplete, and split genomes
if missing_size > 0:
    all_mut_con = pd.concat([wt, new_all_merge, pot_del, full_del, missing])
else:
    all_mut_con = pd.concat([wt, new_all_merge, pot_del, full_del])

if incomplete_size > 0:
    all_mut_con = pd.concat([all_mut_con, incomplete])
else:
    all_mut_con = all_mut_con.copy()

if split_size > 0:
    all_mut_con = pd.concat([all_mut_con, split])
else:
    all_mut_con = all_mut_con.copy()

# Remove duplicates, keeping first instance
all_mut_con = all_mut_con.drop_duplicates(subset=['Genome'], keep='first')

# Save
all_mut_con.to_csv(f"{gene}_all_functions_incomplete.csv", index = False)

# Obtain only genome and mut status
new_all_mut_con = all_mut_con[['Genome', 'Mut_Status']]

# Save
new_all_mut_con.to_csv(f"{gene}_functions_incomplete.csv", index = False)

# Change no function to 1, functional to 0
new_all_mut_con['Mut_Status'] = new_all_mut_con['Mut_Status'].str.replace('No function', '1')
new_all_mut_con['Mut_Status'] = new_all_mut_con['Mut_Status'].str.replace('Functional', '0')

# Save
new_all_mut_con.to_csv(f"{gene}_functions_pres_abs_incomplete.csv", index = False)