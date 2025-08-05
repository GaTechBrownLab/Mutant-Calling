#!/usr/bin/env python3

from Bio import AlignIO
import pandas as pd
import numpy as np
import os
import sys

if __name__ == "__main__":
    base = sys.argv[1]
    genome = sys.argv[2]
    file = sys.argv[3]

# Load the alignment using AlignIO
alignment = AlignIO.read(file, "fasta")

# Set the first sequence as the reference sequence
reference_seq = str(alignment[0].seq)

data = []

# Iterate through each position in the alignment, setting each aa in PAO1 protein as reference
for i in range(len(reference_seq)):
    ref_aa = reference_seq[i]

    # Iterate through each sequence in the alignment
    for record in alignment:
        sequence_id = record.id
        sequence = str(record.seq)
        seq_aa = sequence[i]

        # Check for differences
        if seq_aa != ref_aa:
            data.append({'ID': sequence_id, 'Pos': i + 1, 'Reference_AA': ref_aa, 'Sample_AA': seq_aa})

df = pd.DataFrame(data, columns=['ID', 'Pos', 'Reference_AA', 'Sample_AA'])

#Set deletion as an empty column (this will represent deletions or insertions)
df['Deletion'] = ''

#If the sample amino acid is deleted ("-"), find the starting position of that amino acid and the next location where the next amino acid is not a deletion
#then pull the index of the position that is one before that next amino acid change, which is at the end of the deletion

start_num = None
end_num = None
for index, row in df.iterrows():
    if row['Sample_AA'] == '-' and start_num is None:
        start_num = row['Pos']
    if index < len(df) - 1:
        if row['Sample_AA'] != '-' and start_num is not None or df.loc[index + 1, 'Pos'] != row['Pos'] + 1 and start_num is not None:
            end_num = df.loc[index, 'Pos']
            df.loc[index, 'Deletion'] = f'Deletion {start_num} to {end_num}'
            start_num = None

# Handle the last consecutive sequence if any
if start_num is not None:
    df.loc[index, 'Deletion'] = f'Deletion {start_num} to {df.iloc[-1]["Pos"]}'

#Same code but for insertions
start_num = None
end_num = None
for index, row in df.iterrows():
    if row['Reference_AA'] == '-' and start_num is None:
        start_num = row['Pos']
    elif index < len(df) - 1:
        if row['Reference_AA'] != '-' and start_num is not None or df.loc[index + 1, 'Pos'] != row['Pos'] + 1 and start_num is not None:
            end_num = df.loc[index, 'Pos']
            df.loc[index, 'Deletion'] = f'Insertion {start_num} to {end_num}'
            start_num = None

if start_num is not None:
    df.loc[index, 'Deletion'] = f'Insertion {start_num} to {df.iloc[-1]["Pos"]}'

#Identify specific mutation in correct notations, including deletions and insertions
df['Pos'] = df['Pos'].astype(str)
for index, row in df.iterrows():
    if row['Sample_AA'] != '-' and row['Reference_AA'] != '-':
        df.loc[index,'Mutation'] = df.loc[index,'Reference_AA'] + df.loc[index,'Pos'] + df.loc[index,'Sample_AA']
    elif row['Sample_AA'] == '-' or row['Reference_AA'] == '-':
        df.loc[index,'Mutation'] = df.loc[index,'Deletion']


df['Mutation'].replace('', np.nan, inplace=True)
df.dropna(subset=['Mutation'], inplace=True)

df = df[['ID', 'Mutation']]

#Merge cells based on protein ID
merged_df = df.groupby('ID')['Mutation'].agg(lambda x: ', '.join(map(str, x))).reset_index()

merged_df.to_csv(f"{base}/alignments/csv/{genome}.csv", index = False, header = False)