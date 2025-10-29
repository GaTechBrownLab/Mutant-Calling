#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob, os
import sys

if __name__ == "__main__":
     base = sys.argv[1]
     cd_hit_cluster = sys.argv[2]


# Open panaroo initial cd-hit output
with open(f"{cd_hit_cluster}", "r") as f:
	lines = f.read().strip().split("\n")

# Split into clusters based on lines starting with ">Cluster"
clusters = []
current_cluster = []
current_name = None

for line in lines:
    if line.startswith(">Cluster"):
        # Save previous cluster if it exists
        if current_cluster:
            df = pd.DataFrame(current_cluster, columns=["raw"])
            df["Cluster"] = current_name
            clusters.append(df)
            current_cluster = []
        current_name = line.strip().split()[1]  # e.g., "0", "1", "2"
    else:
        current_cluster.append([line.strip()])

# Add the last cluster
if current_cluster:
    df = pd.DataFrame(current_cluster, columns=["raw"])
    df["Cluster"] = current_name
    clusters.append(df)

# Combine all clusters into a single dataframe (optional)
combined_df = pd.concat(clusters, ignore_index=True)

# Split to get the correct reference
ident = combined_df['raw'].str.split(" ", expand = True)

ident.rename(columns = {0:'Length', 1:'ID', 2:'Status'}, inplace = True)

combined_df['ID'] = ident['ID']
combined_df['Status'] = ident['Status']

length = ident["Length"].str.split("\t", expand = True)

combined_df['Length'] = length[1]

combined_df_sub = combined_df[["ID", "Cluster", "Length", "Status"]]

# Reformat
combined_df_sub['ID'] = combined_df_sub['ID'].str.replace('>','')
combined_df_sub['ID'] = combined_df_sub['ID'].str.replace('...','')
combined_df_sub['Length'] = combined_df_sub['Length'].str.replace('aa,','')
combined_df_sub['Length'] = combined_df_sub['Length'].str.replace(',,','')
combined_df_sub['Status'] = combined_df_sub['Status'].str.replace('*','ref')
combined_df_sub['Status'] = combined_df_sub['Status'].str.replace('at','clustered')

combined_df_sub.to_csv(f"{base}_cd_hit_table.csv", index=False)

combined_df_sub = combined_df_sub[combined_df_sub['Status'].str.contains("ref")]

combined_df_sub = combined_df_sub.drop(combined_df_sub.columns[[1, 2, 3]], axis = 1)

combined_df_sub.to_csv(f"{base}_muts.txt", index=False, header=False)