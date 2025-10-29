#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os, glob
import sys

if __name__ == "__main__":
    gene = sys.argv[1]
    functions = sys.argv[2]
    ani_clusters = sys.argv[3]
    envs = sys.argv[4]

# Read in functions
functions_df = pd.read_csv(functions, sep = ",")

# Read in ani clusters file
ani_clusters_df = pd.read_csv(ani_clusters, sep = ",")

merge = ani_clusters_df.merge(functions_df, how = "inner", on = "Genome")

ani_clusters_df_counts = merge.groupby(['Cluster','Threshold']).size().reset_index(name='Cluster_counts')

ani_clusters_df_counts_cutoff = merge.merge(ani_clusters_df_counts, on = ['Cluster','Threshold'], how = "outer")

filtered_ani_clusters_df = ani_clusters_df_counts_cutoff[ani_clusters_df_counts_cutoff['Cluster_counts'] > 1]


# Created categories
def cluster_status(status_series):
    unique_status = set(status_series)
    if unique_status == {'Functional'}:
        return 'Functional'
    elif unique_status == {'No function'}:
        return 'No function'
    else:
        return 'Mixed'

# Apply the function for each Cluster and Threshold group
filtered_ani_clusters_df['Cluster_Status'] = filtered_ani_clusters_df.groupby(['Cluster','Threshold'])['Mut_Status'].transform(cluster_status)

# Merge environments
envs_df = pd.read_csv(envs, sep = ",")

envs_df_group = envs_df[['Genome', 'Group']]

merge_envs = filtered_ani_clusters_df.merge(envs_df_group, on = "Genome", how = "left")

merge_envs_cluster = merge_envs[["Cluster", "Threshold", "Cluster_Status"]]

merge_envs_cluster = merge_envs_cluster.drop_duplicates()

cluster_envs_counts = merge_envs.groupby(['Cluster','Threshold','Group']).size().reset_index(name='count')

cluster_envs_counts_unique = cluster_envs_counts.groupby(['Cluster','Threshold'])['Group'].nunique().reset_index(name='Unique_Group_Count')

cluster_envs_counts_status = cluster_envs_counts_unique.merge(merge_envs_cluster, on = ['Cluster','Threshold'], how = "left")

cluster_envs_counts_status.to_csv(f"{gene}_ANI_cluster_env_counts_status.csv", index = False)