#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os, glob
import sys

functions_df = pd.read_csv("/data1/I_Irby/Mutant-Calling/results/ouputs_for_graphing/PA1430/PA1430_functions_incomplete.csv", sep = ",")

def dif_clusters(file_path):
    value = file_path.split("/")[-1].replace("_output.csv", "")
    ani_clusters_df = pd.read_csv(file_path, sep = ",")
    
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

    filtered_ani_clusters_df_sub = filtered_ani_clusters_df[['Genome', 'Cluster', 'Cluster_Status']]

    filtered_ani_clusters_df_sub.rename(columns = {'Cluster':f'{value}_Cluster', 'Cluster_Status':f'{value}_Cluster_Status'}, inplace = True)

    return(filtered_ani_clusters_df_sub)

ani_99_5 = dif_clusters("/data1/I_Irby/Mutant-Calling/results/cluster_ani_output/99.50000000_output.csv")
ani_99_5_3 = dif_clusters("/data1/I_Irby/Mutant-Calling/results/cluster_ani_output/99.53500000_output.csv")

ani_99_5_3.to_csv("Test.csv", index = False)
merge_high = ani_99_5.merge(ani_99_5_3, on = "Genome", how = "outer")

merge_high_sub = merge_high[merge_high["99.50000000_Cluster_Status"] != merge_high["99.53500000_Cluster_Status"]]

merge_high_sub.to_csv("Incorrect_clusters.csv", index = False)