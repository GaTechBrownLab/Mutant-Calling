#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os, glob
import sys

if __name__ == "__main__":
    base = sys.argv[1]
    table = sys.argv[2]
    cd_hit_table = sys.argv[3]

table = pd.read_csv(table, names=['ID', 'Mutation'])

#Read in cd-hit final table
cd_hit = pd.read_table(cd_hit_table, sep=",")

merge = cd_hit.merge(table, how='outer', on='ID')

#Save final table with mutations annotated for each cluster
merge.to_csv(f"{base}_final_mut_table.csv", index = False)

#Simpler reference table
merge = merge[merge['Status'].str.contains("ref")]

merge.to_csv(f"{base}_final_mut_table_refs.csv", index = False)