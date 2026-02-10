#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

if __name__ == "__main__":
    input = sys.argv[1]
    subset_to_pull = sys.argv[2]
    data = sys.argv[3]
    outdir = sys.argv[4]

table = pd.read_csv(f"{data}/{input}")

table_sub = table.drop_duplicates(subset=['gene'], keep='first')

mask = (
    table_sub['Non.unique.Gene.name'].isna() |
    ~table_sub['Non.unique.Gene.name'].duplicated(keep='first')
)

table_cleaned = table_sub[mask]

top_20 = table_cleaned.head(int(subset_to_pull))

top_20.to_csv(f"{outdir}/Top_20_genes_adult_CF.csv", index = None)

top_20_genes = top_20[["gene"]]

top_20_genes.to_csv(f"{outdir}/Top_20_genes.txt", sep = "\t", index = None, header = None)