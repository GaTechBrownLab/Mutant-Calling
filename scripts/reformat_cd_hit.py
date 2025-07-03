import pandas as pd
import numpy as np
import glob, os
import sys

if __name__ == "__main__":
    base = sys.argv[1]

data=[]

# Read in split CD-Hit file
for file in glob.glob("x*"):
	if os.path.getsize(file) == 0:
		print(file)
	else:
		f=open(file, 'r')
		results = pd.read_table(f, header=None)
		results['Cluster']=file
		acc = results['Cluster'].str.split("x", expand = True)
		results['Cluster']= acc[2]
		data.append(results)

table = pd.concat([dfi.rename({old: new for new, old in enumerate(dfi.columns)}, axis=1) for dfi in data], ignore_index=True)
table.rename(columns = {0:'Number', 1:'ID', 2:'Cluster'}, inplace = True)

# Reformat CD-Hit table in order to parse
drop = table.drop(table.columns[[0]], axis = 1)

ident = drop['ID'].str.split(" ", expand = True)

drop['ID'] = ident[1]
drop['Length'] = ident[0]
drop['Status'] = ident[2]
if len(ident.columns) > 3:
    drop['Percent_similar'] = ident[3]
else:
    drop['Percent_similar'] = ''

drop['ID'] = drop['ID'].str.replace('>','')
drop['Length'] = drop['Length'].str.replace('aa,','')
drop['Status'] = drop['Status'].str.replace('*','ref')
drop['Status'] = drop['Status'].str.replace('at','clustered')

drop = drop[["ID", "Cluster", "Length", "Status"]]

drop.to_csv(f"{base}/{base}_cd_hit_table.csv", index=False)

drop = drop[drop['Status'].str.contains("ref")]

drop = drop.drop(drop.columns[[1, 2, 3]], axis = 1)

drop.to_csv(f"{base}/{base}_muts.txt", index=False, header=False)

