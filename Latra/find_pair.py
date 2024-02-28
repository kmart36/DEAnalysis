# %%
import pandas as pd
import numpy as np

# %%
#Make this easier to change via Nextflow
path_to_kraken = "/home/kam071/DEAnalysis/Latra/kraken_files/"
path = "/home/kam071/DEAnalysis/Latra/"

# %%
reads = pd.read_csv(path + "latra_reads.csv", index_col=0)

# %%
males = reads[reads.Sex == "M"].sort_values(by=["After Kraken"], ascending=False)
females = reads[reads.Sex == "F"].sort_values(by=["After Kraken"], ascending=False)

# %%
males = males[males['Locality'].isin(females['Locality'])]
females = females[females['Locality'].isin(males['Locality'])]
if males.empty or females.empty:
    print("No males and females in the same locality.")
    exit()
males = males[males['Tissue'] == "antennae"]
females = females[females['Tissue'] == "antennae"]


# %%
sample_m = reads[reads.index.str.contains(males.index[0].split('_')[0])].index
sample_f = reads[reads.index.str.contains(females.index[0].split('_')[0])].index

# %%
left = [path_to_kraken + x + "_R_1.unclassified.out.fq" for x in sample_m] + [path_to_kraken + x + "_R_1.unclassified.out.fq" for x in sample_f]
right = [path_to_kraken + x + "_R_2.unclassified.out.fq" for x in sample_m] + [path_to_kraken + x + "_R_2.unclassified.out.fq" for x in sample_f]
left = ','.join(left)
right = ','.join(right)

# %%
with open(path + 'left_kraken.txt', 'w') as f:
    f.write(left)
with open(path + 'right_kraken.txt', 'w') as f:
    f.write(right)
print(left)
print(right)


