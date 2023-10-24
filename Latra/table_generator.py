# %%
import pandas as pd
import glob
import numpy as np

path = "/home/kam071/DEAnalysis/Latra/"

# %%
# Change these paths so they locate the fastqc_data.txt in each of the fastqc folders generated. 
raws_fastqc = glob.glob(path + 'fastqc_summaries/*_1_fastqc/fastqc_data.txt')
trim_fastqc = glob.glob(path + 'fastqc_summaries/*_1_paired_fastqc/fastqc_data.txt')
krak_fastqc = glob.glob(path + 'fastqc_files/*_R_2_fastqc/*.unclassified.out_fastqc/fastqc_data.txt')

# %%
raws_fastqc = sorted(raws_fastqc, key=str.lower)
trim_fastqc = sorted(trim_fastqc, key=str.lower)
krak_fastqc = sorted(krak_fastqc, key=str.lower)

# %%
# Change this path so it leads to where your fastqc files are stored, but do not change anything after the last /.
ind = glob.glob(path + 'fastqc_summaries/*_1_fastqc')
ind = ['_'.join(x.split('/')[-1].split('_')[0:2]) for x in ind]

# %%
raws = []
trims = []
kraks = []
for i in raws_fastqc:
    raws.append(open(i, "r").readlines(150)[-1].split('\t')[-1].strip())
for j in trim_fastqc:
    trims.append(open(j, "r").readlines(170)[-1].split('\t')[-1].strip())
for k in krak_fastqc:
    kraks.append(open(k, "r").readlines(170)[-1].split('\t')[-1].strip())

# %%
table = pd.DataFrame(data=[raws, trims, kraks], index=["Raw Reads", "After Trimmomatic", "After Kraken"], columns=ind)

# %%
table = table.transpose()

# %%
# Use this line to drop any samples if necessary
table = table.drop(["Pcorr585_BL", "Pcorr585_ant"])

# %%
# Make sure you point to the metadata file that includes tissue and locality.
metadata = pd.read_csv(path + 'latra_metadata.csv')
metadata = metadata.set_index(table.index)

# %%
new_tbl = pd.concat([metadata, table], axis=1)

# %%
new_tbl.to_csv(path + "Latra_reads.csv")