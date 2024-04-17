# %%
import pandas as pd
import numpy as np
import glob

# %%
path_to_kraken = "/home/kam071/DEAnalysis/Latra/kraken_files/"
path = "/home/kam071/DEAnalysis/Latra/"

# %%
reads = pd.read_csv("Copy of Pyralis Injection  - Injection.csv", index_col=0)

# %%
infected = reads[reads.Treatment.str.contains('PBS')].Identification
noninfected = ~reads[reads.Treatment.str.contains('PBS')].Identification

# %%
infected_left = [path_to_kraken + x + '_R1*' for x in infected]
infected_right = [path_to_kraken + x + '_R2*' for x in infected]
non_infe_left = [path_to_kraken + x + '_R1*' for x in noninfected]
non_infe_right = [path_to_kraken + x + '_R2*' for x in noninfected]

# %%
infected_left = ','.join(infected_left)
infected_right = ','.join(infected_right)
non_infe_left = ','.join(non_infe_left)
non_infe_right = ','.join(non_infe_right)

# %%
with open(path + 'left_infected_kraken.txt', 'w') as f:
    f.write(infected_left)
with open(path + 'right_infected_kraken.txt', 'w') as f:
    f.write(infected_right)
with open(path + 'left_noninfected_kraken.txt', 'w') as f:
    f.write(non_infe_left)
with open(path + 'right_non_infected_kraken.txt', 'w') as f:
    f.write(non_infe_right)


