# %%
import pandas as pd
import numpy as np
import glob

# %%
path_to_kraken = "/home/kam071/DEAnalysis/Latra/kraken_files/"
path = "/home/kam071/DEAnalysis/Latra/"

# %%
infected_left = glob.glob(path_to_kraken + '*_infected_1.fq')
infected_right = glob.glob(path_to_kraken + '*_infected_2.fq')
non_infe_left = glob.glob(path_to_kraken + '*_noninfected_1.fq')
non_infe_right = glob.glob(path_to_kraken + '*_noninfected_2.fq')

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


