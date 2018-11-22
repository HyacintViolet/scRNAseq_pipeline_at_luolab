import os
import pandas as pd



# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def/'


nametable = pd.read_table(os.path.join(parent_wd,'gencode.vM4.annotation.tab'), sep="\t")

nametable.rename(columns = {'Unnamed: 0': 'gene_id', 'Unnamed: 1': 'gene_name'})