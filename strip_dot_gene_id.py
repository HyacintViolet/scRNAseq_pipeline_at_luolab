"""
In the R pipeline, the CELL CYCLE CHECK pipeline only takes in ENSEMBL STABLE IDS (ENSMUSG________.__) with a dot-less
format. As a result, the genes must only be renamed after concatenation of the counts matrix. Previous codes dealing
with gene renaming in individual library count's spread-sheet are obsolete.
"""

import os
from collections import Counter
import numpy as np
import pandas as pd


def has_duplicates(list_of_values):
    if len(np.unique(list_of_values)) < len(list_of_values):
        print("The list \'" + list_of_values.name + "\' contains duplicate values.")
    else:
        print("The list \'" + list_of_values.name + "\' does not contain duplicate values.")


# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/dat_gfp/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/dat_gfp/mapping/'
os.chdir(parent_wd)

# Load name table [ENSEMBL STABLE ID, gene name]
nametable = pd.read_table(os.path.join(parent_wd, 'gencode.vM19.annotation.tab'), sep="\t", names=["stable_id",
                                                                                                   "gene_name"])
# nametable = nametable.rename(columns={'Unnamed: 0': 'stable_id', 'Unnamed: 1': 'gene_name'})

# Note that not all the elements in gene_name is unique.
has_duplicates(nametable.stable_id)
has_duplicates(nametable.gene_name)

# Strip dots
ensmusg = nametable.stable_id.str.split(".", n=1, expand=True)
nametable_new = pd.concat([nametable.stable_id, ensmusg[0]], axis=1)
nametable_new.columns = ['stable_id', 'stable_id_dotless']

# Convert nametable(df) to dict. Copied from stackOverflow. Works. Read later.
dictionary = nametable_new.set_index('stable_id')['stable_id_dotless'].T.to_dict()

# Load expression matrix (QC2)
expression_mat = pd.read_csv('counts_stbid_QC1.txt', sep=' ')
print(expression_mat.shape)

expression_mat.id = expression_mat.id.map(dictionary)

expression_mat.to_csv('counts_QC1_dotless.txt', sep=' ', index=False)
print('Finished.')