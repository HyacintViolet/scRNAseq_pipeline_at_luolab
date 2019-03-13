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
parent_wd = '/media/luolab/ZA1BT1ER/yanting/vM19/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM19/mapping/'
os.chdir(parent_wd)

# Load name table [ENSEMBL STABLE ID, gene name]
nametable = pd.read_table(os.path.join(parent_wd, 'gencode.vM19.annotation.tab'), sep="\t", names=["stable_id",
                                                                                                   "gene_name"])
# nametable = nametable.rename(columns={'Unnamed: 0': 'stable_id', 'Unnamed: 1': 'gene_name'})

# Note that not all the elements in gene_name is unique. If not fixed, the FindVariableGenes pipe will return error.
has_duplicates(nametable.stable_id)
has_duplicates(nametable.gene_name)

# Add suffix to duplicate gene_names
names = nametable.gene_name.tolist()
counts = Counter(names)
for s, num in counts.items():
    if num > 1:
        for suffix in range(1, num+1):
            # *** For some unknown reason, the following line doesn't workwithout '_'
            names[names.index(s)] = s + '_' + str(suffix)
nametable_new = pd.concat([nametable.stable_id, pd.Series(names, name='gene_name')], axis=1)

# # Debugging above snippet.
# test_list = ['x','x','x1','x1','x1','y','y','z']
# counts_test = Counter(test_list)
# for s, num in counts_test.items():
#     if num>1:
#         for suf in range(1, num+1):
#             test_list[test_list.index(s)] = s + str(suf)

# Confirm again
has_duplicates(nametable_new.stable_id)
has_duplicates(nametable_new.gene_name)

# Remove dots (version id)
ensmusg = nametable_new.stable_id.str.split(".", n=1, expand=True)
nametable_new.stable_id = ensmusg[0]

# Convert nametable(df) to dict. Copied from stackOverflow. Works. Read later.
dictionary = nametable_new.set_index('stable_id')['gene_name'].T.to_dict()

# Load expression matrix (QC2)
expression_mat = pd.read_csv('counts_stbid_QC1.txt', sep=' ')
print(expression_mat.shape)

expression_mat.id = expression_mat.id.map(dictionary)

expression_mat.to_csv('counts_QC1_renamed.txt', sep=' ', index=False)

# ----------------------------------------------------------------------------------------------------------------------
# The following codes are obsolete due to pipeline adjustment.
# ----------------------------------------------------------------------------------------------------------------------


# # Iteratively enter each data folder
# for folder in os.listdir(data_wd):
#
#     # Change working directory
#     os.chdir(os.path.join(data_wd, folder))
#
#     # Extract names parts
#     match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)
#
#     nameparts = [match.group(1), 'counts', 'suffix.tsv']
#
#     filename = '_'.join(nameparts)
#
#     print('parsing ' + filename + ' ...')
#
#     counts_suffix = pd.read_csv(filename, sep='\t')
#
#     # df.replace is so slow for large dataframe
#     # counts.replace({'gene': dictionary})
#
#     # map() function is much more efficient
#     mapped = counts_suffix['gene'].map(dictionary)
#     counts_suffix.gene = mapped
#
#
#     # Output
#     out_nameparts = [match.group(1), 'suffix', 'renamed.tsv']
#     out_filename = '_'.join(out_nameparts)
#     counts_suffix.to_csv(out_filename, sep="\t", index=False)


