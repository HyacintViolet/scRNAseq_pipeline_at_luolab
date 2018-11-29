"""
In the R pipeline, the CELL CYCLE CHECK pipeline only takes in ENSEMBL STABLE IDS (ENSMUSG________.__) with a dot-less
format. As a result, the genes must only be renamed after concatenation of the counts matrix. Previous codes dealing
with gene renaming in individual library count's spread-sheet are obsolete.
"""

import os
import pandas as pd


# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4/'
os.chdir(parent_wd)

# Load name table [ENSEMBL STABLE ID, gene name]
nametable = pd.read_table(os.path.join(parent_wd, 'gencode.vM4.annotation.tab'), sep="\t")

nametable = nametable.rename(columns={'Unnamed: 0': 'stable_id', 'Unnamed: 1': 'gene_name'})

# Remove dots (version id)
ensmusg = nametable.stable_id.str.split(".", n=1, expand=True)
nametable.stable_id = ensmusg[0]

# Copied from stackOverflow. Works. Read later.
dictionary = nametable.set_index('stable_id')['gene_name'].T.to_dict()

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


