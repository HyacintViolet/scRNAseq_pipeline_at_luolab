import os, re
import pandas as pd



# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def_2/'

nametable = pd.read_table(os.path.join(parent_wd, 'gencode.vM4.annotation.tab'), sep="\t")

nametable = nametable.rename(columns = {'Unnamed: 0': 'gene_id', 'Unnamed: 1': 'gene_name'})

# Copied from stackOverflow. Works. Read later.
dictionary = nametable.set_index('gene_id')['gene_name'].T.to_dict()

# Iteratively enter each data folder
for folder in os.listdir(data_wd):

    # Change working directory
    os.chdir(os.path.join(data_wd, folder))

    # Extract names parts
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)

    nameparts = [match.group(1), 'counts', 'suffix.tsv']

    filename = '_'.join(nameparts)

    print('parsing ' + filename + ' ...')

    counts_suffix = pd.read_csv(filename, sep='\t')

    # df.replace is so slow for large dataframe
    # counts.replace({'gene': dictionary})

    # map() function is much more efficient
    mapped = counts_suffix['gene'].map(dictionary)
    counts_suffix.gene = mapped


    # Output
    out_nameparts = [match.group(1), 'suffix', 'renamed.tsv']
    out_filename = '_'.join(out_nameparts)
    counts_suffix.to_csv(out_filename, sep="\t", index=False)


