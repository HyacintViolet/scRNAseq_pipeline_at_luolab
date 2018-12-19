import os
import re
import pandas as pd


# Generic function to add suffix to the 'cell' column in counts.tsv DataFrame.
# suffix is equivalent to codename
def add_suffix(df, suffix):
    df.cell = df.cell + suffix


# Set up some default parameters, i.e. working directory and filename
parent_wd = '/media/luolab/ZA1BT1ER/linrui/DR_DAT/'
data_wd = '/media/luolab/ZA1BT1ER/linrui/DR_DAT/vM19/'

# Change working directory
os.chdir(parent_wd)

# Import experimental design, which stores library name, index, codename & experiment setup.
exp_design = pd.read_excel('experimental_design.xlsx')

for folder in os.listdir(data_wd):

    # Iteratively enter each library
    library_wd = os.path.join(data_wd, folder)
    os.chdir(library_wd)

    # Extract names from folder
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)

    nameparts = [match.group(1), 'counts.tsv']

    filename = '_'.join(nameparts)

    print('parsing ' + filename + ' ...')

    # Read counts data
    data = pd.read_csv(filename, sep="\t")

    # Find codename w.r.t. lib_prefix to use as suffix
    codename = exp_design.codename[exp_design.lib_prefix == match.group(1)]

    # Use iloc to access by position rather than label (awesome!)
    suffix = '_' + codename.iloc[0]

    # Add suffix
    add_suffix(data, suffix)

    # Output XXXXXXXX_counts_suffix.tsv
    # *** The naming convention should absolutely be simplified in the future ***
    out_nameparts = [match.group(1), 'counts', 'suffix.tsv']
    out_filename = '_'.join(out_nameparts)
    data.to_csv(out_filename, sep="\t", index=False)

    print('finished\n')
# # Example code to import data
# gene_cell_counts = pd.read_csv('YT101001_TKR181000273_HNGCTCCXY_L6_counts.tsv', sep="\t")
