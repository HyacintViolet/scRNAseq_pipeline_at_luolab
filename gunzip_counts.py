# This is a tiny script that iteratively uncompress counts.tsv.gz files in each scRNAseq library folder.
# The "umi_tools counts" function outputs .tsv.gz files, which is cumbersome to work with. Here we iteratively
# uncompress them by invoking 'gunzip' with os.system() in python.

import os
import re

# Set working directory
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'

# Change working directory
os.chdir(data_wd)

# List all folders in wd
folders = os.listdir()

for folder in folders:  # Iterate over wd, go into one folder each time

    os.chdir(os.path.join(data_wd, folder))

    # Grab folder name
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)
    prefix = match.group(1)

    if not os.path.exists(os.path.join(data_wd, folder, '_'.join([prefix, 'counts.tsv']))):
        print('Extracting ' + prefix + 'counts.tsv.gz')
        for filename in os.listdir():  # Iterate over all files in each folder
            if filename.endswith(".tsv.gz"):  # Find counts.tsv.gz
                cmd = 'gunzip -k '+filename  # Construct command line
                os.system(cmd)  # Execute
    else:
        print(prefix + 'counts.tsv already exists.')

