# This is a tiny script that iteratively uncompress counts.tsv.gz files in each scRNAseq library folder.
# The "umi_tools counts" function outputs .tsv.gz files, which is cumbersome to work with. Here we iteratively
# uncompress them by invoking 'gunzip' with os.system() in python.

import os

# Set working directory
parent_wd = '/media/luolab/ZA1BT1ER/yanting/dat_gfp/mapping/'

# Change working directory
os.chdir(parent_wd)

# List all folders in wd
folders = os.listdir()

for folder in folders: # Iterate over wd, go into one folder each time
    print('Extracting'+folder)
    os.chdir(os.path.join(parent_wd, folder))
    for filename in os.listdir():  # Iterate over all files in each folder
        if filename.endswith(".tsv.gz"):  # Find counts.tsv.gz
            cmd = 'gunzip -k '+filename  # Construct command line
            os.system(cmd)  # Execute

