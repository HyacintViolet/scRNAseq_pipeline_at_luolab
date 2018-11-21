import gzip, os
import pandas as pd
from os import system, listdir, chdir, path

# Set some default parameters, i.e. working directory and filename
path_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def'
dummy = 'whitelist.txt'

# Change working directory
chdir(path_wd)

for folder in listdir(path_wd):
    folder_path = path.join(path_wd,folder)
    for filename in listdir():
        if filename.endswith(".tsv.gz"):
            f = gzip.open(filename, mode='rt', encoding=None)
            f_content = f.read()

    data = pd.read_csv(path.join(path_wd, folder, filename), sep="\t", header=None)


