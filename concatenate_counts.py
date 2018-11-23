import os, re
import pandas as pd

# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def/'
os.chdir(parent_wd)

folder_list = os.listdir(data_wd)

# I know I know... Iteration is stupid... Let's fix this later.
for i in range(len(folder_list)):
    if i == 0:

        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', folder_list[0])

        count_filename_parts = [match.group(1), match.group(2), match.group(3), match.group(4), 'suffix', 'renamed.tsv']

        count_filename = '_'.join(count_filename_parts)

        counts = pd.read_csv(os.path.join(data_wd,folder_list[0],count_filename),sep='\t')

    else:
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', folder_list[i])

        count_filename_parts = [match.group(1), match.group(2), match.group(3), match.group(4), 'suffix', 'renamed.tsv']

        count_filename = '_'.join(count_filename_parts)

        this_count = pd.read_csv(os.path.join(data_wd,folder_list[i],count_filename),sep='\t')

        counts = counts.append(this_count)

# Output
out_filename = 'counts.txt'
counts.to_csv(out_filename, sep="\t", index=False)