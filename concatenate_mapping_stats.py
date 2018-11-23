# Same operation as concatenate_counts.py
import os, re
import pandas as pd

# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def/'
os.chdir(parent_wd)

folder_list = os.listdir(data_wd)

# I know I know... Iteration is stupid... Will fix this later.
for i in range(len(folder_list)):
    if i == 0:

        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', folder_list[0])

        mapstats_filename_parts = [match.group(1), match.group(2), match.group(3), match.group(4), 'mapping', 'stats.txt']

        mapstats_filename = '_'.join(mapstats_filename_parts)

        mapping_stats = pd.read_csv(os.path.join(data_wd, folder_list[0], mapstats_filename), sep='\t')

    else:
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', folder_list[i])

        mapstats_filename_parts = [match.group(1), match.group(2), match.group(3), match.group(4), 'mapping', 'stats.txt']

        mapstats_filename = '_'.join(mapstats_filename_parts)

        this_mapstats = pd.read_csv(os.path.join(data_wd, folder_list[i], mapstats_filename), sep='\t')

        mapping_stats = mapping_stats.append(this_mapstats)

# Output
out_filename = 'mapping_stats.txt'
mapping_stats.to_csv(out_filename, sep="\t", index=False)