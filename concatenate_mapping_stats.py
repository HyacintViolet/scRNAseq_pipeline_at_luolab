# Same operation as concatenate_counts.py

# ----------------------------------------------------------------------------------------------------------------------
# Strategy 1: by recursion
# ----------------------------------------------------------------------------------------------------------------------

import os
import pandas as pd

# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/linrui/DR_DAT/'
data_wd = '/media/luolab/ZA1BT1ER/linrui/DR_DAT/vM19/'
os.chdir(parent_wd)
mapping_stats_new = pd.DataFrame(columns=['Nreads', 'Nuniquemap', 'cell'])


def navigate_and_concatenate_mapping_stats(src, mapping_stats):

    for item in os.listdir(src):
        s = os.path.join(src, item)
        if os.path.isdir(s):
            mapping_stats = navigate_and_concatenate_mapping_stats(s, mapping_stats)
        elif os.path.basename(s).endswith('mapping_stats.txt'):
            print(s)
            mapping_stats = mapping_stats.append(pd.read_csv(s, sep='\t'), ignore_index=True)
            print(mapping_stats.shape)

    return mapping_stats


mapping_stats_concatenated = navigate_and_concatenate_mapping_stats(data_wd, mapping_stats_new)

mapping_stats_concatenated.to_csv('mapping_stats_all.txt', sep='\t', index=False)
#
# # --------------------------------------------------------------------------------------------------------------------
# # Strategy 2: by iteration
# # --------------------------------------------------------------------------------------------------------------------
# import os
# import re
# import pandas as pd
#
# # Set up working directories
# parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
# data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def_2/'
# os.chdir(parent_wd)
#
# folder_list = os.listdir(data_wd)
#
# # I know I know... Iteration is stupid... Will fix this later.
# for i in range(len(folder_list)):
#     if i == 0:
#
#         match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder_list[0])
#
#         mapstats_filename_parts = [match.group(1), 'mapping', 'stats.txt']
#
#         mapstats_filename = '_'.join(mapstats_filename_parts)
#
#         mapping_stats = pd.read_csv(os.path.join(data_wd, folder_list[0], mapstats_filename), sep='\t')
#
#     else:
#         match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder_list[i])
#
#         mapstats_filename_parts = [match.group(1), 'mapping', 'stats.txt']
#
#         mapstats_filename = '_'.join(mapstats_filename_parts)
#
#         this_mapstats = pd.read_csv(os.path.join(data_wd, folder_list[i], mapstats_filename), sep='\t')
#
#         mapping_stats = mapping_stats.append(this_mapstats)
#
# # Output
# out_filename = 'mapping_stats.txt'
# mapping_stats.to_csv(out_filename, sep="\t", index=False)
