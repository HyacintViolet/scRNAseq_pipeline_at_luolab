"""
Concatenate suffix added counts (YTXXXXXX_counts_suffix.tsv) from individual library. Store them in a huge spread-sheet.
"""
# ----------------------------------------------------------------------------------------------------------------------
# Strategy 1: by recursion
# ----------------------------------------------------------------------------------------------------------------------

import os
import pandas as pd

# Set up working directories
parent_wd = '/media/luolab/ZA1BT1ER/yanting/vM21/'
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM21/mapping/'
os.chdir(parent_wd)
counts_new = pd.DataFrame(columns=['gene', 'cell', 'count'])


def navigate_and_concatenate_counts(src, counts):

    for item in os.listdir(src):
        s = os.path.join(src, item)
        if os.path.isdir(s):
            counts = navigate_and_concatenate_counts(s, counts)
        elif os.path.basename(s).endswith('counts_suffix.tsv'):
            print(s)
            counts = counts.append(pd.read_csv(s, sep='\t'), ignore_index=True)
            print(counts.shape)

    return counts


counts_concatenated = navigate_and_concatenate_counts(data_wd, counts_new)

# counts_concatenated.to_csv('counts_all.txt', sep='\t', index=False)

# ----------------------------------------------------------------------------------------------------------------------
# Un-dotting gene name, continued from above.
# ----------------------------------------------------------------------------------------------------------------------

# # Run this line if counts are not loaded.
# counts_concatenated = pd.read_csv('counts_all.txt', sep='\t')

# Iteration is too slow
# for index, row in counts_concatenated.iterrows():

ensmusg = counts_concatenated.gene.str.split(".", n=1, expand=True)

counts_concatenated.gene = ensmusg[0]

# Write out dot-less counts spread-sheet.
counts_concatenated.to_csv('counts_stbid.txt.gz', sep='\t', index=False, compression='gzip')


# # --------------------------------------------------------------------------------------------------------------------
# # Strategy 2: by iteration
# # --------------------------------------------------------------------------------------------------------------------
# import os
# import re
# import pandas as pd
#
#
# # Set up working directories
# parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
# data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4/vM4_CaiT/'
# os.chdir(parent_wd)
#
# folder_list = os.listdir(data_wd)
#
# # I know I know... Iteration is stupid... Let's fix this later.
# for i in range(len(folder_list)):
#     if i == 0:
#
#         match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder_list[0])
#
#         count_filename_parts = [match.group(1), 'counts', 'suffix.tsv']
#
#         count_filename = '_'.join(count_filename_parts)
#
#         counts = pd.read_csv(os.path.join(data_wd,folder_list[0],count_filename),sep='\t')
#
#     else:
#         match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder_list[i])
#
#         count_filename_parts = [match.group(1), 'counts', 'suffix.tsv']
#
#         count_filename = '_'.join(count_filename_parts)
#
#         this_count = pd.read_csv(os.path.join(data_wd,folder_list[i],count_filename),sep='\t')
#
#         counts = counts.append(this_count)
#
# # Output
# out_filename = 'counts.txt'
# counts.to_csv(out_filename, sep="\t", index=False)
