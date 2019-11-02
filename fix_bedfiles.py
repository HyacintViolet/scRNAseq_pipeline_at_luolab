import os
import re
# import numpy as np
import pandas as pd
import multiprocessing as mp
import subprocess


def work(cmd):
    lib_name = re.search('(YT[0-9]*)', cmd).group(1)
    # Display progress
    print('Fixing: head -n -... library:' + lib_name)
    return subprocess.call(cmd, shell=True)


parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/'
mapping_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'

num_lines_table = pd.read_csv(parent_dir + 'num_lines.csv', index_col='library')

libs = sorted(os.listdir(mapping_dir))
cmd_all_head = []
for l in libs:
    # Setting up working directory
    wd = os.path.join(mapping_dir, l)

    # Grab folder name prefix
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
    prefix = match.group(1)

    # File to correct
    filename_closest_bed = '_'.join([prefix, 'closest.bed'])
    path_to_input = os.path.join(wd, filename_closest_bed)
    # Number of lines to correct
    sub = num_lines_table.at[prefix, 'sub']

    # Output file
    filename_output = '_'.join([prefix, 'fixed_closest.bed'])
    path_to_output = os.path.join(wd, filename_output)

    if not os.path.exists(path_to_output):
        # Construct command
        cmd_this_head = 'head -n -' + str(sub) + ' ' + path_to_input + ' > ' + path_to_output
        cmd_all_head.append(cmd_this_head)

# Parallel run by Pool
pool = mp.Pool(32)
pool.map(work, cmd_all_head)
print('Correct YT..._closest.bed tail lines: finished.')

# Check line numbers are okay
for l in libs:
    # Setting up input/output directory
    wd = os.path.join(mapping_dir, l)

    # Grab folder name prefix
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
    prefix = match.group(1)

    # File to check
    filename_fixed_bed = '_'.join([prefix, 'fixed_closest.bed'])
    path_to_check = os.path.join(wd, filename_fixed_bed)

    # Construct command
    cmd_count_closest_bed = ['wc', '-l', path_to_check]
    process = subprocess.Popen(
        cmd_count_closest_bed, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    num_line_closest_bed = stdout.decode().split()[0]

    if num_line_closest_bed == num_lines_table.at[prefix, 'num_line_stranded_nonoverlap']:
        print(prefix + ': OK')
    else:
        print(prefix + ': Something is wrong')

# parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'
# gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr.annotation.gtf'
# bed = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr.annotation.fixed.bed'
#
#
# def work(cmd):
#     lib = re.search('(YT[0-9]*)', cmd).group(1)
#     if 'intersect' in cmd:
#         print('Processing: bedtools intersect... library=' + lib)
#     elif 'closest' in cmd:
#         print('Processing: bedtools closest... library=' + lib)
#     return subprocess.call(cmd, shell=True)
#
#
# libs = sorted(os.listdir(parent_dir))
# num_libs = len(libs)
#
# # Get index for dataframe
# idx = []
# for l in libs:
#     match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
#     prefix = match.group(1)
#     idx.append(prefix)
# # Initialize empty dataframe
# df_line_counts = pd.DataFrame(data=None, index=idx,
#                               columns=['num_line_closest', 'num_line_stranded_nonoverlap'], dtype=int)
#
# cmd_all = []
# for s in libs:
#
#     # Setting up input/output directory
#     wd = os.path.join(parent_dir, s)
#
#     # Grab folder name prefix
#     match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', s)
#     prefix = match.group(1)
#     # Display progress
#     print('Parsing lib:' + prefix)
#
#     # File to count
#     closest_bed = '_'.join([prefix, 'closest.bed'])
#     path_to_closest_bed = os.path.join(wd, closest_bed)
#
#     stranded_nonoverlap = '_'.join([prefix, 'stranded_nonoverlap.bam'])
#     path_to_stranded_nonoverlap = os.path.join(wd, stranded_nonoverlap)
#
#     # Initialize some variable to store line num
#     num_line_closest_bed = 0
#     num_line_stranded_nonoverlap = 0
#     # Count _closest.bed line number and assign to DataFrame
#     # Check if input file exists.
#     if os.path.exists(path_to_closest_bed):
#         # Construct command
#         cmd_count_closest_bed = ['wc', '-l', path_to_closest_bed]
#         process = subprocess.Popen(
#             cmd_count_closest_bed, stdout=subprocess.PIPE, stderr=subprocess.PIPE
#         )
#         stdout, stderr = process.communicate()
#         num_line_closest_bed = stdout.decode().split()[0]
#         if num_line_closest_bed == 0:
#             warnings.warn('Library ' + prefix + ' parsed line number: 0')
#         df_line_counts.at[prefix, 'num_line_closest'] = num_line_closest_bed
#     else:
#         warnings.warn(path_to_closest_bed + 'doesn\'t exists')
#
#     # Count _stranded_nonoverlap.bam line number and assign to DataFrame
#     if not os.path.exists(path_to_stranded_nonoverlap):
#         # Construct command
#         cmd_count_closest_bed = ['wc', '-l', path_to_closest_bed]
#         process = subprocess.Popen(
#             cmd_count_closest_bed, stdout=subprocess.PIPE, stderr=subprocess.PIPE
#         )
#         stdout, stderr = process.communicate()
#         num_line_closest_bed = stdout.decode().split()[0]
#         if num_line_closest_bed == 0:
#             warnings.warn('Library ' + prefix + ' parsed line number: 0')
#         df_line_counts.at[prefix, 'num_line_closest'] = num_line_closest_bed
#     else:
#         warnings.warn(path_to_stranded_nonoverlap + 'doesn\'t exists')
