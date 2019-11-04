import os
import re
import shlex
# import numpy as np
import pandas as pd
import multiprocessing as mp
from subprocess import Popen, PIPE


def work(cmd):
    lib_name = re.search('(YT[0-9]*)', cmd).group(1)
    # Display progress
    print('Fixing: head -n -... library:' + lib_name)
    return subprocess.call(cmd, shell=True)


def get_libs(parent_dir):
    # Input path to mapping dir. Output a list of library names.
    libs = sorted(os.listdir(parent_dir))
    return libs


def get_prefix(lib):
    # Input library folder name. Output library prefix.
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', lib)
    prefix = match.group(1)
    return prefix


def get_idx(libs):
    idx = []
    for l in libs:
        i = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l).group(1)
        idx.append(i)
    return idx


def count_files(file_to_count, parent_dir):
    libs = get_libs(parent_dir)
    idx = get_idx(libs)
    counts = pd.DataFrame(data=None, index=idx, columns=[file_to_count], dtype=int)
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

        # File to count
        filename_count = '_'.join([prefix, file_to_count])
        path_to_file_count = os.path.join(wd, filename_count)
        # Display progress
        print('Counting ' + filename_count + ':')

        # Construct command
        if filename_count.endswith('.bam'):
            cmd_expand_file = 'samtools view ' + path_to_file_count
        else:
            cmd_expand_file = "cat " + path_to_file_count
        cmd_count = "wc -l"

        # Open process
        p1 = Popen(shlex.split(cmd_expand_file), stdout=PIPE, stderr=PIPE)
        p2 = Popen(shlex.split(cmd_count), stdin=p1.stdout, stdout=PIPE)
        p1.stdout.close()

        # Wait for process to finish and then read stdout
        while True:
            output = p2.communicate()[0]
            if p2.poll() is not None:
                break
        counts.at[prefix, file_to_count] = int(output.decode())
        # Display progress
        print(str(output.decode()) + 'lines.')
    return counts


def count_file_type(ftype="", parent_dir='/media/luolab/ZA1BT1ER/yanting/vM23/mapping/',
                    num_lines_table=pd.DataFrame(data=None)):
    # ftype: 0 = all; 1 = _closest.bed; 2 = _stranded_nonoverlap; 3 = fixed_closest.bed
    valid_types = {1, "closest.bed", 2, "stranded_nonoverlap.bam", 3, "fixed_closest.bed"}
    if ftype not in valid_types:
        raise ValueError("results: ftype must be one of %r." % valid_types)

    if ftype == 1 or ftype == "closest.bed":
        file_to_count = "closest.bed"
        counts = count_files(file_to_count, parent_dir)
        num_lines_table.update(counts)
    elif ftype == 2 or ftype == "stranded_nonoverlap.bam":
        file_to_count = "stranded_nonoverlap.bam"
        counts = count_files(file_to_count, parent_dir)
        num_lines_table.update(counts)
    elif ftype == 3 or ftype == "fixed_closest.bed":
        file_to_count = "fixed_closest.bed"
        counts = count_files(file_to_count, parent_dir)
        num_lines_table.update(counts)

    return num_lines_table


def fix_bed_tails(parent_dir, num_lines_table):
    libs = get_libs(parent_dir)
    cmd_all_head = []
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
        prefix = match.group(1)

        # File to fix
        filename_closest_bed = '_'.join([prefix, 'closest.bed'])
        path_to_input = os.path.join(wd, filename_closest_bed)
        # Number of lines to correct
        sub = num_lines_table.at[prefix, 'sub']

        # Output file
        filename_output = '_'.join([prefix, 'fixed_closest.bed'])
        path_to_output = os.path.join(wd, filename_output)

        if not os.path.exists(path_to_output) and sub != 0:
            # Construct command
            cmd_this_head = 'head -n -' + str(sub) + ' ' + path_to_input + ' > ' + path_to_output
            cmd_all_head.append(cmd_this_head)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_all_head) is not 0:
        pool.map(work, cmd_all_head)
    print('Correct YT..._closest.bed tail lines: finished.')


def main():
    analysis_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/'
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'

    libs = get_libs(parent_dir)
    path_to_num_lines_table = os.path.join(analysis_dir, "num_lines.csv")

    if os.path.exists(path_to_num_lines_table):
        # Read table if it already exists
        num_lines_table = pd.read_csv(path_to_num_lines_table, index_col='library')
    else:
        # Get index
        idx = get_idx(libs)
        # Initialize empty DataFrame
        num_lines_table = pd.DataFrame(data=None, index=idx, columns=["closest.bed",
                                                                      "stranded_nonoverlap.bam",
                                                                      "sub", "fixed_closest.bed"])
    num_lines_new = count_file_type(ftype="", parent_dir=parent_dir, num_lines_table=num_lines_table)
    num_lines_table.update(num_lines_new)

    num_lines_table.to_csv('/media/luolab/ZA1BT1ER/yanting/vM23/file.csv', index_label='library')

    cmd_all_head = []
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)

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

        if not os.path.exists(path_to_output) and sub != 0:
            # Construct command
            cmd_this_head = 'head -n -' + str(sub) + ' ' + path_to_input + ' > ' + path_to_output
            cmd_all_head.append(cmd_this_head)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_all_head) is not 0:
        pool.map(work, cmd_all_head)
    print('Correct YT..._closest.bed tail lines: finished.')

    # Check line numbers are okay
    for l in libs:
        # Setting up input/output directory
        wd = os.path.join(parent_dir, l)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
        prefix = match.group(1)

        # File to check
        filename_fixed_bed = '_'.join([prefix, 'fixed_closest.bed'])
        path_to_check = os.path.join(wd, filename_fixed_bed)

        # Construct command
        cmd_count_closest_bed = "wc -l " + path_to_check
        process = subprocess.Popen(
            shlex.split(cmd_count_closest_bed), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        while True:
            stdout, stderr = process.communicate()
            if process.poll() is not None:
                break
        num_line_closest_bed = int(stdout.decode().split()[0])
        num_line_stranded_nonoverlap = int(num_lines_table.at[prefix, 'num_line_stranded_nonoverlap'])

        if num_line_closest_bed == num_line_stranded_nonoverlap:
            print(prefix + '_fixed_closest.bed: ' + str(num_line_closest_bed) + ' lines, _stranded_nonoverlap.bam ' +
                  str(num_line_stranded_nonoverlap) + ' lines: OK')
        else:
            print(prefix + '_fixed_closest.bed: ' + str(num_line_closest_bed) + ' lines, _stranded_nonoverlap.bam ' +
                  str(num_line_stranded_nonoverlap) + ' lines: Something is wrong')

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
