import os
import re
# import logging
import multiprocessing as mp
import subprocess


parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'
gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr.annotation.gtf'


def work(cmd):
    lib = re.search('(YT[0-9]*)', cmd).group(1)
    print('Parsing ' + lib)
    return subprocess.call(cmd, shell=True)


cmd_intersect = []
for s in os.listdir(parent_dir):

    # Setting up input/output directory
    wd = os.path.join(parent_dir, s)
    # out_dir = os.path.join(dst, s)

    # Grab folder name prefix
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', s)
    prefix = match.group(1)

    # Output file name
    out_name_intersect = '_'.join([prefix, 'Unassigned.bam'])
    path_to_intersect_out = os.path.join(wd, out_name_intersect)

    # Input file name
    file_A = '_'.join([prefix, 'assigned_sorted.bam'])
    path_to_file_A = os.path.join(wd, file_A)

    # Check if output already exists. If not, construct command.
    if not os.path.exists(path_to_intersect_out):
        # Construct command
        cmd_this_intersect = 'bedtools intersect -v -a ' + path_to_file_A + ' -b ' + \
                             gtf + ' > ' + path_to_intersect_out
        cmd_intersect.append(cmd_this_intersect)


# Parallel run by Pool
pool = mp.Pool(16)
pool.map(work, cmd_intersect)
print('bedtools intersect unassigned: finished.')