import os
import re
# import logging
import multiprocessing as mp
import subprocess


parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'
gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr.annotation.gtf'
bed = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr.annotation.fixed.bed'


def work(cmd):
    lib = re.search('(YT[0-9]*)', cmd).group(1)
    if 'intersect' in cmd:
        print('Processing: bedtools intersect... library:' + lib)
    elif 'closest' in cmd:
        print('Processing: bedtools closest... library:' + lib)
    return subprocess.call(cmd, shell=True)


cmd_all = []
for s in os.listdir(parent_dir):

    # Setting up input/output directory
    wd = os.path.join(parent_dir, s)
    # out_dir = os.path.join(dst, s)

    # Grab folder name prefix
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', s)
    prefix = match.group(1)

    # Input file name
    file_A = '_'.join([prefix, 'assigned_sorted.bam'])
    path_to_file_A = os.path.join(wd, file_A)

    # Output file name
    out_name_intersect = '_'.join([prefix, 'stranded_nonoverlap.bam'])
    path_to_intersect_out = os.path.join(wd, out_name_intersect)

    out_name_closest = '_'.join([prefix, 'closest.bed'])
    path_to_closest_out = os.path.join(wd, out_name_closest)

    # Check if output already exists. If not, construct command.
    if not os.path.exists(path_to_intersect_out):
        # Construct command
        cmd_this_intersect = 'bedtools intersect -v -s -a ' + path_to_file_A + ' -b ' + \
                             gtf + ' > ' + path_to_intersect_out
        cmd_all.append(cmd_this_intersect)
        cmd_this_intersect = []

    if not os.path.exists(path_to_closest_out):
        cmd_this_closest = 'bedtools closest -s -D a -io -t first -a ' + path_to_intersect_out + \
                           ' -b ' + bed + ' > ' + path_to_closest_out
        cmd_all.append(cmd_this_closest)
        cmd_this_closest = []


# Parallel run by Pool
pool = mp.Pool(1)
pool.map(work, cmd_all)
print('bedtools intersect & closest: finished.')
