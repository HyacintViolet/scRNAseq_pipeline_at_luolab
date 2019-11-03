import os
import re
# import logging
import multiprocessing as mp
import subprocess


def work(cmd):
    lib = re.search('(YT[0-9]*)', cmd).group(1)
    # Display progress
    print('Extracting library: ' + lib)
    return subprocess.call(cmd, shell=True)


def awk_extract(libs, parent_dir):
    cmd_awk_all = []
    for l in libs:

        # Setting up input/output directory
        wd = os.path.join(parent_dir, l)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
        prefix = match.group(1)

        # Input file name
        filename_fixed_bed = '_'.join([prefix, 'fixed', 'closest.bed'])
        path_to_input = os.path.join(wd, filename_fixed_bed)

        # Output file name
        filename_output = '_'.join([prefix, 'extracted', 'part.bed'])
        path_to_output = os.path.join(wd, filename_output)

        # Check if output already exists. If not, construct command.
        if not os.path.exists(path_to_output):
            # Construct command
            cmd_this_awk = 'awk \'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{print $1,$2,$3,$4,$14,$15,$16,$23}\' ' + \
                           path_to_input + ' > ' + path_to_output
            cmd_awk_all.append(cmd_this_awk)

    # Parallel run by Pool
    pool = mp.Pool(1)
    pool.map(work, cmd_awk_all)
    print('awk extract desired fields: finished.')


def main():
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'
    libs = sorted(os.listdir(parent_dir))
    awk_extract(libs, parent_dir)


if __name__ == '__main__':
    main()
