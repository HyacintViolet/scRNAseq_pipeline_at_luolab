import os
import re
import pandas as pd
import shlex
# import logging
import multiprocessing as mp
import subprocess


def work(cmd):
    lib = re.search('(YT[0-9]*)', cmd).group(1)
    # Display progress
    print('Extracting library: ' + lib)
    return subprocess.call(cmd, shell=True)


def awk_extract(libs, parent_dir, num_lines_table):
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
            cmd_this_awk = 'awk \'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{print $1,$2,$3,$4,$14,$15,$16,$18,$23}\' ' + \
                           path_to_input + ' > ' + path_to_output
            cmd_awk_all.append(cmd_this_awk)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_awk_all) is not 0:
        pool.map(work, cmd_awk_all)
    print('awk extract desired fields: finished.')

    cmd_awk_extract_paste_all = []
    for l in libs:

        # Setting up input/output directory
        wd = os.path.join(parent_dir, l)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
        prefix = match.group(1)

        # Input file name
        filename_fixed_bed = '_'.join([prefix, 'fixed', 'closest.bed'])
        path_to_input = os.path.join(wd, filename_fixed_bed)
        filename_extracted_part = '_'.join([prefix, 'extracted', 'part.bed'])
        path_to_pasted_file = os.path.join(wd, filename_extracted_part)
        path_to_temp = os.path.join(wd, "tmp")

        # Output file name
        filename_output = '_'.join([prefix, 'extracted.bed'])
        path_to_output = os.path.join(wd, filename_output)

        # Check if output already exists. If not, construct command.
        if not os.path.exists(path_to_output):
            # Construct command
            cmd_this_extract_paste = 'awk \'BEGIN{FS="gene_name";OFS=\"\\t\"}{print $2}\' ' + path_to_input + \
                            ' | awk \'BEGIN{FS=";";OFS="\\t"}{print $1,"' + prefix + '"}\' | tr -d \'"\' | paste ' + \
                            path_to_pasted_file + ' - > ' + path_to_temp + ' && mv ' + path_to_temp + ' ' + \
                            path_to_output
            cmd_awk_extract_paste_all.append(cmd_this_extract_paste)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_awk_extract_paste_all) is not 0:
        pool.map(work, cmd_awk_extract_paste_all)
    print('awk extract gene_name: finished.')

    # Check line numbers are okay
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
        prefix = match.group(1)

        # File to check
        filename_extracted_bed = '_'.join([prefix, 'extracted.bed'])
        path_to_check = os.path.join(wd, filename_extracted_bed)

        # Construct command
        cmd_count_closest_bed = "wc -l " + path_to_check
        process = subprocess.Popen(
            shlex.split(cmd_count_closest_bed), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        while True:
            stdout, stderr = process.communicate()
            if process.poll() is not None:
                break
        num_line_extracted_bed = int(stdout.decode().split()[0])
        num_line_stranded_nonoverlap = int(num_lines_table.at[prefix, 'num_line_stranded_nonoverlap'])

        if num_line_extracted_bed == num_line_stranded_nonoverlap:
            print(prefix + '_extracted.bed: ' + str(num_line_extracted_bed) + ' lines, _stranded_nonoverlap.bam ' +
                  str(num_line_stranded_nonoverlap) + ' lines: OK')
        else:
            print(prefix + '_extracted.bed: ' + str(num_line_extracted_bed) + ' lines, _stranded_nonoverlap.bam ' +
                  str(num_line_stranded_nonoverlap) + ' lines: Something is wrong')
    print('Check .._extracted.bed line number: finished.')


def main():
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'
    libs = sorted(os.listdir(parent_dir))
    num_lines_table = pd.read_csv('/media/luolab/ZA1BT1ER/yanting/vM23/num_lines.csv', index_col='library')
    awk_extract(libs, parent_dir, num_lines_table)


if __name__ == '__main__':
    main()
