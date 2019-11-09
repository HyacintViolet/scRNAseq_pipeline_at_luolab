import os
import re
import shlex
import warnings
import pandas as pd
import multiprocessing as mp
import subprocess


def work(cmd):
    # Display progress
    lib = re.search('(YT[0-9]*)', cmd)
    if lib is not None:
        lib = lib.group(1)
        if '$1,$2,$3,$4' in cmd:
            print('Extracting part from fixed_closest.bed. Library: ' + lib)
        elif 'gene_type' in cmd:
            print('Extracting & pasting to part. Library: ' + lib)
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


def awk_extract_part(libs, parent_dir):
    cmd_awk_extract_part_all = []
    for l in libs:

        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

        # Input file name
        filename_fixed_bed = '_'.join([prefix, 'fixed', 'closest.bed'])
        path_to_input = os.path.join(wd, filename_fixed_bed)

        # Output file name
        filename_output = '_'.join([prefix, 'extracted', 'part.bed'])
        path_to_output = os.path.join(wd, filename_output)

        # Check if output already exists. If not, construct command.
        if not os.path.exists(path_to_output):
            # Construct command
            cmd_awk_extract_part_this = 'awk \'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{print $1,$2,$3,$4,$14,$15,$16,$18,' \
                                        '$23}\' ' + path_to_input + ' > ' + path_to_output
            cmd_awk_extract_part_all.append(cmd_awk_extract_part_this)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_awk_extract_part_all) is not 0:
        pool.map(work, cmd_awk_extract_part_all)
    print('awk extract desired fields: finished.')


def awk_extract_paste(libs, parent_dir):
    cmd_awk_extract_paste_all = []
    for l in libs:

        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

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
            cmd_this_extract_paste = 'awk \'BEGIN{FS="\\t"}{print $22}\' ' + path_to_input + ' | awk \'BEGIN{FS="; ";' \
                                     'OFS="\\t"}{ for(group=1;group<=NF;group++) if ($group~"gene_type"){gene_type=' \
                                     '$group} else if ($group~"gene_name"){gene_name=$group} print gene_name,' \
                                     'gene_type,"' + prefix + '"}\' | sed \'s/gene_type//g\' | sed \'s/gene_name//g\'' \
                                     ' | tr -d \'" \' | paste ' + path_to_pasted_file + ' - > ' + path_to_temp + ' &&' \
                                     ' mv ' + path_to_temp + ' ' + path_to_output
            cmd_awk_extract_paste_all.append(cmd_this_extract_paste)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_awk_extract_paste_all) is not 0:
        pool.map(work, cmd_awk_extract_paste_all)
    print('awk extract gene_name: finished.')


def remove_intermediate(libs, parent_dir, suffix_file_to_remove):
    # Remove intermediate files with 'find ... -delete'
    # Check that file numbers match
    num_libs = len(libs)
    # Count file number with pipe
    cmd_count_extracted_part_1 = 'find ' + parent_dir + ' -name "*' + suffix_file_to_remove + '" -type f'
    cmd_count_extracted_part_2 = 'wc -l'
    # Open process
    p1 = subprocess.Popen(shlex.split(cmd_count_extracted_part_1), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(shlex.split(cmd_count_extracted_part_2), stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    while True:
        output = p2.communicate()[0]
        if p2.poll() is not None:
            break
    num_count = int(output.decode())
    # Command to remove files
    cmd_rm_intermediate = ['find ' + parent_dir + ' -name "*' + suffix_file_to_remove + '" -type f -delete']
    pool = mp.Pool(1)
    if num_count == num_libs and len(cmd_rm_intermediate) is not 0:
        pool.map(work, cmd_rm_intermediate)
    else:
        warnings.warn('Trying to remove ..' + suffix_file_to_remove + ', but file number does not match with library '
                      'number. Abort.')


def check_line(libs, parent_dir, suffix_file_to_check, suffix_file_ground_truth, num_lines_table):
    # Check line numbers are okay
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

        # File to check
        filename_extracted_bed = '_'.join([prefix, suffix_file_to_check])
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
        num_line_fixed = int(num_lines_table.at[prefix, suffix_file_ground_truth])

        if num_line_extracted_bed == num_line_fixed:
            print(prefix + '_extracted.bed: ' + str(num_line_extracted_bed) + ' lines, fixed_closest.bed ' +
                  str(num_line_fixed) + ' lines: OK')
        else:
            print(prefix + '_extracted.bed: ' + str(num_line_extracted_bed) + ' lines, fixed_closest.bed ' +
                  str(num_line_fixed) + ' lines: Something is wrong')
    print('Check .._extracted.bed line number: finished.')


def main():
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'
    libs = get_libs(parent_dir)
    num_lines_table = pd.read_csv('/media/luolab/ZA1BT1ER/yanting/vM23/num_lines_table.csv', index_col='library')

    awk_extract_part(libs, parent_dir)
    awk_extract_paste(libs, parent_dir)

    remove_intermediate(libs, parent_dir, 'extracted_part.bed')
    # Check that line number of final output matches input
    check_line(libs, parent_dir, 'extracted.bed', 'fixed_closest.bed', num_lines_table)


if __name__ == '__main__':
    main()

# Old codes in awk_extract:
# cmd_this_extract_paste = 'awk \'BEGIN{FS="gene_name";OFS=\"\\t\"}{print $2}\' ' + path_to_input + \
#                          ' | awk \'BEGIN{FS=";";OFS="\\t"}{print $1,"' + prefix + '"}\' | tr -d \'" \' | paste ' + \
#                          path_to_pasted_file + ' - > ' + path_to_temp + ' && mv ' + path_to_temp + ' ' + \
#                          path_to_output
