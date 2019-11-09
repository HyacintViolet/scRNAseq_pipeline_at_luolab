import os
import re
import shlex
import warnings
# import numpy as np
import pandas as pd
import multiprocessing as mp
import subprocess


def work(cmd):
    lib_name = re.search('(YT[0-9]*)', cmd)
    # Display progress
    if lib_name is not None:
        lib_name = lib_name.group(1)
        if 'head -n' in cmd:
            print('Trimming gibberish tails. head -n ... Library:' + lib_name)
        elif '$14!=-1' in cmd:
            print('Removing telomeric reads. awk ... $14!=-1 ... Library:' + lib_name)
        elif '$18=="+" && $23>0' in cmd:
            print('Removing reads upstream to 5\'UTRs. Library:' + lib_name)
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
        p1 = subprocess.Popen(shlex.split(cmd_expand_file), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(shlex.split(cmd_count), stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()

        # Wait for process to finish and then read stdout
        while True:
            output = p2.communicate()[0]
            if p2.poll() is not None:
                break
        num_lines = int(output.decode())
        counts.at[prefix, file_to_count] = num_lines
        # Display progress
        print(str(num_lines) + ' lines.')
    return counts


def count_file_type(ftype="", parent_dir='/media/luolab/ZA1BT1ER/yanting/vM23/mapping/',
                    num_lines_table=pd.DataFrame(data=None), do_count=True):
    # ftype: 0 = all; 1 = _closest.bed; 2 = _stranded_nonoverlap; 3 = fixed_closest.bed
    valid_types = {"closest.bed", "stranded_nonoverlap.bam", "fixed_closest.bed"}
    if ftype not in valid_types:
        raise ValueError("results: ftype must be one of %r." % valid_types)

    # Top level control of whether to count or not
    if do_count:
        if ftype == "closest.bed":
            file_to_count = "closest.bed"
            counts = count_files(file_to_count, parent_dir)
            num_lines_table.update(counts)
        elif ftype == "stranded_nonoverlap.bam":
            file_to_count = "stranded_nonoverlap.bam"
            counts = count_files(file_to_count, parent_dir)
            num_lines_table.update(counts)
        elif ftype == "fixed_closest.bed":
            file_to_count = "fixed_closest.bed"
            counts = count_files(file_to_count, parent_dir)
            num_lines_table.update(counts)

    return num_lines_table


def trim_bed_tails(parent_dir, num_lines_table):
    libs = get_libs(parent_dir)
    cmd_all_head = []
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

        # File to fix
        filename_input = '_'.join([prefix, 'closest.bed'])   # 'closest.bed'
        path_to_input = os.path.join(wd, filename_input)
        # Number of lines to correct
        tail_to_trim = num_lines_table.at[prefix, 'tail_to_trim']

        # Output file
        filename_output = '_'.join([prefix, 'temp', 'closest.bed'])   # 'temp_closest.bed'
        path_to_output = os.path.join(wd, filename_output)

        if not os.path.exists(path_to_output) and tail_to_trim != 0:
            # Construct command
            cmd_this_head = 'head -n -' + str(tail_to_trim) + ' ' + path_to_input + ' > ' + path_to_output
            cmd_all_head.append(cmd_this_head)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_all_head) is not 0:
        pool.map(work, cmd_all_head)
    print('Trim YT..._closest.bed tail lines: finished.')


def remove_telomeric_reads(libs, parent_dir, suffix_input, suffix_output):
    cmd_all_awk = []
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

        # File to fix
        filename_input = '_'.join([prefix, suffix_input])  # 'temp_closest.bed'
        path_to_input = os.path.join(wd, filename_input)

        # Output file
        filename_output = '_'.join([prefix, suffix_output])  # 'temp2_closest.bed'
        path_to_output = os.path.join(wd, filename_output)

        if not os.path.exists(path_to_output):
            # Construct command
            cmd_this_awk = 'awk \'BEGIN{FS=\"\\t\"} $14!=-1 {print $0}\' ' + path_to_input + ' > ' + path_to_output
            cmd_all_awk.append(cmd_this_awk)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_all_awk) is not 0:
        pool.map(work, cmd_all_awk)
    print('Remove YT..._closest.bed telomeric reads: finished.')


def remove_5_prime_upstream(libs, parent_dir, suffix_input, suffix_output):
    # We're only interested in reads mapped to downstream regions of 3'UTRs. Remove reads mapped to upstream regions
    # of 5'UTRs.
    cmd_all_awk = []
    for l in libs:
        # Setting up working directory
        wd = os.path.join(parent_dir, l)
        # Grab folder name prefix
        prefix = get_prefix(l)

        # File to fix
        filename_input = '_'.join([prefix, suffix_input])   # 'temp2_closest.bed'
        path_to_input = os.path.join(wd, filename_input)

        # Output file
        filename_output = '_'.join([prefix, suffix_output])   # 'fixed_closest.bed'
        path_to_output = os.path.join(wd, filename_output)

        if not os.path.exists(path_to_output):
            # Construct command
            cmd_this_awk = 'awk \'BEGIN{FS="\\t";OFS="\\t"}{if($18=="+" && $23>0){print $0}else if($18=="-" && $23<0)' \
                           '{print $0}}\' ' + path_to_input + ' > ' + path_to_output
            cmd_all_awk.append(cmd_this_awk)

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_all_awk) is not 0:
        pool.map(work, cmd_all_awk)
    print('Remove YT..._closest.bed telomeric entries: finished.')


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


def main():
    analysis_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/'
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'

    libs = get_libs(parent_dir)
    path_to_num_lines_table = os.path.join(analysis_dir, "num_lines_table.csv")

    if os.path.exists(path_to_num_lines_table):
        # Read table if it already exists
        num_lines_table = pd.read_csv(path_to_num_lines_table, index_col='library')
    else:
        # Get index
        idx = get_idx(libs)
        # Initialize empty DataFrame
        num_lines_table = pd.DataFrame(data=None, index=idx, columns=["closest.bed", "stranded_nonoverlap.bam",
                                                                      "tail_to_trim", "fixed_closest.bed"])

    # Count YT..._closest.bed line numbers
    num_lines_new = count_file_type(ftype="closest.bed", parent_dir=parent_dir, num_lines_table=num_lines_table,
                                    do_count=False)  # Don't count, ad interim.
    num_lines_table.update(num_lines_new)
    num_lines_table.to_csv('/media/luolab/ZA1BT1ER/yanting/vM23/num_lines_table.csv', index_label='library')

    # Count YT..._stranded_nonoverlap.bam line numbers
    num_lines_new = count_file_type(ftype="stranded_nonoverlap.bam", parent_dir=parent_dir,
                                    num_lines_table=num_lines_table, do_count=False)  # Don't count, ad interim.
    num_lines_table.update(num_lines_new)
    num_lines_table.to_csv('/media/luolab/ZA1BT1ER/yanting/vM23/num_lines_table.csv', index_label='library')

    # Trim tail lines: YT..._closest.bed
    num_lines_table['tail_to_trim'] = num_lines_table['closest.bed'] - num_lines_table['stranded_nonoverlap.bam'] + 1
    trim_bed_tails(parent_dir, num_lines_table)

    # Remove some useless reads mapped to the telomeric region of chromosomes
    # remove_telomeric_reads(libs, parent_dir, suffix_input='temp_closest.bed', suffix_output='temp2_closest.bed')  # Don't execute, ad interim.
    # Remove intermediate files ending with _temp_closest.bed
    remove_intermediate(libs, parent_dir, suffix_file_to_remove='temp_closest.bed')

    # Remove reads mapped to upstream regions of 5'UTRs.
    remove_5_prime_upstream(libs, parent_dir, suffix_input='temp2_closest.bed', suffix_output='fixed.closest.bed')
    # Remove intermediate files: .._temp2_closest.bed
    remove_intermediate(libs, parent_dir, suffix_file_to_remove='temp2_closest.bed')

    # Count YT..._fixed_closest.bed line numbers
    num_lines_new = count_file_type(ftype="fixed_closest.bed", parent_dir=parent_dir, num_lines_table=num_lines_table)
    num_lines_table.update(num_lines_new)
    # Save to num_lines_table.csv
    num_lines_table.to_csv('/media/luolab/ZA1BT1ER/yanting/vM23/num_lines_table.csv', index_label='library')


if __name__ == '__main__':
    main()

    # # Check line numbers are okay
    # for l in libs:
    #     # Setting up input/output directory
    #     wd = os.path.join(parent_dir, l)
    #
    #     # Grab folder name prefix
    #     match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l)
    #     prefix = match.group(1)
    #
    #     # File to check
    #     filename_fixed_bed = '_'.join([prefix, 'fixed_closest.bed'])
    #     path_to_check = os.path.join(wd, filename_fixed_bed)
    #
    #     # Construct command
    #     cmd_count_closest_bed = "wc -l " + path_to_check
    #     process = subprocess.Popen(
    #         shlex.split(cmd_count_closest_bed), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    #     )
    #     while True:
    #         stdout, stderr = process.communicate()
    #         if process.poll() is not None:
    #             break
    #     num_line_closest_bed = int(stdout.decode().split()[0])
    #     num_line_stranded_nonoverlap = int(num_lines_table.at[prefix, 'num_line_stranded_nonoverlap'])
    #
    #     if num_line_closest_bed == num_line_stranded_nonoverlap:
    #         print(prefix + '_fixed_closest.bed: ' + str(num_line_closest_bed) + ' lines, _stranded_nonoverlap.bam ' +
    #               str(num_line_stranded_nonoverlap) + ' lines: OK')
    #     else:
    #         print(prefix + '_fixed_closest.bed: ' + str(num_line_closest_bed) + ' lines, _stranded_nonoverlap.bam ' +
    #               str(num_line_stranded_nonoverlap) + ' lines: Something is wrong')

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
