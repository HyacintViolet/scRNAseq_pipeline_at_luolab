import os
import re
import shlex
import warnings
# import numpy as np
import pandas as pd
import multiprocessing as mp
import subprocess


def work(cmd_this):
    index = cmd_this[0][0]
    gene_name = cmd_this[0][1]
    length_of_extension = cmd_this[0][2]
    cmd = cmd_this[1]
    print(' '.join([str(index), 'Parsing', gene_name, '... Extend by', str(length_of_extension)]))
    return subprocess.call(cmd, shell=True)


def extend_gtf_coordinates(parent_dir, filename_extension_profile, filename_gtf):

    path_to_gtf = os.path.join(parent_dir, filename_gtf)

    path_to_extension_profile = os.path.join(parent_dir, filename_extension_profile)
    extension_profile = pd.read_table(path_to_extension_profile, sep="\t")

    cmd_all = []
    for index, row in extension_profile.iterrows():
        gene_id = row['gene_id']
        gene_name = row['gene_name']
        strand_this_gene = row['strand']
        coord_to_extend = int(row['coord_to_extend'])
        coord_after_extend = int(row['coord_after_extend'])
        length_of_extension = int(row['length_of_extension'])

        if length_of_extension is 0:
            continue
        else:
            if strand_this_gene == "+":
                cmd = 'awk \'BEGIN{FS="\\t";OFS="\\t"}{if($9~"' + gene_id + '" && ($3=="gene" || $3=="exon") &&' \
                           ' $5==' + str(coord_to_extend) + '){print $1,$2,$3,$4,' + str(coord_after_extend) + \
                           ',$6,$7,$8,$9}else{print $0}}\' ' + path_to_gtf + ' > tmp && mv tmp ' + path_to_gtf
                cmd_all.append([(index, gene_name, length_of_extension), cmd])
            elif strand_this_gene == "-":
                cmd = 'awk \'BEGIN{FS="\\t";OFS="\\t"}{if($9~"' + gene_id + '" && ($3=="gene" || $3=="exon") &&' \
                           ' $4==' + str(coord_to_extend) + '){print $1,$2,$3,' + str(coord_after_extend) + \
                           ',$5,$6,$7,$8,$9}else{print $0}}\' ' + path_to_gtf + ' > tmp && mv tmp ' + path_to_gtf
                cmd_all.append([(index, gene_name, length_of_extension), cmd])

    # Parallel run by Pool
    pool = mp.Pool(1)
    if len(cmd_all) is not 0:
        pool.map(work, cmd_all)
    print('Extend gtf coordinates: finished.\n')


def main():
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM21/'

    filename_extension_profile = 'extension_profiles.txt'

    filename_gtf = 'gencode.vM21.chr_patch_hapl_scaff.annotation.extended.gtf'

    extend_gtf_coordinates(parent_dir, filename_extension_profile, filename_gtf)


if __name__ == '__main__':
    main()

# def get_libs(parent_dir):
#     # Input path to mapping dir. Output a list of library names.
#     libs = sorted(os.listdir(parent_dir))
#     return libs
#
#
# def get_prefix(lib):
#     # Input library folder name. Output library prefix.
#     match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', lib)
#     prefix = match.group(1)
#     return prefix
#
#
# def get_idx(libs):
#     idx = []
#     for l in libs:
#         i = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', l).group(1)
#         idx.append(i)
#     return idx
#
#
# def count_files(file_to_count, parent_dir):
#     libs = get_libs(parent_dir)
#     idx = get_idx(libs)
#     counts = pd.DataFrame(data=None, index=idx, columns=[file_to_count], dtype=int)
#     for l in libs:
#         # Setting up working directory
#         wd = os.path.join(parent_dir, l)
#         # Grab folder name prefix
#         prefix = get_prefix(l)
#
#         # File to count
#         filename_count = '_'.join([prefix, file_to_count])
#         path_to_file_count = os.path.join(wd, filename_count)
#         # Display progress
#         print('Counting ' + filename_count + ':')
#
#         # Construct command
#         if filename_count.endswith('.bam'):
#             cmd_expand_file = 'samtools view ' + path_to_file_count
#         else:
#             cmd_expand_file = "cat " + path_to_file_count
#         cmd_count = "wc -l"
#
#         # Open process
#         p1 = subprocess.Popen(shlex.split(cmd_expand_file), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         p2 = subprocess.Popen(shlex.split(cmd_count), stdin=p1.stdout, stdout=subprocess.PIPE)
#         p1.stdout.close()
#
#         # Wait for process to finish and then read stdout
#         while True:
#             output = p2.communicate()[0]
#             if p2.poll() is not None:
#                 break
#         num_lines = int(output.decode())
#         counts.at[prefix, file_to_count] = num_lines
#         # Display progress
#         print(str(num_lines) + ' lines.')
#     return counts
#
#
# def count_file_type(ftype="", parent_dir='/media/luolab/ZA1BT1ER/yanting/vM23/mapping/',
#                     num_lines_table=pd.DataFrame(data=None), do_count=True):
#     # ftype: 0 = all; 1 = _closest.bed; 2 = _stranded_nonoverlap; 3 = fixed_closest.bed
#     valid_types = {"closest.bed", "stranded_nonoverlap.bam", "fixed_closest.bed"}
#     if ftype not in valid_types:
#         raise ValueError("results: ftype must be one of %r." % valid_types)
#
#     # Top level control of whether to count or not
#     if do_count:
#         if ftype == "closest.bed":
#             file_to_count = "closest.bed"
#             counts = count_files(file_to_count, parent_dir)
#             num_lines_table.update(counts)
#         elif ftype == "stranded_nonoverlap.bam":
#             file_to_count = "stranded_nonoverlap.bam"
#             counts = count_files(file_to_count, parent_dir)
#             num_lines_table.update(counts)
#         elif ftype == "fixed_closest.bed":
#             file_to_count = "fixed_closest.bed"
#             counts = count_files(file_to_count, parent_dir)
#             num_lines_table.update(counts)
#
#     return num_lines_table
#
#
# def trim_bed_tails(parent_dir, num_lines_table):
#     libs = get_libs(parent_dir)
#     cmd_all_head = []
#     for l in libs:
#         # Setting up working directory
#         wd = os.path.join(parent_dir, l)
#         # Grab folder name prefix
#         prefix = get_prefix(l)
#
#         # File to fix
#         filename_input = '_'.join([prefix, 'closest.bed'])   # 'closest.bed'
#         path_to_input = os.path.join(wd, filename_input)
#         # Number of lines to correct
#         tail_to_trim = num_lines_table.at[prefix, 'tail_to_trim']
#
#         # Output file
#         filename_output = '_'.join([prefix, 'temp', 'closest.bed'])   # 'temp_closest.bed'
#         path_to_output = os.path.join(wd, filename_output)
#
#         if not os.path.exists(path_to_output) and tail_to_trim != 0:
#             # Construct command
#             cmd_this_head = 'head -n -' + str(tail_to_trim) + ' ' + path_to_input + ' > ' + path_to_output
#             cmd_all_head.append(cmd_this_head)
#
#     # Parallel run by Pool
#     pool = mp.Pool(1)
#     if len(cmd_all_head) is not 0:
#         pool.map(work, cmd_all_head)
#     print('Trim YT..._closest.bed tail lines: finished.\n')
#
#
# def remove_telomeric_reads(libs, parent_dir, suffix_input, suffix_output):
#     cmd_all_awk = []
#     for l in libs:
#         # Setting up working directory
#         wd = os.path.join(parent_dir, l)
#         # Grab folder name prefix
#         prefix = get_prefix(l)
#
#         # File to fix
#         filename_input = '_'.join([prefix, suffix_input])  # 'temp_closest.bed'
#         path_to_input = os.path.join(wd, filename_input)
#
#         # Output file
#         filename_output = '_'.join([prefix, suffix_output])  # 'temp2_closest.bed'
#         path_to_output = os.path.join(wd, filename_output)
#
#         if not os.path.exists(path_to_output):
#             # Construct command
#             cmd_this_awk = 'awk \'BEGIN{FS=\"\\t\"} $14!=-1 {print $0}\' ' + path_to_input + ' > ' + path_to_output
#             cmd_all_awk.append(cmd_this_awk)
#
#     # Parallel run by Pool
#     pool = mp.Pool(1)
#     if len(cmd_all_awk) is not 0:
#         pool.map(work, cmd_all_awk)
#     print('Remove YT..._closest.bed telomeric reads: finished.\n')
#
#
# def remove_5_prime_upstream(libs, parent_dir, suffix_input, suffix_output):
#     # We're only interested in reads mapped to downstream regions of 3'UTRs. Remove reads mapped to upstream regions
#     # of 5'UTRs.
#     cmd_all_awk = []
#     for l in libs:
#         # Setting up working directory
#         wd = os.path.join(parent_dir, l)
#         # Grab folder name prefix
#         prefix = get_prefix(l)
#
#         # File to fix
#         filename_input = '_'.join([prefix, suffix_input])   # 'temp2_closest.bed'
#         path_to_input = os.path.join(wd, filename_input)
#
#         # Output file
#         filename_output = '_'.join([prefix, suffix_output])   # 'fixed_closest.bed'
#         path_to_output = os.path.join(wd, filename_output)
#
#         if not os.path.exists(path_to_output):
#             # Construct command
#             cmd_this_awk = 'awk \'BEGIN{FS="\\t";OFS="\\t"}{if($23<0){print $0}}\' ' + path_to_input + ' > ' \
#                            + path_to_output
#             cmd_all_awk.append(cmd_this_awk)
#
#     # Parallel run by Pool
#     pool = mp.Pool(1)
#     if len(cmd_all_awk) is not 0:
#         pool.map(work, cmd_all_awk)
#     print('Remove YT..._closest.bed telomeric entries: finished.\n')
#
#
# def remove_intermediate(libs, parent_dir, suffix_file_to_remove):
#     # Remove intermediate files with 'find ... -delete'
#     # Check that file numbers match
#     num_libs = len(libs)
#     # Count file number with pipe
#     cmd_count_extracted_part_1 = 'find ' + parent_dir + ' -name "*' + suffix_file_to_remove + '" -type f'
#     cmd_count_extracted_part_2 = 'wc -l'
#     # Open process
#     p1 = subprocess.Popen(shlex.split(cmd_count_extracted_part_1), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     p2 = subprocess.Popen(shlex.split(cmd_count_extracted_part_2), stdin=p1.stdout, stdout=subprocess.PIPE)
#     p1.stdout.close()
#     while True:
#         output = p2.communicate()[0]
#         if p2.poll() is not None:
#             break
#     num_count = int(output.decode())
#     # Command to remove files
#     cmd_rm_intermediate = ['find ' + parent_dir + ' -name "*' + suffix_file_to_remove + '" -type f -delete']
#     pool = mp.Pool(1)
#     if num_count == num_libs and len(cmd_rm_intermediate) is not 0:
#         pool.map(work, cmd_rm_intermediate)
#         print('Files ending with ' + suffix_file_to_remove + 'successfully removed.\n')
#     else:
#         warnings.warn('Trying to remove ..' + suffix_file_to_remove + ', but file number does not match with library '
#                       'number. Abort.\n')
