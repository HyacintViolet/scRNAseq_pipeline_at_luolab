import os
import re
import subprocess
import multiprocessing as mp


def get_libs(parent_dir):
    # Input path to mapping dir. Output a list of library names.
    libs = sorted(os.listdir(parent_dir))
    return libs


def get_prefix(lib, option='linrui_ipcr'):
    # Input library folder name. Output library prefix.
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', lib)
    if option is 'linrui_ipcr':
        prefix = '_'.join(match.group(1, 2, 3, 4))
    return prefix


def parse_input_suffix(task, option, suffix_input):
    if task is "mapping":
        if option is "pair_end":
            # TODO: check input has 2 files and they end with r1.fq.gz/r2.fq.gz
            suffix_input_read_1 = suffix_input[0]
            suffix_input_read_2 = suffix_input[1]
            suffix_input_parsed = [suffix_input_read_1, suffix_input_read_2]
    elif task is "index":
        # TODO: check input ends with bam
        suffix_input_parsed = suffix_input
    return suffix_input_parsed


def work(cmd):
    return subprocess.call(cmd, shell=True)


def do_parallel(src_dir=None, dst_dir=None, suffix_input=None, suffix_output=None, task=None, option="pair_end",
                thread=1, genome_index='/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM23/'):



    libs = get_libs(src_dir)
    cmd_all = []
    for lib in libs:
        suffix_input_parsed = parse_input_suffix(task=task, option=None, suffix_input=suffix_input)

        # Setting up input/output directory
        input_dir = os.path.join(src_dir, lib)
        output_dir = os.path.join(dst_dir, lib)
        # Create output directory if it doesn't exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Grab folder name prefix
        prefix = get_prefix(lib)

        # Input file name
        if option is "pair_end":
            file_input_read_1 = '_'.join([prefix, suffix_input_read_1])
            file_input_read_2 = '_'.join([prefix, suffix_input_read_2])
            path_to_input_read_1 = os.path.join(input_dir, file_input_read_1)
            path_to_input_read_2 = os.path.join(input_dir, file_input_read_2)

        # Output file name
        if task is "mapping":
            out_file_name_prefix = os.path.join(output_dir, prefix) + '_'
            file_output_bam = '_'.join([prefix, 'Aligned.sortedByCoord.out.bam'])
            path_to_output_bam = os.path.join(output_dir, file_output_bam)

        # Check if output already exists. If not, construct command.
        if not os.path.exists(path_to_output_bam):
            cmd_this = 'STAR --runThreadN 32 --genomeDir ' + genome_index + \
                       ' --readFilesIn ' + path_to_input_read_1 + ' ' + path_to_input_read_2 + \
                       ' --readFilesCommand zcat --outReadsUnmapped Fastx --outSAMtype BAM ' \
                       'SortedByCoordinate --outFileNamePrefix ' + out_file_name_prefix
            cmd_all.append(cmd_this)

    # Parallel run by Pool
    pool = mp.Pool(1)
    pool.map(work, cmd_all)
    print('Mapping:' + src_dir + 'finished.')
#
# cmd_index = []
#     for out in os.listdir(dst):
#
#         # Setup output directory
#         out_dir = os.path.join(dst, out)
#
#         # Change directory to output folder
#         os.chdir(out_dir)
#
#         # Grab folder name
#         match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
#         prefix = match.group(1)
#
#         # Input file name for samtools index
#         out_name_samtools = '_'.join([prefix, 'assigned_sorted.bam'])
#         index_in = os.path.join(out_dir, out_name_samtools)
#
#         # Construct command
#         cmd_this_index = 'samtools index ' + index_in
#         cmd_index.append(cmd_this_index)
#
#     # Parallel run by Pool
#     pool = mp.Pool(16)
#     pool.map(work, cmd_index)
#     print('samtools index: finished.')