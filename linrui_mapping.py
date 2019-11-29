import os
import re
import multiprocessing as mp
import subprocess


def get_libs(parent_dir):
    # Input path to mapping dir. Output a list of library names.
    libs = sorted(os.listdir(parent_dir))
    return libs


def get_prefix(lib, option='linrui_ipcr'):
    # Input library folder name. Output library prefix.
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', lib)
    if option is 'linrui_ipcr':
        prefix = '_'.join(match.group(1,2,3,4))
    return prefix


def work(cmd):
    return subprocess.call(cmd, shell=True)


def do_parallel(src_dir=None, dst_dir=None, suffix_input=None, suffix_output=None, task="mapping", option="pair_end",
                genome_index='/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM23/'):

    if task is "mapping":
        if option is "pair_end":
            suffix_input_read_1 = suffix_input[0]
            suffix_input_read_2 = suffix_input[1]

        libs = get_libs(src_dir)
        cmd_all = []
        for lib in libs:

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
                file_input_read_1 = '_'.join([prefix,suffix_input_read_1])
                file_input_read_2 = '_'.join([prefix,suffix_input_read_2])
                path_to_input_read_1 = os.path.join(input_dir, file_input_read_1)
                path_to_input_read_2 = os.path.join(input_dir, file_input_read_2)

            # Output file name
            if task is "mapping":
                out_file_name_prefix = os.path.join(output_dir, prefix) + '_'
                file_output_bam = '_'.join([prefix,'Aligned.sortedByCoord.out.bam'])
                path_to_output_bam = os.path.join(output_dir, file_output_bam)

            # Check if output already exists. If not, construct command.
            if not os.path.exists(path_to_output_bam):
                cmd_this = 'STAR --runThreadN 32 --genomeDir ' + genome_index + \
                                   ' --readFilesIn ' + path_to_input_read_1 + ' ' + path_to_input_read_2 + \
                                   ' --readFilesCommand zcat --outSAMmultNmax 1 --outReadsUnmapped Fastx' \
                                   ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + out_file_name_prefix
                cmd_all.append(cmd_this)

        # Parallel run by Pool
        pool = mp.Pool(1)
        pool.map(work, cmd_all)
        print('Mapping:' + src_dir + 'finished.')


def main():

    data_dir = '/media/luolab/ZA1BT1ER/linrui/data/'
    mapping_dir = '/media/luolab/ZA1BT1ER/linrui/mapping/'
    genome_index = '/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM23/'
    suffix_input = ["r1_NlaIII.fq.gz","r2_NlaIII.fq.gz"]

    do_parallel(data_dir, mapping_dir, suffix_input, task="mapping", option="pair_end", genome_index=genome_index)


if __name__ == '__main__':
    main()