import os
import re
import shlex
import warnings
import multiprocessing as mp
import subprocess
import pandas as pd

idle = True


def set_idle_status(status):
    global idle
    idle = status


def proceed_if_idle():
    global idle
    while True:
        if idle is True:
            set_idle_status(False)
            break


def get_libs(dir):
    # Input path to mapping dir. Output a list of library names.
    libs = sorted(os.listdir(dir))
    return libs


def get_prefix(lib):
    # Input library folder name. Output library prefix.
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', lib)
    prefix = match.group(1)
    return prefix


def check_file_exists(path_to_file):
    if not os.path.exists(path_to_file):
        raise SystemExit(' '.join(['Aborting because', os.path.basename(path_to_file), 'was not found. Please prepare '
                                   'the file before running the pipe. See readme.txt for detail.']))


# def args_int2str(arg): TODO: check args in function calls, convert int to string
#     if type(arg) is int:
#         arg = str(arg)
#     return arg


def work(cmd_this):
    lib = cmd_this[0][0]
    task = cmd_this[0][1]
    cmd = cmd_this[1]
    print(' '.join([task, lib, '...\ncommand:', cmd]))
    return subprocess.call(cmd, shell=True)


def parse_input_output(src_dir, dst_dir, lib, task=None, set_cell_number=80, supplementary_file=None):
    in_dir = os.path.join(src_dir, lib)
    out_dir = os.path.join(dst_dir, lib)
    prefix = get_prefix(lib)

    # Coerce args to str
    set_cell_number = str(set_cell_number)

    # Initialize empty dictionary
    input_args = dict()
    output_args = dict()

    if task is None:
        print("Error in parse_input_output: task not specified.")

    elif task is "bedtools_intersect":
        check_file_exists(supplementary_file)
        input_args['input_bam'] = os.path.join(in_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))
        input_args['input_gtf'] = supplementary_file
        output_args['output'] = os.path.join(in_dir, '_'.join([prefix, 'stranded_nonoverlap.bam']))

    elif task is "bedtools_closest":
        check_file_exists(supplementary_file)
        input_args['input_bam'] = os.path.join(in_dir, '_'.join([prefix, 'stranded_nonoverlap.bam']))
        input_args['input_bed'] = supplementary_file
        output_args['output'] = os.path.join(in_dir, '_'.join([prefix, 'closest.bed']))

    return input_args, output_args


def parse_command(input_args, output_args, task=None, num_thread=None, genome_index=None, genome_gtf=None):

    # Coerce args to string
    num_thread = str(num_thread)

    # Initialize empty var.
    cmd = ''

    if task is None:
        print("Error in parse_command: task not specified.")

    elif task is "bedtools_intersect":
        cmd = ' '.join(['bedtools', 'intersect', '-v', '-s', '-a', input_args['input_bam'], '-b',
                        input_args['input_gtf'], '>', output_args['output']])

    elif task is "bedtools_closest":
        cmd = ' '.join(['bedtools', 'closest', '-s', '-D', 'a', '-io', '-t', 'first', '-a',
                        input_args['input_bam'], '-b', input_args['input_bed'], '>', output_args['output']])
        # -D a: signed dist. w.r.t. bam's strand. Desired output: negative. Reversed sign in
        # extract_desired_bed_files.py, awk_extract_part() function.
    return cmd


def do_parallel(src_dir=None, dst_dir=None, task=None, overwrite=True, num_process=1, num_thread=1, genome_index=None,
                genome_gtf=None, supplementary_file=None):

    # Check idle status. If idle is True, proceed and set idle to False.
    proceed_if_idle()

    # Initiate
    libs = get_libs(src_dir)
    cmd_all = []
    for lib in libs:
        prefix = get_prefix(lib)

        # Parse input/output args
        input_args, output_args = parse_input_output(src_dir, dst_dir, lib, task=task,
                                                     supplementary_file=supplementary_file)

        # Construct commands. Store to cmd_all for (parallel) iter run by mp.Pool.
        if overwrite:
            cmd = parse_command(input_args, output_args, task=task, num_thread=num_thread,
                                genome_index=genome_index, genome_gtf=genome_gtf)
            cmd_all.append([(prefix, task), cmd])
        else:
            check_exists = []
            for file, path in output_args.items():
                check_exists.append(os.path.exists(path))
            if any(check_exists):
                print('Some or all of the output file already exists. Skipping since overwrite mode is inactivated.')
            else:
                cmd = parse_command(input_args, output_args, task=task, num_thread=num_thread,
                                    genome_index=genome_index, genome_gtf=genome_gtf)
                cmd_all.append([(prefix, task), cmd])

    # Check process & thread number
    if num_process*num_thread > 47:
        raise SystemExit(' '.join(['Process & thread number exploded.', 'Task:', task]))

    # Parallel run by pool
    pool = mp.Pool(processes=num_process)
    if len(cmd_all) is not 0:
        pool.map(work, cmd_all)
    print(task + ': finished.')
    # Change idle status
    set_idle_status(True)


def remove_intermediate(wd, suffix_file_to_remove):
    # Remove intermediate files with 'find ... -delete'
    # Check that file numbers match
    libs = get_libs(wd)
    num_libs = len(libs)
    # Count file number with pipe
    cmd_count_extracted_part_1 = 'find ' + wd + ' -name "*' + suffix_file_to_remove + '" -type f'
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
    cmd_rm_intermediate = ['find ' + wd + ' -name "*' + suffix_file_to_remove + '" -type f -delete']
    pool = mp.Pool(1)
    if num_count == num_libs and len(cmd_rm_intermediate) is not 0:
        pool.map(work, cmd_rm_intermediate)
        print('Files ending with ' + suffix_file_to_remove + 'successfully removed.\n')
    else:
        warnings.warn('Trying to remove ..' + suffix_file_to_remove + ', but file number does not match with library '
                      'number. Abort.\n')


def main():
    # Source dirs
    src_dir = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting/'
    src_dir2 = '/media/luolab/ZA1BT1ER/yanting/vM21/mapping/'  # For unextended mapping

    # Destination dirs
    dst_dir = '/media/luolab/ZA1BT1ER/yanting/vM21/mapping/'  # For unextended mapping

    # Parent working dir
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM21/'

    # Supplementary files
    gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM21/gencode.vM21.chr_primary.annotation.gtf'
    bed = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM21/gencode.vM21.chr_primary.annotation.fixed.' \
          'sortedByChrom.bed'

    # Create output directory if not exist.barcode_ground_truth
    for out in get_libs(src_dir):
        if not os.path.exists(os.path.join(dst_dir, out)):
            os.makedirs(os.path.join(dst_dir, out))

    # bedtools_intersect
    # do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="bedtools_intersect", supplementary_file=gtf)

    # remove stranded_nonoverlap.bam
    # remove_intermediate(src_dir2, suffix_file_to_remove="stranded_nonoverlap.bam")

    # bedtools_closest
    do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="bedtools_closest", supplementary_file=bed)


if __name__ == '__main__':
    main()

# parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM21/mapping/'
# gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM21/gencode.vM21.chr_primary.annotation.gtf'
# bed = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM21/gencode.vM21.chr_primary.annotation.fixed.bed'
#
#
# def work(cmd):
#     lib = re.search('(YT[0-9]*)', cmd).group(1)
#     if 'intersect' in cmd:
#         print('Processing: bedtools intersect... library:' + lib)
#     elif 'closest' in cmd:
#         print('Processing: bedtools closest... library:' + lib)
#     return subprocess.call(cmd, shell=True)
#
#
# cmd_all = []
# for s in os.listdir(parent_dir):
#
#     # Setting up input/output directory
#     wd = os.path.join(parent_dir, s)
#     # out_dir = os.path.join(dst, s)
#
#     # Grab folder name prefix
#     match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', s)
#     prefix = match.group(1)
#
#     # Input file name
#     file_A = '_'.join([prefix, 'Aligned.sortedByCoord.out.bam'])
#     path_to_file_A = os.path.join(wd, file_A)
#
#     # Output file name
#     out_name_intersect = '_'.join([prefix, 'stranded_nonoverlap.bam'])
#     path_to_intersect_out = os.path.join(wd, out_name_intersect)
#
#     out_name_closest = '_'.join([prefix, 'closest.bed'])
#     path_to_closest_out = os.path.join(wd, out_name_closest)
#
#     # Check if output already exists. If not, construct command.
#     if not os.path.exists(path_to_intersect_out):
#         # Construct command
#         cmd_this_intersect = 'bedtools intersect -v -s -a ' + path_to_file_A + ' -b ' + \
#                              gtf + ' > ' + path_to_intersect_out
#         cmd_all.append(cmd_this_intersect)
#         cmd_this_intersect = []
#
#     if not os.path.exists(path_to_closest_out):
#         cmd_this_closest = 'bedtools closest -s -D a -io -t first -a ' + path_to_intersect_out + \
#                            ' -b ' + bed + ' > ' + path_to_closest_out
#         cmd_all.append(cmd_this_closest)
#         cmd_this_closest = []
#
#
# # Parallel run by Pool
# pool = mp.Pool(1)
# pool.map(work, cmd_all)
# print('bedtools intersect & closest: finished.')
