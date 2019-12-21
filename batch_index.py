import os
import re
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


def parse_input_output(src_dir, dst_dir, lib, task=None, set_cell_number=80):
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

    elif task is "samtools_index":
        input_args['input'] = os.path.join(in_dir, '_'.join([prefix, 'ambiguity.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'ambiguity.bam.bai']))

    return input_args, output_args


def parse_command(input_args, output_args, task=None, num_thread=None, genome_index=None, genome_gtf=None):

    # Coerce args to string
    num_thread = str(num_thread)

    # Initialize empty var.
    cmd = ''

    if task is None:
        print("Error in parse_command: task not specified.")
    elif task is "samtools_index":
        cmd = ' '.join(['samtools', 'index', input_args['input']])
        # Example command:
        # samtools index /path/to/map_result/YT013101_some_file.bam

    return cmd


def get_libs(parent_dir):
    # Input path to mapping dir. Output a list of library names.
    libs = sorted(os.listdir(parent_dir))
    return libs


def get_prefix(lib):
    # Input library folder name. Output library prefix.
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', lib)
    prefix = match.group(1)
    return prefix


def work(cmd_this):
    lib = cmd_this[0][0]
    task = cmd_this[0][1]
    cmd = cmd_this[1]
    print(' '.join([task, lib, '...\ncommand:', cmd]))
    return subprocess.call(cmd, shell=True)


def do_parallel(src_dir=None, dst_dir=None, task=None, overwrite=True, num_process=1, num_thread=1, genome_index=None,
                genome_gtf=None):

    # Check idle status. If idle is True, proceed and set idle to False.
    proceed_if_idle()

    # Initiate
    libs = get_libs(src_dir)
    cmd_all = []
    for lib in libs:
        prefix = get_prefix(lib)

        # Parse input/output args
        input_args, output_args = parse_input_output(src_dir, dst_dir, lib, task=task)

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


def main():

    # Source dirs
    src_dir = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting/'
    src_dir2 = '/media/luolab/ZA1BT1ER/yanting/vM23_extended/mapping/'

    # Destination dirs
    dst_dir = '/media/luolab/ZA1BT1ER/yanting/vM23_extended/mapping/'  # Same as src_dir2

    # Parent working dir
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23_extended/'

    # Create output directory if not exist.barcode_ground_truth
    for out in get_libs(src_dir):
        if not os.path.exists(os.path.join(dst_dir, out)):
            os.makedirs(os.path.join(dst_dir, out))

    # Execute: samtools index
    do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="samtools_index", num_process=32)


if __name__ == '__main__':
    main()
