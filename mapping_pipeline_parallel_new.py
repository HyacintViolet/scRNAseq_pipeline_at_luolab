# It turns out pycharm has a bug that prevents the script from searching into default PATH. Need to manually set
# environment variable PATH.
# Here's how: Run -> Edit Configurations -> Environment Variables -> manually add PATH entry; from terminal: echo PATH,
# copy and paste into value.


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
    elif task is "umitools_whitelist":
        read1_filename = ''
        for file in os.listdir(in_dir):
            if file.endswith('2.fq.gz') or file.endswith('2.clean.fq.gz'):
                read1_filename = file
        input_args['path_to_read1'] = os.path.join(in_dir, read1_filename)
        input_args['set_cell_number'] = set_cell_number
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'whitelist'+set_cell_number+'.txt']))
        output_args['plot_prefix'] = os.path.join(out_dir, '_'.join(['cell_num', set_cell_number]))

    elif task is "wash_whitelist":
        input_args['path_to_whitelist'] = os.path.join(in_dir,
                                                       '_'.join([prefix, 'whitelist'+set_cell_number+'.txt']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'whitelist_washed.txt']))

    elif task is "umitools_extract":
        read1_filename = ''
        read2_filename = ''
        for file in os.listdir(in_dir):
            if file.endswith('2.fq.gz') or file.endswith('2.clean.fq.gz'):
                read1_filename = file
            elif file.endswith('1.fq.gz') or file.endswith('1.clean.fq.gz'):
                read2_filename = file
        input_args['path_to_read1'] = os.path.join(in_dir, read1_filename)
        input_args['path_to_read2'] = os.path.join(in_dir, read2_filename)
        input_args['whitelist'] = os.path.join(out_dir, '_'.join([prefix, 'whitelist'+set_cell_number+'.txt']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'extracted.fq.gz']))

    elif task is "STAR_mapping":
        input_args['extracted'] = os.path.join(in_dir, '_'.join([prefix, 'extracted.fq.gz']))
        output_args['out_prefix'] = os.path.join(out_dir, prefix+'_')
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))

    elif task is "split_bam":
        input_args['input'] = os.path.join(in_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))
        output_args['unmapped_output'] = os.path.join(in_dir, '_'.join([prefix, 'unmapped_sorted.bam']))
        output_args['tmp'] = os.path.join(in_dir, 'tmp')
        output_args['mapped_output'] = input_args['input']

    elif task is "alignment_stats":
        input_args['mapped_input'] = os.path.join(in_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))
        input_args['unmapped_input'] = os.path.join(in_dir, '_'.join([prefix, 'unmapped_sorted.bam']))
        output_args['uniquemapped'] = os.path.join(in_dir, '_'.join([prefix, 'uniquemapped_by_cell.txt']))
        output_args['multimapped'] = os.path.join(in_dir, '_'.join([prefix, 'multimapped_by_cell.txt']))
        output_args['unmapped'] = os.path.join(in_dir, '_'.join([prefix, 'unmapped_by_cell.txt']))

    elif task is "merge_aln_stats":
        input_args['uniquemapped'] = os.path.join(in_dir, '_'.join([prefix, 'uniquemapped_by_cell.txt']))
        input_args['multimapped'] = os.path.join(in_dir, '_'.join([prefix, 'multimapped_by_cell.txt']))
        input_args['unmapped'] = os.path.join(in_dir, '_'.join([prefix, 'unmapped_by_cell.txt']))
        output_args['tmp'] = os.path.join(in_dir, 'tmp')
        output_args['aln_stats'] = os.path.join(in_dir, '_'.join([prefix, 'aln_stats_by_cell.txt']))

    elif task is "featurecounts":
        input_args['mapped'] = os.path.join(in_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'gene_assigned']))

    elif task is "samtools_sort":
        input_args['input'] = os.path.join(in_dir,
                                           '_'.join([prefix, 'Aligned.sortedByCoord.out.bam.featureCounts.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'assigned_sorted.bam']))

    elif task is "samtools_index":
        input_args['input'] = os.path.join(in_dir, '_'.join([prefix, 'assigned_sorted.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'assigned_sorted.bam.bai']))

    elif task is "umitools_count":
        input_args['input'] = os.path.join(in_dir, '_'.join([prefix, 'assigned_sorted.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'counts.tsv.gz']))

    return input_args, output_args


def parse_command(input_args, output_args, task=None, num_thread=None, genome_index=None, genome_gtf=None):

    # Coerce args to string
    num_thread = str(num_thread)

    # Initialize empty var.
    cmd = ''

    if task is None:
        print("Error in parse_command: task not specified.")

    elif task is "umitools_whitelist":
        cmd = ' '.join(['umi_tools', 'whitelist', '--stdin', input_args['path_to_read1'],
                        '--bc-pattern=CCCCCCCCNNNNNNNN', '--set-cell-number=' + input_args['set_cell_number'],
                        '--plot-prefix=' + output_args['plot_prefix'], '-v', '1', '--log2stderr', '>',
                        output_args['output']])
        # Example command:
        # umi_tools whitelist --stdin /path/to/seq_data/YT013101_......_2.fq.gz
        #                     --bc-pattern=CCCCCCCCNNNNNNNN
        #                     --set-cell-number=80 --plot-prefix=/path/to/map_result/cell_num_80
        #                     -v 1
        #                     --log2stderr > /path/to/map_result/YT013101_whitelist80.txt

    elif task is "umitools_extract":
        cmd = ' '.join(['umi_tools', 'extract', '--bc-pattern=CCCCCCCCNNNNNNNN', '--stdin', input_args['path_to_read1'],
                       '--read2-in', input_args['path_to_read2'], '--stdout', output_args['output'], '--read2-stdout',
                        '--filter-cell-barcode', '--error-correct-cell', '--whitelist=' + input_args['whitelist']])
        # Example command:
        # umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN
        #                   --stdin /path/to/seq_data/YT013101_..._2.fq.gz
        #                   --read2-in /path/to/seq_data/YT013101_..._1.fq.gz
        #                   --stdout /path/to/map_result/YT013101_extracted.fq.gz
        #                   --read2-stdout --filter-cell-barcode --error-correct-cell
        #                   --whitelist=/path/to/map_result/YT013101_whitelist_washed.txt

    elif task is "STAR_mapping":
        cmd = ' '.join(['STAR', '--runThreadN', num_thread, '--genomeDir', genome_index,
                        '--readFilesIn', input_args['extracted'], '--readFilesCommand', 'zcat',
                        '--outFilterMultimapNmax', '1', '--outFilterType', 'BySJout',
                        '--outSAMstrandField', 'intronMotif', '--outFilterIntronMotifs', 'RemoveNoncanonical',
                        '--outFilterMismatchNmax', '6', '--outSAMtype', 'BAM', 'SortedByCoordinate',
                        '--outFileNamePrefix', output_args['out_prefix'], '--outSAMunmapped', 'Within'])
        # Example command:
        # STAR --runThreadN 32 --genomeDir /media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM23_extended/
        #      --readFilesIn /path/to/map_result/YT013101_extracted.fq.gz
        #      --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout
        #      --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical
        #      --outFilterMismatchNmax 6 --outSAMtype BAM SortedByCoordinate
        #      --outFileNamePrefix /path/to/map_result/YT013101_
        #      --outSAMunmapped Within

    elif task is "split_bam":
        cmd = ' '.join(['samtools', 'view', '-f4', '-bS', input_args['input'], '-@', num_thread, '|',
                        'samtools', 'sort', '-', '-o', output_args['unmapped_output'], '-@', num_thread, '&&',
                        'samtools', 'view', '-F4', '-bS', input_args['input'], '-@', num_thread, '|',
                        'samtools', 'sort', '-', '-o', output_args['tmp'], '-@', num_thread, '&&',
                        'mv', output_args['tmp'], output_args['mapped_output']])
        # Example command:
        # samtools view -f4 -bS YT013101_Aligned.sortedByCoord.out.bam -@ 32 | samtools sort - -o
        # YT013101_unmapped_sorted.bam -@ 32 && samtools view -F4 -bS YT013101_Aligned.sortedByCoord.out.bam -@ 32 |
        # samtools sort - -o tmp -@ 32 && mv tmp YT013101_Aligned.sortedByCoord.out.bam

    elif task is "alignment_stats":
        cmd = ' '.join(['samtools', 'view', input_args['mapped_input'], '-@', num_thread, '|', 'cut', '-f2', '-d',
                        '\'_\'', '|', 'sort', '|', 'uniq', '-c', '>', output_args['uniquemapped'], '&&',
                        'samtools', 'view', input_args['unmapped_input'], '-@', num_thread, '|', 'grep', '-w',
                        '"uT:A:3"', '|', 'cut', '-f2', '-d', '\'_\'', '|', 'sort', '|', 'uniq', '-c', '>',
                        output_args['multimapped'], '&&',
                        'samtools', 'view', input_args['unmapped_input'], '-@', num_thread, '|', 'grep', '-w',
                        '"uT:A:[^3]"', '|', 'cut', '-f2', '-d', '\'_\'', '|', 'sort', '|', 'uniq', '-c', '>',
                        output_args['unmapped']])
        # Example command:
        # samtools view YT013101_Aligned.sortedByCoord -@ 32 | cut -f2 -d '_' | sort | uniq -c >
        # YT013101_uniquemapped_by_cell.txt &&
        # samtools view YT013101_unmapped_sorted.bam -@ 32 | grep -w "uT:A:3" | cut -f2 -d '_' | sort | uniq -c >
        # YT013101_multimapped_by_cell.txt &&
        # samtools view YT013101_unmapped_sorted.bam -@ 32 | grep -w "uT:A:[^3]" | cut -f2 -d '_' | sort | uniq -c >
        # YT013101_unmapped_by_cell.txt

    elif task is "merge_aln_stats":
        cmd = ' '.join(['join', '-j', '2', input_args['uniquemapped'], input_args['multimapped'], '|', 'join', '-1',
                        '1', '-2', '2', '-', input_args['unmapped'], '>', output_args['aln_stats'], '&&',
                        'echo', '\'Barcode Uniquemapped Multimapped Unmapped\'', '|',
                        'cat', '-', output_args['aln_stats'], '>', output_args['tmp'], '&&',
                        'mv', output_args['tmp'], output_args['aln_stats']])
        # Example command:
        # join -j 2 YT013101_uniquemapped_by_cell.txt YT013101_multimapped_by_cell.txt | join -1 1 -2 2 -
        # YT013101_unmapped_by_cell.txt > YT013101_aln_stats_by_cell.txt &&
        # echo 'Barcode Uniquemapped Multimapped Unmapped' | cat - YT013101_aln_stats_by_cell.txt > tmp &&
        # mv tmp YT013101_aln_stats_by_cell.txt

    elif task is "featurecounts":
        cmd = ' '.join(['featureCounts', '-s', '1', '-a', genome_gtf, '-o', output_args['output'], '-R', 'BAM',
                        input_args['mapped'], '-T', num_thread])
        # Example command:
        # featureCounts -s 1 -a /media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23_extended/gencode...gtf
        #               -o /path/to/map_result/YT013101_gene_assigned
        #               -R BAM /path/to/map_result/YT013101_Aligned.sortedByCoord.out.bam
        #               -T 32

    elif task is "samtools_sort":
        cmd = ' '.join(['samtools', 'sort', input_args['input'], '-o', output_args['output'], '-@', num_thread])
        # Example command:
        # samtools sort /path/to/map_result/YT013101_Aligned.sortedByCoord.out.bam.featureCounts.bam
        #               -o /path/to/map_result/YT013101_assigned_sorted.bam -@ 32

    elif task is "samtools_index":
        cmd = ' '.join(['samtools', 'index', input_args['input']])
        # Example command:
        # samtools index /path/to/map_result/YT013101_assigned_sorted.bam

    elif task is "umitools_count":
        cmd = ' '.join(['umi_tools', 'count', '--per-gene', '--gene-tag=XT', '--per-cell', '-I', input_args['input'],
                        '-S', output_args['output']])
        # Example command:
        # umi_tools count --per-gene --gene-tag=XT --per-cell -I /path/to/map_result/YT013101_assigned_sorted.bam
        #                 -S /path/to/map_result/YT013101_counts.tsv.gz

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


def wash_whitelist(src_dir, dst_dir, parent_dir, task="wash_whitelist", overwrite=True):

    # Check idle status. If idle is True, proceed and set idle to False.
    proceed_if_idle()

    # Read barcode ground truth list
    barcode_ground_truth_raw = pd.read_excel(os.path.join(parent_dir, 'barcode_ground_truth_checklist.xlsx'))
    barcode_ground_truth = barcode_ground_truth_raw['Primer_sequence'].str.extract(
        r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})', expand=False)

    libs = get_libs(src_dir)
    for lib in libs:
        # Display progress
        print('Washing barcode whitelist. Library: ' + lib)

        # Parse input/output args
        input_args, output_args = parse_input_output(src_dir, dst_dir, lib, task=task)

        # Read whitelist80
        whitelist = pd.read_csv(input_args["path_to_whitelist"], sep="\t",
                                names=['cell', 'candidate', 'Nreads', 'Ncandidate'])

        # Remove rows whose 'cell' value is not found in barcode_ground_truth
        whitelist_washed = whitelist[whitelist['cell'].isin(barcode_ground_truth)]
        # Reset index
        whitelist_washed = whitelist_washed.reset_index(drop=True)

        # Write to output
        if overwrite:
            whitelist_washed.to_csv(output_args['output'], sep="\t", index=False, header=False)
        else:
            if not os.path.exists(output_args['output']):
                whitelist_washed.to_csv(output_args['output'], sep="\t", index=False, header=False)
            else:
                print(lib + "whitelist_washed.txt already exists. Skipping.")
    # Change idle status
    set_idle_status(True)


def main():

    # Source dirs
    src_dir = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting/'
    src_dir2 = '/media/luolab/ZA1BT1ER/yanting/vM23_extended/mapping/'
    src_dir3 = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'  # For unextended mapping

    # Destination dirs
    dst_dir = '/media/luolab/ZA1BT1ER/yanting/vM23_extended/mapping/'  # Same as src_dir2
    dst_dir2 = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'  # For unextended mapping

    # Parent working dir
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/'

    # Path to genome annotation and index
    genome_gtf_unextended = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr_patch_hapl_scaff.' \
                            'annotation.gtf'
    genome_gtf_extended = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr_patch_hapl_scaff.' \
                          'annotation.extended.gtf'
    genome_gtf_extended_clean = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr_patch_hapl_' \
                                'scaff.annotation.extended.clean.gtf'
    genome_index = '/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM23_extended/'

    # Create output directory if not exist.barcode_ground_truth
    for out in get_libs(src_dir):
        if not os.path.exists(os.path.join(dst_dir, out)):
            os.makedirs(os.path.join(dst_dir, out))

    # STEP 1: umi_tools whitelist

    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="umitools_whitelist", num_process=32)

    # STEP 2: wash whitelist
    # wash_whitelist(src_dir=src_dir2, dst_dir=dst_dir, parent_dir=parent_dir, task="wash_whitelist", overwrite=True)

    # STEP 3: umi_tools extract
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="umitools_extract", num_process=32)

    # STEP 4: STAR mapping
    # do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="STAR_mapping", genome_index=genome_index, num_thread=32)

    # STEP 5: split aligned bam file into mapped and unmapped
    # do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="split_bam", num_thread=32)

    # STEP 6: extract alignment statistics
    # do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="alignment_stats", num_thread=32)

    # STEP 7: merge alignment statistics
    # do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="merge_aln_stats", num_process=32)

    # STEP 8: featureCounts
    do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="featurecounts", genome_gtf=genome_gtf_extended_clean,
                num_thread=32)

    # STEP 9: samtools sort
    do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="samtools_sort", num_thread=32)

    # STEP 10: samtools index
    do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="samtools_index", num_process=32)

    # STEP 11: umitools count
    do_parallel(src_dir=src_dir2, dst_dir=dst_dir, task="umitools_count", num_process=24, overwrite=False)


if __name__ == '__main__':
    main()
