# It turns out pycharm has a bug that prevents the script from searching into default PATH. Need to manually set
# environment variable PATH.
# Here's how: Run -> Edit Configurations -> Environment Variables -> manually add PATH entry; from terminal: echo PATH,
# copy and paste into value.


import os
import re
import multiprocessing as mp
import subprocess
import pandas as pd


def parse_input_output(src_dir, dst_dir, l, task=None, set_cell_number=80):
    in_dir = os.path.join(src_dir, l)
    out_dir = os.path.join(dst_dir, l)
    prefix = get_prefix(l)
    input_args = dict()
    output_args = dict()

    if task is None:
        print("Error in parse_input_output: task not specified.")
    elif task is "umitools_whitelist":
        for file in os.listdir(in_dir):
            if file.endswith('2.fq.gz') or file.endswith('2.clean.fq.gz'):
                read1_filename = file
        input_args['path_to_read1'] = os.path.join(in_dir, read1_filename)
        input_args['set_cell_number'] = str(set_cell_number)
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'whitelist'+str(set_cell_number)+'.txt']))
        output_args['plot_prefix'] = os.path.join(out_dir, '_'.join(['cell_num', str(set_cell_number)]))

    elif task is "wash_whitelist":
        input_args['path_to_whitelist'] = os.path.join(out_dir,
                                                       '_'.join([prefix, 'whitelist'+str(set_cell_number)+'.txt']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'whitelist_washed.txt']))

    elif task is "umitools_extract":
        for file in os.listdir(in_dir):
            if file.endswith('2.fq.gz') or file.endswith('2.clean.fq.gz'):
                read1_filename = file
            elif file.endswith('1.fq.gz') or file.endswith('1.clean.fq.gz'):
                read2_filename = file
        input_args['path_to_read1'] = os.path.join(in_dir, read1_filename)
        input_args['path_to_read2'] = os.path.join(in_dir, read2_filename)
        input_args['whitelist'] = os.path.join(out_dir, '_'.join([prefix, 'whitelist'+str(set_cell_number)+'.txt']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'extracted.fq.gz']))

    elif task is "STAR_mapping":
        input_args['extracted'] = os.path.join(out_dir, '_'.join([prefix, 'extracted.fq.gz']))
        output_args['out_prefix'] = os.path.join(out_dir, prefix+'_')
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))

    elif task is "featurecounts":
        input_args['mapped'] = os.path.join(out_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'gene_assigned']))

    elif task is "samtools_sort":
        input_args['input'] = os.path.join(out_dir,
                                           '_'.join([prefix, 'Aligned.sortedByCoord.out.bam.featureCounts.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'assigned_sorted.bam']))

    elif task is "samtools_index":
        input_args['input'] = os.path.join(out_dir, '_'.join([prefix, 'assigned_sorted.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'assigned_sorted.bam.bai']))

    elif task is "umitools_count":
        input_args['input'] = os.path.join(out_dir, '_'.join([prefix, 'assigned_sorted.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'counts.tsv.gz']))

    elif task is "nuniquemapped":
        input_args['input'] = os.path.join(out_dir, '_'.join([prefix, 'Aligned.sortedByCoord.out.bam']))
        output_args['output'] = os.path.join(out_dir, '_'.join([prefix, 'Nuniqmapped.txt']))

    return input_args, output_args


def parse_command(input_args, output_args, task=None, genome_index=None, genome_gtf=None):
    if task is None:
        print("Error in parse_command: task not specified.")
    elif task is "umitools_whitelist":
        cmd = 'umi_tools whitelist --stdin ' + input_args['path_to_read1'] + \
              ' --bc-pattern=CCCCCCCCNNNNNNNN ' + '--set-cell-number=' + input_args['set_cell_number'] +\
              ' --plot-prefix=' + output_args['plot_prefix'] + ' -v 1 --log2stderr > ' + output_args['output']
    elif task is "umitools_extract":
        cmd = 'umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN' \
              ' --stdin ' + input_args['path_to_read1'] + \
              ' --read2-in ' + input_args['path_to_read2'] + \
              ' --stdout ' + output_args['output'] + \
              ' --read2-stdout --filter-cell-barcode --error-correct-cell' \
              ' --whitelist=' + input_args['whitelist']
    elif task is "STAR_mapping":
        cmd = 'STAR --runThreadN 32 --genomeDir ' + genome_index + \
              ' --readFilesIn ' + input_args['extracted'] + \
              ' --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout' \
              ' --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical' \
              ' --outFilterMismatchNmax 6 --outSAMtype BAM SortedByCoordinate ' \
              '--outFileNamePrefix ' + output_args['out_prefix'] + ' --outReadsUnmapped Fastx'
    elif task is "featurecounts":
        cmd = 'featureCounts -s 1 -a ' + genome_gtf + ' -o ' + output_args['output'] + \
              ' -R BAM ' + input_args['mapped'] + ' -T 32'
    elif task is "samtools_sort":
        cmd = 'samtools sort ' + input_args['input'] + ' -o ' + output_args['output']
    elif task is "samtools_index":
        cmd = 'samtools index ' + input_args['input']
    elif task is "umitools_count":
        cmd = 'umi_tools count --per-gene --gene-tag=XT --per-cell -I ' + input_args['input'] + \
              ' -S ' + output_args['output']
    elif task is "nuniquemapped":
        cmd = 'samtools view -F4 ' + input_args['input'] + ' | cut -f 2 -d \'_\' |sort|uniq -c > ' + \
              output_args['output']
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
    print(task + lib + '...')
    return subprocess.call(cmd, shell=True)


def do_parallel(src_dir=None, dst_dir=None, task=None, overwrite=True, thread=1, genome_index=None, genome_gtf=None):

    libs = get_libs(src_dir)
    cmd_all = dict()
    for l in libs:
        # Parse input/output args
        input_args, output_args = parse_input_output(src_dir, dst_dir, l, task=task)

        # Construct commands. Store to cmd_all for (parallel) iter run by mp.Pool.
        if overwrite:
            cmd = parse_command(input_args, output_args, task=task, genome_index=genome_index, genome_gtf=genome_gtf)
            cmd_all[(l, task)] = cmd
        else:
            if not os.path.exists(output_args['output']):
                cmd = parse_command(input_args, output_args, task=task, genome_index=genome_index,
                                    genome_gtf=genome_gtf)
                cmd_all[(l, task)] = cmd
            else:
                continue

    # Parallel run by pool
    pool = mp.Pool(processes=thread)
    if len(cmd_all) is not 0:
        pool.map(work, cmd_all)
    print(task + ': finished.')


def wash_whitelist(src_dir, dst_dir, parent_dir, task="wash_whitelist", overwrite=True):

    # Read barcode ground truth list
    barcode_ground_truth_raw = pd.read_excel(os.path.join(parent_dir, 'barcode_ground_truth_checklist.xlsx'))
    barcode_ground_truth = barcode_ground_truth_raw['Primer_sequence'].str.extract(
        r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})', expand=False)

    libs = get_libs(src_dir)
    for l in libs:
        # Display progress
        print('Washing barcode whitelist. Library: ' + l)

        # Parse input/output args
        input_args, output_args = parse_input_output(src_dir, dst_dir, l, task=task)

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
                print(l + "whitelist_washed.txt exists. Skip.")


def main():

    # Source dir: sequencing reads data
    src_dir = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting/'

    # Destination dir: mapping results
    dst_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/mapping/'

    # Parent working dir
    parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM23/'

    # Path to genome annotation and index
    genome_gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM23/gencode.vM23.chr_patch_hapl_scaff.annotation.' \
                 'extended.gtf'
    genome_index = '/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM23/'

    # Create output directory if not exist.barcode_ground_truth
    os.chdir(dst_dir)
    for out in get_libs(src_dir):
        if not os.path.exists(os.path.join(dst_dir, out)):
            os.makedirs(os.path.join(dst_dir, out))

    # STEP 1: umi_tools whitelist
    do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="umitools_whitelist", thread=32)

    # STEP 2: wash whitelist
    wash_whitelist(src_dir=src_dir, dst_dir=dst_dir, parent_dir=parent_dir, task="wash_whitelist", overwrite=True)

    # STEP 3: umi_tools extract
    do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="umitools_extract", thread=32)

    # STEP 4: STAR mapping
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="STAR_mapping", genome_index=genome_index)  # Default thread = 1

    # STEP 5: featureCounts
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="featurecounts", genome_gtf=genome_gtf)  # Default thread = 1

    # STEP 6: samtools sort
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="samtools_sort", thread=16)

    # STEP 7: samtools index
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="samtools_index", thread=16)

    # STEP 8: umitools count
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="umitools_count", thread=16)

    # STEP 9: umitools count
    # do_parallel(src_dir=src_dir, dst_dir=dst_dir, task="nuniquemapped", thread=16)


if __name__ == '__main__':
    main()