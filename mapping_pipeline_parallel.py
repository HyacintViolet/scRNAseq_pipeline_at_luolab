# It turns out pycharm has a bug that prevents the script from searching into default PATH. Need to manually set
# environment variable PATH.
# Here's how: Run -> Edit Configurations -> Environment Variables -> Manually add PATH entry; from terminal: echo PATH,
# copy and paste into value.

# ----------------------------------------------------------------------------------------------------------------------
# THE FULL PARALLEL DOESN'T WORK YET. THE SECOND PARALLEL COMMAND WAS STARTED, BEFORE THE FIRST PARALLEL COMMAND WAS
# FINISHED. NEED TO FIX WAITING ISSUE.
#
# CURRENT SOLUTION: RUN PIPE CHUNK BY CHUNK.
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 1 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

import os
import re
# import logging
import subprocess
import pandas as pd


def wash_whitelist(out_dir, bc_ground_truth, match):
    print('Washing ' + out_dir)
    os.chdir(out_dir)

    # Read whitelist80
    whitelist80 = pd.read_csv('whitelist80.txt', sep="\t", names=['cell', 'candidate', 'Nreads', 'Ncandidate'])

    # Remove rows whose 'cell' value is not found in barcode_ground_truth
    whitelist_washed = whitelist80[whitelist80['cell'].isin(bc_ground_truth)]

    # Reset index
    whitelist_washed = whitelist_washed.reset_index(drop=True)

    # Output
    out_name_parts = [match.group(1), 'whitelist_washed.txt']
    filename = '_'.join(out_name_parts)
    whitelist_washed.to_csv(filename, sep="\t", index=False, header=False)


# ----------------------------------------------------------------------------------------------------------------------

def main():

    # Source dir: sequencing reads data
    src = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting_previous/'
    folder_name_list = os.listdir(src)

    # Destination dir: mapping results
    dst = '/media/luolab/ZA1BT1ER/yanting/vM19/yanting_previous/'

    # Path to genome annotation and index
    genome_anno = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf'
    genome_index = '/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM19/'

    # Parent working dir: to write aggregate results
    parent_wd = '/media/luolab/ZA1BT1ER/yanting/'

    # Read barcode ground truth list
    barcode_ground_truth_raw = pd.read_excel(os.path.join(parent_wd, 'barcode_ground_truth_checklist.xlsx'))
    barcode_ground_truth = barcode_ground_truth_raw['Primer_sequence'].str.extract(
        r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})', expand=False)

    # Create output directory if not exist.
    os.chdir(dst)
    for out_dir in folder_name_list:
        if not os.path.exists(os.path.join(dst, out_dir)):
            os.makedirs(os.path.join(dst, out_dir))

    # ------------------------------------------------------------------------------------------------------------------

    # # logging
    # log = parent_wd + 'pipeline_running.log'
    # logging.basicConfig(filename=log, level=logging.DEBUG)
    #
    # # console handler
    # console = logging.StreamHandler()
    # console.setLevel(logging.ERROR)
    # logging.getLogger("").addHandler(console)

    # ------------------------------------------------------------------------------------------------------------------
    # STEP 1: Identify correct cell barcodes.
    # Command to use: umi_tools whitelist. Script to use: whitelist_wash.py
    # This is the parent directory for all output directories
    processes_whitelist = set()
    max_processes_whietlist = 16
    for s in os.listdir(src):

        # Setting up input/output directory
        input_dir = os.path.join(src, s)
        out_dir = os.path.join(dst, s)

        # Change directory to output folder
        os.chdir(out_dir)

        # Fetch file name. Read 1 suffix: _2.fq.gz, read 2 suffix: _1.fq.gz.
        for item in os.listdir(input_dir):
            if item.endswith('2.fq.gz'):
                read1_file_name = item
            elif item.endswith('1.fq.gz'):
                read2_file_name = item

        # # LOGGING NOT WORKING. Fix later.
        # try:
        #     read1_file_name
        # except NameError:
        #     logging.debug('Library ' + s + ' not parsed. Read 1 Missing.')

        # Construct command and execute
        command_whitelist = 'umi_tools whitelist --stdin ' + os.path.join(input_dir, read1_file_name) +\
                            ' --bc-pattern=CCCCCCCCNNNNNNNN --set-cell-number=80 --plot-prefix=cell_num_80 -v 1' \
                            ' --log2stderr > whitelist80.txt'

        # Parallel running umi_tools whitelist
        processes_whitelist.add(subprocess.Popen(command_whitelist, shell=True))
        if len(processes_whitelist) >= max_processes_whietlist:
            os.wait()
            processes_whitelist.difference_update([p for p in processes_whitelist if p.poll() is not None])

    print('umi_tools whitelist: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 1 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 2 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    for out in os.listdir(dst):

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)

        # Wash whitelist, output whitelist_washed.txt
        out_dir = os.path.join(dst, out)
        wash_whitelist(out_dir, barcode_ground_truth, match)

    print('Wash whitelist: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 2 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 3 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # STEP 2: Extract barcodes and UMIs and add to read names
    # Command to use: umi_tools extract. Because this step takes considerable amount of time but takes only 1 thread
    # to run, we parallel it in 4 batches. (8 library each parallel run, 4 runs)

    processes_extract = set()
    max_processes_extract = 12
    for s in os.listdir(src):

        # Setting up input/output directory
        input_dir = os.path.join(src, s)
        out_dir = os.path.join(dst, s)

        # Change directory to output folder
        os.chdir(out_dir)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', s)
        out_prefix = match.group(1)

        # Fetch file name. Read 1 suffix: _2.fq.gz, read 2 suffix: _1.fq.gz.
        for item in os.listdir(input_dir):
            if item.endswith('2.fq.gz'):
                read1_file_name = item
            elif item.endswith('1.fq.gz'):
                read2_file_name = item

        # Construct output file name
        out_name_extract = '_'.join([out_prefix, 'extracted.fq.gz'])

        # Construct command: umi_tools extract
        command_extract = 'umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN' \
                          ' --stdin ' + os.path.join(input_dir, read1_file_name) + \
                          ' --read2-in ' + os.path.join(src, input_dir, read2_file_name) + \
                          ' --stdout ' + out_name_extract + \
                          ' --read2-stdout --filter-cell-barcode --error-correct-cell' \
                          ' --whitelist=' + out_prefix + '_whitelist_washed.txt'

        # Parallel running umi_tools extract
        processes_extract.add(subprocess.Popen(command_extract, shell=True))
        if len(processes_extract) >= max_processes_extract:
            os.wait()
            processes_extract.difference_update([p for p in processes_extract if p.poll() is not None])

    print('umi_tools extract: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 3 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 4 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # STEP 3: Map reads
    # Command to use: STAR
    # Construct command and execute. DON'T PARALLELIZE.
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Change directory to output folder
        os.chdir(out_dir)

        print('Mapping ' + out_dir)

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        out_prefix = match.group(1)
        out_name_extract = '_'.join([out_prefix, 'extracted.fq.gz'])

        command_star_mapping = 'STAR --runThreadN 32 --genomeDir ' + genome_index + \
                               ' --readFilesIn ' + out_name_extract + \
                               ' --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout' + \
                               ' --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical' + \
                               ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + out_prefix + '_'
        p = subprocess.Popen(command_star_mapping, shell=True)
        p.wait()

        # STEP 4: Assign reads to genes
        # Command to use: featureCounts, samtools sort, samtools index
        # Construct commands for featureCounts. DON'T PARALLELIZE.
        out_name_feature_counts = '_'.join([out_prefix, 'gene_assigned'])
        in_name_aligned = '_'.join([out_prefix, 'Aligned.sortedByCoord.out.bam'])
        command_feature_counts = 'featureCounts -a ' + genome_anno + ' -o ' + out_name_feature_counts + \
                                 ' -R BAM ' + in_name_aligned + ' -T 32'
        p = subprocess.Popen(command_feature_counts, shell=True)
        p.wait()

    print('STAR & featureCounts: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 4 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 5 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # samtools sort parallelized
    processes_sort = set()
    max_processes_sort = 8
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Change directory to output folder
        os.chdir(out_dir)

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        out_prefix = match.group(1)

        # Input file name for samtools sort
        in_name_feature_counted = '_'.join([out_prefix, 'Aligned.sortedByCoord.out.bam.featureCounts.bam'])

        # Output file name
        out_name_samtools = '_'.join([out_prefix, 'assigned_sorted.bam'])

        # Construct commands
        command_samtools_sort = 'samtools sort ' + in_name_feature_counted + ' -o ' + out_name_samtools

        # Parallel running samtools sort
        processes_sort.add(subprocess.Popen(command_samtools_sort, shell=True))
        if len(processes_sort) >= max_processes_sort:
            os.wait()
            processes_sort.difference_update([p for p in processes_sort if p.poll() is not None])

    # print('samtools sort: finished.') # This doesn't work. Alert appears before job finishing.

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 5 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 6 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # samtools index parallelized
    processes_index = set()
    max_processes_index = 16
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Change directory to output folder
        os.chdir(out_dir)

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        out_prefix = match.group(1)

        # Input & output file name for samtools index
        out_name_samtools = '_'.join([out_prefix, 'assigned_sorted.bam'])

        # Construct command
        command_samtools_index = 'samtools index ' + out_name_samtools

        # Parallel running samtools index
        processes_index.add(subprocess.Popen(command_samtools_index, shell=True))
        if len(processes_index) >= max_processes_index:
            os.wait()
            processes_index.difference_update([p for p in processes_index if p.poll() is not None])

    print('samtools index: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 6 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 7 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------
# SOMETHING IS WRONG WITH THIS PARALLELIZATION.

    # STEP 5: Count UMIs per gene per cell
    # Command to use: umi_tools count. Parallelized
    processes_count = set()
    max_processes_count = 16
    for out_dir in os.listdir(dst):

        # Change directory to output folder
        os.chdir(os.path.join(dst, out_dir))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out_dir)
        out_prefix = match.group(1)

        # Input file name for umi_tools count
        out_name_samtools = '_'.join([out_prefix, 'assigned_sorted.bam'])

        # Output file name
        out_name_count = '_'.join([out_prefix, 'counts.tsv.gz'])

        # Construct command
        command_count = 'umi_tools count --per-gene --gene-tag=XT --per-cell -I ' + out_name_samtools + ' -S ' + \
                        out_name_count
        processes_count.add(subprocess.Popen(command_count, shell=True))
        if len(processes_count) >= max_processes_count:
            os.wait()
            processes_count.difference_update([p for p in processes_count if p.poll() is not None])

    print('umi_tools count: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 7 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 8 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    processes_uniqmapped = set()
    max_processes_uniqmapped = 16
    for out_dir in os.listdir(dst):

        # Change directory to output folder
        os.chdir(os.path.join(dst, out_dir))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out_dir)
        out_prefix = match.group(1)

        in_name_aligned = '_'.join([out_prefix, 'Aligned.sortedByCoord.out.bam'])

        # Extract uniquely mapped reads *** this can be paralleled with other commands ***
        out_nuniqmapped = '_'.join([out_prefix, 'Nuniqmapped.txt'])
        command_nuniqmapped = 'samtools view -F4 ' + in_name_aligned + ' | cut -f 2 -d \'_\' |sort|uniq -c > ' +\
                              out_nuniqmapped

        # Nuniqmapped extract parallelized
        processes_uniqmapped.add(subprocess.Popen(command_nuniqmapped, shell=True))
        if len(processes_uniqmapped) >= max_processes_uniqmapped:
            os.wait()
            processes_uniqmapped.difference_update([p for p in processes_uniqmapped if p.poll() is not None])

    print('samtools extract Nuniqmapped: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 8 ENDS
# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
