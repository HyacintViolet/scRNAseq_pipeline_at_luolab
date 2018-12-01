# It turns out pycharm has a bug that prevents the script from searching into default PATH. Need to manually set
# environment variable PATH.
# Here's how: Run -> Edit Configurations -> Environment Variables -> manually add PATH entry; from terminal: echo PATH,
# copy and paste into value.


# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 1 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

import os
import re
# import logging
import multiprocessing as mp
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


def work(cmd):
    return subprocess.call(cmd, shell=True)

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
    for out in folder_name_list:
        if not os.path.exists(os.path.join(dst, out)):
            os.makedirs(os.path.join(dst, out))

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
    # Command to use: umi_tools whitelist.

    # Generate list of strings as commands
    cmd_whitelist = []
    for s in os.listdir(src):

        # Setting up input/output directory
        input_dir = os.path.join(src, s)
        out_dir = os.path.join(dst, s)

        # Fetch file name. Read 1 suffix: _2.fq.gz, read 2 suffix: _1.fq.gz.
        for item in os.listdir(input_dir):
            if item.endswith('2.fq.gz'):
                read1_file_name = item
            elif item.endswith('1.fq.gz'):
                read2_file_name = item

        # Construct commands
        cmd_this_whitelist = 'umi_tools whitelist --stdin ' + os.path.join(input_dir, read1_file_name) +\
                             ' --bc-pattern=CCCCCCCCNNNNNNNN --set-cell-number=80 --plot-prefix=cell_num_80 -v 1' \
                             ' --log2stderr > ' + out_dir + 'whitelist80.txt'
        cmd_whitelist.append(cmd_this_whitelist)

    pool = mp.Pool(12)
    pool.map(work, cmd_whitelist)
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

    # STEP 2: Extract barcode and UMIs and add to read names
    # Command to use: umi_tools extract. Because this step takes considerable amount of time but only runs on 1 thread,
    # it is critical to parallelize this command for a speed boost.

    cmd_extract = []
    for s in os.listdir(src):

        # Setting up input/output directory
        input_dir = os.path.join(src, s)
        out_dir = os.path.join(dst, s)

        # Grab folder name prefix
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', s)
        prefix = match.group(1)

        # Fetch file name. Read 1 suffix: _2.fq.gz, read 2 suffix: _1.fq.gz.
        for item in os.listdir(input_dir):
            if item.endswith('2.fq.gz'):
                read1_file_name = item
            elif item.endswith('1.fq.gz'):
                read2_file_name = item

        # Construct output file name
        out_name_extract = '_'.join([prefix, 'extracted.fq.gz'])
        extract_out = os.path.join(out_dir, out_name_extract)

        # Construct command
        cmd_this_extract = 'umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN' \
                           ' --stdin ' + os.path.join(input_dir, read1_file_name) + \
                           ' --read2-in ' + os.path.join(input_dir, read2_file_name) + \
                           ' --stdout ' + extract_out + \
                           ' --read2-stdout --filter-cell-barcode --error-correct-cell' \
                           ' --whitelist=' + out_dir + prefix + '_whitelist_washed.txt'
        cmd_extract.append(cmd_this_extract)

    # Parallel run by Pool
    pool = mp.Pool(12)
    pool.map(work, cmd_extract)
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
        prefix = match.group(1)
        out_name_extract = '_'.join([prefix, 'extracted.fq.gz'])

        command_star_mapping = 'STAR --runThreadN 32 --genomeDir ' + genome_index + \
                               ' --readFilesIn ' + out_name_extract + \
                               ' --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout' + \
                               ' --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical' + \
                               ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + prefix + '_'
        p = subprocess.Popen(command_star_mapping, shell=True)
        p.wait()

        # STEP 4: Assign reads to genes
        # Command to use: featureCounts, samtools sort, samtools index
        # Construct commands for featureCounts. DON'T PARALLELIZE.
        aligned_out = '_'.join([prefix, 'Aligned.sortedByCoord.out.bam'])
        featurecounts_out = '_'.join([prefix, 'gene_assigned'])
        cmd_featurecounts = 'featureCounts -a ' + genome_anno + ' -o ' + featurecounts_out + \
                            ' -R BAM ' + aligned_out + ' -T 32'
        p = subprocess.Popen(cmd_featurecounts, shell=True)
        p.wait()

    print('STAR & featureCounts: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 4 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 5 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # samtools sort parallelized
    # Generate list of strings as commands.
    cmd_sort = []
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Grab folder name
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        prefix = match.group(1)

        # Input file name for samtools sort
        featurecounts_bam_out = '_'.join([prefix, 'Aligned.sortedByCoord.out.bam.featureCounts.bam'])
        sort_in = os.path.join(out_dir, featurecounts_bam_out)

        # Output file name
        out_name_samtools = '_'.join([prefix, 'assigned_sorted.bam'])
        sort_out = os.path.join(out_dir, out_name_samtools)

        # Construct commands
        cmd_this_sort = 'samtools sort ' + sort_in + ' -o ' + sort_out
        cmd_sort.append(cmd_this_sort)

    # Parallel run by Pool
    pool = mp.Pool(12)
    pool.map(work, cmd_sort)
    print('samtools sort: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 5 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 6 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # samtools index parallelized
    # Generate list of strings as commands.
    cmd_index = []
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Change directory to output folder
        os.chdir(out_dir)

        # Grab folder name
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        prefix = match.group(1)

        # Input file name for samtools index
        out_name_samtools = '_'.join([prefix, 'assigned_sorted.bam'])
        index_in = os.path.join(out_dir, out_name_samtools)

        # Construct command
        cmd_this_index = 'samtools index ' + index_in
        cmd_index.append(cmd_this_index)

    # Parallel run by Pool
    pool = mp.Pool(12)
    pool.map(work, cmd_index)
    print('samtools index: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 6 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 7 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # STEP 5: Count UMIs per gene per cell. Command to use: umi_tools count.
    # Generate list of strings as commands.
    cmd_count = []
    for out in os.listdir(dst)[8:]:

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Grab folder name
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        prefix = match.group(1)

        # Input file name for umi_tools count
        out_name_samtools = '_'.join([prefix, 'assigned_sorted.bam'])
        samtools_in = os.path.join(out_dir, out_name_samtools)

        # Output file name
        out_name_count = '_'.join([prefix, 'counts.tsv.gz'])
        count_out = os.path.join(out_dir, out_name_count)

        # Construct commands
        cmd_this_count = 'umi_tools count --per-gene --gene-tag=XT --per-cell -I ' + samtools_in + ' -S ' + count_out
        cmd_count.append(cmd_this_count)

    # Parallel run by Pool
    pool = mp.Pool(12)
    pool.map(work, cmd_count)
    print('umi_tools count: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 7 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 8 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # Construct commands
    cmd_nuniquemap = []
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Grab folder name
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)
        prefix = match.group(1)

        aligned_out = '_'.join([prefix, 'Aligned.sortedByCoord.out.bam'])
        nuniquemap_in = os.path.join(out_dir, aligned_out)
        out_name_nuniquemap = '_'.join([prefix, 'Nuniqmapped.txt'])
        nuniquemap_out = os.path.join(out_dir, out_name_nuniquemap)

        cmd_this_nuniquemap = 'samtools view -F4 ' + nuniquemap_in + ' | cut -f 2 -d \'_\' |sort|uniq -c > ' + \
                              nuniquemap_out
        cmd_nuniquemap.append(cmd_this_nuniquemap)

    # Parallel run by Pool
    pool = mp.Pool(32)
    pool.map(work, cmd_nuniquemap)
    print('samtools view extract Nuniqmapped: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 8 ENDS
# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
