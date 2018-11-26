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
import subprocess
import pandas as pd


def wash_whitelist(output_folder_path, barcode_ground_truth, match):
    print('Washing ' + output_folder_path)
    os.chdir(output_folder_path)

    # Read whitelist80
    whitelist80 = pd.read_csv('whitelist80.txt', sep="\t", names=['cell', 'candidate', 'Nreads', 'Ncandidate'])

    # Remove rows whose 'cell' value is not found in barcode_ground_truth
    whitelist_washed = whitelist80[whitelist80['cell'].isin(barcode_ground_truth)]

    # Reset index
    whitelist_washed = whitelist_washed.reset_index(drop=True)

    # Output
    out_nameparts = [match.group(1), 'whitelist_washed.txt']
    out_filename = '_'.join(out_nameparts)
    whitelist_washed.to_csv(out_filename, sep="\t", index=False, header=False)


# ----------------------------------------------------------------------------------------------------------------------

def main():

    # Setup input data path
    path_to_seq_data = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting_181023/'
    input_folder_list = os.listdir(path_to_seq_data)

    # Output (mapping) parent path
    path_to_mapping_directory = '/media/luolab/ZA1BT1ER/yanting/vM4_def_2/'

    # Path to genome annotation and index
    path_to_genome_anno = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/gencode.vM4.annotation.gtf'
    path_to_genome_index = '/media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM4_def/'

    # Read ground truth barcode list
    parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
    barcode_ground_truth_raw = pd.read_excel(os.path.join(parent_wd, 'barcode_ground_truth_checklist.xlsx'))
    barcode_ground_truth = barcode_ground_truth_raw['Primer_sequence'].str.extract(
        r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})', expand=False)

    # Create output directory if not exist.
    os.chdir(path_to_mapping_directory)
    for output_folder in input_folder_list:
        output_folder = output_folder + '_vM4_def'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    # ------------------------------------------------------------------------------------------------------------------
    # STEP 1: Identify correct cell barcodes.
    # Command to use: umi_tools whitelist. Script to use: whitelist_wash.py
    # This is the parent directory for all output directories
    output_folder_list = os.listdir(path_to_mapping_directory)
    processes_whitelist = set()
    max_processes_whietlist = 8
    for input_folder in input_folder_list:

        output_folder = input_folder + '_vM4_def'

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        input_file_name_prefix = '_'.join([match.group(1), match.group(2), match.group(3), match.group(4)])

        read1_file_name = '_'.join([input_file_name_prefix, '2.fq.gz'])

        # Construct command and execute
        command_whitelist = 'umi_tools whitelist --stdin ' + os.path.join(path_to_seq_data,
                                                                          input_folder,
                                                                          read1_file_name) +\
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

    for output_folder in output_folder_list:

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)

        # Wash whitelist, output whitelist_washed.txt
        output_folder_path = os.path.join(path_to_mapping_directory, output_folder)
        wash_whitelist(output_folder_path, barcode_ground_truth, match)

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
    max_processes_extract = 8
    for input_folder in input_folder_list:

        output_folder = input_folder + '_vM4_def'

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        input_file_name_prefix = '_'.join([match.group(1), match.group(2), match.group(3), match.group(4)])
        out_file_name_prefix = match.group(1)

        read1_file_name = '_'.join([input_file_name_prefix, '2.fq.gz'])

        read2_file_name = '_'.join([input_file_name_prefix, '1.fq.gz'])

        # Construct output file name
        out_name_extract = '_'.join([out_file_name_prefix, 'extracted.fq.gz'])

        # Construct command: umi_tools extract
        command_extract = 'umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN --stdin ' + os.path.join(path_to_seq_data,
                                                                                                    input_folder,
                                                                                                    read1_file_name) + \
                          ' --read2-in ' + os.path.join(path_to_seq_data, input_folder, read2_file_name) + \
                          ' --stdout ' + out_name_extract + \
                          ' --read2-stdout --filter-cell-barcode --error-correct-cell' \
                          ' --whitelist=' + out_file_name_prefix + '_whitelist_washed.txt'

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
    for output_folder in output_folder_list:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        out_file_name_prefix = match.group(1)
        out_name_extract = '_'.join([out_file_name_prefix, 'extracted.fq.gz'])

        command_star_mapping = 'STAR --runThreadN 42 --genomeDir ' + path_to_genome_index + \
                               ' --readFilesIn ' + out_name_extract + \
                               ' --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout' + \
                               ' --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical' + \
                               ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + out_file_name_prefix
        os.system(command_star_mapping)

        # STEP 4: Assign reads to genes
        # Command to use: featureCounts, samtools sort, samtools index
        # Construct commands for featureCounts. DON'T PARALLELIZE.
        out_name_feature_counts = '_'.join([out_file_name_prefix, 'gene_assigned'])
        in_name_aligned = '_'.join([out_file_name_prefix, 'Aligned.sortedByCoord.out.bam'])
        command_feature_counts = 'featureCounts -a ' + path_to_genome_anno + ' -o ' + out_name_feature_counts + \
                                 ' -R BAM ' + in_name_aligned + ' -T 32'
        os.system(command_feature_counts)

    print('STAR & featureCounts: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 4 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 5 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # samtools sort parallelized
    processes_sort = set()
    max_processes_sort = 8
    for output_folder in output_folder_list:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        out_file_name_prefix = match.group(1)

        # Input file name for samtools sort
        in_name_feature_counted = '_'.join([out_file_name_prefix, 'Aligned.sortedByCoord.out.bam.featureCounts.bam'])

        # Output file name
        out_name_samtools = '_'.join([out_file_name_prefix, 'assigned_sorted.bam'])

        # Construct commands
        command_samtools_sort = 'samtools sort ' + in_name_feature_counted + ' -o ' + out_name_samtools

        # Parallel running samtools sort
        processes_sort.add(subprocess.Popen(command_samtools_sort, shell=True))
        if len(processes_sort) >= max_processes_sort:
            os.wait()
            processes_sort.difference_update([p for p in processes_sort if p.poll() is not None])

    print('samtools sort: finished.')

# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 5 ENDS
# ----------------------------------------------------------------------------------------------------------------------
# CHUNK 6 STARTS HERE
# ----------------------------------------------------------------------------------------------------------------------

    # samtools index parallelized
    processes_index = set()
    max_processes_index = 16
    for output_folder in output_folder_list:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        out_file_name_prefix = match.group(1)

        # Input & output file name for samtools index
        out_name_samtools = '_'.join([out_file_name_prefix, 'assigned_sorted.bam'])

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

    # STEP 5: Count UMIs per gene per cell
    # Command to use: umi_tools count. Parallelized
    processes_count = set()
    max_processes_count = 16
    for output_folder in output_folder_list:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        out_file_name_prefix = match.group(1)

        # Input file name for umi_tools count
        out_name_samtools = '_'.join([out_file_name_prefix, 'assigned_sorted.bam'])

        # Output file name
        out_name_count = '_'.join([out_file_name_prefix, 'counts.tsv.gz'])

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
    for output_folder in output_folder_list:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        out_file_name_prefix = match.group(1)

        in_name_aligned = '_'.join([out_file_name_prefix, 'Aligned.sortedByCoord.out.bam'])

        # Extract uniquely mapped reads *** this can be paralleled with other commands ***
        out_nuniqmapped = '_'.join([out_file_name_prefix, 'Nuniqmapped.txt'])
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
