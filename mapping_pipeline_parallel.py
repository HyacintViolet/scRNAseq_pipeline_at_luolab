import os, re, subprocess
import pandas as pd

# It turns out pycharm has a bug that prevents the script from searching into default PATH. Need to manually set
# environment variable PATH.
# Here's how: Run -> Edit Configurations -> Environment Variables -> Manually add PATH entry; from terminal: echo PATH,
# copy and paste into value.


def wash_whitelist(folder):
    print('Washing ' + folder)
    os.chdir(os.path.join(path_to_mapping_directory, folder))

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

    # This is the parent directory for all output directories
    output_folder_list = os.listdir(path_to_mapping_directory)
    for output_folder in output_folder_list[1:3]:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory,output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        input_file_name_prefix = '_'.join([match.group(1), match.group(2), match.group(3), match.group(4)])

# ----------------------------------------------------------------------------------------------------------------------

    # STEP 1: Identify correct cell barcodes.
    # Command to use: umi_tools whitelist. Script to use: whitelist_wash.py
        read1_file_name = '_'.join([input_file_name_prefix, '2.fq.gz'])

        # Construct command and execute
        command_whitelist = 'umi_tools whitelist --stdin ' + os.path.join(path_to_seq_data, output_folder, read1_file_name) +\
                            ' --bc-pattern=CCCCCCCCNNNNNNNN --set-cell-number=80 --plot-prefix=cell_num_80 -v 1 --log2stderr > whitelist80.txt'
        subprocess.Popen(command_whitelist, shell=True)

        # Wash whitelist, output whitelist_washed.txt
        wash_whitelist(output_folder)

# ----------------------------------------------------------------------------------------------------------------------

    # STEP 2: Extract barcodes and UMIs and add to read names
    # Command to use: umi_tools extract. Because this step takes considerable amount of time but takes only 1 thread
    # to run, we parallel it in 4 batches. (8 library each parallel run, 4 runs)
    # Construct output file name

    for output_folder in output_folder_list[0:9]:

        # Change directory to output folder
        os.chdir(os.path.join(path_to_mapping_directory, output_folder))

        # Grab folder name and construct input file names. For the data in this example, the read1, read2 naming
        # convention is reversed.
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)
        input_file_name_prefix = '_'.join([match.group(1), match.group(2), match.group(3), match.group(4)])
        out_file_name_prefix = match.group(1)

        read1_file_name = '_'.join([input_file_name_prefix, '2.fq.gz'])

        read2_file_name = '_'.join([input_file_name_prefix, '1.fq.gz'])

        out_name_extract = '_'.join([out_file_name_prefix, 'extracted.fq.gz'])

        # Construct command and execute: umi_tools extract
        command_extract = 'umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN --stdin ' + os.path.join(path_to_seq_data, output_folder, read1_file_name) + \
                          ' --read2-in ' + os.path.join(path_to_seq_data, output_folder, read2_file_name) + \
                          ' --stdout ' + out_name_extract + \
                          ' --read2-stdout --filter-cell-barcode --whitelist=' + out_file_name_prefix +'_whitelist_washed.txt'
        subprocess.Popen(command_extract, shell=True)


        # --------------------------------------------------------------------------------------------------------------

        # STEP 3: Map reads
        # Command to use: STAR
        # Construct command and execute
        command_STAR_mapping = 'STAR --runThreadN 32 --genomeDir ' + path_to_genome_index + \
                               ' --readFilesIn ' + out_name_extract + \
                               ' --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterType BySJout' + \
                               ' --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical' + \
                               ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + out_file_name_prefix
        os.system(command_STAR_mapping)

    # ------------------------------------------------------------------------------------------------------------------

        # STEP 4: Assign reads to genes
        # Command to use: featureCounts, samtools sort, samtools index
        # Construct commands for featureCounts
        out_name_featureCounts = '_'.join([out_file_name_prefix, 'gene_assigned'])
        in_name_Aligned = '_'.join([out_file_name_prefix, 'Aligned.sortedByCoord.out.bam'])
        command_featureCounts = 'featureCounts -a ' + path_to_genome_anno + ' -o ' + out_name_featureCounts + ' -R BAM ' +\
                                in_name_Aligned + ' -T 32'
        os.system(command_featureCounts)

        # Construct commands for samtools
        in_name_featureCounted = '_'.join([out_file_name_prefix, 'Aligned.sortedByCoord.out.bam.featureCounts.bam'])
        out_name_samtools = '_'.join([out_file_name_prefix, 'assigned_sorted.bam'])
        command_samtools_sort = 'samtools sort ' + in_name_featureCounted + ' -o ' + out_name_samtools
        os.system(command_samtools_sort)

        command_samtools_index = 'samtools index ' + out_name_samtools
        os.system(command_samtools_index)

    # ------------------------------------------------------------------------------------------------------------------

        # STEP 5: Count UMIs per gene per cell
        # Command to use: umi_tools count
        out_name_count = '_'.join([out_file_name_prefix, 'counts.tsv.gz'])
        command_count = 'umi_tools count --per-gene --gene-tag=XT --per-cell -I ' + out_name_samtools + ' -S ' + \
                        out_name_count
        os.system(command_count)


if __name__ == '__maine__':
    main()