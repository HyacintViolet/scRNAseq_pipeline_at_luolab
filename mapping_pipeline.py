# For some unknown reason the os.syste() command never works.

import os, re
#
# app_path = os.path.join('home','luolab','anaconda3','bin','umi_tools')
# os.environ["PATH"] += os.pathsep + app_path

path_to_seq_data = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting_181023/'
input_folder_list = os.listdir(path_to_seq_data)

path_to_mapping_directory = '/media/luolab/ZA1BT1ER/yanting/vM4_def/'
folder_list = os.listdir(path_to_mapping_directory)

for i in range(len(folder_list)):

    input_folder = input_folder_list[i]
    output_folder = input_folder+'_vM4_def'

    os.chdir(os.path.join(path_to_mapping_directory, output_folder))

    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', output_folder)

    input_file_parts = [match.group(1), match.group(2), match.group(3), match.group(4), '2.fq.gz']
    input_file_name = '_'.join(input_file_parts)

    command_whitelist = 'umi_tools whitelist --stdin ' + os.path.join(path_to_seq_data, input_folder, input_file_name) + ' --bc-pattern=CCCCCCCCNNNNNNNN --set-cell-number=80 --plot-prefix=cell_num_80 -v 1 --log2stderr > whitelist80.txt'
    os.system(command_whitelist)

