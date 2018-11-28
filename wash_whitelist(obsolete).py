import os, re
import pandas as pd

# Setup working directory
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
path_to_mapping_directory = '/media/luolab/ZA1BT1ER/yanting/vM4_def_2/'

# Read ground truth barcode list
barcode_ground_truth_raw = pd.read_excel(os.path.join(parent_wd, 'barcode_ground_truth_checklist.xlsx'))

# Extract ground truth barcodes
barcode_ground_truth = barcode_ground_truth_raw['Primer_sequence'].str.extract(r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})',
                                                                               expand=False)

for folder in os.listdir(path_to_mapping_directory):

    # Enter mapping data directory
    print('Parsing ' + folder)
    os.chdir(os.path.join(path_to_mapping_directory, folder))

    # Grab folder name
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)

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
