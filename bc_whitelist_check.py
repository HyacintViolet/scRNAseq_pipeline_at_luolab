import os, re
import pandas as pd

# Setup working directory
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
os.chdir(parent_wd)
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4_def/'

# Read ground truth barcode list
barcode_ground_truth_raw = pd.read_excel('barcode_ground_truth_checklist.xlsx')

# Extract ground truth barcodes
barcode_ground_truth = barcode_ground_truth_raw['Primer_sequence'].str.extract(r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})', expand=False)

# Iteratively load in whitelist80.txt, check barcodes.
# Output: presence_in_eighty: found/miss; rank_in_whitlist: [1-80]/NaN; Nreads: /d
for i in range(len(os.listdir(data_wd))):
    os.chdir(os.path.join(data_wd,os.listdir(data_wd)[i]))

    print('Parsing '+os.listdir(data_wd)[i])

    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)_vM4_def$', os.listdir(data_wd)[i])

    # Read whitelist80, adjust column display position
    whitelist80 = pd.read_csv('whitelist80.txt', sep="\t", names = ['cell', 'candidate', 'Nreads', 'Ncandidate'])

    # reordered_col = ['cell', 'Nreads', 'candidate', 'Ncandidate']
    # whitelist80 = whitelist80.reindex(columns=reordered_col)

    # Sort by 'Nreads' in descending order
    sorted_whitelist80 = whitelist80.sort_values('Nreads', ascending=False)
    sorted_whitelist80 = sorted_whitelist80.reset_index(drop=True)
    sorted_whitelist80['ranking'] = range(1,len(sorted_whitelist80)+1)

    reindexed_barcode_ground_truth = pd.Series(barcode_ground_truth.index.values, index=barcode_ground_truth)

    ground_truth_dict = reindexed_barcode_ground_truth.to_dict()

    # pandas mapping only works on series.
    presence_in_ground_truth = sorted_whitelist80.cell.map(ground_truth_dict)
    sorted_whitelist80['presence_in_ground_truth'] = presence_in_ground_truth

    output = sorted_whitelist80[['cell', 'presence_in_ground_truth', 'ranking', 'Nreads']]

    # Output
    out_nameparts = [match.group(1), 'barcode_checklist.txt']
    out_filename = '_'.join(out_nameparts)
    output.to_csv(out_filename, sep="\t", index=False)