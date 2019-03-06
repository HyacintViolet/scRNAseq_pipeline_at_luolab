import os
import re
import pandas as pd

# Setup working directory
parent_wd = '/media/luolab/ZA1BT1ER/yanting/'
os.chdir(parent_wd)
data_wd = '/media/luolab/ZA1BT1ER/yanting/vM19/yanting_190301/'

# Read ground truth barcode list
bc_ground_truth_raw = pd.read_excel('barcode_ground_truth_checklist.xlsx')

# Extract ground truth barcodes
bc_ground_truth = bc_ground_truth_raw['Primer_sequence'].str.extract(r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})',
                                                                     expand=False)

# Iteratively load in whitelist80.txt, check barcodes.
# Output: presence_in_ground_truth: [0-47]/NaN; rank_in_whitlist: [1-80]; Nreads: /d
for folder in os.listdir(data_wd):
    os.chdir(os.path.join(data_wd, folder))

    print('Parsing '+folder)

    # Get file id for output name
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)

    # Read whitelist80
    whitelist80 = pd.read_csv('whitelist80.txt', sep="\t", names=['cell', 'candidate', 'Nreads', 'Ncandidate'])

    # reordered_col = ['cell', 'Nreads', 'candidate', 'Ncandidate']
    # whitelist80 = whitelist80.reindex(columns=reordered_col)

    # Sort by 'Nreads' in descending order
    sorted_whitelist80 = whitelist80.sort_values('Nreads', ascending=False)
    sorted_whitelist80 = sorted_whitelist80.reset_index(drop=True)
    sorted_whitelist80['ranking'] = range(1, len(sorted_whitelist80)+1)

    # # Compute Nreads + Ncandidate, write to a new column.
    # # In fact this is not very useful, as the candidate barcodes 1 Hamming Distance away from whitelist barcodes
    # # is automatically parsed into the 'umi_tools counts' function. And for sequencing QC, Nuniquemapped/Nreads is
    # # good enough.
    # sorted_whitelist80['Nreads_plus_candidate'] = 0
    # for j in range(len(sorted_whitelist80)):
    #     row = sorted_whitelist80.Ncandidate.iloc[j]
    #     row = row.split(',')
    #     row = list(map(int, row))
    #     # print(row)
    #     sorted_whitelist80['Nreads_plus_candidate'].loc[j] = sum(row) + sorted_whitelist80['Nreads'].loc[j]

    # Swap column
    reindexed_barcode_ground_truth = pd.Series(bc_ground_truth.index.values, index=bc_ground_truth)

    ground_truth_dict = reindexed_barcode_ground_truth.to_dict()

    # pandas mapping only works on series.
    presence_in_ground_truth = sorted_whitelist80.cell.map(ground_truth_dict)
    sorted_whitelist80['presence_in_ground_truth'] = presence_in_ground_truth

    output = sorted_whitelist80[['cell', 'presence_in_ground_truth', 'ranking', 'Nreads']]

    if len(output.presence_in_ground_truth.dropna()) != 48:
        print('Something is missing. Longer whitelist may be needed.\n')

    # Output
    out_name_parts = [match.group(1), 'barcode_checklist.txt']
    out_filename = '_'.join(out_name_parts)
    output.to_csv(out_filename, sep="\t", index=False)
