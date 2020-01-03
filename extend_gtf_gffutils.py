# Check this useful thread: https://github.com/daler/gffutils/issues/148

import os
import pandas as pd
import gffutils


parent_dir = '/media/luolab/ZA1BT1ER/yanting/vM21/extension/'

filename_extension_profile = 'extension_profiles.txt'
filename_gtf = 'gencode.vM21.chr_patch_hapl_scaff.annotation.gtf'
filename_gtf_new = 'gencode.vM21.chr_patch_hapl_scaff.annotation.extended.gtf'

path_to_gtf = os.path.join(parent_dir, filename_gtf)
path_to_gtf_new = os.path.join(parent_dir, filename_gtf_new)
path_to_extension_profile = os.path.join(parent_dir, filename_extension_profile)
df = pd.read_table(path_to_extension_profile, sep="\t")

counter = 1

# To get a GTF from a db:
with open(path_to_gtf_new, 'w') as fout:
    for feature in gffutils.DataIterator(path_to_gtf):

        gene_id = feature.attributes['gene_id'][0]
        gene_name = feature.attributes['gene_name'][0]
        if not any(df['gene_id'].str.contains(gene_id)):
            fout.write(str(feature) + '\n')
            print(' '.join(["Line", str(counter), gene_name, "no extension"]))
            counter = counter + 1
            continue

        coord_to_extend = int(df[df['gene_id'].str.match(gene_id)]['coord_to_extend'])
        coord_after_extend = int(df[df['gene_id'].str.match(gene_id)]['coord_after_extend'])

        strand = feature.strand
        start = feature.start
        end = feature.end

        if strand is '+':
            if end == coord_to_extend:
                feature.end = coord_after_extend
                fout.write(str(feature) + '\n')
                print(' '.join(["Line", str(counter), gene_name, "extend from", str(coord_to_extend), "to",
                                str(coord_after_extend)]))
            else:
                fout.write(str(feature) + '\n')
                print(' '.join(["Line", str(counter), gene_name, "no extension"]))
        elif strand is '-':
            if start == coord_to_extend:
                feature.start = coord_after_extend
                fout.write(str(feature) + '\n')
                print(' '.join(["Line", str(counter), gene_name, "extend from", str(coord_to_extend), "to",
                                str(coord_after_extend)]))
            else:
                fout.write(str(feature) + '\n')
                print(' '.join(["Line", str(counter), gene_name, "no extension"]))

        counter = counter + 1

# # Doesn't look like we need these anymore
# path_to_db = os.path.join(parent_dir, "vM21_gtf.db")
#
# fn = gffutils.example_filename(path_to_gtf)
# if not os.path.exists(path_to_db):
#     db = gffutils.create_db(fn,
#                             path_to_db,
#                             keep_order=True,
#                             disable_infer_genes=True, disable_infer_transcripts=True)
