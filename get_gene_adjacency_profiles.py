# import gffutils
# gffutiles creates a generator object for accessing gtf data.
# However, this prevents access to previous entries during iteration.

import pandas as pd
import re


def get_next_feature_index_in_strand(gtf, feature_idx):
    feature_idx_next = None
    if gtf.loc[feature_idx, "strand"] is '+':
        feature_idx_next = search_forward(gtf, feature_idx)
    elif gtf.loc[feature_idx, "strand"] is '-':
        feature_idx_next = search_backward(gtf, feature_idx)
    return feature_idx_next


def search_forward(gtf, feature_idx):
    feature_idx_next = None
    for j in range(feature_idx + 1, len(gtf), 1):
        if gtf.loc[j, "seqname"] is gtf.loc[feature_idx, "seqname"] \
                and gtf.loc[j, "strand"] is gtf.loc[feature_idx, "strand"]\
                and gtf.loc[j, "start"] > gtf.loc[feature_idx, "end"]:
            feature_idx_next = j
            break
        else:
            feature_idx_next = None
    return feature_idx_next


def search_backward(gtf, feature_idx):
    feature_idx_next = None
    for j in range(feature_idx - 1, -1, -1):
        if gtf.loc[j, "seqname"] is gtf.loc[feature_idx, "seqname"] \
                and gtf.loc[j, "strand"] is gtf.loc[feature_idx, "strand"]\
                and gtf.loc[j, "end"] < gtf.loc[feature_idx, "start"]:
            feature_idx_next = j
            break
        else:
            feature_idx_next = None
    return feature_idx_next


def extract_attribute(attributes, attr_type=None):
    attribute = None
    if attr_type is "gene_name":
        attribute = re.search('gene_name "([^;]*)"', attributes)
    elif attr_type is "gene_id":
        attribute = re.search('gene_id "([^;]*)"', attributes)
    return attribute.group(1)


def main():
    path_to_gtf = '/media/luolab/ZA1BT1ER/raywang/annotation/Mouse/vM21/' \
                  'gencode.vM21.chr_primary.annotation.gene_only.gtf'
    gtf = pd.read_table(path_to_gtf, header=None)
    gtf.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

    df = pd.DataFrame(data=None, index=None, columns=["gene_name", "gene_id", "strand",
                                                      "gene_name_next", "gene_id_next", "strand_next",
                                                      "distance_to_next"])

    for feature_idx in range(0, len(gtf), 1):

        distance = None

        seqname = gtf.loc[feature_idx, "seqname"]
        start = gtf.loc[feature_idx, "start"]
        end = gtf.loc[feature_idx, "end"]
        strand = gtf.loc[feature_idx, "strand"]
        attributes = gtf.loc[feature_idx, "attribute"]
        gene_name = extract_attribute(attributes, attr_type="gene_name")
        gene_id = extract_attribute(attributes, attr_type="gene_id")

        next_feature_idx = get_next_feature_index_in_strand(gtf, feature_idx)

        if next_feature_idx is not None:
            seqname_next = gtf.loc[next_feature_idx, "seqname"]
            start_next = gtf.loc[next_feature_idx, "start"]
            end_next = gtf.loc[next_feature_idx, "end"]
            strand_next = gtf.loc[next_feature_idx, "strand"]
            attributes_next = gtf.loc[next_feature_idx, "attribute"]
            gene_name_next = extract_attribute(attributes_next, attr_type="gene_name")
            gene_id_next = extract_attribute(attributes_next, attr_type="gene_id")

            if strand is '+':
                distance = start_next - end
            elif strand is '-':
                distance = start - end_next

        if distance is not None:
            adjacent_profile = pd.DataFrame(
                data=[[gene_name, gene_id, strand, gene_name_next, gene_id_next, strand_next, distance]], index=None,
                columns=["gene_name", "gene_id", "strand", "gene_name_next", "gene_id_next", "strand_next",
                         "distance_to_next"])
            df = df.append(adjacent_profile, ignore_index=True)
            print(' '.join(['Parsed line', str(feature_idx),gene_name , 'distance to next feature is', str(distance)]))
        df.to_csv('/media/luolab/ZA1BT1ER/yanting/vM21/adjacent_profiles.txt', sep="\t", index=False, index_label=False)

    # db = gffutils.create_db(path_to_gtf,
    #                         ":memory:",
    #                         keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)
    #
    # df = pd.DataFrame(data=None, index=None, columns=["gene_name", "stable_id", "strand", "distance"])
    # features = db.all_features()

    # features = db.all_features()
    #
    # for feature in features:
    #     gene_name = feature.attributes["gene_name"]
    #     stable_id = feature.attributes["gene_id"]
    #     strand_this = feature.strand
    #     next_feature_in_strand = get_next_feature_in_strand(feature)
    # # Refer to this thread for basic usage: https://bioinformatics.erc.monash.edu/home/kirill/gfftalk/

    # # Abandoned usage of gffutils


if __name__ == '__main__':
    main()
