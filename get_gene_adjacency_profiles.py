import gffutils
import pandas as pd

def main():
    path_to_gtf = 'D:/scRNAseq/yanting/vM21/gencode.vM21.chr_primary.annotation.gene_only.gtf'
    db = gffutils.create_db(path_to_gtf,
                            ":memory:",
                            keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)

    df = pd.DataFrame(data=None, index=None, columns=["gene_name", "stable_id", "strand", "distance"])
    for g in db.all_features():
        g.

if __name__ == '__main__':
    main()
