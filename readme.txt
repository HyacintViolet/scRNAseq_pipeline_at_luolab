scRNAseq reads mapping pipeline: modified based on umi_tools protocol

1. Edit params (i.e. input/output dirs, parallel threads) in mapping_pipeline_parallel.py. Run.

2. gunzip_counts.py

3. add_suffix.py

4. extract_mapping_stats.py

5. concatenate_counts.py

6. concatenate_mapping_stats.py

# R pipeline starts here

7. rename_gene_id.py




Fix annotation pipeline:

1. Prepare fixed gtf and bed annotations:
   a. for gtf: remove non-primary assembly entries (first column is not chrX, i.e. KZ28...)
      awk -F "\t" 'BEGIN{OFS="\t"} $1~"chr"{print $0}' gencode.vM23.chr_hapl_scaffold.annotation.gtf \
      > gencode.vM23.chr.annotation.gtf
   b. for bed: convert gencode_vM23.chr.annotation.gtf with BEDOPS gtf2bed (this also performs sorting, which is
      required by 'bedtools closest' operation)
      - (BEDOPS gtf2bed checks GTF column 9 and requires BOTH gene_id and transcript_id entries, which is not satisfied
      by original gencode gtf files. Generate a fixed file by adding placeholders: transcript_id "";.
      See this thread: www.biostars.org/p/206342)
      awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.vM23.chr.annotation.gtf
      > gencode.vM23.chr.annotation.fixed.gtf
      - Now convert with gtf2bed:
      gtf2bed < gencode.vM23.chr.annotation.fixed.gtf > gencode.vM23.chr.annotation.fixed.bed

2. Run pipe: bedtools_intersect_closest.py

3. Fix bed file errors: fix_bedfiles.py TODO: code to count lines -> num_lines.csv

