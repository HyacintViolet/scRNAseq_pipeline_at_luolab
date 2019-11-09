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
   a. for gtf, remove non-primary assembly entries with the following bash line (first column not chrX, i.e. KZ28...):
      - Run command:
      awk -F "\t" 'BEGIN{OFS="\t"} $1~"chr"{print $0}' gencode.vM23.chr_hapl_scaffold.annotation.gtf \
      > gencode.vM23.chr.annotation.gtf
   b. for bed, convert gencode_vM23.chr.annotation.gtf with BEDOPS gtf2bed (this also performs sorting, which is
      required by 'bedtools closest' command)
      - BEDOPS gtf2bed checks GTF column 9 and requires BOTH gene_id and transcript_id entries, which is not satisfied
      by original gencode gtf files. Generate a fixed file by adding placeholders: transcript_id "";
      See this thread: www.biostars.org/p/206342
      - Run command:
      awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.vM23.chr.annotation.gtf
      > gencode.vM23.chr.annotation.fixed.gtf
   c. Now convert with gtf2bed:
      gtf2bed < gencode.vM23.chr.annotation.fixed.gtf > gencode.vM23.chr.annotation.fixed.bed

2. bedtools_intersect_closest.py # run pipe

3. fix_bedfiles.py # fix bed file errors

4. extract_desired_bed_fields.py

5. Extract fields from individual .._extracted.bed files and cat together, steps as follow:
   - Go to mapping directory
   - Run command:
     find . "*_extracted.bed" -type f -exec cat {} \; | gzip -9c > .../YT_extracted_all.txt.gz

6. Further filtering on final output
   - Run command:
     zcat YT_extracted_all.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"} $9<0 && $9 > -20000 {print $0}' | gzip -9c > \
     YT_extracted_less_than_20k.txt.gz
   - Run command:
     zcat YT_extracted_less_than_20k.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"}{print $7, $10, $11, $9}' | gzip -9c > \
     distance_to_gene.txt.gz
