************************************************************************************************************************
scRNAseq reads mapping pipeline: modified based on umi_tools protocol

1. Edit params (i.e. input/output dirs, parallel threads) in mapping_pipeline_parallel.py. Run.

2. gunzip_counts.py

3. add_suffix.py

*4. extract_mapping_stats.py

5. concatenate_counts.py

*6. concatenate_mapping_stats.py

# R pipeline starts here
    i.  Alignment QC: alignment_QC.R
    ii. Cell QC:      cell_QC.R


- R1 (Renamed) count matrix : pituitary_QC.R -> Count matrix -> 7
- R2 (Stable id) count matrix : seurat_pipeline.R

7. rename_gene_id.py



************************************************************************************************************************
Fix (extend) annotation pipeline:

1. Prepare fixed gtf and bed annotations:
   a. for gtf, remove non-primary assembly entries with the following bash line (first column not chrX, i.e. KZ28...):
      - Run command:
      awk 'BEGIN{FS="\t";OFS="\t"} $1~"chr"{print $0}' gencode.vM21.chr_patch_hapl_scaff.annotation.gtf \
      > gencode.vM21.chr_primary.annotation.gtf
   b. for bed, convert gencode_vM21.chr_primary.annotation.gtf with BEDOPS gtf2bed (this also performs sorting, which is
      required by 'bedtools closest' command)
      - BEDOPS gtf2bed checks GTF column 9 and requires BOTH gene_id and transcript_id entries, which is not satisfied
      by original gencode gtf files. Generate a fixed file by adding placeholders: transcript_id "";
      See this thread: www.biostars.org/p/206342
      - Run command:
      awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.vM21.chr_primary.annotation.gtf
      > gencode.vM21.chr_primary.annotation.fixed.gtf
   c. Now convert with gtf2bed:
      gtf2bed < gencode.vM21.chr_primary.annotation.fixed.gtf > gencode.vM21.chr_primary.annotation.fixed.bed
   d. The chromosome gencode.vM21.chr_primary.annotation.fixed.bed is not in desired order. Sort them so that they're
      compatible to run with STAR aligned.bam files:
      bedtools sort -faidx /media/luolab/ZA1BT1ER/raywang/STAR_index_mm10_vM21/chrName.txt -i gencode.vM21.chr_primary.
      annotation.fixed.bed > gencode.vM21.chr_primary.annotation.fixed.sortedByChrom.bed

2. Run bedtools_intersect_closest.py
    # - input intersect: ..._Aligned.sortedByCoord.out.bam
      - output intersect: ..._stranded_nonoverlap.bam
    # - input closest: ..._stranded_nonoverlap.bam
      - output closest: ..._closest.bed

3. fix_bedfiles.py # fix bed file errors

4. extract_desired_bed_fields.py

5. Extract fields from individual .._extracted.bed files and cat together, steps as follow:
   - Go to mapping directory
   - Run command:
     find . "*_extracted.bed" -type f -exec cat {} \; | gzip -9c > .../YT_extracted_all.txt.gz

6. Further filtering on final output
   - Run command:
     zcat YT_extracted_all.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"} $9>0 && $9 < 10000 {print $0}' | gzip -9c > \
     YT_extracted_less_than_10k.txt.gz
   - Run command:
     zcat YT_extracted_less_than_10k.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"}{if($8=="+"){print $6, $7, $8, $10, $11, $9} \
     else {print $5+1, $7, $8, $10, $11, $9}}' | awk 'BEGIN{FS="\t";OFS="\t"}$5~"protein_coding"{print $0}' | gzip -9c > \
     distance_to_gene_10k_protein_coding.txt.gz


     zcat YT_extracted_all.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"}{print $7, $10, $11, $9}' | gzip -9c > \
     distance_to_gene.txt.gz
     zcat YT_extracted_less_than_20k.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"}{print $7, $10, $11, $9}' | gzip -9c > \
     distance_to_gene_20k.txt.gz
     zcat distance_to_gene_20k.txt.gz | awk 'BEGIN{FS="\t";OFS="\t"}$3=="protein_coding"{print $0}' | gzip -9c > \
     distance_to_gene_20k_protein_coding.txt.gz



7. # Prepare adjacent_profiles.txt  # No longer required.

   Go to R, run extension_length.R, acquire extension_profiles.txt

8. Run extend_gtf_gffutils.py  # extend gtf files based on extension_profiles.txt

************************************************************************************************************************
Extract ambiguity reads

1. Run mapping_pipeline_parallel_new.py:
   do_parallel("unassigned_ambiguity_summary")
   do_parallel("ambiguity_bam_to_bed")

2. Go to mapping directory:
   samtools merge all_ambiguity.bam `find mapping/ -name "*ambiguity.bam"` -@ 32
   This merges all ambiguity reads into one BAM file for easier visualization in the future.

3. Go to mapping directory:
   find mapping/ -name "*ambiguity.bed" -type f -exec cat {} \; | bedtools sort -i - | uniq - > all_ambiguity_uniq_bed

4. Process all_ambiguity_uniq.bed:
   bedtools cluster -d 150 -i all_ambiguity_uniq.bed > all_ambiguity_uniq_clusters.bed

   uniq -f1 -f2 -f3 -c all_ambiguity_uniq_clusters.bed | sed 's/^\s*//' | tr ' ' '\t' >
   all_ambiguity_uniq_clusters_stats.bed

   awk 'BEGIN{FS="\t";OFS="\t"} $1>60 {print $0}' all_ambiguity_uniq_clusters_stats.bed > all_ambiguity_final.bed


************************************************************************************************************************
What's inside ...:
    batch_process.py:
        - samtools index
        - md5-checksum
