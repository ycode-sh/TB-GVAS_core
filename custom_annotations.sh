#!/bin/bash

# This script performs a sequential custom annotation on vcf files using annotion files. 
# It annotates the variants in a vcf file that intersects with lineage, amr_region, gene_annotation, and variable_regions information
# $1= gzipped and tabix-indexed gene_annotation_file, $2 = gene_annotation_header, $3 = lineage_annotation_file, $4 = lineage_annotation_header, 
# $5 = amr_region_file, $6 = amr_region_header, $7 = variable_region_annotation_file, $8 = variable_region_annotation_header $9 = vcf_file

bgzip $1

bgzip $3

bgzip $5

bgzip $7


file_name=$(basename $9 .vcf)
bcftools annotate --annotations  ${3} --header-lines $4  --columns CHROM,from,to,lineage $9 > ${file_name}_la.vcf
bcftools annotate --annotations  ${5} --header-lines $6  --columns CHROM,from,to,related_antibiotics,antibiotics_gene ${file_name}_la.vcf > ${file_name}_ar.vcf
bcftools annotate --annotations  ${7} --header-lines $8  --columns CHROM,from,to,locus_tag,comment ${file_name}_ar.vcf > ${file_name}_vr.vcf
bcftools annotate --annotations  ${1} --header-lines $2  --columns CHROM,from,to,gene_ann_region,gene_annotation ${file_name}_vr.vcf > ${file_name}_c_ann.vcf


rm ${file_name}_la.vcf ${file_name}_ar.vcf ${file_name}_vr.vcf 
#bcftools annotate --annotations $1 --header-lines $2  --columns CHROM,from,to,gene_ann_region,gene_annotation $3 > ${file_name}_ga.vcf