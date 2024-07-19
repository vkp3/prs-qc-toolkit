#!/bin/sh
# Removes multi-allelic variants from VCF
# - Useful for conversion to Plink BED
#   for input into PRSice

# vcf_ 1
vcftools --gzvcf vcf_1.vcf.gz \
--max-alleles 2 \
--keep-INFO-all \
--recode \
--out vcf_1.biallelic
# gzip vcf_1.biallelic.recode.vcf

# vcf_ 2
vcftools --gzvcf vcf_2.vcf.gz \
--max-alleles 2 \
--keep-INFO-all \
--recode \
--out vcf_2.biallelic
# gzip vcf_2.biallelic.recode.vcf
