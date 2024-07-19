#!/bin/sh

# vcf_1
vcftools --gzvcf vcf_1.vcf.gz \
--min-alleles 3 \
--max-alleles 3 \
--recode \
--out multi_allelic_variants_vcf_1

cat multi_allelic_variants_vcf_1.recode.vcf | grep -v '#' | cut -f3 > multi_allelic_variants_vcf_1.recode.rsIDs

# vcf_ 2
vcftools --gzvcf vcf_2.vcf.gz \
--min-alleles 3 \
--max-alleles 3 \
--recode \
--out multi_allelic_variants_vcf_2

cat multi_allelic_variants_vcf_2.recode.vcf | grep -v '#' | cut -f3 > multi_allelic_variants_vcf_2.recode.rsIDs
