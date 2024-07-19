#!/bin/sh

# vcf_1
~/bin/bin/plink2 \
--vcf vcf_1.biallelic.recode.vcf.gz \
--make-bed \
--out vcf_1.biallelic

# vcf_2
~/bin/bin/plink2 \
--vcf vcf_2.biallelic.recode.vcf.gz \
--make-bed \
--out vcf_2.biallelic
