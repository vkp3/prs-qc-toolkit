#!/bin/sh
# Generates polygenic scores (grs) for PRS1, PRS2 x 2 vcfs
plink2=/Users/vkp/bin/bin/plink2
dat=dat/
out=grs/
log=logs/

# prs1 ----
## vcf_1
${plink2} \
  --pfile ${dat}/geno/vcf_1 \
  --score ${dat}/grs_models/prs1.tsv 1 2 3 header list-variants \
  --read-freq ${dat}/effect_alleles.afreq \
  --out ${out}/prs1.vcf_1

## vcf_2
${plink2} \
  --pfile ${dat}/geno/vcf_2 \
  --score ${dat}/grs_models/prs1.tsv 1 2 3 header list-variants \
  --read-freq ${dat}/effect_alleles.afreq \
  --out ${out}/prs1.vcf_2

# prs2 ----
## vcf_1
${plink2} \
  --pfile ${dat}/geno/vcf_1 \
  --score ${dat}/grs_models/prs2.tsv 1 2 3 header list-variants \
  --read-freq ${dat}/effect_alleles.afreq \
  --out ${out}/prs2.vcf_1

## vcf_2
${plink2} \
  --pfile ${dat}/geno/vcf_2 \
  --score ${dat}/grs_models/prs2.tsv 1 2 3 header list-variants \
  --read-freq ${dat}/effect_alleles.afreq \
  --out ${out}/prs2.vcf_2

mv ${out}/*.log ${log}
