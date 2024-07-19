#!/bin/sh
prsice=/Users/vkp/bin/bin/PRSice/PRSice_mac
dat=dat/
out=PRS/prsice/

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bed_file_prefix> <PRS_filename>. Prefix/filenames map to dat/geno/<prefix> and dat/PRS_models/<filename>"
    exit 1
fi

BED_PREFIX=$1
PRS_FILE=$2
PRS_BASE=$(basename "${PRS_FILE}" .tsv)
OUT_PREFIX="${out}/${BED_PREFIX}_${PRS_BASE}"

${prsice} \
--base ${dat}/PRS_models/${PRS_FILE} \
--target ${dat}/geno/${BED_PREFIX} \
--out ${OUT_PREFIX} \
--stat effect_weight \
--beta \
--A1 effect_allele \
--A2 REF \
--chr chr \
--bp position_hg19 \
--snp variant \
--pvalue P \
--no-clump \
--lower 1 \
--upper 1 \
--interval 1 \
--no-regress
