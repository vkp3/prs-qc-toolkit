from prs_qc import evaluate_confidence

confidence_score = evaluate_confidence('dat/geno/vcf_1.biallelic.recode.vcf', 'dat/grs_models/prs1.tsv', 'dat/effect_alleles.afreq')
print(f"Confidence Score: {confidence_score}")
