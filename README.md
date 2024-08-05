# PRS QC Toolkit

This toolkit provides quality control checks for Polygenic Risk Scores (PRS).

It offers a suite of metrics and analyses to assess the reliability and applicability of PRS models to specific genetic datasets. By evaluating factors such as genetic coverage, allele frequency distributions, and genotyping quality, it helps researchers and clinicians gauge the confidence level in PRS results for genotyped samples. The toolkit aims to enhance the accuracy and interpretability of PRS in both research and clinical settings.

Author: Vamsee Pillalamarri (vpillal1@alumni.jh.edu)

## Features
- VCF file validation
- Confidence score calculation using metrics related to PRS models applied to specific genotyping datasets
- Current implemented metrics:
  - Coverage calculation
  - EAF-weighted coverage
  - EAF-ratio
  - Alt VAF and MQ quality assessment
  - LOW QUAL rate
  - Missing genotype rate

**Please see below [Notes](#Notes) section for detailed description of current implemented metrics and notes on future metrics.**

## Installation

You can install the PRS QC Toolkit directly from GitHub:

```bash
pip install git+https://github.com/vkp3/prs-qc-toolkit.git
```

For development, clone the repository and install in editable mode:

```bash
git clone https://github.com/vkp3/prs-qc-toolkit.git
cd prs-qc-toolkit
pip install -e .
```

## Usage

To run the main script:

```bash
python src/prs_qc.py
```

Make sure your directory structure matches the expected layout:

```
prs-qc-toolkit/
├── src/
│   └── prs_qc.py
└── dat/
    ├── geno/
    │   ├── vcf_1.biallelic.recode.vcf // Input VCFs with biallelic records only.
    │   └── vcf_2.biallelic.recode.vcf
    ├── grs_models/
    │   ├── prs1.tsv
    │   └── prs2.tsv
    └── effect_alleles.afreq
```

## Example Input
- VCF files must be gunzipped (Future change will allow gzip) and placed in `./dat/geno`
- PRS models (i.e., discovery GWAS summary statistics) in `./dat/grs_models/` should follow the same column format as example models found in `dat/grs_models`

## Example Output

```
Metrics for vcf_1.biallelic.recode.vcf and prs1.tsv:
Coverage        EAF Ratio       Model EAF       Covered EAF     Missing Rate    LowQual Prop    Alt_VAF_and_MQ
0.9980          1.0005          0.1980          0.1981          0.0000          0.0026          0.9607

Metrics for vcf_1.biallelic.recode.vcf and prs2.tsv:
Coverage        EAF Ratio       Model EAF       Covered EAF     Missing Rate    LowQual Prop    Alt_VAF_and_MQ
0.9968          1.0005          0.1985          0.1986          0.0000          0.0026          0.9607

Metrics for vcf_2.biallelic.recode.vcf and prs1.tsv:
Coverage        EAF Ratio       Model EAF       Covered EAF     Missing Rate    LowQual Prop    Alt_VAF_and_MQ
0.9979          1.0003          0.1980          0.1981          0.0000          0.0010          0.6780

Metrics for vcf_2.biallelic.recode.vcf and prs2.tsv:
Coverage        EAF Ratio       Model EAF       Covered EAF     Missing Rate    LowQual Prop    Alt_VAF_and_MQ
0.9969          1.0004          0.1985          0.1986          0.0000          0.0010          0.6780

  Sample Model  Confidence
0  vcf_1  prs1    0.987165
1  vcf_1  prs2    0.986571
2  vcf_2  prs1    0.933028
3  vcf_2  prs2    0.932502

Scores that should be masked (below confidence threshold):
  Sample Model  Confidence
2  vcf_2  prs1    0.933028
3  vcf_2  prs2    0.932502
```

## Example Usage in Python

This toolkit can be used in your Python scripts. Here's a basic example:

```python
from prs_qc import evaluate_confidence

vcf_file = 'path/to/your/vcf_file.vcf'
model_file = 'path/to/your/model_file.tsv'
eaf_file = 'path/to/your/eaf_file.afreq'

confidence_score = evaluate_confidence(vcf_file, model_file, eaf_file)
print(f"Confidence Score: {confidence_score}")
```

To run the example script:

```bash
python examples/example_usage.py
```

## Contributing

Please feel free to submit a Pull Request!


## Notes
To provide QC for applying a specific PRS to a sample VCF, various heuristic and informed apporaches can be taken. The following describes the thought process.
## Heuristic Approaches

### Implemented QC Metrics

1. **Coverage**:
   - Calculates the proportion of variants in the GRS model that are present in the VCF data.
   - Transformed using a squared function (min(coverage^2, 1)) to emphasize high coverage.

2. **EAF-weighted Coverage**:
   - Calculates the proportion of total "frequency" of all variants in the GRS Model that is covered by the genotyped variants.
   - Potential change would be to use the square root of EAF to balance the influence of rare and common variants.

3. **EAF Ratio**:
   - Compares the mean EAF of covered variants to the mean EAF of all model variants.
   - Transformed using 1 / (1 + abs(log(eaf_ratio))) to penalize deviations from an ideal ratio of 1.

4. **Missing Genotype Rate**:
   - Calculates the proportion of variants with missing genotypes (including partial missingness on one or both alleles).
   - Transformed using an exponential decay function (exp(-5 * missing_rate)) to heavily penalize high missing rates.

5. **Low Quality Variant Proportion**:
   - Calculates the proportion of variants marked as "LowQual" in the FILTER column of the VCF.
   - Transformed using a steeper exponential decay function (exp(-10 * low_qual_prop)) to strongly penalize low quality data.

6. **ALT VAF and MQ**
   - Compute the effect variant allele frequency using read depth info in VCF
   - Use mapping quality (MQ) output by mapping method (e.g. BWA-Mem) to filter poor-quality genotyping calls.

### Future Metrics where Additional Data Needed

To refine the decision-making process, the following additional data would be beneficial:

1. **Percentile-based Risk Estimate**:
   - Requires a population dataset for comparison.

2. **Relative Risk**:
   - Requires prevalence of disease among deciles of PRS compared to general population prevalence.
   - Utilize UK Biobank PRSs and general population prevalence.

3. **Absolute Risk**:
   - Requires incidence of disease from a population similar to the sample's ancestry.
   - Multiply by the relative risk of the score compared to a PRS distribution.

**Citation:** 1. Lewis, A.C.F., Green, R.C., and Vassy, J.L. (2021). Polygenic risk scores in the clinic: Translating risk into action. Human Genetics and Genomics Advances 2, 100047. https://doi.org/10.1016/j.xhgg.2021.100047.


For this, the following would be needed:

1. More samples in a population for comparison:
   - For percentile ranking of PRS.
   - For measures of relative and absolute risk.

2. Ancestry components to ensure matching between discovery GRS models and samples.

3. Calculated risks for other traits (even if not being reported) to inform masking decisions.

## Potential Additional Metrics

With access to more data, the following metrics could be implemented to further refine the confidence heuristic:

1. **Hardy-Weinberg Equilibrium (HWE) Deviation**:
   - Measure the proportion of variants that significantly deviate from HWE.
   - Requires population-level genotype data from

2. **Ancestry-Specific MAF Comparison**:
   - Compare sample MAFs to ancestry-matched reference populations.
   - Requires detailed ancestry information for the sample and reference MAF data for various populations.

3. **Linkage Disequilibrium (LD) Coverage**:
   - Assess whether key LD blocks in the GRS model are well-represented in the sample data.
   - Requires LD information for the variants in the GRS model.

4. **Imputation Quality**:
   - Add imputation quality scores into the confidence metric if variants are imputed.

5. **Technical Batch Effects**:
   - Assess and account for any batch effects in genotyping or sequencing.
   - Requires technical metadata about sample processing and ideally a set of technical control samples.

6. **Phenotype-Informed Metrics**:
   - For validated GRS models, incorporate information about which variants are most predictive of the phenotype.
   - Remove variants with pleiotropic associations to other diseases.
   - Consider variants which are eQTLs in enriched tissues for particular disease, etc.
   - Requires detailed information about variant-phenotype associations (GWAS catalog, etc.)

## Utilizing Family Information to Enhance GRS Confidence

Additional information from family members (probands, siblings, or parents) could significantly increase our confidence in Genetic Risk Scores (GRS).

Here's a plan for incorporating this data:

1. **Family History of Disease**:
   - Collect detailed pedigree for diseases relevant to the GRS.
   - Use this information to adjust the prior probability of high-risk scores.
   - For instance, if there's a strong family history of a disease, we might have more confidence in a high GRS for that condition.
   - If there is a strong family history but low GRS, very interesting to explore why from QC perspective and/or presence of a rare mutation that's not being captured by GRS.

2. **Parental Genotypes**:
   - Compare the sample's genotypes with parental genotypes.
   - Verify Mendelian inheritance patterns to identify potential genotyping errors.
   - Use trio-based phasing to improve haplotype estimation, maybe improve accuracy of the GRS?

3. **Sibling Comparisons**:
   - Compare GRS and underlying genotypes across siblings.
   - Unusually large differences might indicate potential issues with one of the scores.
   - Consistent patterns across siblings could increase confidence in the scores.
   - Take advanage of common allelic background.

4. **Rare Variant Analysis**:
   - Identify rare variants shared among family members from WGS.
   - If these rare variants are known to be associated with the disease of interest, it could increase confidence in seeing a high GRS.

5. **Polygenic Transmission Disequilibrium (pTDT)**:
   - Analyze how risk alleles are transmitted from parents to offspring.
   - Consistent over-transmission of risk alleles in affected families could increase confidence in the GRS?

6. **Family-Based Imputation**:
   - Use family data to improve imputation of missing genotypes.
   - This could increase the number of variants included in the GRS calculation.
