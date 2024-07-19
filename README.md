# PRS QC Toolkit

This toolkit provides quality control checks for Polygenic Risk Scores (PRS).

Author: Vamsee Pillalamarri

## Features
- VCF file validation
- Coverage calculation
- EAF-weighted coverage
- Alt VAF and MQ quality assessment
- Confidence score calculation

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
    │   ├── vcf_1.biallelic.recode.vcf
    │   └── vcf_2.biallelic.recode.vcf
    ├── grs_models/
    │   ├── prs1.tsv
    │   └── prs2.tsv
    └── effect_alleles.afreq
```

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
