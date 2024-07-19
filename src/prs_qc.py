import pandas as pd
import numpy as np
from pysam import VariantFile
import os
import sys

def load_vcf(vcf_file):
    """
    Load and parse a VCF file.
    Currently assumes VCF has only one sample.

    Args:
    vcf_file (str): Path to the Variant Call Format (VCF) file.

    Returns:
    pd.DataFrame: Parsed VCF data containing variant information.

    Raises:
    FileNotFoundError: If the specified VCF file is not found.
    """
    try:
        vcf = VariantFile(vcf_file)
    except FileNotFoundError:
        print(f"Error: VCF file not found: {vcf_file}")
        sys.exit(1)
    
    variant_data = []
    for record in vcf.fetch():
        # Extract relevant information from each VCF record
        alt = record.alts[0] if record.alts else None
        gt = record.samples[0]['GT']
        ad = record.samples[0].get('AD', (None, None))
        dp = record.samples[0].get('DP', None)
        
        variant_data.append({
            'CHROM': record.chrom,
            'POS': record.pos,
            'variant': record.id,
            'REF': record.ref,
            'ALT': alt,
            'QUAL': record.qual,
            'FILTER': ','.join(record.filter.keys()),
            'AC': record.info.get('AC', [None])[0],
            'MQ': record.info.get('MQ', None),
            'GT': '/'.join(map(str, gt)) if gt else None,
            'REF_DP': ad[0] if ad[0] is not None else None,
            'ALT_DP': ad[1] if ad[1] is not None else None,
            'DP': dp
        })
    
    df = pd.DataFrame(variant_data)
    # Calculate Variant Allele Frequencies (VAF) of REF and ALT
    df['REF_VAF'] = df['REF_DP'] / df['DP']
    df['ALT_VAF'] = df['ALT_DP'] / df['DP']
    return df

def load_model(model_file, eaf_file):
    """
    Load and merge the GRS model and Effect Allele Frequency (EAF) data.

    Args:
    model_file (str): Path to the GRS model file.
    eaf_file (str): Path to the Effect Allele Frequency file.

    Returns:
    pd.DataFrame: Merged data containing model and EAF information.

    Raises:
    FileNotFoundError: If either the model file or EAF file is not found.
    """
    try:
        model = pd.read_csv(model_file, sep='\t')
        eaf = pd.read_csv(eaf_file, sep='\t')
    except FileNotFoundError as e:
        print(f"Error: File not found: {e.filename}")
        sys.exit(1)
    
    eaf = eaf.rename(columns={'ALT_FREQS': 'ALT_FREQ'})
    model = pd.merge(model, eaf[['ID', 'REF', 'ALT', 'ALT_FREQ']], 
                     left_on=['variant', 'REF', 'ALT'], 
                     right_on=['ID', 'REF', 'ALT'])
    return model

def calculate_coverage(vcf_data, model_data):
    """
    Calculate the coverage of model variants in the VCF data.

    Args:
    vcf_data (pd.DataFrame): Parsed VCF data.
    model_data (pd.DataFrame): GRS model data.

    Returns:
    float: Proportion of model variants present in the VCF data.
    """
    present_snps = sum(model_data['variant'].isin(vcf_data['variant']))
    total_snps = len(model_data)
    return present_snps / total_snps

def calculate_alt_vaf_mq_quality(vcf_data, min_depth=10, min_alt_vaf_het=0.25, min_alt_vaf_hom=0.75, min_mq=20):
    """
    Calculate the proportion of variants meeting quality criteria based on 
    alternate allele VAF and mapping quality.

    Args:
    vcf_data (pd.DataFrame): Parsed VCF data.
    min_depth (int): Minimum read depth required.
    min_alt_vaf_het (float): Minimum ALT VAF for heterozygous calls. See plots under dat/
    min_alt_vaf_hom (float): Minimum ALT VAF for homozygous calls.
    min_mq (int): Minimum mapping quality required.

    Returns:
    float: Proportion of variants meeting the quality criteria.
    """
    good_quality_variants = vcf_data[
        (vcf_data['DP'] >= min_depth) &
        (vcf_data['MQ'] >= min_mq) &
        (
            ((vcf_data['GT'] == '0/1') & (vcf_data['ALT_VAF'] >= min_alt_vaf_het) & (vcf_data['ALT_VAF'] <= (1 - min_alt_vaf_het))) |
            ((vcf_data['GT'] == '1/1') & (vcf_data['ALT_VAF'] >= min_alt_vaf_hom))
        )
    ]
    return len(good_quality_variants) / len(vcf_data)

def calculate_eaf_coverage(vcf_data, model_data, af_data):
    """
    Calculate EAF-weighted coverage metrics.
    - Simple coverage proportion (number of model variants present in VCF / total number of model variants)
    - Effect Allele Frequency (EAF) for all variants in the model and for the covered variants.
    - EAF ratio (mean EAF of covered variants / mean EAF of all model variants)

    Args:
    vcf_data (pd.DataFrame): Parsed VCF data.
    model_data (pd.DataFrame): GRS model data.
    af_data (pd.DataFrame): Allele frequency data.

    Returns:
    dict: Dictionary containing various EAF coverage metrics.
    """
    model_variants = set(model_data['variant'])
    vcf_variants = set(vcf_data['variant'])
    covered_variants = model_variants.intersection(vcf_variants)
    coverage_prop = len(covered_variants) / len(model_variants)
    
    model_mean_eaf = af_data[af_data['ID'].isin(model_variants)]['ALT_FREQS'].mean()
    covered_mean_eaf = af_data[af_data['ID'].isin(covered_variants)]['ALT_FREQS'].mean()
    eaf_ratio = covered_mean_eaf / model_mean_eaf
    
    return {
        'coverage_prop': coverage_prop,
        'model_mean_eaf': model_mean_eaf,
        'covered_mean_eaf': covered_mean_eaf,
        'eaf_ratio': eaf_ratio,
        'covered_count': len(covered_variants),
        'model_count': len(model_variants)
    }

def evaluate_confidence(vcf_file, model_file, eaf_file):
    """
    - Evaluate the confidence of applying a GRS model to a VCF file.
    - Applies transformations to metrics and computes a confidence score.
    - The EAF ratio is transformed (using 1 / (1 + abs(log(eaf_ratio))))
    to create a score that's sensitive to deviations from the ideal ratio of 1.
    - Missing rate and 'LOW QUAL' rates - exponential decay penalizes high missing rates

    Args:
    vcf_file (str): Path to the VCF file.
    model_file (str): Path to the GRS model file.
    eaf_file (str): Path to the Effect Allele Frequency file.

    Returns:
    float: Confidence score for applying the GRS model to the VCF data.
    """
    vcf_data = load_vcf(vcf_file)
    model_data = load_model(model_file, eaf_file)
    af_data = pd.read_csv(eaf_file, sep='\t')

    # Metrics & transformations
    coverage = calculate_coverage(vcf_data, model_data)
    eaf_coverage_metrics = calculate_eaf_coverage(vcf_data, model_data, af_data)
    eaf_coverage = eaf_coverage_metrics['coverage_prop']
    eaf_ratio = eaf_coverage_metrics['eaf_ratio']
    eaf_ratio_score = 1 / (1 + abs(np.log(eaf_ratio)))

    missing_rate = vcf_data['GT'].str.contains(r'\.', na=False).mean()
    low_qual_prop = (vcf_data['FILTER'] == 'LowQual').mean()
    alt_vaf_mq_quality = calculate_alt_vaf_mq_quality(vcf_data)
    
    # Calculate confidence score using weighted sum of metrics
    confidence_score = (
        min(coverage**2, 1) * 0.25 +
        eaf_ratio_score * 0.25 +
        alt_vaf_mq_quality * 0.2 +
        np.exp(-15 * missing_rate) * 0.15 +
        np.exp(-10 * low_qual_prop) * 0.15
    )
    
    # Print a table of detailed metrics per PRS x VCF combination
    print(f"\nMetrics for {os.path.basename(vcf_file)} and {os.path.basename(model_file)}:")
    print(f"{'Coverage':<15} {'EAF Ratio':<15} {'Model EAF':<15} {'Covered EAF':<15} {'Missing Rate':<15} {'LowQual Prop':<15} {'Alt_VAF_and_MQ':<15}")
    print(f"{eaf_coverage_metrics['coverage_prop']:<15.4f} {eaf_coverage_metrics['eaf_ratio']:<15.4f} {eaf_coverage_metrics['model_mean_eaf']:<15.4f} {eaf_coverage_metrics['covered_mean_eaf']:<15.4f} {missing_rate:<15.4f} {low_qual_prop:<15.4f} {alt_vaf_mq_quality:<15.4f}")
    
    return confidence_score

def main():
    """
    Main function to run the PRS confidence evaluation.
    Runs evaluate_confidence() on each PRS x VCF combination,
    which computes confidence in applying a particular PRS model
    on a particular sample using VCF's genotyping.
    """
    samples = ['vcf_1', 'vcf_2']
    models = ['prs1', 'prs2']

    # Using relative paths
    geno_dir = os.path.join('dat', 'geno')
    model_dir = os.path.join('dat', 'grs_models')
    eaf_file = os.path.join('dat', 'effect_alleles.afreq')

#    script_dir = os.path.dirname(os.path.abspath(__file__))
#    geno_dir = os.path.join(script_dir, 'dat', 'geno')
#    model_dir = os.path.join(script_dir, 'dat', 'grs_models')

    # Iterate over samples x models, run compute confidence scores
    # TODO: Could use pandas' MultiIndex for more efficient iteration (like expand.grid() in R)
    results = []
    for sample in samples:
        for model in models:
            vcf_file = os.path.join(geno_dir, f"{sample}.biallelic.recode.vcf")
            model_file = os.path.join(model_dir, f"{model}.tsv")
            # eaf_file = os.path.join(script_dir, 'dat', 'effect_alleles.afreq')

            # Check if all required files exist
            if not all(os.path.exists(f) for f in [vcf_file, model_file, eaf_file]):
                print(f"Error: One or more required files not found for {sample} and {model}")
                continue

            # Compute confidence scores
            confidence = evaluate_confidence(vcf_file, model_file, eaf_file)
            results.append({
                'Sample': sample,
                'Model': model,
                'Confidence': confidence
            })

    # Print confidence scores for all GRS x VCF combinations
    results_df = pd.DataFrame(results).sort_values('Confidence', ascending=False)
    print(results_df)

    # Identify and output masked GRS x VCF combinations
    threshold = 0.95
    masked_res = results_df[results_df['Confidence'] < threshold]
    
    if not masked_res.empty:
        print("\nScores that should be masked (below confidence threshold):")
        print(masked_res)

if __name__ == "__main__":
    main()
