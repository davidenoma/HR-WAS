import pandas as pd
import numpy as np
import argparse


# Load genotype dataset with quality control
def load_genotype_data(input_file):
    """Load genotype data and apply quality control."""
    genotype = pd.read_csv(input_file, sep='\t', comment='#', header=None)
    genotype.columns = ['CHR', 'LOC', 'ID'] + [f'IND-{i}' for i in range(genotype.shape[1] - 3)]

    # Apply quality control filters
    genotype = genotype.dropna()
    genotype = genotype[(genotype.iloc[:, 3:].apply(lambda x: x.isnull().mean(),
                                                    axis=1) < 0.1)]  # Remove variants with >10% missingness

    return genotype


# Generate formatted genotype file
def generate_genotype_file(genotype, output_file):
    """Create a tissue-specific input genotype file with SNP IDs as columns."""
    formatted_data = genotype.copy()
    formatted_data = formatted_data.drop(columns=['ID'])

    # Convert genotype calls to numeric format (assuming 0,1,2 encoding)
    formatted_data.iloc[:, 2:] = formatted_data.iloc[:, 2:].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    # Save the formatted file
    formatted_data.to_csv(output_file, index=False)
    return formatted_data


# Load and process HAR-related datasets
def load_har_data(genotype_file, hars_file):
    """Load genotype and HAR annotations"""
    genotype = pd.read_csv(genotype_file, sep='\t', header=None, names=['CHR', 'ID', 'CM', 'LOC', 'A1', 'A2'])
    hars = pd.read_csv(hars_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'har_id'])  # HARs data
    return genotype, hars


# Expand HARs to ensure at least 10 SNPs
def expand_hars(hars, genotype):
    """Expand HAR regions to left and right until they contain at least 10 SNPs."""
    expanded_hars = []
    for _, row in hars.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        har_id = row['har_id']

        # Ensure at least 10 SNPs in the region
        snps_in_region = genotype[(genotype['CHR'] == chrom) & (genotype['LOC'].between(start, end))]

        while len(snps_in_region) < 10:
            start -= 500  # Expand left
            end += 500  # Expand right
            snps_in_region = genotype[(genotype['CHR'] == chrom) & (genotype['LOC'].between(start, end))]

        expanded_hars.append([chrom, start, end, har_id])

    return pd.DataFrame(expanded_hars, columns=['chrom', 'start', 'end', 'har_id'])


# Generate SNP annotation file
def generate_snp_annotation_file(genotype, output_file):
    """Generate a SNP annotation file with top 50K regulatory variants."""
    # Simulate filtering for top regulatory variants (sorting by CHR and LOC as a placeholder)
    top_snps = genotype.sort_values(by=['CHR', 'LOC']).head(50000)

    snp_annotation_df = top_snps[['ID', 'CHR', 'LOC', 'A1', 'A2']]
    snp_annotation_df.columns = ['SNP', 'varID', 'chr', 'pos', 'ref', 'effect']

    snp_annotation_df.to_csv(output_file, index=False)
    return snp_annotation_df


# Main workflow
def main():
    parser = argparse.ArgumentParser(description="Process genotype and HAR data.")
    parser.add_argument('--genotype', required=True, help="Input genotype file")
    parser.add_argument('--hars', required=True, help="Input HAR regions file")
    parser.add_argument('--output_genotype', required=True, help="Output formatted genotype file")
    parser.add_argument('--output_snp_annotation', required=True, help="Output SNP annotation file")

    args = parser.parse_args()

    print("Loading genotype data...")
    genotype = load_genotype_data(args.genotype)

    print("Generating formatted genotype file...")
    formatted_genotype = generate_genotype_file(genotype, args.output_genotype)
    print(f"Formatted genotype file saved as {args.output_genotype}")

    print("Loading HAR data...")
    genotype, hars = load_har_data(args.genotype, args.hars)

    print("Expanding HARs to ensure a minimum of 10 SNPs")
    expanded_hars = expand_hars(hars, genotype)
    expanded_hars.to_csv('expanded_HARs.bed', sep='\t', index=False, header=False)

    print("Generating SNP annotation file")
    snp_annotation_df = generate_snp_annotation_file(genotype, args.output_snp_annotation)
    print(f"SNP annotation file saved as {args.output_snp_annotation}")


if __name__ == "__main__":
    main()
