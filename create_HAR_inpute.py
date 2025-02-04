import pandas as pd
import argparse


# Load BIM file, HAR regions, and genotype data efficiently using chunks

def load_har_bim_data(bim_file, hars_file):
    """Load SNP metadata from BIM file and HAR region annotations."""
    bim = pd.read_csv(bim_file, sep='\t', header=None, names=['CHR', 'SNP_ID', 'CM', 'LOC', 'A1', 'A2'])
    hars = pd.read_csv(hars_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'har_id'])
    return bim, hars


# Process genotype file efficiently and filter relevant SNPs

def process_genotype_file(geno_file, bim, hars, output_file, chunksize=100000):
    """Process genotype data efficiently and extract relevant SNPs for HARs."""
    with open(output_file, 'w') as out_f:
        header_written = False

        for chunk in pd.read_csv(geno_file, sep=',', chunksize=chunksize):
            for _, har in hars.iterrows():
                chrom = har['chrom']
                start = har['start']
                end = har['end']
                har_id = har['har_id']

                # Filter SNPs that fall within the HAR region
                snps_in_region = bim[(bim['CHR'] == chrom) & (bim['LOC'].between(start, end))]
                while len(snps_in_region) < 10:
                    start -= 500  # Expand left
                    end += 500  # Expand right
                    snps_in_region = bim[(bim['CHR'] == chrom) & (bim['LOC'].between(start, end))]

                snp_ids = set(snps_in_region['SNP_ID'])
                print(snp_ids)
                chunk_filtered = chunk[chunk['ID'].isin(snp_ids)]

                if chunk_filtered.empty:
                    continue

                # Rename SNP rows with HAR-SNP format
                chunk_filtered.insert(0, 'HAR_SNP', har_id + '_' + chunk_filtered['ID'].astype(str))
                chunk_filtered = chunk_filtered.drop(columns=['ID'])

                # Write header only once
                if not header_written:
                    chunk_filtered.to_csv(out_f, index=False, header=True, mode='w')
                    header_written = True
                else:
                    chunk_filtered.to_csv(out_f, index=False, header=False, mode='a')

    print(f"Processed genotype file and saved output to {output_file}")


# Main workflow
def main():
    parser = argparse.ArgumentParser(description="Process genotype, HAR, and SNP data.")
    parser.add_argument('--bim', required=True, help="Input BIM file (for SNP metadata)")
    parser.add_argument('--hars', required=True, help="Input HAR regions file")
    parser.add_argument('--genotype', required=True, help="Input genotype file (large CSV)")
    parser.add_argument('--output', required=True, help="Output file for processed data")

    args = parser.parse_args()

    print("Loading HAR and SNP metadata...")
    bim, hars = load_har_bim_data(args.bim, args.hars)

    print("Processing genotype file efficiently...")
    process_genotype_file(args.genotype, bim, hars, args.output)


if __name__ == "__main__":
    main()
