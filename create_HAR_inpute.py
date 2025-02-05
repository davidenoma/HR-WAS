import pandas as pd
import argparse


# Load BIM file and HAR regions
def load_har_bim_data(bim_file, hars_file):
    """Load SNP metadata from BIM file and HAR region annotations."""
    bim = pd.read_csv(bim_file, sep='\t', header=None, names=['CHR', 'SNP_ID', 'CM', 'LOC', 'A1', 'A2'])
    hars = pd.read_csv(hars_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'har_id'])
    return bim, hars


# Process genotype file efficiently and filter relevant SNPs
def process_genotype_file(geno_file, bim, hars, output_file, chunksize=500000):
    """Process genotype data efficiently and extract relevant SNPs for HARs."""
    total_hars = len(hars)
    print(f"Total HAR regions: {total_hars}")

    with open(output_file, 'w') as out_f:
        header_written = False

        for chunk_idx, chunk in enumerate(pd.read_csv(geno_file, sep=',', chunksize=chunksize)):
            print(f"Processing genotype chunk {chunk_idx + 1}...")
            chunk = chunk.set_index('LOC')  # Optimize lookup by indexing LOC

            for har_idx, (_, har) in enumerate(hars.iterrows(), start=1):
                if har_idx % 100 == 0:
                    print(f"Processed {har_idx}/{total_hars} HARs...")

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

                snp_locs = list(snps_in_region['LOC'])
                chunk_filtered = chunk.loc[chunk.index.intersection(snp_locs)]
                if chunk_filtered.empty:
                    continue
                print(snp_locs)
                                # Assign sequential numbers to SNPs within each HAR set
                snp_numbering = {loc: i + 1 for i, loc in enumerate(snp_locs)}
                chunk_filtered.insert(0, 'HAR_SNP', har_id + '_' + chunk_filtered.index.map(snp_numbering).astype(str))

                # Keep LOC and all original GTEX individual genotype columns
                original_individual_columns = [col for col in chunk_filtered.columns if col.startswith("GTEX")]
                columns_to_keep = ['HAR_SNP', 'LOC'] + original_individual_columns
                chunk_filtered = chunk_filtered.reset_index()[columns_to_keep]
                print(chunk_filtered)

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
