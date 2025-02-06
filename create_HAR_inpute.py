import pandas as pd
import argparse
import os
import multiprocessing


# Load BIM file and HAR regions
def load_har_bim_data(bim_file, hars_file):
    """Load SNP metadata from BIM file and HAR region annotations."""
    bim = pd.read_csv(bim_file, sep='\t', header=None, names=['CHR', 'SNP_ID', 'CM', 'LOC', 'A1', 'A2'])
    hars = pd.read_csv(hars_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'har_id'])
    return bim, hars


# Process genotype file directly without splitting into chromosome-specific files
def process_genotype_file(geno_file, bim, hars, chunksize=500000):
    """Process genotype data efficiently from a single genotype CSV file."""
    output_file = "processed_har_genotypes.csv"  # Fixed output file name
    print(f"Processing genotype file {geno_file}...")
    total_hars = len(hars)

    with open(output_file, 'w') as out_f:
        header_written = False

        for chunk in pd.read_csv(geno_file, sep=',', chunksize=chunksize):
            chunk = chunk.set_index('LOC', drop=False)  # Optimize lookup by indexing LOC

            for har_idx, (_, har) in enumerate(hars.iterrows(), start=1):
                if har_idx % 100 == 0:
                    print(f"Processed {har_idx}/{total_hars} HARs...")

                start, end, har_id = har['start'], har['end'], har['har_id']
                snps_in_region = bim[(bim['LOC'].between(start, end))]

                while len(snps_in_region) < 10:
                    start -= 500  # Expand left
                    end += 500  # Expand right
                    snps_in_region = bim[(bim['LOC'].between(start, end))]

                if snps_in_region.empty:
                    continue

                snp_locs = list(snps_in_region['LOC'])
                snp_ids = list(snps_in_region['SNP_ID'])
                # print(f"HAR: {har_id}, SNPs: {snp_ids}")  # Debugging print

                chunk_filtered = chunk.loc[chunk.index.intersection(snp_locs)]
                if chunk_filtered.empty:
                    continue
                print(chunk_filtered)

                snp_numbering = {loc: i + 1 for i, loc in enumerate(snp_locs)}
                chunk_filtered.insert(0, 'HAR_SNP', har_id + '_' + chunk_filtered['LOC'].map(snp_numbering).astype(str))
                original_individual_columns = [col for col in chunk_filtered.columns if col.startswith("GTEX")]
                columns_to_keep = ['HAR_SNP', 'CHR', 'LOC'] + original_individual_columns
                chunk_filtered = chunk_filtered[columns_to_keep]

                if not header_written:
                    chunk_filtered.to_csv(out_f, index=False, header=True, mode='w')
                    header_written = True
                else:
                    chunk_filtered.to_csv(out_f, index=False, header=False, mode='a')

    print(f"Processed genotype file and saved output to {output_file}")


# Main workflow
def main():
    parser = argparse.ArgumentParser(description="Process genotype data from a single file.")
    parser.add_argument('--bim', required=True, help="Input BIM file (for SNP metadata)")
    parser.add_argument('--hars', required=True, help="Input HAR regions file")
    parser.add_argument('--genotype', required=True, help="Input genotype file (CSV format)")

    args = parser.parse_args()
    bim, hars = load_har_bim_data(args.bim, args.hars)

    process_genotype_file(args.genotype, bim, hars)


if __name__ == "__main__":
    main()
