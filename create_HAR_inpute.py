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


# Process genotype file for a single chromosome
def process_genotype_file_for_chr(chromosome, geno_file, bim, hars, output_file, chunksize=500000):
    """Process genotype data efficiently for a specific chromosome."""
    print(f"Processing chromosome {chromosome} from {geno_file}...")
    total_hars = len(hars)

    with open(output_file, 'w') as out_f:
        header_written = False

        for chunk in pd.read_csv(geno_file, sep=',', chunksize=chunksize):
            chunk = chunk[chunk['CHR'] == chromosome].set_index('LOC', drop=False)  # Filter by chromosome

            for har_idx, (_, har) in enumerate(hars.iterrows(), start=1):
                if har['chrom'] != chromosome:
                    continue

                if har_idx % 100 == 0:
                    print(f"Processed {har_idx}/{total_hars} HARs on chromosome {chromosome}...")

                start, end, har_id = har['start'], har['end'], har['har_id']
                snps_in_region = bim[(bim['CHR'] == chromosome) & (bim['LOC'].between(start, end))]

                while len(snps_in_region) < 10:
                    start -= 500  # Expand left
                    end += 500  # Expand right
                    snps_in_region = bim[(bim['CHR'] == chromosome) & (bim['LOC'].between(start, end))]

                if snps_in_region.empty:
                    continue

                snp_locs = list(snps_in_region['LOC'])
                chunk_filtered = chunk.loc[chunk.index.intersection(snp_locs)]
                if chunk_filtered.empty:
                    continue

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

    print(f"Processed chromosome {chromosome} and saved output to {output_file}")


# Run processing in parallel for all chromosomes
def process_all_chromosomes(geno_file, bim, hars, output_dir, num_processes=22):
    """Process all chromosomes in parallel."""
    os.makedirs(output_dir, exist_ok=True)
    chromosomes = hars['chrom'].unique()
    tasks = [(chr_num, geno_file, bim, hars, os.path.join(output_dir, f'chr{chr_num}_processed.csv')) for chr_num in
             chromosomes]

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_genotype_file_for_chr, tasks)

    print("Finished processing all chromosomes.")


# Main workflow
def main():
    parser = argparse.ArgumentParser(description="Process genotype data per chromosome in parallel.")
    parser.add_argument('--bim', required=True, help="Input BIM file (for SNP metadata)")
    parser.add_argument('--hars', required=True, help="Input HAR regions file")
    parser.add_argument('--genotype', required=True, help="Input genotype file (CSV format)")
    parser.add_argument('--output_dir', required=True, help="Output directory for processed files")

    args = parser.parse_args()
    bim, hars = load_har_bim_data(args.bim, args.hars)
    process_all_chromosomes(args.genotype, bim, hars, args.output_dir)


if __name__ == "__main__":
    main()
