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


# Split genotype file by chromosome
def split_genotype_by_chr(geno_file, output_dir, chunksize=500000):
    """Split genotype data into separate files for each chromosome."""
    os.makedirs(output_dir, exist_ok=True)
    print("Splitting genotype data by chromosome...")

    for chunk in pd.read_csv(geno_file, sep=',', chunksize=chunksize):
        for chrom in chunk['CHR'].unique():
            chr_chunk = chunk[chunk['CHR'] == chrom]
            chr_file = os.path.join(output_dir, f'chr{chrom}_geno.csv')
            mode = 'w' if not os.path.exists(chr_file) else 'a'
            header = not os.path.exists(chr_file)
            chr_chunk.to_csv(chr_file, index=False, mode=mode, header=header)
            print(f"Written chromosome {chrom} data to {chr_file}")


# Process genotype file efficiently and filter relevant SNPs per chromosome
def process_genotype_file_per_chr(chromosome, geno_file, bim, hars, output_file, chunksize=500000):
    """Process genotype data efficiently for a specific chromosome and extract relevant SNPs for HARs."""
    total_hars = len(hars)
    print(f"Processing chromosome {chromosome} with {total_hars} HAR regions...")

    with open(output_file, 'w') as out_f:
        header_written = False

        for chunk in pd.read_csv(geno_file, sep=',', chunksize=chunksize):
            chunk = chunk.set_index('LOC')  # Optimize lookup by indexing LOC

            for har_idx, (_, har) in enumerate(hars.iterrows(), start=1):
                if har_idx % 100 == 0:
                    print(f"Processed {har_idx}/{total_hars} HARs for chromosome {chromosome}...")

                if har['chrom'] != chromosome:
                    continue

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

                # Assign sequential numbers to SNPs within each HAR set
                snp_numbering = {loc: i + 1 for i, loc in enumerate(snp_locs)}
                chunk_filtered.insert(0, 'HAR_SNP', har_id + '_' + chunk_filtered.index.map(snp_numbering).astype(str))

                # Keep LOC and all original GTEX individual genotype columns
                original_individual_columns = [col for col in chunk_filtered.columns if col.startswith("GTEX")]
                columns_to_keep = ['HAR_SNP', 'LOC'] + original_individual_columns
                chunk_filtered = chunk_filtered.reset_index()[columns_to_keep]

                # Write header only once
                if not header_written:
                    chunk_filtered.to_csv(out_f, index=False, header=True, mode='w')
                    header_written = True
                else:
                    chunk_filtered.to_csv(out_f, index=False, header=False, mode='a')

    print(f"Processed genotype file for chromosome {chromosome} and saved output to {output_file}")


# Merge processed chromosome files into one final output
def merge_chromosome_outputs(output_dir, final_output_file):
    """Merge all chromosome-specific processed files into a final output file."""
    print("Merging chromosome outputs into final output file...")
    all_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if
                 f.startswith("chr") and f.endswith("_processed.csv")]
    merged_df = pd.concat([pd.read_csv(f) for f in all_files], ignore_index=True)
    merged_df.to_csv(final_output_file, index=False)
    print(f"Final merged file saved to {final_output_file}")


# Main workflow
def main():
    parser = argparse.ArgumentParser(description="Process genotype, HAR, and SNP data in parallel by chromosome.")
    parser.add_argument('--bim', required=True, help="Input BIM file (for SNP metadata)")
    parser.add_argument('--hars', required=True, help="Input HAR regions file")
    parser.add_argument('--genotype', required=True, help="Input genotype file (large CSV)")
    parser.add_argument('--output_dir', required=True, help="Output directory for chromosome-wise processing")
    parser.add_argument('--final_output', required=True, help="Final merged output file")

    args = parser.parse_args()

    print("Loadindg HAR and SNP metadata...")
    bim, hars = load_har_bim_data(args.bim, args.hars)

    print("Splitting genotype file by chromosome...")
    split_genotype_by_chr(args.genotype, args.output_dir)

    # Process each chromosome in parallel
    pool = multiprocessing.Pool(processes=22)  # Run in parallel for 22 chromosomes
    tasks = []
    for chrom in range(1, 23):
        chr_geno_file = os.path.join(args.output_dir, f'chr{chrom}_geno.csv')
        chr_output_file = os.path.join(args.output_dir, f'chr{chrom}_processed.csv')
        if os.path.exists(chr_geno_file):
            tasks.append((chrom, chr_geno_file, bim, hars, chr_output_file))

    pool.starmap(process_genotype_file_per_chr, tasks)
    pool.close()
    pool.join()

    print("Merging chromosome outputs...")
    merge_chromosome_outputs(args.output_dir, args.final_output)


if __name__ == "__main__":
    main()
