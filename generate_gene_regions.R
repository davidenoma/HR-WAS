library(rtracklayer)

# Select gene regions from GFF file

# Load the GFF file
gff_file <- "/Users/davidenoma/Desktop/PhD._BMB/LONG_LAB/Projects/EDAS/EDAS/AFR_SKAT/Homo_sapiens.GRCh38.110.gff3"
gff <- readGFF(gff_file)

# Filter for gene features
gene_features <- subset(gff, type == "gene" & seqid %in% 1:22)
print(nrow(gene_features))
# Assuming gene_features is your data frame
gene_features <- gene_features[complete.cases(gene_features$Name), ]

print(nrow(gene_features))
# Assuming gene_features is your data frame



#   Obtain start and stop positions
flank_size <- 10000  # 500kb
#   add 500kb to the start and stop positions respectively
# Create a data frame to store gene regions
gene_regions <- data.frame(GeneName = character(0), Start = numeric(0), Stop = numeric(0))

# flank_size <- 500000  # 500kbxC

gene_regions <- data.frame(
  GeneName = character(),
  Start = numeric(),
  Stop = numeric(),
  Chromosome = character(),
  stringsAsFactors = FALSE
)
gene_name = gene_features$Name
gene_start = gene_features$start
gene_end = gene_features$end
gene_chromosome = gene_features$seqid
# Assuming gene_start is a vector of start positions and flank_size is a scalar
region_start <- as.integer(pmax(1, gene_start - flank_size))

region_end <- integer(length(gene_end))  # Initialize region_end as an integer vector

for (i in seq_along(gene_end)) {
  region_end[i] <- gene_end[i] + flank_size
}
gene_regions <- rbind(gene_regions, data.frame(GeneName = gene_name, Start = region_start, Stop = region_end,Chromosome = gene_chromosome))
# Save gene regions to a CSV file
write.csv(gene_regions, "gene_regions_500kb.csv", row.names = FALSE)

