library(vroom)
library(ggplot2)

# Path to your SAM file
sam_file <- "/Users/haleymacdonald/Downloads/SA039-A138856-R61-C64_canonical.sam"

# Read SAM file into a data frame
sam_data <- vroom(sam_file, col_names=FALSE, delim = "\t")

# Extract chromosome positions (assuming the positions are in the 4th column, adjust if needed)
positions <- as.numeric(sam_data$X4)

# Extract chromosome names (assuming the names are in the 3rd column, adjust if needed)
chromosomes <- sam_data$X3

calculate_gc <- function(seq) {
  gc_count <- sum(unlist(strsplit(seq, NULL)) %in% c("G", "C"))
  seq_length <- nchar(seq)
  gc_percentage <- gc_count / seq_length
  return(gc_percentage)
}

seq <- sam_data$X10

gc_content <- sapply(seq, calculate_gc)

mapq <- sam_data$X5

# Create a data frame with positions and chromosomes
data <- data.frame(Position = positions, Chromosome = chromosomes, seq, gc_content , mapq)

# Include both standard chromosomes and contigs (adjust as per your data)
#chrom_order <- unique(sam_data$X3)
chrom_order <- c(paste0("chr", c(1:22, "X", "Y", "M")))

# Convert 'Chromosome' to factor with natural sorting, including contigs
data$Chromosome <- factor(data$Chromosome, levels = chrom_order)

#canonical_chromosomes <- c(paste0(1:22), "X", "Y") #for hg19
canonical_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY") # for hg38

data_canon <- subset(data, (Chromosome %in% canonical_chromosomes))
data_contig <- subset(data, !(Chromosome %in% canonical_chromosomes))

data_contig_highGC <- subset(data_contig, gc_content>0.85)



data_canon$gc_content_group <- ifelse(data_canon$gc_content > 0.60, "GC > 0.60", "GC <= 0.60")

# Now create the plot
ggplot(data_canon, aes(x = Position, color = mapq == 60)) +
  geom_histogram(binwidth = 10000000) +
  facet_grid(gc_content_group ~ Chromosome, scales = "free_x", space = "free") +
  labs(title = "NTC (SA039-A138856-R9-C9) Chromosome Position Histogram", x = "Position", y = "Count") + 
  theme(axis.text.x = element_blank())  