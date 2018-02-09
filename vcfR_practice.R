pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

# Read vcf, dna, and gff files
library(vcfR)
vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep = "\t", quote = "")

# Create a chromosome object from the files
chrom <- create.chromR(name = 'Supercontig', vcf = vcf, seq = dna, ann = gff)
plot(chrom)

# use masker to filter data that we do not have confidence in, based on plotting
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9, max_MQ = 60.1)
# min_Qual = 1 because most of our data has quality near 0, not a good metric to filter on
# DP clusters around 500, data with DP > 700 is likely CNVs
# MQ is peaked at 60 w/out much spread, narrow threshold is appropriate
plot(chrom)

# now that we have only high quality variants, we can process the data
chrom <- proc.chromR(chrom, verbose = TRUE)
plot(chrom) # the plot now shows variant count per window

# composite plot of variant, sequence, and annotation data
chromoqc(chrom, dp.alpha = 20)
# zoom in on a feature of interest
chromoqc(chrom, xlim=c(5e+05, 6e+05))
