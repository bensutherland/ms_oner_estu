# Plot FST values across chromosomes
# Initialized 2025-05-01
# Ben J. G. Sutherland (SBIO)

# Requires: 
#   - simple_pop_stats repository
#   - VCF file

# note: all code and directories listed are within the simple_pop_stats repository
# note: all output will go into simple_pop_stats/03_results


#### 00. Front Matter ####
# Clear space

# Source functions by sourcing simple_pop_stats_start.R ( https://github.com/bensutherland/simple_pop_stats )

# Load additional libraries
library("fastman")

# User set variables
vcf.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_EStu_limited_w_tags_MAF0.05_5w50kb.vcf.gz"

run_perloc_stats     <- TRUE


#### 01. Load genotypes and prepare datasets ####
if(run_perloc_stats==TRUE){
  
  # Read in data
  my_vcf <- read.vcfR(file = vcf.FN)
  
  # Convert to genind
  obj <- vcfR2genind(x = my_vcf)
  
  # Set population IDs
  ind.vec <- indNames(obj)
  pop(obj) <- gsub(pattern = "\\_.*", replacement = "", x = ind.vec)
  table(pop(obj))
  
  # Run per locus stats
  per_locus_stats(data = obj)
  
}

#read.vcfR()


# Start from FST vals
fst.df <- read.table("../results_from_thelio/per_locus_stats_EStu_only/per_locus_stats_2025-03-26.txt"
                     , header = T)
head(fst.df)
fst.df <- as.data.frame(fst.df)
head(fst.df)
fst.df <- separate(data = fst.df, col = "mname", into = c("chr", "pos"), sep = "_", remove = T)
fst.df$chr <- gsub(pattern = "a", replacement = ".1", x = fst.df$chr)
fst.df$chr <- gsub(pattern = "b", replacement = ".2", x = fst.df$chr)
fst.df$chr <- as.numeric(fst.df$chr)

fst.df$pos <- as.numeric(fst.df$pos)

table(fst.df$chr)
boxplot(fst.df$Fst)

#pdf(file = paste0("03_results/DAPC_loadings_on_chr_", id, "_LD1.pdf"), width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = fst.df, chr = "chr", bp = "pos"
        , p = "Fst", logp = F
)

pegas_fst.df <- fst.df 

# Read in data
fst.df <- read.table(file = "03_results/paula_vs_driftwood.weir.fst", header = T)

head(fst.df)
str(fst.df)

fst.df$CHROM <- gsub(pattern = "a", replacement = ".1", x = fst.df$CHROM)
fst.df$CHROM <- gsub(pattern = "b", replacement = ".2", x = fst.df$CHROM)
fst.df$CHROM <- as.numeric(fst.df$CHROM)
str(fst.df)


par(mfrow=c(1,2))
boxplot(fst.df$WEIR_AND_COCKERHAM_FST)
boxplot(pegas_fst.df$Fst)

# Plot 

#pdf(file = paste0("03_results/DAPC_loadings_on_chr_", id, "_LD1.pdf"), width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = fst.df, chr = "CHROM", bp = "POS"
        , p = "WEIR_AND_COCKERHAM_FST", logp = F
)
#dev.off()




