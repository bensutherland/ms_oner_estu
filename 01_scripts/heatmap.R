# Generate heatmap based on a specific region of a VCF file
# Initialized 2025-04-15
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

# User set variables
vcf.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18_sel_region.vcf.gz"

subset_to_variants <- -1 # how many variants to limit the VCF file (debugging)


#### 01. Load genotypes and prepare datasets ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN, nrows = subset_to_variants)
my_vcf


#### 02. Prepare genotypes for heatmap
geno.df <- extract.gt(x = my_vcf, element = "GT", as.numeric = F)

geno.df[geno.df=="0/0"] <- 0
geno.df[geno.df=="0/1"] <- 1
geno.df[geno.df=="1/1"] <- 2

mode(geno.df) = "numeric"

geno.df[1:5,1:5]



#head(geno.df)
geno.df<- t(geno.df)
geno.df[1:50,1:5]
dim(geno.df)
#str(geno.df)

pdf(file = "03_results/heatmap.pdf", width = 20, height = 15)
heatmap(x = geno.df, Rowv = NULL, Colv = NA)
dev.off()
