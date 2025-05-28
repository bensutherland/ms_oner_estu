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

# Source additional functions (depends that ms_oner_estu is on same level as current simple_pop_stats)
source("../ms_oner_estu/01_scripts/sample_rename.R")


# User set variables
vcf.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18_sel_region.vcf.gz"
subset_to_variants <- -1 # how many variants to limit the VCF file (debugging) (note: -1 will input all)


#### 01. Load genotypes and prepare datasets ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN, nrows = subset_to_variants)
my_vcf

# Load colour file
colours.df <- read.table(file = "00_archive/colours.csv", header = T, sep = ",")


#### 02. Prepare genotypes for heatmap
geno.df <- extract.gt(x = my_vcf, element = "GT", as.numeric = F)

# Convert genotypes to allelic dosage
geno.df[geno.df=="0/0"] <- 0
geno.df[geno.df=="0/1"] <- 1
geno.df[geno.df=="1/1"] <- 2

mode(geno.df) = "numeric"

geno.df[1:5,1:5]

# Transpose (rows = samples; cols = SNPs)
geno.df<- t(geno.df)
geno.df[1:5,1:5]
dim(geno.df)


#### 03. Prepare samples for heatmap (colours) ####
ind.vec <- rownames(geno.df)

# Use custom function to update the names of the samples
sample_rename(data.vec = ind.vec)
#cbind(ind.vec, data.vec_renamed) # to confirm all OK

sample_colours.df <- as.data.frame(data.vec_renamed)
colnames(sample_colours.df) <- "collection"
sample_colours.df$index <- seq(1:nrow(sample_colours.df))
head(sample_colours.df) # this is in order
sample_colours.df <- merge(x = sample_colours.df, y = colours.df, by = "collection", all.x = T, sort = F)
sample_colours.df <- sample_colours.df[order(sample_colours.df$index), ]


## How to see order of samples (not actually required due to inherent behaviour of heatmap)
# test <- heatmap(x = geno.df, Rowv = NULL, Colv = NA)
# test$rowInd # this is the order of the samples based on the clustering, can use to reorder colours

# Plot heatmap
heatmap(x = geno.df
        , Rowv = NULL
        , Colv = NA # do not cluster by cols 
        , RowSideColors = sample_colours.df$colour
        , col = c("lightgrey", "yellow", "blue")
        )
# lightgrey = 0/0; yellow = 0/1; blue = 1/1


# heatmap(x = geno.df
#         , Rowv = NA # do not cluster by rows
#         , Colv = NA # do not cluster by cols 
#         , RowSideColors = sample_colours.df$colour
#         , col = c("lightgrey", "yellow", "blue")
# )


## Alternate colour options
# RColorBrewer::brewer.pal(n = 3, name = "Set1")

# Save out
pdf(file = "03_results/heatmap_w_sample_cols.pdf", width = 20, height = 15)
heatmap(x = geno.df
        , Rowv = NULL
        , Colv = NA # do not cluster by cols 
        , RowSideColors = sample_colours.df$colour
        , col = c("grey90", "yellow", "darkblue")
)
dev.off()
