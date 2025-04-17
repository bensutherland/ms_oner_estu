# Characterize the region of interest on the specific chromosome BCF
# Initialized 2025-04-16
# Ben J. G. Sutherland (SBIO)

# Requires: 
#   - simple_pop_stats repository
#   - VCF file

# note: all code and directories listed are within the simple_pop_stats repository
# note: all output will go into simple_pop_stats/03_results


#### 00. Front Matter ####
# Clear space

# Source functions by sourcing simple_pop_stats_start.R ( https://github.com/bensutherland/simple_pop_stats )
# Source additional functions (depends that ms_oner_estu is on same level as current simple_pop_stats)
source("../ms_oner_estu/01_scripts/post_dapc.R")
source("../ms_oner_estu/01_scripts/post_pca.R")
source("../ms_oner_estu/01_scripts/sample_rename.R")

# Load additional libraries
library("fastman")

# User set variables
vcf.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_chr18.vcf.gz"

run_perloc_stats     <- TRUE
subset_to_variants <- -1 # how many variants to limit the VCF file (debugging)


#### 01. Load genotypes and prepare datasets ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN, nrows = subset_to_variants)
my_vcf


##### Inspect the region in Stuart region only #####
dataset.id <- "stuart"

## Create region-specific VCF
retain_samples_stuart.vec <- colnames(my_vcf@gt)[grep(pattern = "Bivouac|Paula|Dust|Driftwood|Felix|Takla|Pinchi|Middle|Kuzkwa|Tachie"
                                                       , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_stuart.vec)
my_vcf_stuart <- my_vcf[, c("FORMAT", retain_samples_stuart.vec)]
my_vcf_stuart

# Clean space
rm(my_vcf)

## Convert to genind
obj_stuart <- vcfR2genind(x = my_vcf_stuart)
obj_stuart

## Set populations
pop.vec <- indNames(obj_stuart)
sample_rename(data.vec = pop.vec)
#cbind(pop.vec, data.vec_renamed)
pop(obj_stuart) <- data.vec_renamed
table(pop(obj_stuart))

## Drop monomorphs 
drop_loci(df = obj_stuart, drop_monomorphic = T)
obj_stuart <- obj_filt

## Drop low MAF
maf_filt(data = obj_stuart, maf = 0.01)
obj_stuart <- obj_maf_filt

## DAPC (with number DF based on PCA)
num_df <- 1
dapc_from_genind(data = obj_stuart
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = num_df         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

post_dapc(id = dataset.id) # post DAPC processing (incl. Manhattan plot with DAPC DF loadings)

## Per locus statistics
per_locus_stats(data = obj_stuart)

dim(per_loc_stats.df)
head(per_loc_stats.df)
per_loc_stats_chr.df <- per_loc_stats.df


# Then plot average FST across the chr w/ fastman
per_loc_stats_chr.df <- separate(data = per_loc_stats_chr.df, col = "mname", into = c("chr", "pos"), sep = "_", remove = F)


per_loc_stats_chr.df$chr.num <- per_loc_stats_chr.df$chr # prepare for numeric

per_loc_stats_chr.df$chr.num <- gsub(pattern = "a", replacement = ".1", x = per_loc_stats_chr.df$chr.num)
per_loc_stats_chr.df$chr.num <- gsub(pattern = "b", replacement = ".2", x = per_loc_stats_chr.df$chr.num)

per_loc_stats_chr.df$chr.num <- as.numeric(per_loc_stats_chr.df$chr.num)

per_loc_stats_chr.df$pos <- as.numeric(per_loc_stats_chr.df$pos)

per_loc_stats_chr.df <- per_loc_stats_chr.df[order(per_loc_stats_chr.df$pos), ]

str(per_loc_stats_chr.df)
head(per_loc_stats_chr.df)
tail(per_loc_stats_chr.df)

# e.g.,
pdf(file = paste("03_results/per_loc_FST_on_chr_stuart.pdf"), width = 13, height = 5.5)
fastman(m = per_loc_stats_chr.df, chr = "chr.num", bp = "pos"
        , p = "Fst", logp = F
        , ylim = c(0, 1.2)
)
dev.off()

# Read in the per loc stats for plotting

par(mar=c(5,4,4,2) + 0.1)
plot(x = per_loc_stats_chr.df$pos, y = per_loc_stats_chr.df$Fst, cex = 0.5, pch = 16
     , xlab = "", ylab = "")
abline(h = 0.5, lty = 2)
abline(v = 56293719, lty = 3)
abline(v = 57982347, lty = 3)

# fastman(m = per_loc_stats_chr.df, chr = "chr.num", bp = "pos"
#         , p = "Fit", logp = F
#         , ylim = c(0, 1.2)
# )
