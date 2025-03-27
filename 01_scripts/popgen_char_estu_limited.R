# Analysis of EStu limited dataset
# Initialized 2025-03-26
# Ben J. G. Sutherland (SBIO)

# Requires: 
#   - simple_pop_stats repository
#   - VCF file for EStu only best collections, MAF and LD filtered

# Requires that 01_scripts/popgen_char_2025-02-27.R has already been run, specifically for colour file

# note: all code and directories listed are within the simple_pop_stats repository
# note: all output will go into simple_pop_stats/03_results

#### 00. Front Matter ####
# Clear space

# Source functions by sourcing simple_pop_stats_start.R ( https://github.com/bensutherland/simple_pop_stats )

# Load additional libraries
library("fastman")

# User set variables
vcf_estu.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_estu_limited_w_tags_MAF0.01_5w50kb.vcf" # EStu-sp MAF and LD for popgen

#### 01. Load genotypes and prepare datasets ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf_estu.FN)
my_vcf

## Convert VCF obj to genind
# EStu limited
obj_estu <- vcfR2genind(x = my_vcf)
obj_estu

## Assign pop vectors to each
# EStu
pop.vec <- indNames(obj_estu)
pop.vec <- gsub("[0-9]+", replacement = "", x = pop.vec) # Remove digits (there is not always an underscore)
pop.vec <- gsub("\\_.*", replacement = "", x = pop.vec)  # Remove everything from first underscore onwards
pop.vec <- gsub(pattern = "R$", replacement = "", x = pop.vec, ignore.case = F) # remove river character
table(pop.vec) # confirm it worked
pop(obj_estu) <- pop.vec
table(pop(obj_estu))


# Read in colour file to ensure that even if not recreating, will have it in the enviro
present_pops.df <- read.delim2(file = "00_archive/colours.csv", header = T, sep = ",")


#### 02. Population genetic analysis ####
#obj_estu.bck <- obj_estu

# # Create smaller obj for debugging
# keep <- locNames(obj_estu)[1:10000]
# obj.temp <- obj_estu[, loc = keep]

# FST with bootstraps # note: 1000 bootstraps applied
calculate_FST(format = "genind", dat = obj_estu, separated = FALSE, bootstrap = TRUE)
pairwise_wc_fst_booted.df <- pairwise_wc_fst.df 
pairwise_wc_fst_hfstat.list

# FST without bootstraps
calculate_FST(format = "genind", dat = obj_estu, separated = FALSE, bootstrap = FALSE)

# DAPC
# The goal of the following is to determine if there is any evidence of population separation among the EStu collections specifically
dapc_from_genind(data = obj_estu
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 1         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 


# Retain outputs
dataset.id <- "estu_limited"
dapc_var_contrib_rename.FN <- paste0("03_results/dapc_variance_contrib_", dataset.id, ".csv")
file.copy(from = "03_results/dapc_variance_contrib.csv", to = dapc_var_contrib_rename.FN
          , overwrite = T
)

file.copy(from = "03_results/sample_DAPC.pdf", to = paste0("03_results/sample_DAPC_", dataset.id,".pdf"), overwrite = T)
file.copy(from = "03_results/DAPC_loadings.pdf", to = paste0("03_results/DAPC_loadings_", dataset.id, ".pdf"), overwrite = T)

## Inspect loadings
dapc_var.df <- read.delim2(file = dapc_var_contrib_rename.FN, header = T, sep = ",")
dapc_var.df <- as.data.frame(dapc_var.df)
head(dapc_var.df)
tail(dapc_var.df)
str(dapc_var.df)
dapc_var.df$LD1 <- as.numeric(dapc_var.df$LD1)

# Keep only one representative allele (remove redundancy)
nrow(dapc_var.df)
length(grep(pattern = "\\.0$", x = dapc_var.df$mname, perl = T)) # keep only one allele (redundancy)

dapc_var.df <- dapc_var.df[grep(pattern = "\\.0$", x = dapc_var.df$mname, perl = T, invert = T), ] 
nrow(dapc_var.df)

# Obtain chromosome and positional info
dapc_var.df <- separate(data = dapc_var.df, col = "mname", into = c("chr", "pos"), sep = "_", remove = F)
table(dapc_var.df$chr)
dapc_var.df$chr.num <- dapc_var.df$chr # prepare for numeric

dapc_var.df$chr.num <- gsub(pattern = "a", replacement = ".1", x = dapc_var.df$chr.num)
dapc_var.df$chr.num <- gsub(pattern = "b", replacement = ".2", x = dapc_var.df$chr.num)

dapc_var.df$chr.num <- as.numeric(dapc_var.df$chr.num)

dapc_var.df$pos <- as.numeric(dapc_var.df$pos)

str(dapc_var.df)

# Order 
dapc_var.df <- dapc_var.df[order(dapc_var.df$chr.num), ]
head(dapc_var.df)
tail(dapc_var.df)

# Plot 
pdf(file = paste0("03_results/DAPC_loadings_on_chr_", dataset.id, "_LD1.pdf"), width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = dapc_var.df, chr = "chr.num", bp = "pos"
        , p = "LD1", logp = F
        , ylim = c(0, max(dapc_var.df$LD1)+ max(dapc_var.df$LD1)*0.5)
)
dev.off()


# End 

