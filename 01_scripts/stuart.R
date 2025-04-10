# Analysis of Stuart wgrs genotypes, MAF > 0.01 filter
# Initialized 2025-04-10
# Ben J. G. Sutherland (SBIO)

# Requires: 
#   - simple_pop_stats repository
#   - VCF file for Stuart region, MAF filtered (0.01)
#   - colour file script has been run

# note: all code and directories listed are within the simple_pop_stats repository
# note: all output will go into simple_pop_stats/03_results


#### 00. Front Matter ####
# Clear space

# Source functions by sourcing simple_pop_stats_start.R ( https://github.com/bensutherland/simple_pop_stats )

# Load additional libraries
library("fastman")

# User set variables
vcf.FN <- "02_input_data/stuart_chr18_w_tags_maf0.01.vcf.gz"

max_missing <- 0.3
run_pca     <- FALSE
run_dendro  <- FALSE
bootstraps  <- 10000 # set number bootstraps to run for dendro

subset_to_variants <- -1 # how many variants to limit the VCF file (debugging) (# -1 is default)


#### 01. Load genotypes and prepare datasets ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN, nrows = subset_to_variants)
my_vcf

## Convert VCF objects to genind objects
obj <- vcfR2genind(x = my_vcf)
obj

## Edit and assign pop vector
pop.vec <- indNames(obj)
pop.vec
pop.vec <- gsub("[0-9]+", replacement = "", x = pop.vec) # Remove digits (there is not always an underscore)
pop.vec <- gsub("\\_.*", replacement = "", x = pop.vec)  # Remove everything from first underscore onwards

# Note: there remains a lot of ambiguity here as to what the actual populations should be, so for now we'll just leave it

table(pop.vec)
pop(obj) <- pop.vec
table(pop(obj))

# Read in colour file
present_pops.df <- read.delim2(file = "00_archive/colours.csv", header = T, sep = ",") # will work with all datasets

## Drop monomorphs 
# NOTE: not required, as MAF > 0.01 filter has already been applied to the subset dataset


#### 02. Overview by PCA ####
if(run_pca==TRUE){
  
  pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T
                  , plot_ellipse = T, colour_file = "00_archive/colours.csv"
  )
  
}else{
  
  print("Not running PCA, as per user inputs")
  
}


#### 03. DAPC ####
dapc_from_genind(data = obj
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 1         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 8, dapc.height = 4.5             # PDF filesize for scatterplot
) 

# Retain outputs
dataset.id <- "stuart_chr18"
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


# Per locus FST next 
#### HERE ####



##### PLOT RESULTS OF GEMMA ANALYSIS ####
dataset.id <- "early_summer_with_taseko"
gemma_output.df <- read.delim2(file = "03_results/output/gwas.assoc.txt")
gemma_output.df <- as.data.frame(gemma_output.df)

gemma_output.df <- separate(data = gemma_output.df, col = "rs", into = c("chr", "pos"), sep = "__"
                            , remove = F)
gemma_output.df$chr.num <- gemma_output.df$chr # prepare for numeric

gemma_output.df$chr.num <- gsub(pattern = "a", replacement = ".1", x = gemma_output.df$chr.num)
gemma_output.df$chr.num <- gsub(pattern = "b", replacement = ".2", x = gemma_output.df$chr.num)

gemma_output.df$chr.num <- as.numeric(gemma_output.df$chr.num)

gemma_output.df$pos <- as.numeric(gemma_output.df$pos)

gemma_output.df$p_wald <- as.numeric(gemma_output.df$p_wald)

str(gemma_output.df)

# Order 
gemma_output.df <- gemma_output.df[order(gemma_output.df$chr.num, gemma_output.df$pos), ]
head(gemma_output.df)
tail(gemma_output.df)

# Plot 
pdf(file = paste0("03_results/Manhattan_plot", dataset.id, ".pdf"), width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = gemma_output.df, chr = "chr.num", bp = "pos"
        , p = "p_wald", logp = T
        , maxP = NULL
        #, ylim = c(0, max(dapc_var.df$LD1)+ max(dapc_var.df$LD1)*0.5)
)
dev.off()

#### 04. Save output ####
# optional. 
