# Analysis of MFR and UFR wgrs genotypes
# Initialized 2025-02-27
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
#vcf.FN              <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15_w_tags_MAF0.05_5w50kb.vcf" # MAF and LD filtered
vcf.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP1000_minGQ20_miss0.15.vcf" # no MAF or LD filter
max_missing <- 0.3
run_pca     <- FALSE
run_dendro  <- FALSE
bootstraps  <- 10000 # set number bootstraps to run for dendro

#### 01. Load genotypes and prepare datasets ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)
my_vcf

## Create region-specific VCF objects
# Stuart (EStu and Summer run-timings)
retain_samples_stuart.vec <- colnames(my_vcf@gt)[grep(pattern = "Bivouac|Paula|Dust|Driftwood|Felix|Takla|Pinchi|Middle|Kuzkwa|Tachie"
                                                      , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_stuart.vec)
my_vcf_stuart <- my_vcf[, c("FORMAT", retain_samples_stuart.vec)]
my_vcf_stuart

# ES vs. S (No EStu)
retain_samples_summer_v_early_summer.vec <- colnames(my_vcf@gt)[grep(pattern = "Bowron|Kuzkwa|Middle|Nadina|Pinchi|Stellako|Tachie"
                                                      , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_summer_v_early_summer.vec)
my_vcf_ses <- my_vcf[, c("FORMAT", retain_samples_summer_v_early_summer.vec)]
my_vcf_ses
# note: may need to drop a few outlier samples from Nadina and Stellako, as they may be considered early-summer

# EStu (for within metapop, collections w/ at least 7 inds)
retain_samples_estu.vec <- colnames(my_vcf@gt)[grep(pattern = "Driftwood|Dust|Bivouac|Paula"
                                                      , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_estu.vec)
my_vcf_estu <- my_vcf[, c("FORMAT", retain_samples_estu.vec)]
my_vcf_estu


## Convert VCF obj to genind
# All samples
obj <- vcfR2genind(x = my_vcf)
obj

# Stuart
obj_stuart <- vcfR2genind(x = my_vcf_stuart)
obj_stuart

# Early summer vs. summer
obj_ses <- vcfR2genind(x = my_vcf_ses)
obj_ses

# EStu limited
obj_estu <- vcfR2genind(x = my_vcf_estu)
obj_estu


## Assign pop vectors to each
# All samples
pop.vec <- indNames(obj)
pop.vec <- gsub("[0-9]+", replacement = "", x = pop.vec) # Remove digits (there is not always an underscore)
pop.vec <- gsub("\\_.*", replacement = "", x = pop.vec)  # Remove everything from first underscore onwards
pop.vec <- gsub(pattern = "R$", replacement = "", x = pop.vec, ignore.case = F) # remove river character
table(pop.vec) # confirm it worked
pop(obj) <- pop.vec
table(pop(obj))

# Stuart
pop.vec <- indNames(obj_stuart)
pop.vec <- gsub("[0-9]+", replacement = "", x = pop.vec) # Remove digits (there is not always an underscore)
pop.vec <- gsub("\\_.*", replacement = "", x = pop.vec)  # Remove everything from first underscore onwards
pop.vec <- gsub(pattern = "R$", replacement = "", x = pop.vec, ignore.case = F) # remove river character
table(pop.vec) # confirm it worked
pop(obj_stuart) <- pop.vec
table(pop(obj_stuart))

# S vs. ES
pop.vec <- indNames(obj_ses)
pop.vec <- gsub("[0-9]+", replacement = "", x = pop.vec) # Remove digits (there is not always an underscore)
pop.vec <- gsub("\\_.*", replacement = "", x = pop.vec)  # Remove everything from first underscore onwards
pop.vec <- gsub(pattern = "R$", replacement = "", x = pop.vec, ignore.case = F) # remove river character
table(pop.vec) # confirm it worked
pop(obj_ses) <- pop.vec
table(pop(obj_ses))

# EStu
pop.vec <- indNames(obj_estu)
pop.vec <- gsub("[0-9]+", replacement = "", x = pop.vec) # Remove digits (there is not always an underscore)
pop.vec <- gsub("\\_.*", replacement = "", x = pop.vec)  # Remove everything from first underscore onwards
pop.vec <- gsub(pattern = "R$", replacement = "", x = pop.vec, ignore.case = F) # remove river character
table(pop.vec) # confirm it worked
pop(obj_estu) <- pop.vec
table(pop(obj_estu))


## Prepare colour file
# Only create a new colour file, based on the alphanumeric ordered unique families
if(file.exists(x = "00_archive/colours.csv")){
  
  print("Using existing colour file in 00_archive, not rebuilding.")
  
}else{
  
  print("Creating new colour file based on the alphanumerically-sorted populations in all sample dataset")
  
  # Create a colours file
  present_pops.df <- sort(unique(pop(obj)))
  present_pops.df <- as.data.frame(present_pops.df)
  
  # using https://sashamaps.net/docs/resources/20-colors/
  unique_cols <- c("red", "darkgreen", "yellow", "blue", "orange", "purple", "cyan"
                   , "magenta", "limegreen", "pink", "darkcyan"
                   , "lavender", "sienna", "darkred", "aquamarine", "darkkhaki"
                   , "navajowhite", "navy", "darkgrey", "black"
  )
  
  present_pops.df$colour <- unique_cols[1:nrow(present_pops.df)]
  present_pops.df
  
  # Rename columns to necessary colnames for the colour file
  colnames(present_pops.df) <- c("collection", "colour")
  
  # Write out
  write.csv(x = present_pops.df, file = "00_archive/colours.csv", quote = F, row.names = F)
  
}

# Read in to ensure that even if not recreating, will have it in the enviro
present_pops.df <- read.delim2(file = "00_archive/colours.csv", header = T, sep = ",")
# This colour file will work with all datasets


## Drop monomorphs 
# All sample
drop_loci(df = obj, drop_monomorphic = T)
obj <- obj_filt

# Stuart
drop_loci(df = obj_stuart, drop_monomorphic = T)
obj_stuart <- obj_filt

# Summer vs. Early Summer
drop_loci(df = obj_ses, drop_monomorphic = T)
obj_ses <- obj_filt

# EStu limited
drop_loci(df = obj_estu, drop_monomorphic = T)
obj_estu <- obj_filt


#### 02. All sample analysis ####
## PCA
if(run_pca==TRUE){
  
  pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = F
                  , plot_ellipse = T, colour_file = "00_archive/colours.csv"
  )
  
}else{
  
  print("Not running PCA, as per user inputs")
  
}

## Dendrogram
if(run_dendro==TRUE){

  make_tree(bootstrap = T, boot_obj = obj, nboots = bootstraps
          , dist_metric = "edwards.dist", separated = FALSE
          )
}else{
  
  print("Not running dendrogram, as per user inputs")
  
}


#### 03. All loci, DAPC in EStu vs. Stuart-S ####
dapc_from_genind(data = obj_stuart
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 2         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

# Retain outputs
dataset.id <- "stuart"
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
dapc_var.df$LD2 <- as.numeric(dapc_var.df$LD2)

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

pdf(file = paste("03_results/DAPC_loadings_on_chr_", dataset.id, "_LD2.pdf"), width = 13, height = 5.5)
fastman(m = dapc_var.df, chr = "chr.num", bp = "pos"
        , p = "LD2", logp = F
        , ylim = c(0, max(dapc_var.df$LD2)+ max(dapc_var.df$LD2)*0.5)
)
dev.off()



#### 03. All loci, DAPC in Summer vs. Early-summer (no EStu) ####
dapc_from_genind(data = obj_ses
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 2         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

# Retain outputs
dataset.id <- "ses"
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
dapc_var.df$LD2 <- as.numeric(dapc_var.df$LD2)

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

pdf(file = paste("03_results/DAPC_loadings_on_chr_", dataset.id, "_LD2.pdf"), width = 13, height = 5.5)
fastman(m = dapc_var.df, chr = "chr.num", bp = "pos"
        , p = "LD2", logp = F
        , ylim = c(0, max(dapc_var.df$LD2)+ max(dapc_var.df$LD2)*0.5)
)
dev.off()


#### 05. DAPC in EStu limited ####
# The goal of the following is to determine if there is any evidence of population separation among the EStu collections specifically
dapc_from_genind(data = obj_estu
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 2         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

# Retain outputs
dataset.id <- "estu"
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
dapc_var.df$LD2 <- as.numeric(dapc_var.df$LD2)

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

pdf(file = paste("03_results/DAPC_loadings_on_chr_", dataset.id, "_LD2.pdf"), width = 13, height = 5.5)
fastman(m = dapc_var.df, chr = "chr.num", bp = "pos"
        , p = "LD2", logp = F
        , ylim = c(0, max(dapc_var.df$LD2)+ max(dapc_var.df$LD2)*0.5)
)
dev.off()

# FST per-locus
per_locus_stats(data = obj_estu)
head(per_loc_stats.df)

# Prepare for plotting
# Obtain chromosome and positional info
per_loc_stats.df <- separate(data = per_loc_stats.df, col = "mname", into = c("chr", "pos"), sep = "_", remove = F)
table(per_loc_stats.df$chr)

# Prepare chr as numeric
per_loc_stats.df$chr.num <- per_loc_stats.df$chr # prepare for numeric
per_loc_stats.df$chr.num <- gsub(pattern = "a", replacement = ".1", x = per_loc_stats.df$chr.num) # replace a w/ .1
per_loc_stats.df$chr.num <- gsub(pattern = "b", replacement = ".2", x = per_loc_stats.df$chr.num) # replace b w/ .2
per_loc_stats.df$chr.num <- as.numeric(per_loc_stats.df$chr.num) # make numeric
per_loc_stats.df$pos <- as.numeric(per_loc_stats.df$pos) # make numeric
head(per_loc_stats.df)
str(per_loc_stats.df)

# Order 
per_loc_stats.df <- per_loc_stats.df[order(per_loc_stats.df$chr.num), ]
head(per_loc_stats.df)
tail(per_loc_stats.df)

# Plot  
pdf(file = "03_results/FST_chr_EStu_only.pdf", width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = per_loc_stats.df, chr = "chr.num", bp = "pos"
        , p = "Fst", logp = F
        , ylim = c(-0.2, max(per_loc_stats.df$Fst)+ max(per_loc_stats.df$Fst)*0.5)
)
dev.off()


## Average FST
calculate_FST(format = "genind", dat = obj_estu, separated = F) ### Note: should work on bootstrap option ###



