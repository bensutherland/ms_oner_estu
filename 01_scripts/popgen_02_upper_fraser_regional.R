# Analysis of UFR wgrs genotypes (at regional level), MAF and LD-filtered
# Initialized 2025-04-16
# Ben J. G. Sutherland (SBIO)

# Requires: 
#   - simple_pop_stats repository
#   - script 01_scripts/popgen_01_upper_fraser_all.R has already been run and RData is available in working directory

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


#### 01. Load data from previous ####
## Choose one of the following
#load(file = "03_results/popgen_all_upper_fr.RData") 
load(file = "03_results/popgen_upper_fr_input_no_analysis.RData")

my_vcf # source data

# In case need to see pop options for subsetting below
pops <- colnames(my_vcf@gt)[grep(pattern = "FORMAT", x = colnames(my_vcf@gt), invert = T)]
sample_rename(data.vec = pops) # gives a table output with renamed names
colnames(my_vcf@gt)[grep(pattern = "FORMAT", x = colnames(my_vcf@gt), invert = T)] # Gives full names


#### 03. Regional analyses setup ####
##### Nechako headwaters #####
## Create region-specific VCF
retain_samples_nechako.vec <- colnames(my_vcf@gt)[grep(pattern = "Bivouac|Paula|Dust|Driftwood|Felix|Takla|Pinchi|Middle|Kuzkwa|Tachie|Nadina|Stellako"
                                                      , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_nechako.vec)
my_vcf_nechako <- my_vcf[, c("FORMAT", retain_samples_nechako.vec)]
my_vcf_nechako

## Convert to genind
obj_nechako <- vcfR2genind(x = my_vcf_nechako)
obj_nechako

## Set populations
pop.vec <- indNames(obj_nechako)
sample_rename(data.vec = pop.vec)
#cbind(pop.vec, data.vec_renamed)
pop(obj_nechako) <- data.vec_renamed
table(pop(obj_nechako))

## Drop monomorphs 
drop_loci(df = obj_nechako, drop_monomorphic = T)
obj_nechako <- obj_filt

## Drop low MAF
maf_filt(data = obj_nechako, maf = 0.01)
obj_nechako <- obj_maf_filt

## PCA
pca_from_genind(data = obj_nechako, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T
                , plot_ellipse = T, colour_file = "00_archive/colours.csv"
)

pc1_v_pc2.plot_nechako <- pc1_v_pc2.plot
pc3_v_pc4.plot_nechako <- pc3_v_pc4.plot

# Save out
###TODO: use post_pca function
pdf(file = paste0(result.path, "pca_samples_PC1_v_PC2_nechako.pdf"), width = 11.5, height = 7.5)
print(pc1_v_pc2.plot_nechako)
dev.off()

pdf(file = paste0(result.path, "pca_samples_PC3_v_PC4_nechako.pdf"), width = 11.5, height = 7.5)
print(pc3_v_pc4.plot_nechako)
dev.off()


## DAPC (with number DF based on PCA)
dataset.id <- "nechako"
num_df <- 1
dapc_from_genind(data = obj_nechako
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = num_df         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

post_dapc(id = dataset.id) # post DAPC processing (incl. Manhattan plot with DAPC DF loadings)


##### Chilcotin headwaters #####
## Create region-specific VCF
retain_samples_chilcotin.vec <- colnames(my_vcf@gt)[grep(pattern = "Chilko|Taseko"
                                                       , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_chilcotin.vec)
my_vcf_chilcotin <- my_vcf[, c("FORMAT", retain_samples_chilcotin.vec)]
my_vcf_chilcotin

## Convert to genind
obj_chilcotin <- vcfR2genind(x = my_vcf_chilcotin)
obj_chilcotin

## Set populations
pop.vec <- indNames(obj_chilcotin)
sample_rename(data.vec = pop.vec)
#cbind(pop.vec, data.vec_renamed)
pop(obj_chilcotin) <- data.vec_renamed
table(pop(obj_chilcotin))

## Drop monomorphs 
drop_loci(df = obj_chilcotin, drop_monomorphic = T)
obj_chilcotin <- obj_filt

## Drop low MAF
maf_filt(data = obj_chilcotin, maf = 0.01)
obj_chilcotin <- obj_maf_filt

## PCA
pca_from_genind(data = obj_chilcotin, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T
                , plot_ellipse = T, colour_file = "00_archive/colours.csv"
)

pc1_v_pc2.plot_chilcotin <- pc1_v_pc2.plot
pc3_v_pc4.plot_chilcotin <- pc3_v_pc4.plot

# Save out
###TODO: use post_pca function
pdf(file = paste0(result.path, "pca_samples_PC1_v_PC2_chilcotin.pdf"), width = 11.5, height = 7.5)
print(pc1_v_pc2.plot_chilcotin)
dev.off()

pdf(file = paste0(result.path, "pca_samples_PC3_v_PC4_chilcotin.pdf"), width = 11.5, height = 7.5)
print(pc3_v_pc4.plot_chilcotin)
dev.off()




## DAPC (with number DF based on PCA)
dataset.id <- "chilcotin"
num_df <- 1
dapc_from_genind(data = obj_chilcotin
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = num_df         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

post_dapc(id = dataset.id) # post DAPC processing (incl. Manhattan plot with DAPC DF loadings)


##### Quesnel headwaters #####
dataset.id <- "quesnel"

retain_samples_quesnel.vec <- colnames(my_vcf@gt)[grep(pattern = "Horsefly|BlueLead|McKinley|Mitchell|horsefly|Wasko|Quesnel"
                                                         , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_quesnel.vec)
my_vcf_quesnel <- my_vcf[, c("FORMAT", retain_samples_quesnel.vec)]
my_vcf_quesnel

## Convert to genind
obj_quesnel <- vcfR2genind(x = my_vcf_quesnel)
obj_quesnel

## Set populations
pop.vec <- indNames(obj_quesnel)
sample_rename(data.vec = pop.vec)
#cbind(pop.vec, data.vec_renamed)
pop(obj_quesnel) <- data.vec_renamed
table(pop(obj_quesnel))

## Drop monomorphs 
drop_loci(df = obj_quesnel, drop_monomorphic = T)
obj_quesnel <- obj_filt

## Drop low MAF
maf_filt(data = obj_quesnel, maf = 0.01)
obj_quesnel <- obj_maf_filt

## PCA
pca_from_genind(data = obj_quesnel, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T
                , plot_ellipse = T, colour_file = "00_archive/colours.csv"
)

pc1_v_pc2.plot_quesnel <- pc1_v_pc2.plot
pc3_v_pc4.plot_quesnel <- pc3_v_pc4.plot

post_pca(id = dataset.id)


## DAPC (with number DF based on PCA)
num_df <- 1
dapc_from_genind(data = obj_quesnel
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = num_df         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

post_dapc(id = dataset.id) # post DAPC processing (incl. Manhattan plot with DAPC DF loadings)


##### Stuart #####
dataset.id <- "stuart"

retain_samples_stuart.vec <- colnames(my_vcf@gt)[grep(pattern = "Bivouac|Driftwood|Dust|Paula"
                                                       , x = colnames(my_vcf@gt), perl = T)]
length(retain_samples_stuart.vec)
my_vcf_stuart <- my_vcf[, c("FORMAT", retain_samples_stuart.vec)]
my_vcf_stuart

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

## PCA
pca_from_genind(data = obj_stuart, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T
                , plot_ellipse = T, colour_file = "00_archive/colours.csv"
)

pc1_v_pc2.plot_stuart <- pc1_v_pc2.plot
pc3_v_pc4.plot_stuart <- pc3_v_pc4.plot

post_pca(id = dataset.id)


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


#### Save output ####
save.image(file = "03_results/popgen_regional_upper_fr.RData")
