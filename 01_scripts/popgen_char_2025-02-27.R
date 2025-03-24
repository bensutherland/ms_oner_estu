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
max_missing         <- 0.3
run_pca    <- FALSE
run_dendro <- FALSE
bootstraps <- 10000 # set number bootstraps to run for dendro

#### 01. Load genotypes ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN)
my_vcf

# Convert to genind for simple_pop_stats functions
obj <- vcfR2genind(x = my_vcf)
obj
indNames(x = obj) # indiv names
pop(obj) # pop has not yet been assigned


#### 02. Annotate samples ####
### Create popmap
generate_popmap(df = obj)
# now manually annotate the popmap

# Use the popmap to annotate the obj with population attribute
annotate_from_popmap(df = obj, popmap.FN = "00_archive/my_data_ind-to-pop_annot.txt", convert_to_alt_ID = F)
obj_annot

table(pop(obj_annot))


#### 03. Assign colours ####
# Only create a new colour file, based on the alphanumeric ordered unique families
if(file.exists(x = "00_archive/colours.csv")){
  
  print("Using existing colour file in 00_archive, not rebuilding.")
  
}else{
  
  print("Creating new colour file based on the alphanumerically-sorted populations")
  
  # Create a colours file
  present_pops.df <- sort(unique(pop(obj_annot)))
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


#### 04. All sample analysis ####
# Drop monomorphic loci (if present)
drop_loci(df = obj_annot, drop_monomorphic = T)
obj_annot <- obj_filt


## PCA
if(run_pca==TRUE){
  
  pca_from_genind(data = obj_annot, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = F
                  , plot_ellipse = T, colour_file = "00_archive/colours.csv"
  )
  
}else{
  
  print("Not running PCA, as per user inputs")
  
}

## Dendrogram
if(run_dendro==TRUE){

make_tree(bootstrap = T, boot_obj = obj_annot, nboots = bootstraps
          , dist_metric = "edwards.dist", separated = FALSE
          )
}else{
  
  print("Not running dendrogram, as per user inputs")
  
}


## DAPC, supervised, all samples
dapc_from_genind(data = obj_annot
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 2         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

# Retain outputs
dapc_var_contrib_all.FN <- "03_results/dapc_variance_contrib_all_samples.csv"
file.copy(from = "03_results/dapc_variance_contrib.csv", to = dapc_var_contrib_all.FN
          , overwrite = T
          )

file.copy(from = "03_results/sample_DAPC.pdf", to = "03_results/sample_DAPC_all_samples.pdf", overwrite = T)
file.copy(from = "03_results/DAPC_loadings.pdf", to = "03_results/DAPC_loadings_all_samples.pdf", overwrite = T)

## Inspect loadings
dapc_var_all.df <- read.delim2(file = dapc_var_contrib_all.FN, header = T, sep = ",")
dapc_var_all.df <- as.data.frame(dapc_var_all.df)
head(dapc_var_all.df)
tail(dapc_var_all.df)
str(dapc_var_all.df)
dapc_var_all.df$LD1 <- as.numeric(dapc_var_all.df$LD1)
dapc_var_all.df$LD2 <- as.numeric(dapc_var_all.df$LD2)

# Keep only one representative allele (remove redundancy)
nrow(dapc_var_all.df) # 919400 records
length(grep(pattern = "\\.0$", x = dapc_var_all.df$mname, perl = T)) # keep only one allele (redundancy)

dapc_var_all.df <- dapc_var_all.df[grep(pattern = "\\.1", x = dapc_var_all.df$mname, perl = T, invert = T), ] 
nrow(dapc_var_all.df) # 459700 records

# Obtain chromosome and positional info
dapc_var_all.df <- separate(data = dapc_var_all.df, col = "mname", into = c("chr", "pos"), sep = "_", remove = F)
table(dapc_var_all.df$chr)
dapc_var_all.df$chr.num <- dapc_var_all.df$chr # prepare for numeric

dapc_var_all.df$chr.num <- gsub(pattern = "a", replacement = ".1", x = dapc_var_all.df$chr.num)
dapc_var_all.df$chr.num <- gsub(pattern = "b", replacement = ".2", x = dapc_var_all.df$chr.num)

dapc_var_all.df$chr.num <- as.numeric(dapc_var_all.df$chr.num)

dapc_var_all.df$pos <- as.numeric(dapc_var_all.df$pos)
dapc_var_all.df$LD1 <- as.numeric(dapc_var_all.df$LD1)
str(dapc_var_all.df)

# Order 
dapc_var_all.df <- dapc_var_all.df[order(dapc_var_all.df$chr.num), ]
head(dapc_var_all.df)
tail(dapc_var_all.df)

# Plot 
pdf(file = "03_results/DAPC_loadings_on_chr_all_samples_LD1.pdf", width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = dapc_var_all.df, chr = "chr.num", bp = "pos"
        , p = "LD1", logp = F
        , ylim = c(0, max(dapc_var_all.df$LD1)+ max(dapc_var_all.df$LD1)*0.5)
        )
dev.off()

pdf(file = "03_results/DAPC_loadings_on_chr_all_samples_LD2.pdf", width = 13, height = 5.5)
fastman(m = dapc_var_all.df, chr = "chr.num", bp = "pos"
        , p = "LD2", logp = F
        , ylim = c(0, max(dapc_var_all.df$LD2)+ max(dapc_var_all.df$LD2)*0.5)
)
dev.off()


#### 05. Separate populations ####
obj_annot.list <- seppop(x = obj_annot) # separate pops

#### 04. Stuart analysis (Estu and Summer) ####
# Inspect Stuart system only, incl the EStu and Summer Stuart
obj.Stu <- repool(obj_annot.list$Bivouac
               , obj_annot.list$Paula
               , obj_annot.list$Dust
               , obj_annot.list$Driftwood
               , obj_annot.list$Felix
               , obj_annot.list$Takla
               , obj_annot.list$Pinchi
               , obj_annot.list$Middle
               , obj_annot.list$Kuzkwa
               , obj_annot.list$Tachie
        )


#### 05. EStu-only analysis  ####
obj.EStu <- repool(obj_annot.list$Bivouac
                   , obj_annot.list$Paula
                   , obj_annot.list$Dust
                   , obj_annot.list$Driftwood
                   , obj_annot.list$Felix
                   , obj_annot.list$Takla
                   )
obj.EStu
table(pop(obj.EStu))

# Drop monomorphs
drop_loci(df = obj.EStu, drop_monomorphic = T) # 303 monomorphs dropped
obj.Estu <- obj_filt # rename output back to original
obj.Estu # 101,800 loci

# PCA, Estu only
pca_from_genind(data = obj.Estu
                , PCs_ret = 4
                , plot_eigen = T
                , plot_allele_loadings = F
                , retain_pca_obj = F
                , plot_ellipse = T
                , colour_file = "00_archive/colours.csv"
)


# DAPC, supervised
dapc_from_genind(data = obj.Estu
                 , plot_allele_loadings = TRUE  # plot and export the discrim. fn. locus variance contributions? 
                 , colour_file = "00_archive/colours.csv"           # use custom colours 
                 , n.pca = 10, n.da = 1         # number PC axes to use and discriminant functions to retain
                 , scree.da = TRUE              # plot scree plot for DF
                 , scree.pca = TRUE, posi.pca = "topright"     # plot PCA scree plot
                 , dapc.width = 7, dapc.height = 5             # PDF filesize for scatterplot
) 

# Inspect loadings
dapc_var_all.df <- read.delim2(file = "03_results/dapc_variance_contrib.csv", header = T, sep = ",")
dapc_var_all.df <- as.data.frame(dapc_var_all.df)
head(dapc_var_all.df)
tail(dapc_var_all.df)

# Keep only one representative allele (remove redundancy)
nrow(dapc_var_all.df) # 203600 records
dapc_var_all.df <- dapc_var_all.df[grep(pattern = "\\.1", x = dapc_var_all.df$mname, perl = T, invert = T), ] 
nrow(dapc_var_all.df) # 101800 records

# Obtain chromosome and positional info
dapc_var_all.df <- separate(data = dapc_var_all.df, col = "mname", into = c("chr", "pos"), sep = "_", remove = F)
table(dapc_var_all.df$chr)

# Prepare chr as numeric
dapc_var_all.df$chr.num <- dapc_var_all.df$chr # prepare for numeric
dapc_var_all.df$chr.num <- gsub(pattern = "a", replacement = ".1", x = dapc_var_all.df$chr.num) # replace a w/ .1
dapc_var_all.df$chr.num <- gsub(pattern = "b", replacement = ".2", x = dapc_var_all.df$chr.num) # replace b w/ .2
dapc_var_all.df$chr.num <- as.numeric(dapc_var_all.df$chr.num) # make numeric
dapc_var_all.df$pos <- as.numeric(dapc_var_all.df$pos) # make numeric
dapc_var_all.df$LD1 <- as.numeric(dapc_var_all.df$LD1) # make numeric
str(dapc_var_all.df)

# Order 
dapc_var_all.df <- dapc_var_all.df[order(dapc_var_all.df$chr.num), ]
head(dapc_var_all.df)
tail(dapc_var_all.df)

# Plot 
library("fastman")
pdf(file = "03_results/DAPC_loadings_on_chr_EStu_only.pdf", width = 13, height = 5.5)
par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
fastman(m = dapc_var_all.df, chr = "chr.num", bp = "pos"
        , p = "LD1", logp = F
        , ylim = c(0, max(dapc_var_all.df$LD1)+ max(dapc_var_all.df$LD1)*0.5)
)
dev.off()

# NOTE: save output into a subfolder or else will be overwritten


# FST per-locus
per_locus_stats(data = obj.Estu)
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

# There are some interesting values with highly significant FST scores to consider


# Optional: 
# Isolate to only the most distant sampling locations with the best sample size 
# (Paula off Trembleur, Dust midway up Takla, Driftwood at top of Takla)






