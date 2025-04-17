# Analysis of all UFR wgrs genotypes, MAF and LD-filtered
# Initialized 2025-04-09
# Ben J. G. Sutherland (SBIO)

# Requires: 
#   - simple_pop_stats repository
#   - VCF file for all samples, MAF and LD filtered

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
vcf.FN <- "02_input_data/Oner.BiSNP.MM0.9.MAR0.01.MMD8-100.LCI.chr_retained_w_tags_MAF0.05_5w50kb.vcf.gz" # MAF and LD filtered

all_sample_analysis <- TRUE
run_pca     <- TRUE
run_dendro  <- TRUE
run_FST     <- TRUE
bootstraps  <- 10000 # set number bootstraps to run for dendro

subset_to_variants <- -1 # how many variants to limit the VCF file (debugging)


#### 01. Load genotypes, rename pops, create colour file ####
# Load input VCF
my_vcf <- read.vcfR(file = vcf.FN, nrows = subset_to_variants)
my_vcf

if(all_sample_analysis==TRUE){
  
  ## Convert VCF objects to genind objects
  # All samples
  obj <- vcfR2genind(x = my_vcf)
  obj
  
  
  ## Assign pop vectors to each
  # All samples
  pop.vec <- indNames(obj)
  
  # Use custom function to update the names of the samples
  sample_rename(data.vec = pop.vec)
  
  # Inspect (optional)
  head(cbind(pop.vec, data.vec_renamed))
  tail(cbind(pop.vec, data.vec_renamed))
  
  pop(obj) <- data.vec_renamed
  table(pop(obj))
  
  
  ## Prepare colour file
  # Only create a new colour file, based on the alphanumeric ordered unique families
  if(file.exists(x = "00_archive/colours.csv")){
    
    print("Using existing colour file in 00_archive, not rebuilding.")
    
  }else{
    
    print("Creating new colour file based on the alphanumerically-sorted populations in all sample dataset")
    
    # Create a colours file
    present_pops.df <- sort(unique(pop(obj)))
    present_pops.df <- as.data.frame(present_pops.df)
    
    unique_cols <- c("cornflowerblue", "purple", "darkgreen", "darkkhaki", "blue4"
                     , "blue", "cyan1", "aquamarine", "coral", "pink"
                     , "red", "yellow", "navajowhite", "deepskyblue", "darkgray"
                     , "black", "limegreen", "brown4", "cadetblue4", "sienna"
                     , "darkslateblue"
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
  
  
  #### 02. All sample analysis ####
  ## Drop monomorphs 
  # All sample
  drop_loci(df = obj, drop_monomorphic = T)
  obj <- obj_filt
  
  ## MAF filter
  # not required, already conducted
  
  ## PCA
  if(run_pca==TRUE){
    
    pca_from_genind(data = obj, PCs_ret = 4, plot_eigen = T, plot_allele_loadings = F, retain_pca_obj = T
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
  
  
  ## Average FST
  if(run_FST ==TRUE){
    
    # Average option  
    calculate_FST(format = "genind", dat = obj, separated = F, bootstrap = F) ### Note: should work on bootstrap option ###
    
    # Bootstrap option (95% CI)
    calculate_FST(format = "genind", dat = obj, separated = F, bootstrap = T) ### Note: should work on bootstrap option ###
    
  }
  
}else{
  
  print("Not performing all sample analysis")
  
}


#### 03. Save output ####
if(all_sample_analysis==TRUE){
  
  output.FN <- "03_results/popgen_all_upper_fr.RData"
  
}else{
  
  output.FN <- "03_results/popgen_upper_fr_input_no_analysis.RData"
  
}
  
  

save.image(file = output.FN)

