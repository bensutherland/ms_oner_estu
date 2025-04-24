# Calculate a 95% confidence interval on a numeric vector in R
# B. Sutherland (SBIO)
# requires that vcftools has been run to produce FST values and result files are in 03_results

# source simple_pop_stats

# Identify FST filenames
input.FN <- list.files(path = "03_results/", pattern = ".fst")

# Set up collector df
result.df <- matrix(data = NA, nrow = 1, ncol = 4)
result.df <- as.data.frame(result.df)
colnames(result.df) <- c("contrast", "avg", "ll", "ul")

# Loop to identify 95% CI per contrast
fst.df <- NULL; result.mod <- NULL; all_result.df <- NULL
for(i in 1:length(input.FN)){
  
  # Reporting
  print(paste0("Working on contrast: ", input.FN[i]))
  
  # Read in data
  fst.df <- read.table(file = paste0("03_results/", input.FN[i]), header = T)
  
  # Run one-sample t-test to obtain 95% CI on vector
  result.mod <- t.test(fst.df$WEIR_AND_COCKERHAM_FST)
  
  # Reporting
  print(result.mod$conf.int)
  print(result.mod$estimate)
  
  # Save results
  result.df[,"contrast"] <- input.FN[i]            # contrast name
  result.df[,"avg"]      <- result.mod$estimate    # mean
  result.df[,"ll"]       <- result.mod$conf.int[1] # lower limit of 95% CI
  result.df[,"ul"]       <- result.mod$conf.int[2] # upper limit of 95% CI
  
  # Retain in larger list
  all_result.df <- rbind(all_result.df, result.df)
  
}


