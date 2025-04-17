# Post DAPC processing, plotting on chr
# Will save out a plot into 03_results

post_dapc <- function(id = dataset.id){
  
      # Reporting
      print(paste0("Considering dataset ", id, " with ", num_df, " DFs"))
  
  
      # Retain outputs
      print("Retaining outputs")
      dapc_var_contrib_rename.FN <- paste0("03_results/dapc_variance_contrib_", id, ".csv")
      file.copy(from = "03_results/dapc_variance_contrib.csv", to = dapc_var_contrib_rename.FN
                , overwrite = T
      )
            
      file.copy(from = "03_results/sample_DAPC.pdf", to = paste0("03_results/sample_DAPC_", id,".pdf"), overwrite = T)
      file.copy(from = "03_results/DAPC_loadings.pdf", to = paste0("03_results/DAPC_loadings_", id, ".pdf"), overwrite = T)
      
      ## Inspect loadings
      print("Inspecting loadings")
      dapc_var.df <- read.delim2(file = dapc_var_contrib_rename.FN, header = T, sep = ",")
      dapc_var.df <- as.data.frame(dapc_var.df)
      # head(dapc_var.df)
      # tail(dapc_var.df)
      # str(dapc_var.df)
      dapc_var.df$LD1 <- as.numeric(dapc_var.df$LD1)
      
      if(num_df > 1){
        
        dapc_var.df$LD2 <- as.numeric(dapc_var.df$LD2)
        
      }
      
      # Keep only one representative allele (remove redundancy)
      paste0("Considering ", nrow(dapc_var.df), " alleles")
      length(grep(pattern = "\\.0$", x = dapc_var.df$mname, perl = T)) # keep only one allele (redundancy)
      
      dapc_var.df <- dapc_var.df[grep(pattern = "\\.0$", x = dapc_var.df$mname, perl = T, invert = T), ] 
      nrow(dapc_var.df)
      paste0("Considering ", nrow(dapc_var.df), " loci")
      
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
      # head(dapc_var.df)
      # tail(dapc_var.df)
      
      # Plot 
      print("Plotting loadings on chromosomes")
      pdf(file = paste0("03_results/DAPC_loadings_on_chr_", id, "_LD1.pdf"), width = 13, height = 5.5)
      par(mfrow=c(1,1), mar = c(5,6,4,2)+0.1, mgp = c(5.5, 2.5, 0))
      fastman(m = dapc_var.df, chr = "chr.num", bp = "pos"
              , p = "LD1", logp = F
              , ylim = c(0, max(dapc_var.df$LD1)+ max(dapc_var.df$LD1)*0.5)
      )
      dev.off()
      
      if(num_df > 1){
        
        pdf(file = paste("03_results/DAPC_loadings_on_chr_", id, "_LD2.pdf"), width = 13, height = 5.5)
        fastman(m = dapc_var.df, chr = "chr.num", bp = "pos"
                , p = "LD2", logp = F
                , ylim = c(0, max(dapc_var.df$LD2)+ max(dapc_var.df$LD2)*0.5)
        )
        dev.off()
        
      }


}