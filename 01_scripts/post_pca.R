# Post PCA processing save files
# Will save out a plot into 03_results

post_pca <- function(id = dataset.id){
  
  # Reporting
  print(paste0("Considering dataset ", id))
  
  # Retain outputs
  print("Retaining outputs")
  
  # PCA scores
  file.copy(from = "03_results/pca_scores_per_sample.txt", to = paste0("03_results/pca_scores_per_sample_", id, ".txt")
            , overwrite = T
  )
  
  # PCA plots
  file.copy(from = "03_results/pca_samples_PC1_v_PC2.pdf", to = paste0("03_results/pca_samples_PC1_v_PC2_", id, ".pdf")
            , overwrite = T
  )
  
  file.copy(from = "03_results/pca_samples_PC3_v_PC4.pdf", to = paste0("03_results/pca_samples_PC3_v_PC4_", id, ".pdf")
            , overwrite = T
  )
  
  # PCA eigenvalues
  file.copy(from = "03_results/pca_eigenvalues.pdf", to = paste0("03_results/pca_eigenvalues_", id, ".pdf")
            , overwrite = T
  )
  
  
}