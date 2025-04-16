# Sample rename
# Will save out a plot into 03_results

sample_rename <- function(data.vec = pop.vec){
  

  data.vec <- gsub("[0-9]+", replacement = "", x = data.vec) # Remove digits (there is not always an underscore)
  data.vec <- gsub("\\_.*", replacement = "", x = data.vec)  # Remove everything from first underscore onwards
  
  data.vec <- gsub(pattern = "M$|RL$|RM$|RU$|U$", replacement = "", x = data.vec, ignore.case = F) # remove characters separating Horsefly
  data.vec <- gsub(pattern = "^Uh", replacement = "H", x = data.vec, ignore.case = F) # remove characters separating Horsefly
  data.vec <- gsub(pattern = "R$", replacement = "", x = data.vec, ignore.case = F) # remove river character
  data.vec <- gsub(pattern = "L$|LN$", replacement = "", x = data.vec, ignore.case = F) # remove characters separating Chilko
  print(table(data.vec)) # confirm it worked
  
  assign(x = "data.vec_renamed", value = data.vec, envir = .GlobalEnv)
  
  print("Output provided as data.vec_renamed")
  
    
}

  