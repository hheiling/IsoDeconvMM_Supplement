
# Extract results from the Random Simulation
library(stringr)

header_nas = "/nas/longleaf/home/hheiling"

for(s in 1:10){
  if(s <= 9){
    sim = str_c("0",s)
  }else{
    sim = as.character(s)
  }
  
  files = list.files(path = sprintf("%s/Blueprint/RandomSim_Results/RandSim_%s/Summary_Probs/", header_nas, sim), 
                     full.names = T)
  # print(files)
  
  # Note: There are 5 results files for each mixture file (in order to reduce computation time)
  
  probs = list()
  mix_names = unique(str_sub(basename(files), start = 1, end = 7))
  
  for(i in 1:length(mix_names)){
    files_i = files[which(str_detect(files, mix_names[i]))]
    for(j in 1:length(files_i)){
      # Load SimSummary
      load(files_i[j])
      if(j == 1){
        p_mat = SimSummary[[1]]$p_mat
      }else{
        p_mat = rbind(p_mat, SimSummary[[1]]$p_mat)
      }
      
    }
    
    probs[[i]] = p_mat 
    
  }
  
  names(probs) = mix_names
  
  save(probs, file = sprintf("%s/Blueprint/RandomSim_Results/Probs_RandSim_%s.RData",header_nas,sim))
}

