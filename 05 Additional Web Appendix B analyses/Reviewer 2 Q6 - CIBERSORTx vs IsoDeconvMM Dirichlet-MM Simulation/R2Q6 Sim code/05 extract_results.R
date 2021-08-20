
# Extract results from the Random Simulation
library(stringr)

header_nas = "/nas/longleaf/home/hheiling"

files = list.files(path = sprintf("%s/deconvolution/Simulation_Results/Ref2Q6/Summary_Probs/", header_nas), 
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

save(probs, file = sprintf("%s/deconvolution/Simulation_Results/Ref2Q6/Probs_Ref2Q6.RData",header_nas))
