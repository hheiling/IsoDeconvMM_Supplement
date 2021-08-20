

# Goals:
# Using the file 'blue_key.RData' created in code file 
# "Blueprint_Simulation/data/Files_of_Interest_RandomSim.R",
# Put relevant files into RandSim_Samples folder 

# Run code in Terminal tab of RStudio
## In Terminal, navigate to Longleaf/Blueprint folder in personal computer
## R (to run things in R)
## run the code below

library(stringr)


# All pure reference files prefix
prefix_files = "Pure_ExonSetCts"
# Blueprint materials output prefix
prefix_mat = "Blueprint_Materials"
prefix_out = "RandSim_Samples"

CT_folders = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")

for(j in 1:length(CT_folders)){
  prefix_CT = str_c(prefix_out,CT_folders[j],sep = "/")
  if(!dir.exists(prefix_CT)){dir.create(prefix_CT, recursive = T)}
}

# Load the blue_key data.frame
load(sprintf("%s/blue_key_randomsim.RData", prefix_mat))

# saveRDS(blue_key, file = sprintf("%s/blue_key_randomsim.rds", prefix_mat))
# blue_key = readRDS(sprintf("%s/blue_key_randomsim.rds", prefix_mat))

for(j in 1:length(CT_folders)){
  
  use_files = blue_key[,(j+1)]
  
  all_files = list.files(path = str_c(prefix_files, CT_folders[j], sep="/"))
  
  prefix_CT = str_c(prefix_files, CT_folders[j], sep = "/")
  
  for(f in 1:length(use_files)){
    use_f = all_files[which(str_detect(all_files, use_files[f]))]
    system(sprintf("cp %s/%s %s/%s", prefix_CT, use_f, prefix_out, CT_folders[j]),
           intern = T)
  }

}

# Checks

for(j in 1:length(CT_folders)){
  # Should have 113 files in each prefix_out subdirectory
  print(length(list.files(path = str_c(prefix_out, CT_folders[j], sep = "/"))))

}

##########################################################################################################

print(gc())

q(save='no')

##########################################################################################################
