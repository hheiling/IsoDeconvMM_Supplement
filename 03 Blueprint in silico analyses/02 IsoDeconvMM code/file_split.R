

# Goals:
# Using the files specified in 'blue_mix_key' and 'blue_pure' from 
# "Blueprint_Simulation/data/Files_of_Interest.R",
# split the into Fit_Samples2 and MixCreation_Samples2 folders (mixture files created using pure
# reference files from the same individual)

library(stringr)

# Maneuver to Longleaf/Blueprint folder in personal compuater

# All pure reference files prefix
prefix_files = "Pure_ExonSetCts"
# Blueprint materials output prefix
prefix_mat = "Blueprint_Materials"
# Pure cell type files prefix
prefix_pure = "Fit_Samples2"
# Simulated mixture cell type files prefix
prefix_mix = "MixCreation_Samples2"

CT_folders = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")

for(j in 1:length(CT_folders)){
  prefix_pure_CT = str_c(prefix_pure,CT_folders[j],sep = "/")
  if(!dir.exists(prefix_pure_CT)){dir.create(prefix_pure_CT, recursive = T)}
  prefix_mix_CT = str_c(prefix_mix,CT_folders[j],sep = "/")
  if(!dir.exists(prefix_mix_CT)){dir.create(prefix_mix_CT, recursive = T)}
}

# Load the blue_mix_key data.frame
load(sprintf("%s/blue_mix_key.RData", prefix_mat))

# Load the blue_pure data.frame
load(sprintf("%s/blue_pure.RData", prefix_mat))

for(j in 1:length(CT_folders)){
  
  mix_filenames = blue_mix_key[,(j+1)]
  pure_filenames = blue_pure[,(j+1)]
  
  all_files = list.files(path = str_c(prefix_files, CT_folders[j], sep="/"))
  
  prefix_CT = str_c(prefix_files, CT_folders[j], sep = "/")
  
  for(f in 1:length(mix_filenames)){
    mix_file = all_files[which(str_detect(all_files, mix_filenames[f]))]
    system(sprintf("cp %s/%s %s/%s", prefix_CT, mix_file, prefix_mix, CT_folders[j]))
  }
  
  
  for(f in 1:length(pure_filenames)){
    pure_file = all_files[which(str_detect(all_files, pure_filenames[f]))]
    system(sprintf("cp %s/%s %s/%s", prefix_CT, pure_file, prefix_pure, CT_folders[j]))
  }
  
}

# Checks

for(j in 1:length(CT_folders)){
  # Should have 100 files in each prefix_mix subdirectory
  print(length(list.files(path = str_c(prefix_mix, CT_folders[j], sep = "/"))))
  # Should have 35 files in each prefix_pure subdirectory
  print(length(list.files(path = str_c(prefix_pure, CT_folders[j], sep = "/"))))
}

##########################################################################################################

print(gc())

q(save='no')

##########################################################################################################
