

# Creating mixture count files from Blueprint pure cell type count files

library(stringr)

###################################################################################################

# Define mix_creation function

mix_creation = function(p_vec, CT1_file, CT2_file, CT3_file, total_cts, 
                        out_directory, mix_name){
  
  # Number cell types
  J = 3
  
  # All files are single columns of counts with rownames indicating exon set
  # Assume given full path to files
  CT1_cts = read.table(CT1_file, as.is = T)
  CT2_cts = read.table(CT2_file, as.is = T)
  CT3_cts = read.table(CT3_file, as.is = T)
  
  # Checks
  if(!identical(rownames(CT1_cts),rownames(CT2_cts)) | !identical(rownames(CT2_cts),rownames(CT3_cts))){
    stop("rows of pure count files not in same order \n")
  }
  if(length(p_vec) != J){
    stop("Should have 3 cell types \n")
  }
  
  cts_mat = cbind(CT1_cts, CT2_cts, CT3_cts)
  
  cts_ratio = total_cts / colSums(cts_mat)
  
  mix_parts = cts_mat * matrix(cts_ratio, nrow = nrow(cts_mat), ncol = J, byrow = T) *
    matrix(p_vec, nrow = nrow(cts_mat), ncol = J, byrow = T)
  
  mix_cts = rowSums(round(mix_parts))
  
  mix_df = data.frame(mix_name = mix_cts)
  rownames(mix_df) = rownames(CT1_cts) # rownames should be same for all cell types, just picked one
  
  write.table(mix_df, file = sprintf("%s/%s", out_directory, mix_name))
  
}


###################################################################################################

# Define variables

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Simulated output prefix
prefix_out = str_c(header_pine,"Mixture_Samples2")
if(!dir.exists(prefix_out)){dir.create(prefix_out, recursive = T)}
# Materials prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Pure cell type files prefix (to be used for algorithm fit)
prefix_pure = str_c(header_pine,"MixCreation_Samples2")


# Load blue_mix_key data.frame that links file names with individuals
# Note: this was created in "Blueprint_Simulation/data/Files_of_Interest.R"
load(sprintf("%s/blue_mix_key.RData", prefix_mat))
## Each row has the following information for a single individual:
## DONOR_ID, name of pure reference sample for CT1, name of "" for CT2, name of "" for CT3

# Datasets for cell types
CT = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")

# Number of mixture files to create
M = 100

# Set seed for random selection of pure samples
set.seed(8280)
seeds = sample(1000:9999, size = M+1, replace = F)

# Load p_combos matrix created in original 'mix_creation.R' file
load(sprintf("%s/Blueprint_ProbCombos.RData",prefix_mat))
## Code from the 'mix_creation.R' file to create the p_combos matrix
# library(DirichletReg)
# # Determine probability combinations
# set.seed(seeds[M+1])
# p_combos1 = rdirichlet(n = M*1.5, alpha = rep(2, times = 3))
# ## Remove probability combinations where a cell type prob is < 0.05
# indic = ifelse(p_combos1 < 0.05, 1, 0)
# p_combos2 = p_combos1[which(rowSums(indic) == 0),]
# p_combos = p_combos2[sample(1:nrow(p_combos2), size = M, replace = F),]
# rownames(p_combos) = str_c("pc_",1:M)
# 
# save(p_combos, file = sprintf("%s/Blueprint_ProbCombos.RData",prefix_mat))

CT1_files = list.files(path = str_c(prefix_pure, CT[1], sep = "/"), full.names = T)
CT2_files = list.files(path = str_c(prefix_pure, CT[2], sep = "/"), full.names = T)
CT3_files = list.files(path = str_c(prefix_pure, CT[3], sep = "/"), full.names = T)


###################################################################################################

# Run mix_creation function

for(i in 1:M){
  
  # Each row contains all sample information for a single individual
  
  # Find which sample names are associated with donor_id
  CT1_name = blue_mix_key[i,2]
  CT2_name = blue_mix_key[i,3]
  CT3_name = blue_mix_key[i,4]
  
  # Pick pure cell type samples, one for each cell type
  CT1_file = CT1_files[which(str_detect(CT1_files, CT1_name))]
  CT2_file = CT2_files[which(str_detect(CT2_files, CT2_name))]
  CT3_file = CT3_files[which(str_detect(CT3_files, CT3_name))]
  
  if(i <= 9){
    mix_name = str_c("Mix_","00",i,"_counts.txt")
  }else if(i <= 99){
    mix_name = str_c("Mix_","0",i,"_counts.txt")
  }else{
    mix_name = str_c("Mix_",i,"_counts.txt")
  }
  
  set.seed(seeds[i])
  mix_creation(p_vec = p_combos[i,], CT1_file, CT2_file, CT3_file,
               total_cts = rnorm(n = 1, mean = 1.2*10^7, sd = 2.5*10^6), # mean and sd from mix_experiment.R in Github repo deconvolution subdirectory Blueprint_Simulation/code
               out_directory = prefix_out, mix_name = mix_name)
  
}

q(save = "no")

###################################################################################################