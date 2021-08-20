

# Creating mixture count files from Blueprint pure cell type count files
## Use pure reference samples in 'RandSim_Samples' (partitioned by '02 file_split_randsim.R')

library(stringr)
library(DirichletReg)

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

header_nas = "/nas/longleaf/home/hheiling/"
header_pine = "/pine/scr/h/h/hheiling/"

# Pure cell type files prefix - to partition into use for mixture creation and pure sample fit
## Has 113 samples within each cell type (belonging to 113 individuals)
prefix_samples = str_c(header_pine,"RandSim_Samples")
# Materials prefix
prefix_mat = str_c(header_pine,"Blueprint/Blueprint_Materials")
# Folder to save p_combos files
prefix_pcombos = str_c(prefix_mat,"/RandSim_Probs")
if(!dir.exists(prefix_pcombos)){dir.create(prefix_pcombos, recursive = T)}



# Load blue_key data.frame that links file names with individuals
# Note: this was created in "01 Files_of_Interest.R"
load(sprintf("%s/blue_key_randomsim.RData", prefix_mat))
## Each row has the following information for a single individual:
## DONOR_ID, name of pure reference sample for CT1, name of "" for CT2, name of "" for CT3

# Datasets for cell types
CT = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")

# Number of mixture files to create in each simulation set-up
M = 100

# Total number of simulations
S = 10

CT1_files = list.files(path = str_c(prefix_samples, CT[1], sep = "/"), full.names = T)
CT2_files = list.files(path = str_c(prefix_samples, CT[2], sep = "/"), full.names = T)
CT3_files = list.files(path = str_c(prefix_samples, CT[3], sep = "/"), full.names = T)


###################################################################################################

for(s in 1:S){
  
  # set seed 
  set.seed(4783 + 5*s)
  mix_seeds = sample(10000:99999, size = M, replace = F)
  
  if(s <= 9){
    label_s = str_c("0",s)
  }else{
    label_s = as.character(s)
  }
  
  # Simulated output prefix
  prefix_mix = str_c(header_pine,"Mixture_Samples_RandSim_",label_s)
  if(!dir.exists(prefix_mix)){dir.create(prefix_mix, recursive = T)}
  prefix_pure = str_c(header_pine,"Pure_Samples_RandSim_",label_s)
  if(!dir.exists(prefix_pure)){dir.create(prefix_pure, recursive = T)}
  for(j in 1:3){
    prefix_CT = str_c(prefix_pure, CT[j], sep="/")
    if(!dir.exists(prefix_CT)){dir.create(prefix_CT, recursive = T)}
  }
  
  # Create the p_combos matrix
  library(DirichletReg)
  # Determine probability combinations
  p_combos1 = rdirichlet(n = M*1.5, alpha = rep(2, times = 3))
  ## Remove probability combinations where a cell type prob is < 0.05
  indic = ifelse(p_combos1 < 0.05, 1, 0)
  p_combos2 = p_combos1[which(rowSums(indic) == 0),]
  p_combos = p_combos2[sample(1:nrow(p_combos2), size = M, replace = F),]
  rownames(p_combos) = str_c("pc_",1:M)
  # save p_combos
  save(p_combos, file = sprintf("%s/p_combos_%s.RData",prefix_pcombos,label_s))
  
  # Partition the 113 pure reference samples in each cell type to 100 samples
  ## to use for mixture creation, and 5 samples to use for pure reference sample in
  ## the algorithm fit
  ## Assign mix/pure designation using row index value in blue_key
  use_idx = sample(1:nrow(blue_key), size = 105, replace = F)
  mix_idx = use_idx[1:100]
  pure_idx = use_idx[101:105]
  
  # Save the pure reference samples for use in the algorithm fit
  CT1_pure = blue_key[pure_idx,2]
  CT2_pure = blue_key[pure_idx,3]
  CT3_pure = blue_key[pure_idx,4]
  
  for(k in 1:length(pure_idx)){
    # Find which sample names are associated with donor_id
    CT1_name = CT1_pure[k]
    CT2_name = CT2_pure[k]
    CT3_name = CT3_pure[k]
    
    # Identify the pure sample files
    CT1_file = CT1_files[which(str_detect(CT1_files, CT1_name))]
    CT2_file = CT2_files[which(str_detect(CT2_files, CT2_name))]
    CT3_file = CT3_files[which(str_detect(CT3_files, CT3_name))]
    
    system(sprintf("cp %s %s/%s", CT1_file, 
                   prefix_pure, CT[1]))
    system(sprintf("cp %s %s/%s", CT2_file, 
                   prefix_pure, CT[2]))
    system(sprintf("cp %s %s/%s", CT3_file, 
                   prefix_pure, CT[3]))
    
  }
  
  # Run mix_creation function
  
  for(i in 1:M){
    
    # Each row contains all sample information for a single individual
    indiv = blue_key[mix_idx[i],]
    
    # Find which sample names are associated with donor_id
    CT1_name = as.character(indiv[2])
    CT2_name = as.character(indiv[3])
    CT3_name = as.character(indiv[4])
    
    # Identify the sample file names based on the above sample names
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
    
    set.seed(mix_seeds[i])
    mix_creation(p_vec = p_combos[i,], CT1_file, CT2_file, CT3_file,
                 total_cts = rnorm(n = 1, mean = 1.2*10^7, sd = 2.5*10^6), # mean and sd from mix_experiment.R in Github repo deconvolution subdirectory Blueprint_Simulation/code
                 out_directory = prefix_mix, mix_name = mix_name)
    
  }
  
  print(sprintf("End sim set-up for sim %s",label_s))
  
}



q(save = "no")

###################################################################################################