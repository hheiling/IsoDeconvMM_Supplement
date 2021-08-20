

# Creation of mixture files for R2Q6 simulation

library(DirichletReg)
library(stringr)

# Number cell types
J = 3

header_nas = "/nas/longleaf/home/hheiling/"
header_pine = "/pine/scr/h/h/hheiling/"
# Simulated output prefix
prefix_sim_out = str_c(header_nas,"deconvolution/Simulated_Output")
# Sourced code location prefix
prefix_code = str_c(header_nas,"deconvolution/code")
# Simulated pure cell type files prefix
prefix_pure = str_c(header_pine,"Simulated_Files/Ciber_Sim_Pure_Sample_Counts")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Simulated_Files/Ciber_Sim_Mixture_Sample_Counts")
if(!dir.exists(prefix_mix)){dir.create(prefix_mix, recursive = T)}

# Source mixture sample creation function
source(sprintf("%s/sim_functions.R", prefix_code))


# Find pure reference samples
pure_names_full = list.files(path = prefix_pure,
                             pattern = "counts.txt$", full.names = TRUE)

cat("Number pure files found: ", length(pure_names_full), "\n")

n = length(pure_names_full) / J

CT1_files = pure_names_full[str_detect(pure_names_full, "CT1_")]
CT2_files = pure_names_full[str_detect(pure_names_full, "CT2_")]
CT3_files = pure_names_full[str_detect(pure_names_full, "CT3_")]

# Use replicates 1 through 10 for pure reference samples to use in CIBERSORTx algorithm
## Will eventually chose 5 to use
CT1_files_ref = CT1_files[1:(n/2)]
CT2_files_ref = CT2_files[1:(n/2)]
CT3_files_ref = CT3_files[1:(n/2)]

CT_ref_files = c(CT1_files_ref, CT2_files_ref, CT3_files_ref)

# Check correct selection of files 01 - 10
print(basename(CT_ref_files))

pure_ref_files = cbind(CT_ref_files, c(rep("CT1",times=n/2),rep("CT2",times=n/2),rep("CT3",times=n/2)))

cat("dim pure_ref_files: ",dim(pure_ref_files),"\n")

# Use replicates 11 through 20 for mixture file creation
CT1_files_mixSim = CT1_files[(n/2 + 1):n]
CT2_files_mixSim = CT2_files[(n/2 + 1):n]
CT3_files_mixSim = CT3_files[(n/2 + 1):n]

set_mixSim = list()
set_mixSim[["CT1"]] = CT1_files_mixSim
set_mixSim[["CT2"]] = CT2_files_mixSim
set_mixSim[["CT3"]] = CT3_files_mixSim

# Check files chosen correctly
set_mixSim

# Create probability combinations: p_combos

# Number of mixture files to create
M = 100

# Set seed for random selection of pure samples
set.seed(1269) 
seeds = sample(1000:9999, size = M+1, replace = F)

# Determine probability combinations
set.seed(seeds[M+1])
p_combos1 = rdirichlet(n = M*1.5, alpha = rep(2, times = 3))
## Remove probability combinations where a cell type prob is < 0.05
indic = ifelse(p_combos1 < 0.05, 1, 0)
p_combos2 = p_combos1[which(rowSums(indic) == 0),]
p_combos = p_combos2[sample(1:nrow(p_combos2), size = M, replace = F),]
colnames(p_combos) = c("CT1","CT2","CT3")
rownames(p_combos) = str_c("pc_", c(str_c("0",1:9),10:M))

# Save probability combinations

save(p_combos, file = sprintf("%s/Ciber_Sim_MixProbCombos100.RData",prefix_sim_out))
pc_labels = rownames(p_combos)

# Specify total read counts for mixture replicates
set.seed(seeds[M+1])
total_cts_mix = rnorm(n = nrow(p_combos), mean = 7*10^6, sd = 10^6)
names(total_cts_mix) = pc_labels

for(batch in 1:M){
  
  # if(batch %in% 1:9){
  #   batch_label = str_c("0",batch)
  # }else{
  #   batch_label = as.character(batch)
  # }
  
  if(batch <= 9){
    batch_label = str_c("00",batch)
  }else if(batch <= 99){
    batch_label = str_c("0",batch)
  }else{
    batch_label = as.character(batch)
  }
  
  current_seed = seeds[batch]
  
  # Create mixture samples
  mix_label = str_c("Mix_ProbCombo_",batch_label,"_counts")
  
  # print(mix_label)
  
  # Run simulations
  
  cellTypes = c("CT1","CT2","CT3")
  
  pc_nums = batch
  
  # Probability combination label
  pc = pc_labels[pc_nums]
  
  # Create mixture files
  ## Call mix_creation2 from 'sim_functions.R' file
  mix_creation2(set_mixSim = set_mixSim, out_folder = prefix_mix,
                file_labels = mix_label, total_cts = total_cts_mix[pc_nums],
                probs = matrix(p_combos[pc_nums,], nrow = 1), seed = current_seed)
  
  print(sprintf("Finished mix file %s",mix_label))
  
}

q(save = 'no')

########################################################################################