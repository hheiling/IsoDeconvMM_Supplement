
# This simulation a derivative of the sim_blue11A_longleaf.R simulation
# Fit of Blueprint mixture samples using 10 initial points and 5 pure samples per cell type
# Same clusters and seeds as the 10A simulation

# Fit IsoDeconvMM model fit to simulated data

# Load libraries needed to run the IsoDeconvMM functions:
library(stringr)
library(alabama)
library(ICSNP)
library(gtools)
library(MASS)

# Arrays 1-10^4
## 10 simulation sets
## Per simulation set, 100 mixuture samples
## Per mixture sample, fit 10 genes/clusters at a time, leading to 10 arrays per mixture sample
## 10*100*10 = 10^4
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# which of the 10 simulation sets is being used?
if((array_val %% 1000) != 0){ ## if array_val not a multiple of 1000:
  sim_set = (array_val %/% 1000)+1
}else{ ## if array_val a multiple of 1000
  sim_set = array_val / 1000
}

if(sim_set <= 9){
  sim_label = str_c("0",sim_set)
}else{
  sim_label = as.character(sim_set)
}

# num_mix: number of mixture files per simulation set
num_mix = 100
# total_clusts: total number of clusters to fit
total_clusts = 100
# clust_size: number of clusters to fit per array submission
clust_size = 10
# num_arrays: number of arrays needed to perform IsoDeconvMM on all clusters of interest for a single mixture file
num_arrays = total_clusts/clust_size

# batch should be between 1 and 1000
batch = array_val - (sim_set-1)*(num_mix*num_arrays)
# batch_samp: designates which mixture file is being run
batch_samp = rep(1:num_mix, each = num_arrays)[batch]

if(batch_samp <= 9){
  batch_label = str_c("00",batch_samp)
}else if(batch_samp <= 99 ){
  batch_label = str_c("0",batch_samp)
}else{
  batch_label = as.character(batch_samp)
}
# clust_batch: designates which partition of the clusters is being run
clust_batch = rep(1:num_arrays, times = num_mix)[batch]

# Number cell types
J = 3

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/"

# Blueprint materials output prefix
prefix_mat = str_c(header_pine,"Blueprint/Blueprint_Materials")
# Pure cell type files prefix
prefix_pure = str_c(header_pine,"Pure_Samples_RandSim_",sim_label)
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Mixture_Samples_RandSim_",sim_label)

# Simulation results prefix
prefix_pre = str_c(header_nas,"RandomSim_Results/")
if(!dir.exists(prefix_pre)){dir.create(prefix_pre, recursive = T)}
prefix_results = str_c(prefix_pre,"RandSim_",sim_label)
if(!dir.exists(prefix_results)){dir.create(prefix_results, recursive = T)}
prefix_results1 = str_c(prefix_results,"/Results_Abbrev")
if(!dir.exists(prefix_results1)){dir.create(prefix_results1, recursive = T)}
prefix_results2 = str_c(prefix_results,"/Summary_Probs")
if(!dir.exists(prefix_results2)){dir.create(prefix_results2, recursive = T)}
# prefix_results4 = str_c(prefix_results,"/Full_Results")
# if(!dir.exists(prefix_results4)){dir.create(prefix_results4, recursive = T)}

# Source IsoDeconvMM functions
IsoDeconv_code = list.files(path = "/nas/longleaf/home/hheiling/deconvolution/code/IsoDeconvMM", 
                            full.names = T)
for(i in 1:length(IsoDeconv_code)){
  source(IsoDeconv_code[i])
}

# Use the 100 transcript clusters used in the simulations used for the paper
## character vector 'clusts_11A_best'
load(sprintf("%s/Best Clusters Paper 11A.RData", prefix_mat))
discrim_clusts = clusts_11A_best
length(discrim_clusts)

# Since prob estimates occur separately for each gene, can subset the number of genes
# to run at one time (here, run clust_size at a time). Will concatenate all p estimates later.
select_idx_first = (clust_batch-1)*clust_size+1
(select_idx = (select_idx_first):(select_idx_first+(clust_size-1)))

clusts_to_use = discrim_clusts[select_idx]

# Datasets for cell types
CT = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")
## Assume ordering of datasets in CT correspond to cell types CT1, CT2, and CT3
cellTypes = c("CT1","CT2","CT3")

# Find pure reference files
CT1_files = list.files(path = str_c(prefix_pure, CT[1], sep = "/"), full.names = T)
CT2_files = list.files(path = str_c(prefix_pure, CT[2], sep = "/"), full.names = T)
CT3_files = list.files(path = str_c(prefix_pure, CT[3], sep = "/"), full.names = T)

n = 5

pure_ref_files= cbind(c(CT1_files, CT2_files, CT3_files), 
                      c(rep("CT1",times=n),rep("CT2",times=n),rep("CT3",times=n)))

cat("dim pure_ref_files: ",dim(pure_ref_files),"\n")

# Find fragment length files
## Assume distributions from paired-end fragment lengths, even though not true for CT2 and CT3
fraglens_file = str_c(prefix_mat, "Combined_FragDist.txt", sep="/")

print(basename(fraglens_file))

# initPts:
# Extreme cases (two cell types 0.10, one cell type 0.80)
# Moderate cases (two cell types 0.25, one cell type 0.50), (two cell types 0.40, one cell type 0.20)
# Equal case (all three cell types = 1/3)

initPts = matrix(c(0.10,0.10,0.80,
                   0.10,0.80,0.10,
                   0.80,0.10,0.10,
                   0.25,0.25,0.50,
                   0.25,0.50,0.25,
                   0.50,0.25,0.25,
                   0.20,0.40,0.40,
                   0.40,0.20,0.40,
                   0.40,0.40,0.20,
                   1/3,1/3,1/3), ncol = 3, byrow = T)

colnames(initPts) = c("CT1","CT2","CT3")
print(initPts)

# Run simulations

SimResults = list()
SimSummary = list()
# SimFull = list()
  
  mix_files = list.files(path = prefix_mix, full.names = T)
  
  mix_names = str_remove(basename(mix_files), "_counts.txt")
  
  # Fit IsoDeconvMM algorithm
  SimResults_Full = IsoDeconvMM(directory = NULL, mix_files = mix_files[batch_samp],
                                pure_ref_files = pure_ref_files,
                                fraglens_files = fraglens_file,
                                bedFile = sprintf("%s/gencode.v15.nonoverlap.exon.bed", prefix_mat),
                                knownIsoforms = sprintf("%s/gencode.v15.nonoverlap.exon.knownIsoforms.RData", prefix_mat),
                                discrim_clusts = clusts_to_use,
                                readLen = 100, lmax = 1200, eLenMin = 1,
                                mix_names = mix_names[batch_samp], initPts = initPts,
                                optim_options = optimControl(simple.Init = FALSE))
  
  # Save all results for reference
  # SimFull = SimResults_Full
  
  # Save subset of results
  for(m in mix_names){
    
    results = SimResults_Full[[m]]
    clust_names = names(results)
    
    for(clust in clust_names){
      
      SimResults[[m]][[clust]][["mix"]][["gamma.est"]] = results[[clust]][["mix"]][["gamma.est"]]
      SimResults[[m]][[clust]][["mix"]][["tau.est"]] = results[[clust]][["mix"]][["tau.est"]]
      SimResults[[m]][[clust]][["mix"]][["p.est"]] = results[[clust]][["mix"]][["p.est"]]
      SimResults[[m]][[clust]][["l_tilde"]] = results[[clust]][["l_tilde"]]
      SimResults[[m]][[clust]][["alpha.est"]] = results[[clust]][["alpha.est"]]
      SimResults[[m]][[clust]][["beta.est"]] = results[[clust]][["beta.est"]]
      SimResults[[m]][[clust]][[cellTypes[1]]][["tau.hat"]] = results[[clust]][[cellTypes[1]]][["tau.hat"]]
      SimResults[[m]][[clust]][[cellTypes[2]]][["tau.hat"]] = results[[clust]][[cellTypes[2]]][["tau.hat"]]
      SimResults[[m]][[clust]][[cellTypes[3]]][["tau.hat"]] = results[[clust]][[cellTypes[3]]][["tau.hat"]]
      SimResults[[m]][[clust]][[cellTypes[1]]][["gamma.hat"]] = results[[clust]][[cellTypes[1]]][["gamma.hat"]]
      SimResults[[m]][[clust]][[cellTypes[2]]][["gamma.hat"]] = results[[clust]][[cellTypes[2]]][["gamma.hat"]]
      SimResults[[m]][[clust]][[cellTypes[3]]][["gamma.hat"]] = results[[clust]][[cellTypes[3]]][["gamma.hat"]]
      SimResults[[m]][[clust]][["CellType_Order"]] = results[[clust]][["CellType_Order"]]
      SimResults[[m]][[clust]][["WARN"]] = results[[clust]][["WARN"]]
      
      
    }
    
  }
  
  # Find probability estimates for mixtures
  ## Only interested in p_mat from this result, will concatenate results for all gene batches later
  SimSummary = Summarize_Report(SimResults, plots_options = plotsControl(plots = FALSE))
  SimSummary[["Time"]] = results[["Time"]]
  

save(SimResults, file = sprintf("%s/Mix_%s_GeneBatch_%s_Sim_Results.RData", prefix_results1, batch_label, as.character(clust_batch)))
save(SimSummary, file = sprintf("%s/Mix_%s_GeneBatch_%s_Sim_Summary.RData", prefix_results2, batch_label, as.character(clust_batch)))
# save(SimFull, file = sprintf("%s/Mix_%s_GeneBatch_%s_FullSimRes.RData", prefix_results4, batch_label, as.character(clust_batch)))

##########################################################################################################

print(gc())

q(save='no')

##########################################################################################################
