# Purpose of Sim 11: 
## Redo Blueprint simulations with different cluster selection based on fold change values and 
## Wilcoxon rank sum p-values (see Parallel Selection section of Blueprint_IsoDu_Results.Rmd for 
## details) comparing a single cell type to the two other cell types collectively

# Sim11C: Run fit using 20 pure samples per cell type, truth as initial point


# Fit IsoDeconvMM model fit to simulated data

# Load IsoDeconvMM library:
# library(IsoDeconvMM)
library(stringr)
library(alabama)
library(ICSNP)
library(gtools)
library(MASS)

# Arrays 1-300
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
batch <- array_val

batch_samp = rep(1:100, each = 3)[batch]

if(batch_samp <= 9){
  batch_label = str_c("00",batch_samp)
}else if(batch_samp <= 99 ){
  batch_label = str_c("0",batch_samp)
}else{
  batch_label = as.character(batch_samp)
}

clust_batch = rep(1:3, times = 100)[batch]

# Number cell types
J = 3

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Blueprint materials output prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Pure cell type files prefix
prefix_pure = str_c(header_pine,"Fit_Samples2")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Mixture_Samples2")

# Simulation results prefix
prefix_results = str_c(header_nas,"Simulation_Results/Fit11C")
if(!dir.exists(prefix_results)){dir.create(prefix_results, recursive = T)}
prefix_results1 = str_c(prefix_results,"/Results_Abbrev")
if(!dir.exists(prefix_results1)){dir.create(prefix_results1, recursive = T)}
prefix_results2 = str_c(prefix_results,"/Summary_Probs")
if(!dir.exists(prefix_results2)){dir.create(prefix_results2, recursive = T)}
prefix_results3 = str_c(prefix_results,"/EffLen")
if(!dir.exists(prefix_results3)){dir.create(prefix_results3, recursive = T)}
prefix_results4 = str_c(prefix_results,"/Full_Results")
if(!dir.exists(prefix_results4)){dir.create(prefix_results4, recursive = T)}

# Source IsoDeconvMM functions
IsoDeconv_code = list.files(path = "/nas/longleaf/home/hheiling/deconvolution/code/IsoDeconvMM", 
                            full.names = T)
for(i in 1:length(IsoDeconv_code)){
  source(IsoDeconv_code[i])
}

# Find genes with differential isoform usage
## load nTE_discrimE data.frame (from Blueprint_IsoDu_Results.Rmd)
load(sprintf("%s/gencode.v15.nTE.discrimE.RData", prefix_mat))
length(unique(nTE_discrimE$clustID)) # 130
discrim_clusts = unique(nTE_discrimE$clustID)

# Since prob estimates occur separately for each gene, can subset the number of genes
# to run at one time (here, run 1/3 at a time). Will concatenate all p estimates later.
# Pick 1/3 of the clusters to run at a time
select_idx_first = (clust_batch-1)*45+1
if(clust_batch == 1){ # clusters 1-45
  (select_idx = (select_idx_first):45)
}else if(clust_batch == 2){ # clusters 46-90
  (select_idx = (select_idx_first):90)
}else if(clust_batch == 3){ # clusters 91-130
  (select_idx = (select_idx_first):length(discrim_clusts))
}

clusts_to_use = discrim_clusts[select_idx]

# Datasets for cell types
CT = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")
## Assume ordering of datasets in CT correspond to cell types CT1, CT2, and CT3
cellTypes = c("CT1","CT2","CT3")

# Find pure reference files
CT1_files = list.files(path = str_c(prefix_pure, CT[1], sep = "/"), full.names = T)
CT2_files = list.files(path = str_c(prefix_pure, CT[2], sep = "/"), full.names = T)
CT3_files = list.files(path = str_c(prefix_pure, CT[3], sep = "/"), full.names = T)

set.seed(8500)
CT_ref_files = c(sample(CT1_files, size = 20, replace = F),
                 sample(CT2_files, size = 20, replace = F),
                 sample(CT3_files, size = 20, replace = F))

pure_ref_files_all = cbind(CT_ref_files, c(rep("CT1",times=20),rep("CT2",times=20),rep("CT3",times=20)))

# Select all 20 pure samples from each cell type
pure_ref_files = pure_ref_files_all

cat("dim pure_ref_files: ",dim(pure_ref_files),"\n")

# Find fragment length files
## Assume distributions from paired-end fragment lengths, even though not true for CT2 and CT3
fraglens_file = str_c(prefix_mat, "Combined_FragDist.txt", sep="/")

print(basename(fraglens_file))

# Load probability combinations: p_combos
load(file = sprintf("%s/Blueprint_ProbCombos.RData", prefix_mat))

cat("Num prob combos: ", nrow(p_combos), "\n")

colnames(p_combos) = c("CT1","CT2","CT3")

# initPts:
# Only the truth
initPts = matrix(p_combos[batch_samp,], nrow = 1)
colnames(initPts) = c("CT1","CT2","CT3")
print(initPts)

# Run simulations

SimResults = list()
SimSummary = list()
EffLen_Mat = list()
SimFull = list()

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
SimFull = SimResults_Full

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
    
    EffLen_Mat[[m]][[clust]][["X"]] = results[[clust]][["X"]]
    EffLen_Mat[[m]][[clust]][["X.prime"]] = results[[clust]][["X.prime"]]
    
  }
  
}

# Find probability estimates for mixtures
## Only interested in p_mat from this result, will concatenate results for all gene batches later
SimSummary = Summarize_Report(SimResults, plots_options = plotsControl(plots = FALSE))


save(SimResults, file = sprintf("%s/Mix_%s_GeneBatch_%s_Sim_Results.RData", prefix_results1, batch_label, as.character(clust_batch)))
save(SimSummary, file = sprintf("%s/Mix_%s_GeneBatch_%s_Sim_Summary.RData", prefix_results2, batch_label, as.character(clust_batch)))

##########################################################################################################

print(gc())

q(save='no')

##########################################################################################################

