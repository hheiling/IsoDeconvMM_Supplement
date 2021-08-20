# Use Blueprint alpha values from 11A series in Dirichlet-NB model during simulation
# In this simulation, fit using genes with differential isoform usage

# Note: Mixture files created in "sim1_fit2_mixtures.R"

# Fit IsoDeconvMM model fit to simulated data

library(stringr)
library(alabama)
library(ICSNP)
library(gtools)
library(MASS)

# Arrays 1-250
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
batch <- array_val

batch_samp = rep(1:50, each = 5)[batch]

if(batch_samp %in% 1:9){
  batch_label = str_c("0",batch_samp)
}else{
  batch_label = as.character(batch_samp)
}

gene_batch = rep(1:5, times = 50)[batch]

# Number cell types
J = 3

header_nas = "/nas/longleaf/home/hheiling/"
header_pine = "/pine/scr/h/h/hheiling/"
# Simulated output prefix
prefix_sim_out = str_c(header_nas,"deconvolution/Simulated_Output")
# Human materials output prefix
prefix_HM = str_c(header_nas,"deconvolution/Human_Materials")
# Sourced code location prefix
prefix_code = str_c(header_nas,"deconvolution/code")
# Simulated pure cell type files prefix
prefix_pure = str_c(header_pine,"Simulated_Files/Sim1_Fit2_Pure_Sample_Counts")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Simulated_Files/Sim1_Fit2_Mixture_Sample_Counts")

# Simulation results prefix
prefix_results = str_c(header_nas,"deconvolution/Simulation_Results/Sim1Fit2")
if(!dir.exists(prefix_results)){dir.create(prefix_results, recursive = T)}
prefix_results1 = str_c(prefix_results,"/Results_Abbrev")
if(!dir.exists(prefix_results1)){dir.create(prefix_results1, recursive = T)}
prefix_results2 = str_c(prefix_results,"/Summary_Probs")
if(!dir.exists(prefix_results2)){dir.create(prefix_results2, recursive = T)}

# Source IsoDeconvMM functions
IsoDeconv_code = list.files(path = sprintf("%s/IsoDeconvMM", prefix_code), full.names = T)
for(i in 1:length(IsoDeconv_code)){
  source(IsoDeconv_code[i])
}

# load genes_df data.frame
load(sprintf("%s/genes_simulated_w_diffExp_diffUsg.RData", prefix_sim_out))
# Find 100 genes with differential isoform usage
genesDiffUsg = genes_df$diffUsg

# Since prob estimates occur separately for each gene, can subset the number of genes
# to run at one time. Will concatenate all p estimates later.
if(gene_batch == 1){
  genes_to_use = genesDiffUsg[1:20]
}else if(gene_batch == 2){
  genes_to_use = genesDiffUsg[21:40] 
}else if(gene_batch == 3){
  genes_to_use = genesDiffUsg[41:60]
}else if(gene_batch == 4){
  genes_to_use = genesDiffUsg[61:80] 
}else if(gene_batch == 5){
  genes_to_use = genesDiffUsg[81:100] 
}

# Find pure reference samples
pure_names_full = list.files(path = prefix_pure,
                             pattern = "counts.txt$", full.names = TRUE)

cat("Number pure files found: ", length(pure_names_full), "\n")

n = 5

CT1_files = pure_names_full[str_detect(pure_names_full, "CT1_")]
CT2_files = pure_names_full[str_detect(pure_names_full, "CT2_")]
CT3_files = pure_names_full[str_detect(pure_names_full, "CT3_")]

# Use replicates 1 through 5 for IsoDeconvMM fit algorithm
CT1_files_ref = CT1_files[1:n]
CT2_files_ref = CT2_files[1:n]
CT3_files_ref = CT3_files[1:n]

CT_ref_files = c(CT1_files_ref, CT2_files_ref, CT3_files_ref)

# Check correct selection of files 01 - 05
print(basename(CT_ref_files))

pure_ref_files = cbind(CT_ref_files, c(rep("CT1",times=n),rep("CT2",times=n),rep("CT3",times=n)))

cat("dim pure_ref_files: ",dim(pure_ref_files),"\n")

# Specify file for fragment length distribution
## For this simulation, use one fragment length distribution file for each mixture file
fraglens_file = list.files(path = prefix_sim_out, pattern = "fraglens", full.names = T)
print(basename(fraglens_file))

# Load probability combinations: p_combos
load(file = sprintf("%s/Sim1_Fit2_ProbCombos.RData", prefix_sim_out))

cat("Num prob combos: ",nrow(p_combos), "\n")

colnames(p_combos) = c("CT1","CT2","CT3")
rownames(p_combos) = str_c("pc_", c(str_c("0",1:9),10:nrow(p_combos)))

pc_labels = rownames(p_combos)

# initPts: the 10 generic initial points used in the Blueprint analysis

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

colnames(initPts) = colnames(p_combos)
print(initPts)

# Create mixture samples
mix_label = str_c("Mix_ProbCombo_",batch_label,"_counts")

print(mix_label)

# Run simulations

SimResults = list()
SimSummary = list()

cellTypes = c("CT1","CT2","CT3")
  
  pc_nums = batch_samp
  
  # Probability combination label
  pc = pc_labels[pc_nums]
  
  mix_files = list.files(path = prefix_mix, pattern = mix_label,
                         full.names = T)
  
  mix_names = str_remove(mix_label, "_counts")
  
  # Fit IsoDeconvMM algorithm
  SimResults_Full = IsoDeconvMM(directory = NULL, mix_files = mix_files,
                                pure_ref_files = pure_ref_files,
                                fraglens_files = fraglens_file,
                                bedFile = sprintf("%s/Homo_sapiens.GRCh37.66.nonoverlap.exon.bed", prefix_HM),
                                knownIsoforms = sprintf("%s/Homo_sapiens.GRCh37.66.nonoverlap.exon.knownIsoforms.RData", prefix_HM),
                                discrim_genes = genes_to_use,
                                readLen = 75, lmax = 600, eLenMin = 1,
                                mix_names = mix_names, initPts = initPts,
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
      
    }
    
  }
  
  # Find probability estimates for mixtures
  ## Only interested in p_mat from this result, will concatenate results for all gene batches later
  SimSummary = Summarize_Report(SimResults, plots_options = plotsControl(plots = FALSE))
  

save(SimResults, file = sprintf("%s/Mix_%s_GeneBatch_%s_Sim_Results.RData", prefix_results1, batch_label, as.character(gene_batch)))
save(SimSummary, file = sprintf("%s/Mix_%s_GeneBatch_%s_Sim_Summary.RData", prefix_results2, batch_label, as.character(gene_batch)))

##########################################################################################################

print(gc())

q(save = 'no')

##########################################################################################################
