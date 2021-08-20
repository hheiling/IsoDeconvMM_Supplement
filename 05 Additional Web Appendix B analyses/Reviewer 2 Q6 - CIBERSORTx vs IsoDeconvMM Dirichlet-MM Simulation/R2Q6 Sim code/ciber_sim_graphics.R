
# Create graphics to analyze the CIBERSORTx simulation results in folder
# 'CIBERSORTx Sim results'
## Scatterplot
## Correlations and SSE

library(stringr)
library(ggplot2)
library(reshape2)
library(knitr)

# True mixture probabilities
## p_combos matrix object
load("~/Longleaf/deconvolution/Simulated_Output/Ciber_Sim_MixProbCombos100.RData")

# Source helper functions
source("~/GitHub/deconvolution/Blueprint_Simulation/Graphics11/graphics_code.R")

path_results = "~/GitHub/deconvolution/Simulation_Ideas/Reviewer 2 Q6 Simulation/R2Q6 Sim Results/"
files = list.files(path = path_results, pattern = ".txt")
print(files)

################################################################################################
# Scatter plot - CIBERSORTx
################################################################################################

# Results using 100 clusters
probs_100 = read.table(sprintf("%s/CIBERSORTx_Sim_Mix100_Results.txt", path_results), header = T)

df_100 = p_dfC(p_results = probs_100, CT_cols = 2:4, p_combos = p_combos)
df_100$GeneNum = "100 Genes"

# Results using 50 clusters in signature matrix, 100 clusters in mixture matrix
probs_50 = read.table(sprintf("%s/CIBERSORTx_Sim_Mix050_Results.txt", path_results), header = T)

df_50 = p_dfC(p_results = probs_50, CT_cols = 2:4, p_combos = p_combos)
df_50$GeneNum = "50 Genes"

# Results using 25 clusters in signature matrix, 100 clusters in mixture matrix
probs_25 = read.table(sprintf("%s/CIBERSORTx_Sim_Mix025_Results.txt", path_results), header = T)
df_25 = p_dfC(p_results = probs_25, CT_cols = 2:4, p_combos = p_combos)
df_25$GeneNum = "25 Genes"

# Results using 10 clusters in signature matrix, 100 clusters in mixture matrix
probs_10 = read.table(sprintf("%s/CIBERSORTx_Sim_Mix010_Results.txt", path_results), header = T)
df_10 = p_dfC(p_results = probs_10, CT_cols = 2:4, p_combos = p_combos)
df_10$GeneNum = "10 Genes"

df_cx_all = rbind(df_100, df_50, df_25, df_10)

df_cx_all$GeneNum = factor(df_cx_all$GeneNum, 
                           levels = c("10 Genes","25 Genes","50 Genes","100 Genes"))

df_cx_all = df_cx_all[which(df_cx_all$GeneNum %in% c("10 Genes","100 Genes")),]


ggplot(data = df_cx_all) + geom_point(mapping = aes(x = True_p, y = p_est), shape = 1) +
  facet_grid(GeneNum ~ CellType) + geom_abline(slope = 1, intercept = 0) +
  xlab("True Proportion") + ylab("Proportion Estimate") +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0,1)) +
  ggtitle("CIBERSORTx Method Results") +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_text(angle = 270))
ggsave(filename = sprintf("%s/CIBERSORTx_Sim_Scatter.jpeg",path_results), # .eps
       width = 7, height = 5, units = "in")

x_ciber = corr_SSE(df_list = list(df_100, df_50, df_25, df_10), 
                   col_labels = c("100 Genes","50 Genes","25 Genes","10 Genes"), 
                   CT = c("CT1","CT2","CT3"))

################################################################################################
# Scatter plot - IsoDeconvMM
################################################################################################

# Results using 100 clusters
## load probs list results
load("~/Longleaf/deconvolution/Simulation_Results/R2Q6_IsoDeconvMM/Probs_Ref2Q6.RData")
probs_100 = probs

clust_names_100 = rownames(probs_100[[1]])
set.seed(2021)
clust_names_50 = sample(clust_names_100, 50, replace = F)
clust_names_25 = sample(clust_names_50, 25, replace = F)
clust_names_10 = sample(clust_names_25, 10, replace = F)

df_100 = p_dfA(pList = probs_100, mix_labels = mix_names, clust_names = clust_names_100,
                   p_combos = p_combos, out_remove = F)
df_100$GeneNum = "100 Genes"

df_50 = p_dfA(pList = probs_100, mix_labels = mix_names, clust_names = clust_names_50,
                  p_combos = p_combos, out_remove = F)
df_50$GeneNum = "50 Genes"

df_25 = p_dfA(pList = probs_100, mix_labels = mix_names, clust_names = clust_names_25,
                  p_combos = p_combos, out_remove = F)
df_25$GeneNum = "25 Genes"

df_10 = p_dfA(pList = probs_100, mix_labels = mix_names, clust_names = clust_names_10,
                  p_combos = p_combos, out_remove = F)
df_10$GeneNum = "10 Genes"

df_cx_all = rbind(df_100, df_50, df_25, df_10)

df_cx_all$GeneNum = factor(df_cx_all$GeneNum, 
                           levels = c("10 Genes","25 Genes","50 Genes","100 Genes"))

df_cx_all = df_cx_all[which(df_cx_all$GeneNum %in% c("10 Genes","100 Genes")),]


ggplot(data = df_cx_all) + geom_point(mapping = aes(x = True_p, y = p_est), shape = 1) +
  facet_grid(GeneNum ~ CellType) + geom_abline(slope = 1, intercept = 0) +
  xlab("True Proportion") + ylab("Proportion Estimate") +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0,1)) +
  ggtitle("IsoDeconvMM Method Results") +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_text(angle = 270))
ggsave(filename = sprintf("%s/CIBERSORTx_Sim_Scatter_IsoDeconvMM.jpeg",path_results), # .eps
       width = 7, height = 5, units = "in")

x_iso = corr_SSE(df_list = list(df_100, df_50, df_25, df_10), 
                   col_labels = c("100 Genes","50 Genes","25 Genes","10 Genes"), 
                   CT = c("CT1","CT2","CT3"))

################################################################################################
# Correlations and SSE
################################################################################################


df_corr_a = bar_df(table = x_ciber$corr, 
                   col_labels = c("100 Genes","50 Genes", "25 Genes",
                                  "10 Genes"), 
                   var_name = "Genes", value_name = "Correlation")
df_corr_a$FitType = "CIBERSORTx"
df_corr_b = bar_df(table = x_iso$corr, 
                   col_labels = c("100 Genes","50 Genes", "25 Genes",
                                  "10 Genes"), 
                   var_name = "Genes", value_name = "Correlation")
df_corr_b$FitType = "IsoDeconvMM"

df_corr = rbind(df_corr_a, df_corr_b)
df_corr = df_corr[which(df_corr$Genes %in% c("10 Genes","100 Genes")),]

ggplot(data = df_corr, aes(x = CellType, y = Correlation, fill = FitType)) +
  facet_grid(. ~ Genes) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("black","gray40","gray60","gray80")) +
  coord_cartesian(ylim = c(0,1)) + 
  theme_grey(base_size = 18) + theme(legend.position = "bottom") +
  labs(fill = "Total Genes")
ggsave(filename = sprintf("%s/CIBERSORTx_Sim_Corr.jpeg",path_results), # .eps
       width = 6, height = 4, units = "in")


df_SSE_a = bar_df(table = x_ciber$SSE, 
                   col_labels = c("100 Genes","50 Genes", "25 Genes",
                                  "10 Genes"), 
                   var_name = "Genes", value_name = "SSE")
df_SSE_a$FitType = "CIBERSORTx"
df_SSE_b = bar_df(table = x_iso$SSE, 
                   col_labels = c("100 Genes","50 Genes", "25 Genes",
                                  "10 Genes"), 
                   var_name = "Genes", value_name = "SSE")
df_SSE_b$FitType = "IsoDeconvMM"

df_SSE = rbind(df_SSE_a, df_SSE_b)
df_SSE = df_SSE[which(df_SSE$Genes %in% c("10 Genes","100 Genes")),]


ggplot(data = df_SSE, aes(x = CellType, y = SSE, fill = FitType)) +
  facet_grid(. ~ Genes) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("black","gray40","gray60","gray80")) +
  theme_grey(base_size = 18) + theme(legend.position = "bottom") +
  labs(fill = "Total Genes")
ggsave(filename = sprintf("%s/CIBERSORTx_Sim_SSE.jpeg",path_results), # .eps
       width = 6, height = 4, units = "in")

################################################################################################
################################################################################################
