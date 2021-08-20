# Examination of fragment length distributions

library(plotrix)

fragFiles1 = list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Fragment_Lengths/EGAD00001002671/",
                        full.names = T)
length(fragFiles1)

# par(mfrow = c(3,3), mar = c(3,4,3,2))
# for(i in 1:9){
#   fDist = read.table(fragFiles1[i], as.is = T)
#   weighted.hist(x = fDist[which(fDist[,1]>1000),2], w = fDist[which(fDist[,1]>1000),1])
# }

df_fragLen = read.table(fragFiles1[1], as.is = T)
colnames(df_fragLen) = c("Freq","Len")
df_fragLen$Weight = df_fragLen$Freq / sum(df_fragLen$Freq)
df_fragLen$SampleID = rep(1, times = nrow(df_fragLen))

for(i in 2:50){
  fDist = read.table(fragFiles1[i], as.is = T)
  fDist[,3] = fDist[,1] / sum(fDist[,1])
  fDist[,4] = rep(i, times = nrow(fDist))
  colnames(fDist) = c("Freq","Len","Weight","SampleID")
  
  df_fragLen = rbind(df_fragLen, fDist)
}

library(ggplot2)
ggplot(data = df_fragLen, mapping = aes(x = Len, weight = Weight, group = SampleID)) + geom_density() +
  coord_cartesian(xlim = c(100, 3*10^4)) + # coord_cartesian(xlim = c(100, 3*10^4))
  ggtitle("FragLen Density Plots for 50 Samples") + xlab("Length of Fragment")

###########################################################################################################

# Combine multiple fragment length files together to create one overall fragment length distribution
# file to be used in the Blueprint model fit

set.seed(2020)
files_to_use = sample(fragFiles1, size = 50, replace = F)

df1 = read.table(files_to_use[1], as.is = T)
colnames(df1) = c("Freq","Len")
df2 = read.table(files_to_use[2], as.is = T)
colnames(df2) = c("Freq","Len")

df_merge = merge(df1, df2, by = "Len", all = T)

df_combo = data.frame(Len = df_merge$Len, Freq = rowSums(df_merge[,c("Freq.x","Freq.y")]))

for(i in 3:50){
  fDist = read.table(files_to_use[i], as.is = T)
  colnames(fDist) = c("Freq","Len")
  
  df_merge = merge(df_combo, fDist, by = "Len", all = T)
  df_merge$Freq.x = ifelse(is.na(df_merge$Freq.x), 0, df_merge$Freq.x)
  df_merge$Freq.y = ifelse(is.na(df_merge$Freq.y), 0, df_merge$Freq.y)
  
  df_combo = data.frame(Len = df_merge$Len, Freq = rowSums(df_merge[,c("Freq.x","Freq.y")]))
}

fragLenDist = data.frame(df_combo$Freq, df_combo$Len)
colnames(fragLenDist) = NULL

write.table(fragLenDist, 
            file = "C:/Users/hheiling/Documents/GitHub/deconvolution/Blueprint_Simulation/Materials/Combined_FragDist.txt",
            row.names = F, col.names = F)

head(read.table("C:/Users/hheiling/Documents/GitHub/deconvolution/Blueprint_Simulation/Materials/Combined_FragDist.txt"))


###########################################################################################################
