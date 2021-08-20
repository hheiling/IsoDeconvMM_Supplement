
# Goals:
# Using the dat.rds data information, determine (a) which individuals had pure cell type
# samples for all three cell types and (b) Find the associated EGAF sample id values for these
# individuals. 
# Use the above information to organize the pure samples such that they can easily create mixture
# files composed of pure samples from a single individual and the files used for mixture creation
# and the files used for the algorithm fit are easily distinguishable.

blue_dat0 = readRDS("C:/Users/hheiling/Documents/GitHub/deconvolution/Blueprint_Simulation/data/dat.rds")
dim(blue_dat0) # 614 = 212 + 197 + 205

# First, filter out any files that are specified in the blue_dat0 data.frame but are not present
# in the available samples.

CT1_avail = str_sub(list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Pure_ExonSetCts/EGAD00001002671/"), start = 1, end = 15)
CT2_avail = str_sub(list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Pure_ExonSetCts/EGAD00001002674/"), start = 1, end = 15)
CT3_avail = str_sub(list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Pure_ExonSetCts/EGAD00001002675/"), start = 1, end = 15)

length(CT1_avail)
length(CT2_avail)
length(CT3_avail)

blue_dat = blue_dat0[which(blue_dat0$EGAF %in% union(CT1_avail, union(CT2_avail, CT3_avail))),]
dim(blue_dat) # 573 = 197 + 176 + 200

table(blue_dat$LIBRARY_LAYOUT, blue_dat$DATASET_ID)

# Restrict individuals to those who have samples from all three cell types

table(table(blue_dat$DONOR_ID))

donor_tab = table(blue_dat$DONOR_ID)
donors3 = names(donor_tab[which(donor_tab == 3)])

blueA = blue_dat[which(blue_dat$DONOR_ID %in% donors3),]
dim(blueA)
length(unique(blueA$DONOR_ID))

# Make sure that those in CT EGAD00001002671 all use paired-end reads and all other 
# cell types use single-end reads

table(blueA$LIBRARY_LAYOUT, blueA$DATASET_ID)

lib_lay = blueA$LIBRARY_LAYOUT
datID = blueA$DATASET_ID

rows1 = which(lib_lay == "PAIRED" & datID == "EGAD00001002671")
rows2 = which(lib_lay == "SINGLE" & datID == "EGAD00001002674")
rows3 = which(lib_lay == "SINGLE" & datID == "EGAD00001002675")
blueB = blueA[union(rows1, union(rows2, rows3)),]
dim(blueB)

# After filtering based on paired-end and single-end read requirements, double check
# that all donors have 3 samples
donor_tabB = table(blueB$DONOR_ID)
donors3B = names(donor_tabB[which(donor_tabB == 3)])
blueC = blueB[which(blueB$DONOR_ID %in% donors3B),]
length(unique(blueC$DONOR_ID))
dim(blueC)

# Pick out DONOR_IDs to use for mixture sample creation and pure sample fit
set.seed(3905)
donors_all = unique(blueC$DONOR_ID)
length(donors_all)
donors_mix = sample(donors_all, size = 100, replace = F)
donors_pure1 = donors_all[-which(donors_all %in% donors_mix)]
length(donors_mix)
length(donors_pure1)
blue_mix = blueC[which(blueC$DONOR_ID %in% donors_mix),]
blue_pure1 = blueC[which(blueC$DONOR_ID %in% donors_pure1),]
dim(blue_mix)
dim(blue_pure1)

# Create a data.frame that connects each DONOR_ID with the three EGAF codes (for each cell type)
library(stringr)

## files for mixture creation
blue_mix1 = blue_mix[which(str_detect(blue_mix$DATASET_ID, "2671")),c("DONOR_ID","EGAF")]
blue_mix2 = blue_mix[which(str_detect(blue_mix$DATASET_ID, "2674")),c("DONOR_ID","EGAF")]
blue_mix3 = blue_mix[which(str_detect(blue_mix$DATASET_ID, "2675")),c("DONOR_ID","EGAF")]
blue_mix_key = merge(blue_mix1, blue_mix2, by = "DONOR_ID")
blue_mix_key = merge(blue_mix_key, blue_mix3, by = "DONOR_ID")
colnames(blue_mix_key) = c("DONOR_ID","EGAD00001002671","EGAD00001002674","EGAD00001002675")
dim(blue_mix_key)
# save(blue_mix_key, file = "Blueprint_Simulation/Materials/blue_mix_key.RData")

## files for algorithm fit
## Note: In this case, not important to have samples across cell types come from same individuals
blue_pure1A = blue_pure1[which(str_detect(blue_pure1$DATASET_ID, "2671")), c("DONOR_ID","EGAF")]
blue_pure1B = blue_pure1[which(str_detect(blue_pure1$DATASET_ID, "2674")), c("DONOR_ID","EGAF")]
blue_pure1C = blue_pure1[which(str_detect(blue_pure1$DATASET_ID, "2675")), c("DONOR_ID","EGAF")]
blue_pure2 = merge(merge(blue_pure1A, blue_pure1B, by = "DONOR_ID"), blue_pure1C, by = "DONOR_ID")
colnames(blue_pure2) = c("DONOR_ID","EGAD00001002671","EGAD00001002674","EGAD00001002675")
dim(blue_pure2)
## Add some additional samples to this pure file list
## Start with full blue_dat data set (after filtering out files that are not actually present in data set)
## Remove samples that don't have desired LIBRARY_LAYOUT (based on DATASET_ID)
lib_lay = blue_dat$LIBRARY_LAYOUT
datID = blue_dat$DATASET_ID
rows1 = which(lib_lay == "PAIRED" & datID == "EGAD00001002671")
rows2 = which(lib_lay == "SINGLE" & datID == "EGAD00001002674")
rows3 = which(lib_lay == "SINGLE" & datID == "EGAD00001002675")
blueD = blue_dat[union(rows1, union(rows2, rows3)),]
dim(blueD) # Should match 599 (209 + 194 + 196)
## Ignore samples already associated with the 119 samples already used in either blue_mix_key
##    or blue_pure
blueD2 = blueD[-which(blueD$DONOR_ID %in% donors_all),]
dim(blueD2)
length(unique(blueD2$DONOR_ID))
length(intersect(unique(blueD2$DONOR_ID), donors_all))
## Select from remaining samples
CT1_extra = blueD2[which(str_detect(blueD2$DATASET_ID, "2671")),]
CT2_extra = blueD2[which(str_detect(blueD2$DATASET_ID, "2674")),]
CT3_extra = blueD2[which(str_detect(blueD2$DATASET_ID, "2675")),]
num_extra = 22
set.seed(5475)
pure_extra = cbind(str_c("other",1:num_extra), sample(CT1_extra$EGAF,num_extra), 
                   sample(CT2_extra$EGAF,num_extra), sample(CT3_extra$EGAF,num_extra))
colnames(pure_extra) = c("DONOR_ID","EGAD00001002671","EGAD00001002674","EGAD00001002675")
blue_pure = rbind(blue_pure2, pure_extra)
dim(blue_pure)
# save(blue_pure, file = "Blueprint_Simulation/Materials/blue_pure.RData")

###################################################################################################

# The End