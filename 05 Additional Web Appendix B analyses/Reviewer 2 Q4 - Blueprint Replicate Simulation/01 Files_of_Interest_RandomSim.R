
# Goals:
# Using the dat.rds data information, determine (a) which individuals had pure cell type
# samples for all three cell types and (b) Find the associated EGAF sample id values for these
# individuals. 
# Use the above information to organize the pure samples such that they can easily create mixture
# files composed of pure samples from a single individual and the files used for mixture creation
# and the files used for the algorithm fit are easily distinguishable.

# Purpose:
# For use in the 'random simulation' to address a reviewer comment

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

# Check that all IDs have three separate samples from each of the three cell types
table(blueC$DONOR_ID, blueC$DATASET_ID)

library(stringr)
# Save all DONOR_IDs that met the above criteria
## Create a data.frame that connectrs each DONOR_ID with the three EGAF codes (for each cell type)
blue1 = blueC[which(str_detect(blueC$DATASET_ID, "2671")), c("DONOR_ID","EGAF")]
blue2 = blueC[which(str_detect(blueC$DATASET_ID, "2674")), c("DONOR_ID","EGAF")]
blue3 = blueC[which(str_detect(blueC$DATASET_ID, "2675")), c("DONOR_ID","EGAF")]
blue_key = merge(blue1, blue2, by = "DONOR_ID")
blue_key = merge(blue_key, blue3, by = "DONOR_ID")
colnames(blue_key) = c("DONOR_ID","EGAD00001002671","EGAD00001002674","EGAD00001002675")
dim(blue_key)

save(blue_key, file = "Blueprint_Simulation/Materials/blue_key_randomsim.RData")


###################################################################################################

# The End