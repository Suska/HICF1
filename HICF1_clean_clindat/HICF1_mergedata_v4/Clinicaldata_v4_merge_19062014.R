# Dr. Susanne Weller, 19/06/2014

#IMPORT CLINICAL DATA
clindat = read.csv("clinicaldata_pad_14052014.txt", sep="\t", header=TRUE, strip.white = TRUE, na.strings=c("", "Early death"))

# IMPORT CLONAL DATA
clones = read.table("subclonesv2_GC.txt", sep="\t", header=TRUE, na.strings=c("Withdrew from follow-up collection", "Missing"))
colnames(clones$sample_id) <- clones$ID
names(clones)[1] <- "ID"
names(clones)[2] <- "clones"
names(clones)[7] <- "CNAs"
clones = subset(clones, select=c("ID", "clones", "CNAs"))

clinclone= merge(clindat, clones, by="ID")

#IMPORT GENETIC FACTORS
genfact = read.table("genfactor_pad_14052014.txt", sep="\t", header=TRUE)

#IMPORT BIRCH, ATM and BINET
genBIRCH = read.table("genetics_MRD_v2_for_model.txt", sep="\t", header=TRUE)
genBIRCH$ID <- genBIRCH$sample_id
genBIRCH <- genBIRCH[c(2:5, 40:41)]

#IMPORT BLOOD DATA
blood = read.table("HICF1_nows_blood_03062014.txt", sep="\t", header=TRUE, na.strings=c("Withdrew from follow-up collection", "Missing"))

# MERGE ALL DATA FRAMES
genclin <-merge(clinclone, genfact, by="ID")
genclinv2 <-merge(genclin, genBIRCH, by="ID")
genclinv3 <- merge(genclinv2, blood, by="ID")

save(genclinv3, file="genclinv3.Rda")
