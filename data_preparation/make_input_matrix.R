
## this script takes pre-processed MSK IMPACT genomic and clinical data to
## create input matrices for MHN. additionally, a metadata table is created
## which records the chosen sample per eligible patient (as there sometimes
## is more than one), patient age at which the sample was sequenced (contains
## NAs) and estimated tumor mutational burden per sample. Estimated TMB,
## using the definition from Lengel et al. 2023, is the number of nonsynonymous
## mutations detected in a sample divided by the genomic extent covered by the
## MSK-IMPACT version used to assay the sample, in Mb. 

## this script was written in R 4.3.1

##### CONTENTS

################################################################################
### 01: DATA AND PACKAGE LOADING
### ============================================================================
### the following input files are required: per-sample metadata from GENIE13.1,
### the mutation file from GENIE13.1 and the panel details file from GENIE13.1
################################################################################

################################################################################
### 02: SAMPLE SELECTION
### ============================================================================
### for each patient that has a primary tumor sample annotated with the cancer
### type in question, if there are multiple primary tumor samples, choose one
### according to the following priorities:
### 1. avoid sampling time NAs
### 2. prefer samples at lower sampling age (initial diagnosis)
### 3. prefer newer version of IMPACT assay
### note that when running this step, there may be minor differences to the 
### original sample selection provided with the 2023 manuscript due to dataset 
### versions or negligible differences in the sample selection procedure
################################################################################

################################################################################
### 03: MUTATION EVENT CALLING
### ============================================================================
### firstly, to avoid including events that were not consistently assayed
### across the different versions of MSK-IMPACT, panel details are used to 
### pre-select genes that are assayed in all MSK-IMPACT versions.
### individual mutations are then filtered by binary functionality/pathogenicity
### annotation. then, any sample which has at least one functional variant
### in a gene of interest will get the respective event assigned.
################################################################################

################################################################################
### 04: EVENT SELECTION 
### ============================================================================
### mutations across samples are subsetted for the x most commonly mutated 
### genes and output files are saved
################################################################################



##### 01: DATA AND PACKAGE LOADING
################################################################################

# clear wd
rm(list=ls())

# GENIE metadata
sampleData <- read.delim("data_preparation/data_clinical_sample.txt", comment.char="#")

# load annotated mutation file
muts <- read.delim("data_preparation/data_mutations_annotated_eventFilter.txt")

# load panel details for gene lists and TMB calculation
panelProbes <- read.delim("data_preparation/genomic_information.txt")

################################################################################


##### 02: SAMPLE SELECTION
################################################################################

# set a new name for the dataset to be produced
projName <- "LUAD_n12"

# filter MSK samples by primary tumors of cancer type x
desiredSamples <- sampleData[which(grepl("-MSK-", sampleData$SAMPLE_ID, fixed = T)), ] 
desiredSamples <- desiredSamples[which(desiredSamples$ONCOTREE_CODE %in% c("LUAD")), ]
desiredSamples <- desiredSamples[which(desiredSamples$SAMPLE_TYPE == "Primary"), ]

# convert seq age to numeric - this produces true NAs 
desiredSamples$AGE_AT_SEQ_REPORT <- as.numeric(desiredSamples$AGE_AT_SEQ_REPORT)

# init output
selectionRecord <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(selectionRecord) <- c("patientID", "primID", "nPrim", "ageAtSequencing")

# iterate over selected patients
for (uPat in unique(desiredSamples$PATIENT_ID)) {
  
  # get samples
  allPat <- desiredSamples[which(desiredSamples$PATIENT_ID == uPat), ]
  allPrim <- allPat[which(allPat$SAMPLE_TYPE == "Primary"), ]
  
  # initialise new output row 
  nPrim <- 0; primID <- NA; ageAtSequencing <- FALSE;
  
  # carry out sample selection rationale as described
  if (nrow(allPrim) > 0) {
    
    # sort samples according to the following priorities
    # 1. avoid sampling time NAs
    # 2. prefer samples at lower sampling age (initial diagnosis)
    # 3. prefer newer version of IMPACT assay
    allPrim <- allPrim[order(allPrim$SEQ_ASSAY_ID, decreasing = T, na.last = T), ]
    allPrim <- allPrim[order(allPrim$AGE_AT_SEQ_REPORT, decreasing = F, na.last = T), ]
    
    # choose best
    primID <- allPrim$SAMPLE_ID[1]; nPrim <- nrow(allPrim)
    ageAtSequencing <- sampleData[which(sampleData$SAMPLE_ID == primID), "AGE_AT_SEQ_REPORT"]
    
    
  }
  
  # add to output
  npr <- data.frame(matrix(data = c(uPat, primID, nPrim, ageAtSequencing), nrow = 1, ncol = ncol(selectionRecord)))
  colnames(npr) <- colnames(selectionRecord)
  
  selectionRecord <- rbind(selectionRecord, npr)
  
}

# convert to numeric
selectionRecord$nPrim <- as.numeric(selectionRecord$nPrim); selectionRecord$ageAtSequencing <- as.numeric(selectionRecord$ageAtSequencing)

## remove patients without valid sample (happens if samples are of unspecified type only)
selectionRecord <- selectionRecord[which(selectionRecord$nPrim > 0), ]

# get all included samples
desiredSamples <- sampleData[which(sampleData$SAMPLE_ID %in% selectionRecord$primID), ] 

# for all msk assays, sum up extent of all covered exonic regions and return in Mb units
coveredMb <- vector(length = 6); names(coveredMb) <- c("MSK-IMPACT-HEME-400", "MSK-IMPACT-HEME-468", "MSK-IMPACT341", "MSK-IMPACT410", "MSK-IMPACT468", "MSK-IMPACT505")
for (p in names(coveredMb)) {
  currPanel <- panelProbes[which(panelProbes$SEQ_ASSAY_ID == p), ]
  coveredMb[p] <- (sum(currPanel$End_Position - currPanel$Start_Position)) / 1000000
}


# calculate and add TMB
for (ssr in 1:nrow(selectionRecord)) {
  
  currSID <- selectionRecord[ssr, "primID"]
  currPanel <- sampleData[which(sampleData$SAMPLE_ID == currSID), "SEQ_ASSAY_ID"]
  selectionRecord[ssr, "GENIE_TMB"] <- nrow(muts[which(muts$Tumor_Sample_Barcode == currSID), ]) / coveredMb[currPanel]
  
}

################################################################################


##### 03: MUTATION EVENT CALLING
################################################################################

# first, check mutation validity: is gene measured consistently across arrays

# get all MSK panels
mskPanels <- panelProbes[which(grepl("MSK", panelProbes$SEQ_ASSAY_ID)), ]
table(mskPanels$SEQ_ASSAY_ID)

# record per gene how many positions are assayed in each panel
allGenes <- sort(unique(mskPanels$Hugo_Symbol))[which(sort(unique(mskPanels$Hugo_Symbol)) != "")]
allPanels <- sort(unique(mskPanels$SEQ_ASSAY_ID))

genePerPanel <- data.frame(matrix(nrow = length(allGenes), ncol = length(allPanels), data = NA))
rownames(genePerPanel) <- allGenes; colnames(genePerPanel) <- allPanels

for (p in colnames(genePerPanel)) {
  currPanel <- mskPanels[which(mskPanels$SEQ_ASSAY_ID == p), ]
  for (g in rownames(genePerPanel)) {
    currGene <- currPanel[which(currPanel$Hugo_Symbol == g), ]
    genePerPanel[g, p] <- nrow(currGene)
  }
}

# binarise info: is gene assayed at all in panel?
binGenePerPanel <- genePerPanel
# Iterate over the input table
for (i in 1:nrow(binGenePerPanel)) {
  for (j in 1:ncol(binGenePerPanel)) {
    # Set TRUE if the value is non-zero, FALSE if it's zero
    binGenePerPanel[i, j] <- genePerPanel[i, j] != 0
  }
}

# now, get those genes which are present in all IMPACT panels
validGenesPerPanel <- binGenePerPanel[which(rowSums(binGenePerPanel[, c("MSK-IMPACT341", "MSK-IMPACT410", "MSK-IMPACT468", "MSK-IMPACT505")]) == 4), ]
validGenes <- rownames(validGenesPerPanel)

# apply functionality filter and subset for selected samples
mutsInSamples <- muts[which(muts$eventFilter), ]
mutsInSamples <- mutsInSamples[which(mutsInSamples$Tumor_Sample_Barcode %in% desiredSamples$SAMPLE_ID), ]

# init output
mutEvents <- data.frame(matrix(nrow = length(desiredSamples$SAMPLE_ID), ncol = 0))
rownames(mutEvents) <- desiredSamples$SAMPLE_ID

# per gene present in the samples, create a name, check the samples that have
# a mutation and add to output
mutEvNames <- c()
for (uMut in validGenes) {
  
  newName <- paste0(uMut); mutEvNames <- c(mutEvNames, newName)
  currMuts <- mutsInSamples[which(mutsInSamples$Hugo_Symbol == uMut), ]
  currBinary <- vector(mode = "numeric", length = length(desiredSamples$SAMPLE_ID))
  names(currBinary) <- desiredSamples$SAMPLE_ID
  currBinary[unique(currMuts$Tumor_Sample_Barcode)] <- 1
  
  mutEvents <- cbind(mutEvents, currBinary)
  
}

colnames(mutEvents) <- mutEvNames

# sort cols by freq
mutEvents <- mutEvents[, order(colSums(mutEvents), decreasing = T)]

plot(colSums(mutEvents) / nrow(mutEvents))

################################################################################


##### 04: EVENT SELECTION 
################################################################################

# how many most commonly mutated genes to include?
mcMuts <- 12

# select those events
select <- mutEvents[, (0:mcMuts)]

# check per-sample event counts
table(rowSums(select))

# check that no empty events are in there
table(colSums(select))

# check frequencies
colSums(select) / nrow(select)

# save outputs
write.csv(select, paste0("data/", projName, ".csv"), row.names = F, quote = F)
# write.csv(selectionRecord, paste0("data_preparation/", projName, "_selectionRecord.csv"), row.names = F, quote = F)


################################################################################




