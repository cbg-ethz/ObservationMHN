
## this script uses MHN input matrices to assess whether events are enriched
## in samples with few or many events ("event load enrichment")

## this script was written in R 4.3.1 and requires the 
## "circlize" and "ComplexHeatmap" packages to be installed for plotting



##### DATA LOADING
################################################################################

# clear wd
rm(list = ls())

# load packages
library(circlize)
library(ComplexHeatmap)

# load a MHN input matrix
inputProject <- "LUAD_n12"
input <- read.csv(paste0("data/", inputProject, ".csv"))

################################################################################



##### EVENT LOAD ENRICHMENT ANALYSIS 
################################################################################

# get event proportions over all other events to get a 
# "probability vector" for drawing events
obsEventProps <- colSums(input) / sum(colSums(input))

# get possible event loads
eventLoadRange <- 0:max(rowSums(input))

# get event load proportions to get a "probability vector" for knowing how 
# many events to draw
obsLoadProps <- table(factor(rowSums(input), levels = eventLoadRange)) / nrow(input)

# for nPerm permutations, draw as many samples as present in the input,
# with samples having event loads that follow the original distribution
# and event frequencies also following th original distribution. then, record
# for every event x event load combination how many times it was observed

nPerm <- 10000

permResults <- data.frame(matrix(nrow = length(eventLoadRange), ncol = ncol(input), data = ""))
rownames(permResults) <- paste0("L", eventLoadRange); colnames(permResults) <- colnames(input)

prTemplate <- permResults

for (p in 1:nPerm) {
  
  # get empty copy of results
  currPR <- prTemplate
  
  # draw sample event loads
  currELs <- table(factor(sample(eventLoadRange, nrow(input), prob = obsLoadProps, replace = T), levels = eventLoadRange))
  
  # for each event load, draw as many events as would be present in total and count
  currPR["L0", ] <- 0
  for (l in 1:max(eventLoadRange)) {
    nEventsToDraw <- as.numeric(currELs[as.character(l)] * l)
    currLoadEvents <- table(factor(sample(names(obsEventProps), nEventsToDraw, replace = T, prob = obsEventProps), levels = names(obsEventProps)))
    currPR[paste0("L", l), ] <- as.numeric(currLoadEvents)
  }
  
  permResults[] <- paste(as.matrix(permResults), as.matrix(currPR), sep = ",")
  
}

# find means and SDs of the permutated distribution
permMeans <- prTemplate
permSDs <- prTemplate
for (r in 1:nrow(permResults)) {
  for (c in 1:ncol(permResults)) {
    currDist <- as.numeric(strsplit(permResults[r, c], ",")[[1]])
    permMeans[r, c] <- mean(currDist, na.rm = T)
    permSDs[r, c] <- sd(currDist, na.rm = T)
  }
}

# get observed distribution
obsEventCountsPerLoad <- prTemplate

# for each event load, draw as many events as would be present in total and count
obsEventCountsPerLoad["L0", ] <- 0
for (l in 1:max(eventLoadRange)) {
  currLoadSamples <- input[which(rowSums(input) == l), ]
  currLoadEvents <- colSums(currLoadSamples)
  obsEventCountsPerLoad[paste0("L", l), ] <- as.numeric(currLoadEvents)
}

# convert to numeric
obsEventCountsPerLoad <- sapply(obsEventCountsPerLoad, as.numeric)
permMeans <- sapply(permMeans, as.numeric)
permSDs <- sapply(permSDs, as.numeric)

# find z-scores of observed distribution
zScoreEventsPerLoad <- (as.numeric(obsEventCountsPerLoad) - as.numeric(permMeans)) / permSDs
zScoreEventsPerLoad <- zScoreEventsPerLoad[2:nrow(zScoreEventsPerLoad), ]
rownames(zScoreEventsPerLoad) <- paste0("L", 1:nrow(zScoreEventsPerLoad))

# make color schemes
cfEnrichment = colorRamp2(c(floor(min(as.matrix(zScoreEventsPerLoad))), -2, 0, 2, ceiling(max(as.matrix(zScoreEventsPerLoad)))), c("#0072B2", "#D6F7FF", "white", "#FFECC2", "#E69F00"))

# make legends
lgdEnrichment = Legend(col_fun = cfEnrichment, title = "Deviation to random permutations in SD units", at = c(floor(min(as.matrix(zScoreEventsPerLoad))), 0,  ceiling(max(as.matrix(zScoreEventsPerLoad)))), direction = "vertical", legend_height = unit(5, "cm"), grid_width = unit(1.5, "cm"), border = "black", title_gp = gpar(cex = 1.5, fontface = "bold"))
packedLegend = packLegend(lgdEnrichment, direction = "vertical", column_gap = unit(2, "cm"))

# set up heatmap to plot
HM <- Heatmap(as.matrix(zScoreEventsPerLoad), col = cfEnrichment, name = paste0("Event load enrichments"),
              cluster_rows = F, cluster_columns = F,
              column_dend_reorder = T,
              column_dend_side = "bottom",
              show_heatmap_legend = F,
              rect_gp = gpar(col = "lightgrey", lwd = 2),
              row_names_side = "left",
              row_names_gp = gpar(fontface = "bold", cex = 1.25),
              column_names_side = "top",
              column_names_rot = 45,
              column_names_gp = gpar(fontface = "bold", cex = 1.25),
              row_gap = unit(4, "mm"), column_gap = unit(4, "mm"),
              border = c("black"),
              row_title = NULL,
              column_title = NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.1f", as.matrix(zScoreEventsPerLoad)[i, j]), x, y, gp = gpar(fontsize = 12))}
              )


if (FALSE) {
  # save as pdf
  pdf(paste0("data_preparation/", inputProject, "_enrichmentPlot.pdf"), w = 6, h = 4)
  
  draw(HM, padding = unit(c(0.5, 0.5, 0.5, 1.2), "cm"))
  # draw(packedLegend, x = unit(14, "cm"), y = unit(17, "cm"), just = c("right", "bottom"))
  
  dev.off()
}

################################################################################

