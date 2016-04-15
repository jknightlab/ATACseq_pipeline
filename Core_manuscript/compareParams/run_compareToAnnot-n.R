##### Comparing the Results of Peak Calling Algorithms With Annotations
### JHendry, 2016/03/10
##
##
##
## Built from compareENCODE-v8-byPeak.R
##
## Input Args:
## $1 --> <outputDir>, containing "intersectBeds" & "sampleBeds"
## $2 --> <annotDir>, containing annotation .bed's & key file
## $3 --> <runStyle>, either "-r" or "-d"

args <- commandArgs(T)
outputDir <- args[1]
annotDir <- args[2]
runStyle <- args[3]


### Source Basci ATAC Analysis Functions
source("run_basicATACroutines.R")

### Load Intersect, Sample and Annotation BEDs
# Intersect BEDs
intrsctDir <- file.path(outputDir, "intersectBeds")
intrsctBeds <- list.files(intrsctDir, ".bed$", full.names = T)
intrsctBedOrd <- do.call(c, mapply(function(sampTypes) { grep(sampTypes, basename(intrsctBeds)) }, sampTypes = c("fresh", "frozen", "fixed3d", "fixed7d"), SIMPLIFY = F))

if (length(intrsctBedOrd) < length(intrsctBeds)) { # i.e., the samples aren't all fresh/frozen/fixed3d/fixed7d
	# Load files into a list WITHOUT ordering
	intrsctList <- lapply(intrsctBeds, function(x) { 
			if (file.info(x)$size != 0) { read.table(x) } 
			else { data.frame(chr = 0, start = 0, end = 0, samp = 0, sampChr = 0, sampStart = 0, sampEnd = 0, annot = 0, annotChr = 0, annotStart = 0, annotEnd = 0) }
			})
	names(intrsctList) <- sapply(basename(intrsctBeds), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })	
	intrsctList <- lapply(intrsctList, function(X) { if (!is.null(X)) { colnames(X) <- c("chr", "start", "end", "samp", "sampChr", "sampStart", "sampEnd", "annot", "annotChr", "annotStart", "annotEnd") ; X } })
	
	# This is defined for later when calling color scheme
	allSamps <- F

	} else {
	# Order the files first
	intrsctBeds <- intrsctBeds[intrsctBedOrd]
	
	# Load Ordered Intersect Files Into List
	intrsctList <- lapply(intrsctBeds, function(x) { 
			if (file.info(x)$size != 0) { read.table(x) } 
			else { data.frame(chr = 0, start = 0, end = 0, samp = 0, sampChr = 0, sampStart = 0, sampEnd = 0, annot = 0, annotChr = 0, annotStart = 0, annotEnd = 0) }
			})
	names(intrsctList) <- sapply(basename(intrsctBeds), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })	
	intrsctList <- lapply(intrsctList, function(X) { if (!is.null(X)) { colnames(X) <- c("chr", "start", "end", "samp", "sampChr", "sampStart", "sampEnd", "annot", "annotChr", "annotStart", "annotEnd") ; X } })
	
	# This is defined for later when calling color scheme
	allSamps <- T
}

# Sample BEDs
sampDir <- file.path(outputDir, "sampleBeds")
sampBeds <- list.files(sampDir, ".bed$", full.names = T)
sampList <- lapply(sampBeds, read.table)
names(sampList) <- sapply(basename(sampBeds), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })	
sampList <- lapply(sampList, function(X) { colnames(X) <- c("chr", "start", "end") ; X })

# Annotation BEDs
annotBeds <- list.files(annotDir, ".bed$", full.names = T)
annotList <- lapply(annotBeds, read.table)
names(annotList) <- sapply(basename(annotBeds), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })	
annotList <- lapply(annotList, function(X) { colnames(X) <- c("chr", "start", "end") ; X })
# Switch Names via annotKey
annotKey <- read.table(file.path(annotDir, "annotKey.txt"))
names(annotList) <- annotKey[match(names(annotList), annotKey[ , 1]), 2]


### Functions
# <filterIntrsct>
# this function is producing NA outputs for BED files with NA rows
# why do some BED files have NA rows?
filterIntrsctBed <- function(	intrsctBed, 
								sampFracThresh = 0.5, annotFracThresh = 0.5, numBp = 0, 
								uniqSamp = F, uniqAnnot = F) {
	# (1) Check For Approriate Formatting of intrsctBed
	colNames <- c("chr", "start", "end", "samp", "sampChr", "sampStart", "sampEnd", "annot", "annotChr", "annotStart", "annotEnd")
	if (any(!names(intrsctBed) %in% colNames)) { stop("Incorrect column names.") }
	
	# (2) Calculate Sizes	
	intrsctBed <- within(intrsctBed, {
		intrsctSize <- end - start
		sampSize <- sampEnd - sampStart
		annotSize <- annotEnd - annotStart
		sampFrac <- intrsctSize/sampSize
		annotFrac <- intrsctSize/annotSize
	})
	
	# Here, I need to account for the strange phenomenon of intersects of size 0
	# - this causes NaNs in sampFrac &/or annotFrac
	intrsctBed <- intrsctBed[with(intrsctBed, intrsctSize > 0 & sampSize > 0 & annotSize > 0), ]
	
	# (3) Implement Filter -- causes NAs
	intrsctBedFilt <- intrsctBed[with(intrsctBed, intrsctSize >= numBp & sampFrac >= sampFracThresh & annotFrac >= annotFracThresh), ]
	
	# (4) Implement Unique Peaks Only
	#	- Note, there are issues with removing by unique sample/annotID
	#	- You end up with ONE of the two or more intersections
	#	- !duplicated(sampID) will remove the second entry if there is more than one
	#	- no underlying biological significance to this, the more 3' entry is removed
	#	- will effect total numBps calculated in intersection (see summarizeBed)
	#	---> perhaps, to keep things sensical, if this option is choosen, automatically should keep
	#	only annotation or sample peaks instead of the intersection.
	if (uniqSamp) {
		intrsctBedFilt <- within(intrsctBedFilt, { sampID <- paste(sampChr, sampStart, sampEnd, sep = "-") })
		intrsctBedFilt <- intrsctBedFilt[with(intrsctBedFilt, !duplicated(sampID)), ]
		intrsctBedFilt <- within(intrsctBedFilt, { chr <- sampChr ; start <- sampStart ; end <- sampEnd })
	} 
	
	if (uniqAnnot) {
		intrsctBedFilt <- within(intrsctBedFilt, { annotID <- paste(annotChr, annotStart, annotEnd, sep = "-") })
	 	intrsctBedFilt <- intrsctBedFilt[with(intrsctBedFilt, !duplicated(annotID)), ]
		intrsctBedFilt <- within(intrsctBedFilt, { chr <- annotChr ; start <- annotStart ; end <- annotEnd })
	}
	intrsctBedFilt
}
 
# <normStatBy>
#  NB: output vector matches order of dataDf
normStatBy <- function(dataDf, normDf, groupCol, stat) {
	normVecMatches <- match(dataDf[ , groupCol], normDf[ , groupCol])
	normVec <- dataDf[ , stat]/normDf[normVecMatches, stat]		
	normVec
}

# <calcFScore>
calcFscore <- function(Sensitivity, Specificity, beta) {
	(1 + beta**2)*(Sensitivity*Specificity)/((Specificity*beta**2) + Sensitivity)
}


### Convert List of .beds Into Data Frames Summarizing .bed Statistics
# (A) intrsctSum - Unfiltered, used to calculate bpSpec/Sens
#		- Note, despite extra columns, 
#		- summarizeBed only works with chr, start, and end
#		- which describe the intersection
intrsctSum <- do.call(rbind, lapply(intrsctList, summarizeBed))

# (B) sampIntrsctSum - Filtered, used to calculate peak specificity (pkSpec)
#	- list of the SAMPLE peaks that satisfy the intersect conditions (filter conditions)
#	- this list is UNIQUE, i.e., sample peaks that intersect with more than one annotation
#	- are only listed one time.
sampIntrsctList <- lapply(intrsctList, function(X) { filterIntrsctBed(X, uniqSamp = F, sampFracThresh = 0.5, annotFracThresh = 0, numBp = 1) })
sampIntrsctSum <- do.call(rbind, lapply(sampIntrsctList, summarizeBed))

# (C) annotIntrsctSum - Filtered, used to calculate peak sensitivity (pkSens)
annotIntrsctList <- lapply(intrsctList, function(X) { filterIntrsctBed(X, uniqAnnot = F, sampFracThresh = 0.5, annotFracThresh = 0, numBp = 1) })
annotIntrsctSum <- do.call(rbind, lapply(annotIntrsctList, summarizeBed))

# (D) sampSum and annotSum - Unfiltered
#	- list of ALL sample and annotation peaks
#	- NOT just those involved in intersection
#	- Used to normalize when calculating bp and pk Spec/Sens
sampSum <- do.call(rbind, lapply(sampList, summarizeBed))
annotSum <- do.call(rbind, lapply(annotList, summarizeBed))

## Output of do.call(rbind, ) is a matrix
# - Want to convert to a data frame so can include "character" columns
# - reformatToDf takes rownames() of matrix, and makes them column of dataframe
intrsctSumDf <- reformatToDf(intrsctSum)
sampIntrsctSumDf <- reformatToDf(sampIntrsctSum)
annotIntrsctSumDf <- reformatToDf(annotIntrsctSum)
sampSumDf <- reformatToDf(sampSum, "samp")
annotSumDf <- reformatToDf(annotSum, "annot")


### If running a range of algorithm paramters...
# - extract those paramters and add them as columns to the data frame
# - allows for making parameter matrices later
if (runStyle == "-r") {
	if (length(grep("fseq", sampBeds)) == length(sampBeds)) {
		cat("All input files contain 'fseq' in name.\n")
		cat("Generating parameter matrices using F-seq parameters:\n")
		cat("Parameter 1: 'l', Parameter 2: 't'\n")
		intrsctSumDf <- getParams(intrsctSumDf, paramCol = "name")
		sampIntrsctSumDf <- getParams(sampIntrsctSumDf, paramCol = "name")
		annotIntrsctSumDf <- getParams(annotIntrsctSumDf, paramCol = "name")	
		algo <- "fseq"
	} else if (length(grep("macs", sampBeds)) == length(sampBeds)) {
		cat("All input files contain 'macs' in name.\n")
		cat("Generating parameter matrices using MACS2 parameters.\n")
		cat("Parameter 1: 'q', Parameter 2: 'e'\n")
		intrsctSumDf <- getParams(intrsctSumDf, paramCol = "name", paramOne = "q", paramTwo = "e")
		sampIntrsctSumDf <- getParams(sampIntrsctSumDf, paramCol = "name", paramOne = "q", paramTwo = "e")
		annotIntrsctSumDf <- getParams(annotIntrsctSumDf, paramCol = "name", paramOne = "q", paramTwo = "e")
		algo <- "macs"
	}

} 


### For Intersect .bed, generate individual columns for sample and annotation file intersected
#	- These are used when calculating specificity and sensitivity because...
#	- Sensitivity are calculated on a per-annotation basis
#	- Therefore need to know WHAT annotation was intersected
#	- Specificity is calculated on a per-sample basis
#	- Therefore need to know WHAT sample was intersected
intrsctSumDf <- within(intrsctSumDf, {
			annot <- sapply(strsplit(name, split = "-"), function(x) { x[2] })
			samp <- sapply(strsplit(name, split = "-"), function(x) { x[1] })
			rm(name)
	})
sampIntrsctSumDf <- within(sampIntrsctSumDf, {
			annot <- sapply(strsplit(name, split = "-"), function(x) { x[2] })
			samp <- sapply(strsplit(name, split = "-"), function(x) { x[1] })
			rm(name)
	})
annotIntrsctSumDf <- within(annotIntrsctSumDf, {
			annot <- sapply(strsplit(name, split = "-"), function(x) { x[2] })
			samp <- sapply(strsplit(name, split = "-"), function(x) { x[1] })
			rm(name)
	})


### Recall, intrsctSumDf contains information about the COMPLETE set of intersections
# - not filtered in any way
# - PER intersection, means sample and annotation peaks may be duplicated
# - We will appened all our sensitivity/specificity statistics to this data frame.

# (1) Order it in a a nice way
intrsctSumDf <- intrsctSumDf[with(intrsctSumDf, order(annot)), c(seq(ncol(intrsctSumDf) - 1, ncol(intrsctSumDf)), seq(1, ncol(intrsctSumDf) - 2))]
sampIntrsctSumDf <- sampIntrsctSumDf[with(sampIntrsctSumDf, order(annot)), c(seq(ncol(sampIntrsctSumDf) - 1, ncol(sampIntrsctSumDf)), seq(1, ncol(sampIntrsctSumDf) - 2))]
annotIntrsctSumDf <- annotIntrsctSumDf[with(annotIntrsctSumDf, order(annot)), c(seq(ncol(annotIntrsctSumDf) - 1, ncol(annotIntrsctSumDf)), seq(1, ncol(annotIntrsctSumDf) - 2))]

# (2) By base pair specificity and sensitivity (bpSens, bpSpec)
#	- the numBps is extract from intrsctSumDf, since it contains all unfiltered intersections
#	- normalize to totals from annotSumDf and sampSumDf, which contain all annotation and sample peaks
intrsctSumDf[ , "bpSens"] <- normStatBy(intrsctSumDf, annotSumDf, "annot", "numBps")
intrsctSumDf[ , "bpSpec"] <- normStatBy(intrsctSumDf, sampSumDf, "samp", "numBps")


# (3) By peak specificity and sensitivity (pkSens, pkSpec)
#	- in this case, a peak is considered intersected if it meets certain criterion
#	- criterion are set by intrsctBedFilter
#		* sampFrac > 0.5 , annotFrac > 0, numBps > 0
#	- i.e., at least half of the sample peak has to overlap with the annotation
#	- Intuitively good because if more of the peak is NOT annotation than IS annotation, we say it is not matching the annotation
#	- If the annotation is less than half the size of the sample peak, this is impossible
#		* good because prevents making sample peaks too large
#	- For this, have to be very careful about normalization due to duplicate peaks possible
# (3A) Peak Sensitivity
#	- Take all UNIQUE annotation peaks that met filter conditions (annotIntrsctSumDf)
#	- Normalize to ALL annotation peaks (annotSumDf)
#	- Count the NUMBER of peaks, not number of base pairs (numPks)
intrsctSumDf[ , "pkSens"] <- normStatBy(annotIntrsctSumDf, annotSumDf, "annot", "numPks")

# (3B) Peak Specificity
#	- Take all UNIQUE sample peaks that met filter conditions (sampIntrsctSumDf)
#	- Normalize to ALL sample peaks (sampSumDf)
#	- Count NUMBER of peaks, not number of base pairs (numPks)
intrsctSumDf[ , "pkSpec"] <- normStatBy(sampIntrsctSumDf, sampSumDf, "samp", "numPks")

# (4) Calculate F-scores
intrsctSumDf[  , "bpF"] <- with(intrsctSumDf, calcFscore(bpSens, bpSpec, beta = 1))
intrsctSumDf[  , "bpF2"] <- with(intrsctSumDf, calcFscore(bpSens, bpSpec, beta = 2))
intrsctSumDf[  , "bpFp5"] <- with(intrsctSumDf, calcFscore(bpSens, bpSpec, beta = 0.5))
intrsctSumDf[  , "pkF"] <- with(intrsctSumDf, calcFscore(pkSens, pkSpec, beta = 1))
intrsctSumDf[  , "pkF2"] <- with(intrsctSumDf, calcFscore(pkSens, pkSpec, beta = 2))
intrsctSumDf[  , "pkFp5"] <- with(intrsctSumDf, calcFscore(pkSens, pkSpec, beta = 0.5))


### Common Sense Check
# Specificity, Sensitivity, F, should all be [0-1]
senseCheck <- function(dataFrame, cols = c("bpSens", "bpSpec", "pkSens", "pkSpec", "bpF", "pkF")) {
	for (col in cols) {
	cat(paste("All values in", col, "between 0 - 1?"))
		valuedCol <- dataFrame[!(is.nan(dataFrame[ , col]) | is.na(dataFrame[ , col])), col]
		if (all(valuedCol >= 0 & valuedCol <= 1)) {
			cat(paste("    True.\n"))
		} else {
			cat(paste("    FALSE!!!\n"))
			noWrong <- sum(!(valuedCol >= 0 & valuedCol <= 1))
			cat(paste("    Found", noWrong, "entries outside of 0 - 1.\n"))
		}
	}
}
senseCheck(intrsctSumDf)


### numPks and % Tables
# Generate with number of peaks calculated both
# (1) Allowing multiple intersections with one annotation
numPkTbl <- do.call(rbind, lapply(split(intrsctSumDf, intrsctSumDf$samp), function(X) {
					numPks <- X$numPks
					names(numPks) <- X$annot
					numPks
					})
			)
perPkTbl <- round(100*sweep(numPkTbl, 1, rowSums(numPkTbl), "/"), 5)
perPkTbl <- cbind(perPkTbl, rowSums(numPkTbl))
colnames(perPkTbl)[ncol(perPkTbl)] <- "TotalPks"

# (2) Limiting scoring to only one intersection per annotation
# --> count once will work for BOTH the "open" and "closed" annotation types
#	- i.e. can only lose points for intersecting with repressed region the first time
annotNumPkTbl <- do.call(rbind, lapply(split(annotIntrsctSumDf, annotIntrsctSumDf$samp), function(X) {
					numPks <- X$numPks
					names(numPks) <- X$annot
					numPks
					})
			)
annotPerPkTbl <- round(100*sweep(annotNumPkTbl, 1, rowSums(annotNumPkTbl), "/"), 5)
annotPerPkTbl <- cbind(annotPerPkTbl, rowSums(annotNumPkTbl))
colnames(annotPerPkTbl)[ncol(annotPerPkTbl)] <- "TotalPks"


### Write Output to .csv's
write.csv(annotSumDf, file.path(outputDir, "annotSummaryData.csv"))
write.csv(sampSumDf, file.path(outputDir, "sampSummaryData.csv"))
write.csv(intrsctSumDf, file.path(outputDir, "intrsctSummaryData.csv"))
write.csv(sampIntrsctSumDf, file.path(outputDir, "sampIntrsctSummaryData.csv"))
write.csv(annotIntrsctSumDf, file.path(outputDir, "annotIntrsctSummaryData.csv"))

write.csv(numPkTbl, file.path(outputDir, "numPkTbl.csv"))
write.csv(annotNumPkTbl, file.path(outputDir, "annotNumPkTbl.csv"))
write.csv(perPkTbl, file.path(outputDir, "perPkTbl.csv"))
write.csv(annotPerPkTbl, file.path(outputDir, "annotPerPkTbl.csv"))


### Plots
# (A) Plot Ordered Specificity v. Sensitivity
# --- & Plots of Annotation Cateogry by Significance Threshold Parameter
library(lattice)

if (runStyle == "-r") {
	if (algo == "fseq") {
		# Generate For All Annotations
		# --- Also generate for both peak and bp sensitivity
		pdf(file.path(outputDir, "paramOrdSpec-Fseq.pdf"), height = 8, width = 8)
		for (annotTP in unique(intrsctSumDf$annot)) {
			
			# Select Annotation And Order By Specificity
			intrsctSumDf.TP <- intrsctSumDf[with(intrsctSumDf, annot == annotTP), ]
			intrsctSumDf.TP <- intrsctSumDf.TP[with(intrsctSumDf.TP, order(bpSpec)), ]
			
			# Find Best Sensitivity Within Specificity Ranges
			bpSpecCuts <- seq(1, 0.1, by = -0.1)
			bestSensDf <- do.call(rbind, {
								sumSpecRangeList <- lapply(bpSpecCuts, function(x) {
															indxs <- with(intrsctSumDf.TP, bpSpec > (x - 0.1) & bpSpec < x)
															limIntrsct <- intrsctSumDf.TP[indxs, ]	
															limIntrsct[with(limIntrsct, which.max(bpSens)), c("bpSens", "bpSpec", "numPks", "meanPkWidth", "samp", "annot")]
															})
								names(sumSpecRangeList) <- paste(bpSpecCuts, bpSpecCuts - 0.1, sep = "-")
								sumSpecRangeList
							})
			write.csv(bestSensDf, file.path(outputDir, paste("bestSens", annotTP, ".csv", sep = "")))
			
			# Plot it
			with(intrsctSumDf.TP,{
				# Build Window Layout
				mLinePlot <- matrix(1, ncol = 10, nrow = 5)
				mBarPlots <- matrix(c(rep(2, 10), rep(3, 10)), byrow = T, ncol = 10)
				mLayout <- rbind(mLinePlot, mBarPlots)
				layout(mLayout)
				
				# Generate Spec/Sens Line Plot
				op <- par("mar")
				par("mar" = c(1.1, op[2:4]))
				plot(	bpSpec,
						main = paste("F-seq Parameters vs.", annotTP),
						ylab = "Specificity & Sensitivity",
						ylim = c(0, 1), yaxp = c(0, 1, 10),
						xaxt = "n",
						col = "firebrick",
						pch = 16, 
						cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
						type = "o",
						las = 1
						)
				lines(	bpSens,
						col = "black",
						pch = 16, 
						cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
						type = "o"
						)
				#lines(	pkSpec,
				#		col = "firebrick",
				#		pch = 15, 
				#		cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
				#		type = "o"
				#		)
				#lines(	pkSens,
				#		col = "black",
				#		pch = 15, 
				#		cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
				#		type = "o"
				#		)
				abline(h = seq(0, 1, by = 0.1), lty = "dotted", col = "grey")
				legend(	"topleft",
						c("Specificity", "Sensitivity"),
						pch = 16,
						lwd = 1,
						cex = 1.2,
						bg = "white",
						col = c("firebrick", "black")
						)
						
				# Add Parameter Barcharts Below
				par("mar" = c(1.1, op[2], 1.1, op[4]))	
				barplot(t, ylab = "Threshold (-t)",
						axes = F,
						col = ifelse(samp %in% bestSensDf$samp, "steelblue", "grey"),
						las = 1
						)
				axis(2, at = sort(unique(t)), las = 1)
				box()
				barplot(l, ylab = "Length (-l)",
						axes = F,
						col = ifelse(samp %in% bestSensDf$samp, "darkorange", "grey"),
						las = 1
						)
				axis(2, at = sort(unique(l)), las = 1)
				box()
				par("mar" = op)
			})
		}
		dev.off()
		
		### Improve & Expand This Analysis
		# Stacked Barplot of Number of Peaks / Annotation Cateogry
		pdf(file.path(outputDir, "annotStackedBar-Fseq.pdf"), height = 6, width = 10)
	print(
			barchart(	numPks ~ paramOne | paramTwo, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(	
			barchart(	numPks ~ paramTwo | paramOne, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(
			barchart(	numPks ~ paramOne | paramTwo, 
				data = annotIntrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(	
			barchart(	numPks ~ paramTwo | paramOne, 
				data = annotIntrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(
			barchart(	bpSpec ~ paramOne | paramTwo, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(	
			barchart(	bpSpec ~ paramTwo | paramOne, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		dev.off()
			
	} else if (algo == "macs") {	
		# Generate For All Annotations
		pdf(file.path(outputDir, "paramOrdSpec-Macs.pdf"), height = 8, width = 8)
		for (annotTP in unique(intrsctSumDf$annot)) {
		
			intrsctSumDf.TP <- intrsctSumDf[with(intrsctSumDf, annot == annotTP), ]
			intrsctSumDf.TP <- intrsctSumDf.TP[with(intrsctSumDf.TP, order(bpSpec)), ]

			# Find Best Sensitivity Within Specificity Ranges
			bpSpecCuts <- seq(1, 0.1, by = -0.1)
			bestSensDf <- do.call(rbind, {
								sumSpecRangeList <- lapply(bpSpecCuts, function(x) {
															indxs <- with(intrsctSumDf.TP, bpSpec > (x - 0.1) & bpSpec < x)
															limIntrsct <- intrsctSumDf.TP[indxs, ]	
															limIntrsct[with(limIntrsct, which.max(bpSens)), c("bpSens", "bpSpec", "numPks", "meanPkWidth", "samp", "annot")]
															})
								names(sumSpecRangeList) <- paste(bpSpecCuts, bpSpecCuts - 0.1, sep = "-")
								sumSpecRangeList
							})
			write.csv(bestSensDf, file.path(outputDir, paste("bestSens", annotTP, ".csv", sep = "")))

			with(intrsctSumDf.TP,{
				# Build Window Layout
				mLinePlot <- matrix(1, ncol = 10, nrow = 5)
				mBarPlots <- matrix(c(rep(2, 10), rep(3, 10)), byrow = T, ncol = 10)
				mLayout <- rbind(mLinePlot, mBarPlots)
				layout(mLayout)
				
				# Generate Spec/Sens Line Plot
				op <- par("mar")
				par("mar" = c(1.1, op[2:4]))
				plot(	bpSpec,
						main = paste("MACS2 Parameters vs.", annotTP),
						ylab = "Specificity & Sensitivity",
						ylim = c(0, 1), yaxp = c(0, 1, 10),
						xaxt = "n",
						col = "firebrick",
						pch = 16, 
						cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
						type = "o",
						las = 1
						)
				lines(	bpSens,
						col = "black",
						pch = 16, 
						cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
						type = "o"
						)
				#lines(	pkSpec,
				#		col = "firebrick",
				#		pch = 15, 
				#		cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
				#		type = "o"
				#		)
				#lines(	pkSens,
				#		col = "black",
				#		pch = 15, 
				#		cex = ifelse(samp %in% bestSensDf$samp, 1.4, 0.8), 
				#		type = "o"
				#		)
				abline(h = seq(0, 1, by = 0.1), lty = "dotted", col = "grey")
				legend(	"topleft",
						c("Specificity", "Sensitivity"),
						pch = 16,
						lwd = 1,
						cex = 1.2,
						bg = "white",				
						col = c("firebrick", "black")
						)
						
				# Add Parameter Barcharts Below
				par("mar" = c(1.1, op[2], 1.1, op[4]))	
				barplot(log10(q), ylab = "FDR (-q)",
						axes = F,
						col = ifelse(samp %in% bestSensDf$samp, "steelblue", "grey"),
						las = 1
						)
				axis(2, at = log10(sort(unique(q))), labels = round(log10(sort(unique(q))), 2), las = 1)
				box()
				barplot(e, ylab = "Extsize (-e)",
						axes = F,
						col = ifelse(samp %in% bestSensDf$samp, "darkorange", "grey"),
						las = 1
						)
				axis(2, at = sort(unique(e)), las = 1)
				box()
				par("mar" = op)
			})
		}		
		dev.off()
		
		### Improve & Expand This Analysis
		# Stacked Barplot of Number of Peaks / Annotation Cateogry
		pdf(file.path(outputDir, "annotStackedBar-Macs.pdf"), height = 6, width = 10)
		print(
			barchart(	numPks ~ paramOne | paramTwo, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(	
			barchart(	numPks ~ paramTwo | paramOne, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(
			barchart(	numPks ~ paramOne | paramTwo, 
				data = annotIntrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(	
			barchart(	numPks ~ paramTwo | paramOne, 
				data = annotIntrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
	print(
			barchart(	bpSpec ~ paramOne | paramTwo, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		print(	
			barchart(	bpSpec ~ paramTwo | paramOne, 
				data = intrsctSumDf,
				groups = annot,
				stack = TRUE,
				auto.key = list(space = "right")
				)
			)
		dev.off()			
	}
}


### (B) Lattice Bar Chart / Statistic ~ Sample | Encode Region
# Prepare Colors & Order Samples for Lattice Barchart
intrsctSampNames <- unique(intrsctSumDf$samp)
if (allSamps) {
	colSet <- getColorSet(intrsctSampNames, colScheme)
} else {
	# Apply a generic rainbow color pallete
	colSet <- rainbow(length(intrsctSampNames))
	names(colSet) <- intrsctSampNames
}
# Prepare Panel Settings
intrsctSumDf$samp <- ordered(intrsctSumDf$samp, levels = unique(intrsctSumDf$samp))  # for lattice, samples need to be a factor vector
sampSettings <- list(
	plot.polygon = list(col = colSet, border = "transparent"),
	strip.background = 	list(col = "grey")
)
# Make Barchart
pdf(file.path(outputDir, "trellisBarcharts.pdf"), width = 8, height = 8)
numericCols <- names(intrsctSumDf)[sapply(intrsctSumDf, is.numeric)]
for (numericCol in numericCols) {
	print(
	barchart(	intrsctSumDf[  , numericCol] ~ intrsctSumDf[ , "samp"] | intrsctSumDf[ , "annot"], 
				scales = list(x = list(rot = 45)),
				xlab = "ATAC-seq Samples",
				ylab = numericCol,
				par.settings = sampSettings
				)
			)
	}	
dev.off()

	
### (C) If runStyle = "-r" ; Display Statistics as Parameter Matrix
if (runStyle == "-r") {
	# Run genParamMat and seeParamMat for All Annotation Overlaps & Numeric Stats.
	numericCols <- names(intrsctSumDf)[sapply(intrsctSumDf, is.numeric)]
	for (curAnnot in unique(intrsctSumDf$annot)) {
			intrsctSumDfAnnot <- subset(intrsctSumDf, annot == curAnnot)
			# Make an output directory for all the matrices if it doesn't exist
			paramMatDir <- file.path(outputDir, paste(curAnnot, "-paramMatrices", sep = ""))
			mkDirCmd <- paste(	"if [ ! -d", paramMatDir, " ]; then",
								"mkdir", paramMatDir, ";",
								"fi"
								)
			system(mkDirCmd)
			if (algo == "fseq") {
				# Make one pdf containing all matrices, plus a .csv for every matrix
				pdf(file.path(paramMatDir, paste(curAnnot, "-paramMatrices.pdf", sep = "")), height = 8, width = 8)
					lapply(numericCols, function(numericCol) {
						numericMat <- genParamMat(intrsctSumDfAnnot, summaryStat = numericCol, paramCol = "samp")
						write.csv(numericMat, file.path(paramMatDir, paste(curAnnot, numericCol, "-paramMatrices.csv", sep = "")))
						pltLogA <- ifelse(max(numericMat, na.rm = T) <= 1, F, T)
						seeParamMat(numericMat, pltLog = pltLogA, c(numericCol, "F-seq", paste("\n in", curAnnot, "Intersection")), mkPdf = F)	
					})
				dev.off()			
			} else if (algo == "macs") {
				# Make one pdf containing all matrices, plus a .csv for every matrix
				pdf(file.path(paramMatDir, paste(curAnnot, "-paramMatrices.pdf", sep = "")), height = 8, width = 8)
					lapply(numericCols, function(numericCol) {
						numericMat <- genParamMat(intrsctSumDfAnnot, summaryStat = numericCol, paramCol = "samp", paramOne = "q", paramTwo = "e")
						write.csv(numericMat, file.path(paramMatDir, paste(curAnnot, numericCol, "-paramMatrices.csv", sep = "")))
						pltLogA <- ifelse(max(numericMat, na.rm = T) <= 1, F, T)
						seeParamMat(numericMat, pltLog = pltLogA, c(numericCol, "MACS2", paste("\n in", curAnnot, "Intersection")), mkPdf = F)	
					})
				dev.off()
		}
	}	 
}





