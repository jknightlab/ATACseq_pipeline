##### Analyze Properties of Peaks called on a Set of BED Files
## JHendry, 2016/03/09



### Load Basic ATAC-seq Data Processing Functions
source("run_basicATACroutines.R")

### Parse Command Line Argument(s)
args <- commandArgs(T)

cat("Running analyzePeaks.R:\n")
cat("Input argument(s):", args, "\n")
if (length(args) == 1) {
	cat("Running as a single directory containing .bed files.\n")
	setwd(args)
	bedFiles <- list.files(pattern = ".bed")
	outputPrefix <- args
} else {
	cat("Running as a set of", length(args), ".bed files.\n")
	bedFiles <- args
	outputPrefix <- paste(sapply(basename(bedFiles), function(x) { strsplit(x, ".", fixed = T)[[1]][1] }), collapse = "-")
}

# Re-order Such That: Fresh, Frozen, Fixed3d, Fixed7d
bedFileOrder <- do.call(c, mapply(function(sampTypes) { grep(sampTypes, bedFiles) }, sampTypes = c("fresh", "frozen", "fixed3d", "fixed7d")))

if (length(bedFileOrder) < length(bedFiles)) { # i.e., the samples aren't all fresh/frozen/fixed3d/fixed7d
	# Load files into a list WITHOUT ordering
	sampList <- lapply(bedFiles, read.table)
	names(sampList) <- sapply(basename(bedFiles), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })
	sampList <- lapply(sampList, function(X) { colnames(X) <- c("chr", "start", "end") ; X })	
	
	# Apply a generic rainbow color pallete
	colSet <- rainbow(length(bedFiles))
	names(colSet) <- names(sampList)
	} else {
	# Load Ordered Sample Files Into List
	bedFiles <- bedFiles[bedFileOrder]
	sampList <- lapply(bedFiles, read.table)
	names(sampList) <- sapply(basename(bedFiles), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })
	sampList <- lapply(sampList, function(X) { colnames(X) <- c("chr", "start", "end") ; X })

	# Get Color Set Matching Algo/Samples
	colSet <- getColorSet(names(sampList), colScheme)
}

### Define a Function To:
# (A) Calculate Total Number of Peaks
#	* Whole Genome & Per Chromosome
# (B) Calculate Peak Widths
# (C) Calculate Distance Between Peaks
analyzePeaks <- function(sampList, mkPlot = T) {	
	
	### Step (1): Compute Per Peak Statistics For Each Sample
	# Peak Widths
	sampList.pkStats <- lapply(sampList, function(X) {
		within(X, {
			pkWidth <- end - start
			log.pkWidth <- log10(pkWidth)
			})
	})
	
	# <distToNxt> computes distance from 3' of one peak
	# to 5' of the following peak.
	distToNxt <- function(X) {
		Xnow <- X[-nrow(X), ]
		Xnxt <- X[-1, ]
		names(Xnxt) <- paste("nxt", names(X), sep = ".")
		Xboth <- cbind(Xnow, Xnxt)
		distNxt <- with(Xboth, ifelse(chr == nxt.chr, nxt.start - end, NA))
		distNxt <- c(distNxt, NA) # last peak doesn't have a distance
		within(X, distNxt <- distNxt)
		}
	sampList.pkStats <- lapply(sampList.pkStats, distToNxt)	
	
	### Step (2): Compute Summary Data Frame For All Samples
	sampSummary <- do.call(rbind, lapply(sampList, summarizeBed))
	sampSummaryDf <- reformatToDf(sampSummary)	
	
	### Part 2: Generate Plots
	if (mkPlot) {
		### (1) Calculate Total Number of Peaks
		pdf(paste(outputPrefix, "_numPeaks.pdf", sep = ""), width = 8, height = 8)
		op <- par("mar")
		par("mar" = op + c(5, 0, 0, 0))		
		barplot(sampSummary[ , "numPks"], 
				main = "Number of Peak Calls Across Samples",
				ylab = "Peak Calls",
				col = colSet,
				las = 2)
		par("mar" = op)
		# Number of Peaks // Chr	
		chrOrder <- paste("chr", c(seq(1, 22), "X"), sep = "")	
		numPksByChr <- do.call(	rbind, lapply(sampList, function(X) { 
												chrCounts <- with(X, tapply(chr, chr, length))
												chrCountsOrd <- chrCounts[match(chrOrder, names(chrCounts))]
												if (any(is.na(chrCountsOrd))) {
													missingChr <- which(is.na(chrCountsOrd))
													chrCountsOrd[missingChr] <- 0
													names(chrCountsOrd)[missingChr] <- chrOrder[missingChr]
													}
												chrCountsOrd
												})
								)
		barplot(	numPksByChr[ , rev(colnames(numPksByChr))],
					horiz = T,
					col = colSet,
					las = 1,
					main = "Number of Consensus Peak Calls per Chromosome",
					xlab = "Number of Peaks"
					)
		legend("bottomright", rownames(numPksByChr), col = colSet, pch = 15, pt.cex = 2, bty = "n")
		dev.off()
		
		### (2) Calculate Peak Widths
		logPkWidthList <- lapply(sampList.pkStats, function(X) { X$log.pkWidth })
		density.logPkWidthList <- lapply(logPkWidthList, density)
		pdf(paste(outputPrefix, "_widthPeaks.pdf", sep = ""), width = 8, height = 8)
		op <- par("mar")
		par("mar" = op + c(5, 0, 0, 0))
		boxplot(logPkWidthList,
				main = "Width Distribution of Peaks",
				ylab = "Peak Width (bp)",
				col = colSet,
				las = 2
		)
		par("mar" = op)
		plot(	x = NULL, y = NULL, 
				xlim = c(0, 5), ylim = c(0, 1.8), 
				main = "Distribution of Peak Width Across Samples",
				xlab = "Peak Width [log(bp)]",
				ylab = "Density")
		lapply(rev(names(density.logPkWidthList)), function(sampName) {
			sampDens <- density.logPkWidthList[[sampName]]
			sampCol <- colSet[match(sampName, names(colSet))]
			lines(	sampDens, 
					lwd = 1.5, 
					col = rgb(t(col2rgb(sampCol)/255), alpha = 0.5), 
					type = "h")
		})
		legend("topright", names(density.logPkWidthList), col = colSet, pch = 15, pt.cex = 2, bty = "n")
		dev.off()

		### (3) Calculate Distance Between Peaks
		bwPeakList <- lapply(sampList.pkStats, function(X) { log10(X$distNxt) })
		bwPeakDensityList <- lapply(bwPeakList, density, na.rm = T)
		pdf(paste(outputPrefix, "_densityPeaks.pdf", sep = ""), width = 8, height = 8)
		op <- par("mar")
		par("mar" = op + c(5, 0, 0, 0))
		boxplot(bwPeakList,
				main = "Distribution of Between-peak Distances Across Samples",
				ylab = "Distance Between Peaks (bp)",
				col = colSet,
				las = 2
		)
		par("mar" = op)
		plot(	x = NULL, y = NULL, 
				xlim = c(0, 8), ylim = c(0, 0.6), 
				main = "Distribution of Between-peak Distances Across Samples",
				xlab = "Between-peak Distance [log(bp)]",
				ylab = "Density")
		lapply(rev(names(bwPeakDensityList)), function(sampName) {
			sampDens <- bwPeakDensityList[[sampName]]
			sampCol <- colSet[match(sampName, names(colSet))]
			lines(	sampDens, 
					lwd = 1.5,
					col = rgb(t(col2rgb(sampCol)/255), alpha = 0.5), 
					type = "h")
		})
		legend("topright", names(bwPeakDensityList), col = colSet, lwd = 2, bty = "n")
		dev.off()
	}
	
	write.csv(sampSummaryDf, paste(outputPrefix, "_sampleSummary.csv", sep = ""), row.names = F)
}


### Run analyze Peaks
analyzePeaks(sampList)





