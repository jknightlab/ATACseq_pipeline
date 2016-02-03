##### Analyze Peak Calling Results Between Samples -- Not Looking At Overlap
## JHendry, 2016/01/23
##
##
## Execute as:
## 	- Rscript analyzePeaks-cli.R <samp1.bed> <samp2.bed> ... <sampN.bed>
##								OR
##	- Rscript analyzePeaks-cli.R <samp-dir>
## Note that an arbitrary number of samples can be fed to the script.
## These files are all combined into a single list, upon which
## all subsequent analyses are done.

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
bedFiles <- bedFiles[bedFileOrder]

# Load Sample Files Into List
sampList <- lapply(bedFiles, read.table)
names(sampList) <- sapply(basename(bedFiles), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })
sampList <- lapply(sampList, function(X) { colnames(X) <- c("chr", "start", "end") ; X })

# Determine an appropriate color scheme, including gradients for repeated algorithm/sample instances
colScheme <- data.frame("Algorithm" = rep(c("fseq", "macs"), each = 4),
						"Samples" = rep(c("fresh", "frozen", "fixed3d", "fixed7d"), 2),
						"Colors" = c("darkorange", "orangered", "orangered3", "tomato4", "lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"),
						stringsAsFactors = F)
						
getSampColor <- function(name) { 
						algo <- substr(name, 1, 4)
						algoColScheme <- colScheme[colScheme$Algorithm == algo, ]
						sample <- strsplit(name, "_", fixed = T)[[1]][2]
						algoColScheme[match(sample, algoColScheme$Samples), "Colors"]
					}
colSet <- sapply(names(sampList), getSampColor)
# Linear gradient across RGB for duplicate algorithm/sample instances
calcColorGrad <- function(color, N) {
	asRgb <- t(col2rgb(color)/255)
	oneStep <- asRgb/(N + 1)
	mapply(function(Rstep, Gstep, Bstep) { rgb(asRgb - c(Rstep, Gstep, Bstep)) }, 
				Rstep = seq(0, oneStep[1]*(N - 1), by = oneStep[1]), 
				Gstep = seq(0, oneStep[2]*(N - 1), by = oneStep[2]),
				Bstep = seq(0, oneStep[3]*(N - 1), by = oneStep[3])
				) 
	}
if (any(duplicated(colSet))) {
	# check if any are duplicated
	# get name(s) of duplicated color
	# replace duplicated value(s) with gradient of that color
	dupColors <- unique(colSet[which(duplicated(colSet))])
	for (dupColor in dupColors) {
		numDup <- length(which(colSet == dupColor))
		colSet[which(colSet == dupColor)] <- calcColorGrad(dupColors, numDup)
		}
}


### Define Functions
# <analyzePeaks> calculates the number of peaks and mean peak length
# Input: list() of .bed files delineating peak boundries
# Output: Summary Matrix of Num.Peaks, MeanPkSize, and barplot/boxplot visualizations

analyzePeaks <- function(sampList, mkPlot = T) {
	sampList <- lapply(sampList, function(X) {	within(X, { lenBp <- end - start ; logLenBp <- log10(lenBp) }) })
	sampList <- lapply(sampList, function(X) {
		# Compute "Density"
		# (a) distance to next peak
		# (b) sliding window
		distToNxt <- rep(NA, nrow(X)) 
		for (i in 1:(nrow(X) - 1)) {
			if (X[i, "chr"] == X[i + 1, "chr"]) {
				distToNxt[i] <- X[i + 1, "start"] - X[i, "end"]
			} else {
				distToNxt[i] <- NA
			}
		}
		within(X, distNxt <- distToNxt)	
	})
	sampListSummary <- lapply(sampList, function(X) { c(nrow(X), mean(X$lenBp)) })
	sampMat <- do.call(rbind, sampListSummary)
	colnames(sampMat) <- c("Num.Peaks", "Mean.Bp")
	if (mkPlot) {
		### Total Number of Peaks
		pdf(paste(outputPrefix, "_numPeaks.pdf", sep = ""), width = 8, height = 8)
		op <- par("mar")
		par("mar" = op + c(5, 0, 0, 0))		
		barplot(sampMat[ , "Num.Peaks"], 
				main = "Number of Peak Calls Across Samples",
				ylab = "Peak Calls",
				col = colSet,
				las = 2)
		par("mar" = op)
		# Number of Peaks // Chr	
		chrOrder <- paste("chr", c(seq(1, 22), "X"), sep = "")	
		pkNumByChr <- do.call(	rbind, lapply(sampList, function(X) { 
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
		barplot(	pkNumByChr[ , rev(colnames(pkNumByChr))],
					horiz = T,
					col = colSet,
					las = 1,
					main = "Number of Consensus Peak Calls per Chromosome",
					xlab = "Number of Peaks"
					)
		legend("bottomright", rownames(pkNumByChr), col = colSet, pch = 15, pt.cex = 2, bty = "n")
		dev.off()
		
		### Plot Peak Size
		lengthList <- lapply(sampList, function(X) { X$logLenBp })
		LenDensityList <- lapply(lengthList, density)
		pdf(paste(outputPrefix, "_widthPeaks.pdf", sep = ""), width = 8, height = 8)
		op <- par("mar")
		par("mar" = op + c(5, 0, 0, 0))
		boxplot(lengthList,
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
		lapply(rev(names(LenDensityList)), function(sampName) {
			sampDens <- LenDensityList[[sampName]]
			sampCol <- colSet[match(sampName, names(colSet))]
			lines(	sampDens, 
					lwd = 1.5, 
					col = rgb(t(col2rgb(sampCol)/255), alpha = 0.5), 
					type = "h")
		})
		legend("topright", names(LenDensityList), col = colSet, pch = 15, pt.cex = 2, bty = "n")
		dev.off()
		
		### Plot Between Peak Distannces
		bwPeakList <- lapply(sampList, function(X) { log10(X$distNxt) })
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
	sampMat
}

### Run Functions
analyzePeaks(sampList)







