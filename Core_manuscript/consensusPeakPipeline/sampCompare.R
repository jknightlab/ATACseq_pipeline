##### Compare Overlap Between Peaks Between Two Samples
## JHendry, 2016/01/23
##
## Add:
##	- How similar are the sizes of overlapping peaks
##	- How similar are the boundry positions of overlapping peaks
##	- Total and by chromosome number of peaks
## Issues:
##	- Nothing really feels clean or informative enough
##	- Add Venn Diagram Function
##	- What about duplicate peaks in intersect?
	# A peak in sample 1 can appear in the intersect more than once
	# if it intersects with 2 or more peaks from sample two. When
	# computing the number of peaks in sample 1 in the intersect,
	# must remove these duplicates.
	# If correct, the total number of peaks in sample 1
	# and sample 2 should be the same as the sum of those
	# in the intersect with those not in the intersect.


args <- commandArgs(T)

### Load Data From Both Samples & Intersect
samp1 <- read.table(file = args[1], header = F)
samp2 <- read.table(file = args[2], header = F)
intrsctAll <- read.table(file = args[3], header = F)

# Grab Sample Names
samp1Name <- strsplit(basename(args[1]), ".", fixed = T)[[1]][1]
samp2Name <- strsplit(basename(args[2]), ".", fixed = T)[[1]][1]

# Develop appropriate color scheme
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
colSamp1 <- getSampColor(samp1Name)
colSamp2 <- getSampColor(samp2Name)	
# If both colors are the same, make colSamp2 darker
calcColorGrad <- function(color, N) {
	asRgb <- t(col2rgb(color)/255)
	oneStep <- asRgb/(N + 1)
	mapply(function(Rstep, Gstep, Bstep) { rgb(asRgb - c(Rstep, Gstep, Bstep)) }, 
				Rstep = seq(0, oneStep[1]*(N - 1), by = oneStep[1]), 
				Gstep = seq(0, oneStep[2]*(N - 1), by = oneStep[2]),
				Bstep = seq(0, oneStep[3]*(N - 1), by = oneStep[3])
				) 
	}
if (colSamp1 == colSamp2) { colSamp2 <- calcColorGrad(colSamp2, 2)[2] }										
# Intersect is defined as grey
colIntrsct <- "darkgrey"
colSet <- c(colSamp1, colIntrsct, colSamp2)

# Name Columns
names(samp1) <- c("chr", "start", "end")
names(samp2) <- c("chr", "start", "end")
names(intrsctAll) <- c("chr", "start", "end", "samp1", "samp1Chr", "samp1Start", "samp1End", "samp2", "samp2Chr", "samp2Start", "samp2End")

# Compute on intrsctAll
#	- lengths of peaks for intrsct, samp1 and samp2
intrsctAll <- within(intrsctAll, { 
			intrsctLen <- end - start
			samp1Len <- samp1End - samp1Start
			samp2Len <- samp2End - samp2Start
			samp1Frac <- intrsctLen/samp1Len
			samp2Frac <- intrsctLen/samp2Len
			samp1 <- paste(samp1Chr, samp1Start, samp1End, sep = "-")
			samp2 <- paste(samp2Chr, samp2Start, samp2End, sep = "-")
	})


### Define Functions
# <calcOverlapByPercent> determines the number of overlapping peaks between two samples
# that satisfy a given reciprical percent overlap. E.g., if the percent overlap requested is 50% (0.5)
# then the overlapping region must be greater than 50% of the size of both samples peaks.
# Input:
#	- Sample 1 and Sample 2 .bed files, delineating peak boundries
#	- Intersect .bed file, dileneating regions where peaks from Sample 1 and 2 overlap
#	- perOverlap --> vector of percentages to check for number of overlaps
# Output:
#	- if mkPlot = T, barplot of overlaps
#		- would like to add boxplot of sizes for each reciprocal cut off (vs. total)
#	- Matrix of counts of sample 1, sample 2 and intersect for each percentage
calcOverlapByPercent <- function(	samp1, samp2, intrsct, samp1Name, samp2Name,
									perOverlap = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99), 
									mkPlot = T) {
	samp1NumPk <- nrow(samp1)
	samp2NumPk <- nrow(samp2)
	
	perOverlapMat <- do.call(cbind, lapply(perOverlap, function(y) {
		intrsctFrac <- intrsct[with(intrsct, samp1Frac > y & samp2Frac > y), ]
		numOfIntersections <- nrow(intrsctFrac)
		pkDupSamp1 <- with(intrsctFrac, samp1[duplicated(samp1)])
		pkDupSamp2 <- with(intrsctFrac, samp2[duplicated(samp2)])
	
		OneToOne <- subset(intrsctFrac, !(samp1 %in% pkDupSamp1 | samp2 %in% pkDupSamp2))
		NotOneToOne <- subset(intrsctFrac, samp1 %in% pkDupSamp1 | samp2 %in% pkDupSamp2)
	
		numOneToOne <- nrow(OneToOne)
		numSamp1NotOneToOne <- length(unique(NotOneToOne$samp1))
		numSamp2NotOneToOne <- length(unique(NotOneToOne$samp2))
	
		samp1NotIntrsct <- samp1NumPk - numOneToOne - numSamp1NotOneToOne
		samp2NotIntrsct <- samp2NumPk - numOneToOne - numSamp2NotOneToOne
	
		if ((samp1NumPk + samp2NumPk) != (2*numOneToOne + numSamp1NotOneToOne + samp1NotIntrsct + numSamp2NotOneToOne + samp2NotIntrsct)) {
			stop("Error: Total number of peaks not equal to intersect.")
		}
		c(samp1NotIntrsct, samp2NotIntrsct, numSamp1NotOneToOne, numSamp2NotOneToOne, numOneToOne)
	}))
	rownames(perOverlapMat) <- c(samp1Name, samp2Name, paste("Intersect n:n", samp1Name), paste("Intersect n:n", samp2Name), "Intersect 1:1")
	colnames(perOverlapMat) <- paste(100*perOverlap, "%", sep = "")

	# Generate Plot
	if (mkPlot) { 
				op <- par("mar")
				pdf(paste(samp1Name, samp2Name, "OverlapByPercent.pdf", sep = "-"), width = 8, height = 6)				
				par("mar" = op + c(0, 0, 0, 6))
				par("xpd" = T)
				barplot(	perOverlapMat[c(1, 3, 5, 4, 2), ],
							ylab = "Number of Peaks",
							xlab = "Extent of Overlap",
							main = paste(samp1Name, "vs.", samp2Name),
							col = rep(colSet, each = 2)[-3],
							density = c(300, 25, 300, 25, 300)
							) 
				legend(	"right", inset = c(-0.25, 0),
						rownames(perOverlapMat)[c(1, 5, 2)],
						pch = 15,
						pt.cex = 1.5,
						col = colSet,
						bty = "n"
						)
				dev.off()
				par("xpd" = F)
				par("mar" = op)
				}
	perOverlapMat
}


# <calcOverlapByBp> similar to above, but interesect is conditioned on a certain number of bps
# instead of on a reciprocal percentage. I.e., 100bp means the intersect must contain at least
# 100 bps.
calcOverlapByBp <- function(samp1, samp2, intrsct, samp1Name, samp2Name,
							bpOverlap = c(1, 10, 50, 100, 200, 500, 1000), 
							mkPlot = T) {
	samp1NumPk <- nrow(samp1)
	samp2NumPk <- nrow(samp2)
	bpOverlapMat <- do.call(cbind, lapply(bpOverlap, function(y) {
		intrsctFrac <- intrsct[with(intrsct, intrsctLen > y), ]
		numOfIntersections <- nrow(intrsctFrac)
		pkDupSamp1 <- with(intrsctFrac, samp1[duplicated(samp1)])
		pkDupSamp2 <- with(intrsctFrac, samp2[duplicated(samp2)])
	
		OneToOne <- subset(intrsctFrac, !(samp1 %in% pkDupSamp1 | samp2 %in% pkDupSamp2))
		NotOneToOne <- subset(intrsctFrac, samp1 %in% pkDupSamp1 | samp2 %in% pkDupSamp2)
	
		numOneToOne <- nrow(OneToOne)
		numSamp1NotOneToOne <- length(unique(NotOneToOne$samp1))
		numSamp2NotOneToOne <- length(unique(NotOneToOne$samp2))
	
		samp1NotIntrsct <- samp1NumPk - numOneToOne - numSamp1NotOneToOne
		samp2NotIntrsct <- samp2NumPk - numOneToOne - numSamp2NotOneToOne
		
		if ((samp1NumPk + samp2NumPk) != (2*numOneToOne + numSamp1NotOneToOne + samp1NotIntrsct + numSamp2NotOneToOne + samp2NotIntrsct)) {
			stop("Error: Total number of peaks not equal to intersect.")
		}
		c(samp1NotIntrsct, samp2NotIntrsct, numSamp1NotOneToOne, numSamp2NotOneToOne, numOneToOne)
	}))
	rownames(bpOverlapMat) <- c(samp1Name, samp2Name, paste("Intersect n:n", samp1Name), paste("Intersect n:n", samp2Name), "Intersect 1:1")
	colnames(bpOverlapMat) <- paste(bpOverlap, "bp", sep = "")
	
	# Generate Plot
	if (mkPlot) {
				op <- par("mar")
				pdf(paste(samp1Name, samp2Name, "OverlapByBp.pdf", sep = "-"), width = 8, height = 6)				
				par("mar" = op + c(0, 0, 0, 6))
				par("xpd" = T)
				barplot(	bpOverlapMat[c(1, 3, 5, 4, 2), ],
							ylab = "Number of Peaks",
							xlab = "Extent of Overlap",
							main = paste(samp1Name, "vs.", samp2Name),
							col = rep(colSet, each = 2)[-3],
							density = c(300, 25, 300, 25, 300)
							) 
				legend(	"right", inset = c(-0.25, 0),
						rownames(bpOverlapMat)[c(1, 5, 2)],
						pch = 15,
						pt.cex = 1.5,
						col = colSet,
						bty = "n"
						)
				dev.off()
				par("xpd" = F)
				par("mar" = op)
				}
	bpOverlapMat
}


# <calcHighQualOverlap>
# To deal with the confounding and insatisfactory nature of using either
# % reciprocal overlap or bp individually as a standard for an overlap.
# Instead, define high quality overlap as >Nbp (tentatively, N=20) and
# > 50% overlap.

calcHighQualityOverlap <- function(	samp1, samp2, intrsct, samp1Name, samp2Name,
									perOverlap = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99),
									bpOverlap = 20,
									mkPlot = T) {
	samp1NumPk <- nrow(samp1)
	samp2NumPk <- nrow(samp2)
	
	hqOverlapMat <- do.call(cbind, lapply(perOverlap, function(y) {
		intrsctFrac <- intrsct[with(intrsct, samp1Frac > y & samp2Frac > y & intrsctLen > bpOverlap), ]
		numOfIntersections <- nrow(intrsctFrac)
		pkDupSamp1 <- with(intrsctFrac, samp1[duplicated(samp1)])
		pkDupSamp2 <- with(intrsctFrac, samp2[duplicated(samp2)])
	
		OneToOne <- subset(intrsctFrac, !(samp1 %in% pkDupSamp1 | samp2 %in% pkDupSamp2))
		NotOneToOne <- subset(intrsctFrac, samp1 %in% pkDupSamp1 | samp2 %in% pkDupSamp2)
	
		numOneToOne <- nrow(OneToOne)
		numSamp1NotOneToOne <- length(unique(NotOneToOne$samp1))
		numSamp2NotOneToOne <- length(unique(NotOneToOne$samp2))
	
		samp1NotIntrsct <- samp1NumPk - numOneToOne - numSamp1NotOneToOne
		samp2NotIntrsct <- samp2NumPk - numOneToOne - numSamp2NotOneToOne
	
		if ((samp1NumPk + samp2NumPk) != (2*numOneToOne + numSamp1NotOneToOne + samp1NotIntrsct + numSamp2NotOneToOne + samp2NotIntrsct)) {
			stop("Error: Total number of peaks not equal to intersect.")
		}
		c(samp1NotIntrsct, samp2NotIntrsct, numSamp1NotOneToOne, numSamp2NotOneToOne, numOneToOne)
	}))
	rownames(hqOverlapMat) <- c(samp1Name, samp2Name, paste("Intersect n:n", samp1Name), paste("Intersect n:n", samp2Name), "Intersect 1:1")
	colnames(hqOverlapMat) <- paste(100*perOverlap, "%", sep = "")

	# Generate Plot
	if (mkPlot) { 
				op <- par("mar")
				pdf(paste(samp1Name, samp2Name, "OverlapHq.pdf", sep = "-"), width = 8, height = 6)				
				par("mar" = op + c(0, 0, 0, 6))
				par("xpd" = T)
				barplot(	hqOverlapMat[c(1, 3, 5, 4, 2), ],
							ylab = "Number of Peaks",
							xlab = "Extent of Overlap",
							main = paste(samp1Name, " vs. ", samp2Name, "\nHigh-Quality Overlaps >", bpOverlap, "bp", sep = ""),
							col = rep(colSet, each = 2)[-3],
							density = c(300, 25, 300, 25, 300)
							) 
				legend(	"right", inset = c(-0.25, 0),
						rownames(hqOverlapMat)[c(1, 5, 2)],
						pch = 15,
						pt.cex = 1.5,
						col = colSet,
						bty = "n"
						)
				dev.off()
				par("xpd" = F)
				par("mar" = op)
				}
	hqOverlapMat
}


# <toPercentages> Convert the contents of a matrix to percentages of either row or column totals
toPercentages <- function(dataMat, totalsBy = "columns", roundTo = 4) {
	if (totalsBy == "columns") {
		Totals <- colSums(dataMat)
		dataMatCmbd <- rbind(dataMat, Totals)
		dataMatFracs <- sweep(dataMatCmbd, 2, Totals, "/")
		dataMatPers <- round(100*dataMatFracs, roundTo)
		}
	if (totalsBy == "rows") {
		Totals <- rowSums(dataMat)
		dataMatCmbd <- cbind(dataMat, Totals)
		dataMatFracs <- sweep(dataMatCmbd, 1, Totals, "/")
		dataMatPers <- round(100*dataMatFracs, roundTo)
	}
	dataMatPers
}


### Run Functions
bpOverlapMat <- calcOverlapByBp(samp1, samp2, intrsctAll, samp1Name, samp2Name)
bpOverlapMatPercents <- toPercentages(bpOverlapMat)
PercentOverlapMat <- calcOverlapByPercent(samp1, samp2, intrsctAll, samp1Name, samp2Name)
PercentOverlapMatPercents <- toPercentages(PercentOverlapMat)
hqOverlapMat <- calcHighQualityOverlap(samp1, samp2, intrsctAll, samp1Name, samp2Name)
hqOverlapMatPercents <- toPercentages(hqOverlapMat)

write.table(bpOverlapMat, file = paste(samp1Name, samp2Name, "OverlapByBp-absTable.txt", sep = "-"), sep = "\t")
write.table(PercentOverlapMat, file = paste(samp1Name, samp2Name, "OverlapByPercent-absTable.txt", sep = "-"), sep = "\t")
write.table(hqOverlapMat, file = paste(samp1Name, samp2Name, "OverlapHq-absTable.txt", sep = "-"), sep = "\t")
write.table(bpOverlapMatPercents, file = paste(samp1Name, samp2Name, "OverlapByBp-perTable.txt", sep = "-"), sep = "\t")
write.table(PercentOverlapMatPercents, file = paste(samp1Name, samp2Name, "OverlapByPercent-perTable.txt", sep = "-"), sep = "\t")
write.table(hqOverlapMatPercents, file = paste(samp1Name, samp2Name, "OverlapHq-perTable.txt", sep = "-"), sep = "\t")


