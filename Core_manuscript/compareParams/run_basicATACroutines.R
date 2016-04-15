 ### Repository of Useful Functions
# Requisites for inclusion:
# Used in more than one script


### Basic Environmentals
# --- Below is All For Plot Colouring --- #
# Basic idea is, given a set of bed file names, return a set of colours
# colours are determine based on algorithm name and sample name, from the
# color scheme "colScheme" data frame.
# Assumptions: characters 1-4 give algorithm name "macs", "fseq"
# Sample name is after "_"
colScheme <- data.frame(	"Algorithm" = rep(c("fseq", "macs"), each = 4),
							"Samples" = rep(c("fresh", "frozen", "fixed3d", "fixed7d"), 2),
							"Colors" = c("darkorange", "orangered", "orangered3", "tomato4", "lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"),
							stringsAsFactors = F)
							
getSampColor <- function(name) { 
							algo <- substr(name, 1, 4)
							algoColScheme <- colScheme[colScheme$Algorithm == algo, ]
							sample <- strsplit(name, "_", fixed = T)[[1]][2]
							algoColScheme[match(sample, algoColScheme$Samples), "Colors"]
						}

calcColorGrad <- function(color, N) {
		asRgb <- t(col2rgb(color)/255)
		oneStep <- asRgb/(N + 1)
		mapply(function(Rstep, Gstep, Bstep) { rgb(asRgb - c(Rstep, Gstep, Bstep)) }, 
					Rstep = seq(0, oneStep[1]*(N - 1), by = oneStep[1]), 
					Gstep = seq(0, oneStep[2]*(N - 1), by = oneStep[2]),
					Bstep = seq(0, oneStep[3]*(N - 1), by = oneStep[3])
					) 
		}
		
getColorSet <- function(fileNames, colScheme) {
	colSet <- sapply(fileNames, getSampColor)
	if (any(duplicated(colSet))) {
		# check if any are duplicated
		# get name(s) of duplicated color
		# replace duplicated value(s) with gradient of that color
		dupColors <- unique(colSet[which(duplicated(colSet))])
		for (dupColor in dupColors) {
			numDup <- length(which(colSet == dupColor))
			colSet[which(colSet == dupColor)] <- calcColorGrad(dupColor, numDup)
			}
	}
	colSet
}	
		

### Basic Functions
### ======== Summary Statistics of BED files ====== #
# <summarizeBed>
summarizeBed <- function(X) {
	numPks <- nrow(X)
	meanPkWidth <- mean(with(X, end - start))
	medPkWidth <- median(with(X, end - start))
	numBps <- sum(with(X, as.numeric(end - start))) # see .Machine$integer.max
	numPks <- ifelse(numBps == 0, 0, numPks)
	c("numPks" = numPks, "meanPkWidth" = meanPkWidth, "medPkWidth" = medPkWidth, "numBps" = numBps)
}

# <reformatToDf>
reformatToDf <- function(df, nameCol = "name") {
		df <- as.data.frame(df)
		df[ , nameCol] <- rownames(df) ; rownames(df) <- 1:nrow(df)
		df[ , c(ncol(df), seq(1, ncol(df) - 1))] # put name in front
	}	

### ======== Generating Parameter Matrices ======== ###
# <getParams>
getParams <- function(X, paramCol = NULL, paramOne = "t", paramTwo = "l") {
	paramOneExpr <- paste(paramOne, "[0-9]*", sep = "")
	paramTwoExpr <- paste(paramTwo, "[0-9]*", sep = "")
	paramStrings <- X[ , paramCol]
	X[ , "paramOne"] <- regmatches(paramStrings, regexpr(paramOneExpr, paramStrings))
	X[ , paramOne] <- as.numeric(substr(X[ , "paramOne"], 2, nchar(X[, "paramOne"])))
	X[ , "paramTwo"] <- regmatches(paramStrings, regexpr(paramTwoExpr, paramStrings))
	X[ , paramTwo] <- as.numeric(substr(X[ , "paramTwo"], 2, nchar(X[, "paramTwo"]))) 
	X
}

# <genParamMat> Version 1: Get Parameters From Rownames
# genParamMat <- function(X, summaryStat = "numPeaks", paramOne = "t", paramTwo = "l") {
#	# paramOne is y-axis
#	# paramTwo is x-axis
#	summaryStatDf <- data.frame(	"summaryStat" = X[ , summaryStat],
#									"paramOneVal" = NA,
#									"paramTwoVal" = NA
#									)
#	paramOneExpr <- paste(paramOne, "[0-9]*", sep = "")
#	paramTwoExpr <- paste(paramTwo, "[0-9]*", sep = "")
#	summaryStatDf <- within(summaryStatDf, {
#							paramOneVal <- regmatches(rownames(summaryStatDf), regexpr(paramOneExpr, rownames(summaryStatDf)))
#							paramTwoVal <- regmatches(rownames(summaryStatDf), regexpr(paramTwoExpr, rownames(summaryStatDf)))
#						})
#
#	summaryStatDfOrd <- summaryStatDf[with(summaryStatDf, order(paramOneVal, paramTwoVal)), ]
#	summaryStatMat <- with(summaryStatDfOrd, 	matrix(	summaryStat, 
#												nrow = length(unique(paramTwoVal)), # because sorted FIRST on paramOne, and populated rows first
#												dimnames = list(unique(paramTwoVal), unique(paramOneVal))
#												))
#	summaryStatMat
#}


# <genParamMat> Version 2: From a Define Parameter Column (Instead of From Row Names)
# Parameter values from each sample are parsed from a define column, as opposed to from
# the rownames. In essence, works from a data frame instead of from a matrix.
# --> This is a much more flexible implementation, as can take any column as
# parameter column.
genParamMat <- function(X, summaryStat = "numPks", paramCol = "name", paramOne = "t", paramTwo = "l") {
		# paramOne is y-axis
		# paramTwo is x-axis
		summaryStatDf <- data.frame(	"summaryStat" = X[ , summaryStat],
										"paramSet" = X[ , paramCol],
										"paramOneVal" = NA,
										"paramTwoVal" = NA,
										stringsAsFactors = F
										)
		paramOneExpr <- paste(paramOne, "[0-9]*", sep = "")
		paramTwoExpr <- paste(paramTwo, "[0-9]*", sep = "")
		summaryStatDf <- within(summaryStatDf, {
								paramOneVal <- regmatches(paramSet, regexpr(paramOneExpr, paramSet))
								paramTwoVal <- regmatches(paramSet, regexpr(paramTwoExpr, paramSet))
							})
						
		summaryStatDfOrd <- summaryStatDf[with(summaryStatDf, order(paramOneVal, paramTwoVal)), ]
		summaryStatMat <- with(summaryStatDfOrd, 	matrix(	summaryStat, 
														nrow = length(unique(paramTwoVal)), 
														dimnames = list(unique(paramTwoVal), unique(paramOneVal))
												))
		summaryStatMat
}

seeParamMat <- function(X, pltLog = T, title = c("No. Peaks", "F-seq", ""), xlabVal = "Length Value (l)", ylabVal = "Threshold Value (t)", fileName= "file", mkPdf = T) {
		library(lattice)
		lattice.options(axis.padding = list(factor = 0.5))
		if (pltLog) {
			X <- log10(X)
		}
		colGrad <- colorRampPalette(c("darkorange", "white", "blue"))
		if (mkPdf) { pdf(paste("paramSpace-", fileName, ".pdf", sep = ""), width = 8, height = 8) }
		print(
			levelplot(	X,
					col.regions = colGrad(nrow(X)*ncol(X)),
					main = paste(title[1], ifelse(pltLog, "(log10)", ""), "Across", title[2], "Parameter Space", title[3]),
					xlab = xlabVal,
					ylab = ylabVal,
					border = T,
					scales = list(tck = 0)				
					)
			)
		if (mkPdf) { dev.off() }
}
