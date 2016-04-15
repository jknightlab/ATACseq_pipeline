##### Analyzing The Effect of Varying Two Algorithm Parameters on Peak Calling
## JHendry, 2016/03/10


## Note: all output files MUST be prefixed with "paramSpace"
## * "paramSpace" files are subsequently moved into output directory by run_paramSpace-v2.sh

### Source Basic ATAC-seq Data Processing Functions
source("run_basicATACroutines.R")

args <- commandArgs(T)

### Process Command Line Input
## Two Input Options
# (1) Directory Containing .bed Files
# (2) Range values
# In this step, "bedFiles" (a vector of .bed file names) is generated
# for either input options, later to be used to load data
if (length(args) == 1) {
	setwd(args)
	bedFiles <- list.files(pattern = ".bed")
	cat("Running as a single directory containing", length(bedFiles), ".bed files.\n")
} else {
	bedFiles <- args
	cat("Running as a set of", length(bedFiles), ".bed files.\n")
}


### Load & Clean Data
sampList <- lapply(bedFiles, read.table)
names(sampList) <- sapply(basename(bedFiles), function(x) { strsplit(x, ".", fixed = T)[[1]][1] })
sampList <- lapply(sampList, function(X) { colnames(X) <- c("chr", "start", "end") ; X })


### Calculate Statistics on Each .bed File
pkSummary <- do.call(rbind, lapply(sampList, summarizeBed))
pkSummaryDf <- reformatToDf(pkSummary)
write.csv(pkSummaryDf, "paramSpace-samplesSummary.csv")


### Attempt To Determine Which Type of Peak Calling Algorithm Was Used
# Outcome of this influences parameters fed to genParamMat() and seeParamMat
if (length(grep("fseq", bedFiles)) == length(bedFiles)) {
	cat("All input files contain 'fseq' in name.\n")
	cat("Generating parameter matrices using F-seq parameters:\n")
	cat("Parameter 1: 'l', Parameter 2: 't'\n")
	algo <- "fseq"
} else if (length(grep("macs", bedFiles)) == length(bedFiles)) {
	cat("All input files contain 'macs' in name.\n")
	cat("Generating parameter matrices using MACS2 parameters.\n")
	cat("Parameter 1: 'q', Parameter 2: 'e'\n")
	algo <- "macs"
}


### Generate Parameter Space Heatmap Visualizations and CSV files
if (algo == "fseq") {
	# Total Number of Peaks & Number of Basepairs
	numPksMat <- genParamMat(pkSummaryDf, summaryStat = "numPks")
	seeParamMat(numPksMat, fileName = "numPks-Fseq")
	write.csv(numPksMat, file = "paramSpace-numPks-Fseq.csv")
	numBpsMat <- genParamMat(pkSummaryDf, summaryStat = "numBps")
	seeParamMat(numBpsMat, title = c("No. Bps", "F-seq"), fileName = "numBps-Fseq")
	write.csv(numBpsMat, file = "paramSpace-numBps-Fseq.csv")
	
	# Mean and Median Peak Width
	meanPkWidthMat <- genParamMat(pkSummaryDf, summaryStat = "meanPkWidth")
	seeParamMat(meanPkWidthMat, title = c("Avg. Peak Width", "F-seq"), fileName = "meanWidth-Fseq")
	write.csv(meanPkWidthMat, file = "paramSpace-meanWidth-Fseq.csv")
	medPkWidthMat <- genParamMat(pkSummaryDf, summaryStat = "medPkWidth")
	seeParamMat(medPkWidthMat, title = c("Median Peak Width", "F-seq"), fileName = "medWidth-Fseq")
	write.csv(medPkWidthMat, file = "paramSpace-medWidth-Fseq.csv")
} else if (algo == "macs") {
	# Total Number of Peaks	
	numPksMat <- genParamMat(pkSummaryDf, summaryStat = "numPks", paramOne = "q", paramTwo = "e")
	seeParamMat(numPksMat, title = c("No. Peaks", "MACS2"), xlabVal = "Extension Size (--extsize)", ylabVal = "FDR (--qvalue)", fileName = "numPks-Macs")
	write.csv(numPksMat, file = "paramSpace-numPks-Macs.txt")
	numBpsMat <- genParamMat(pkSummaryDf, summaryStat = "numBps", paramOne = "q", paramTwo = "e")
	seeParamMat(numBpsMat, title = c("No. Bps", "MACS2"), xlabVal = "Extension Size (--extsize)", ylabVal = "FDR (--qvalue)", fileName = "numBps-Macs")
	write.csv(numBpsMat, file = "paramSpace-numBps-Macs.csv")
	
	# Mean and Median Peak Width
	meanPkWidthMat <- genParamMat(pkSummaryDf, summaryStat = "meanPkWidth", paramOne = "q", paramTwo = "e")
	seeParamMat(meanPkWidthMat, title = c("Avg. Peak Width", "MACS2"), xlabVal = "Extension Size (--extsize)", ylabVal = "FDR (--qvalue)", fileName = "meanWidth-Macs")
	write.csv(meanPkWidthMat, file = "paramSpace-meanWidth-Macs.csv")
	medPkWidthMat <- genParamMat(pkSummaryDf, summaryStat = "medPkWidth", paramOne = "q", paramTwo = "e")
	seeParamMat(medPkWidthMat, title = c("Avg. Peak Width", "MACS2"), xlabVal = "Extension Size (--extsize)", ylabVal = "FDR (--qvalue)", fileName = "medWidth-Macs")
	write.csv(medPkWidthMat, file = "paramSpace-medWidth-Macs.csv")
}






