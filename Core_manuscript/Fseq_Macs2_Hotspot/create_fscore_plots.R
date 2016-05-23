args <- commandArgs(T)
input_file <- args[1]
peakcaller <- args[2]
datatype <- args[3]

temp <- read.table(input_file, header=TRUE)


figurename <- paste (peakcaller, datatype, "QC.pdf", sep=".")
pdf(figurename, width = 10, height = 10)

maintitle <- paste ("FDR for ", peakcaller, datatype, sep=" ")
par(mar=c(15.1,4.1,4.1,2.1))
barplot (temp$called_nonTSS/temp$all_bases, col="darkslategrey", ylim=c(0,1), names=temp$parameters, main=maintitle, las=2, cex.main=3, cex.axis=1.2, cex.names=1)

maintitle <- paste ("Specificity for ", peakcaller, datatype, sep=" ")
par(mar=c(15.1,4.1,4.1,2.1))
barplot (temp$called_TSS/temp$all_bases, col="darkslategrey", ylim=c(0,1), names=temp$parameters, main=maintitle, las=2, cex.main=3, cex.axis=1.2, cex.names=1)

maintitle <- paste ("Sensitivity for ", peakcaller, datatype, sep=" ")
par(mar=c(15.1,4.1,4.1,2.1))
barplot (temp$called_TSS/temp$all_TSS, col="darkslategrey", ylim=c(0,1), names=temp$parameters, main=maintitle, las=2, cex.main=3, cex.axis=1.2, cex.names=1)

sn <- temp$called_TSS/temp$all_TSS
sp <- temp$called_TSS/temp$all_bases
beta <- 0.5

Fscore <- (1+beta*beta)*((sn*sp)/(sn+beta*beta*sp))

maintitle <- paste ("FScore for ", peakcaller, datatype, sep=" ")
par(mar=c(15.1,4.1,4.1,2.1))
barplot (Fscore, col="darkslategrey", ylim=c(0,1), names=temp$parameters, main=maintitle, las=2, cex.main=3, cex.axis=1.2, cex.names=1)

fdr <- temp$called_nonTSS/temp$all_bases
tpr <- temp$called_TSS/temp$all_TSS

maintitle <- paste ("FDR/TPR for ", peakcaller, datatype, sep=" ")
par(mar=c(15.1,4.1,4.1,2.1))
barplot (fdr/tpr, col="darkslategrey", ylim=c(0,max(fdr/tpr)), names=temp$parameters, main=maintitle, las=2, cex.main=3, cex.axis=1.2, cex.names=1)

dev.off()












