data <- read.table("reports_summary.txt", header=TRUE, sep='\t')

# ============ BOXPLOT ============
par(mfrow=c(3,3))

boxplot(as.integer(data[1,seq(1,3)]), as.integer(data[1,seq(4,6)]), main="Number of mapped reads after filtering", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[2,seq(1,3)]), as.integer(data[2,seq(4,6)]), main="Total number of peaks", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[3,seq(1,3)]), as.integer(data[3,seq(4,6)]), main="Number of peaks per chromosome", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[4,seq(1,3)]), as.integer(data[4,seq(4,6)]), main="Number of peaks normalized by total mapped reads", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[5,seq(1,3)]), as.integer(data[5,seq(4,6)]), main="Number of peaks normalized by reads mapped per chr", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[6,seq(1,3)]), as.integer(data[6,seq(4,6)]), main="Reads in peaks normalized by total mapped reads", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[7,seq(1,3)]), c(3.07623, 3.50766, 3.92653), main="Signal to noise ratio", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[8,seq(1,3)]), as.integer(data[8,seq(4,6)]), main="Average peak width", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))

boxplot(as.integer(data[9,seq(1,3)]), as.integer(data[9,seq(4,6)]), main="Overlap between replicates", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"))



