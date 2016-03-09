par(mfrow=c(3,3))

boxplot(as.numeric(data[1,seq(1,3)]), as.numeric(data[1,seq(4,6)]), main="Number of mapped reads after filtering", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(1000000, 60000000), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[3,seq(1,3)]), as.numeric(data[3,seq(4,6)]), main="Total number of peaks", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(40000, 180000), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[4,seq(1,3)]), as.numeric(data[4,seq(4,6)]), main="Number of peaks per chromosome", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(2000, 10000), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[5,seq(1,3)]), as.numeric(data[5,seq(4,6)]), main="Number of peaks normalized by total mapped reads", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(50, 250), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[7,seq(1,3)]), as.numeric(data[7,seq(4,6)]), main="Reads in peaks normalized by total mapped reads", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(200, 15000), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[8,seq(1,3)]), as.numeric(data[8,seq(4,6)]), main="Signal to noise ratio", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(2, 6), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[9,seq(1,3)]), as.numeric(data[9,seq(4,6)]), main="Average peak width", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(150, 300), cex.lab=1.5, cex.axis=3, cex.main=2.5)

boxplot(as.numeric(data[10,seq(1,3)]), as.numeric(data[10,seq(4,6)]), main="Overlap between replicates", names=c("dnase", "atac"), col=c("mediumpurple", "aquamarine4"), ylim=c(50, 70), cex.lab=1.5, cex.axis=3, cex.main=2.5)

