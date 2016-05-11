macs_default <- read.table("macs.default.bed", sep='\t', header=FALSE)
macs_chosen <- read.table("macs.chosen.bed", sep='\t', header=FALSE)
fseq_chosen <- read.table("fseq_chosen.bed", sep='\t', header=FALSE)
fseq_fedault <- read.table("fseq_default.bed", sep='\t', header=FALSE)

f_def <- density(as.numeric(fseq_fedault[,3]-fseq_fedault[,2]))
f_cho <- density(as.numeric(fseq_chosen[,3]-fseq_chosen[,2]))
macs_def <- density(as.numeric(macs_default[,3]-macs_default[,2]))
macs_cho <- density(as.numeric(macs_chosen[,3]-macs_chosen[,2]))

plot(f_cho, main="Peak width distribution", cex.axis=2, cex.lab=2, xlab="", ylab="", cex.main=3)

polygon(f_def, border="black", col=adjustcolor("deepskyblue",alpha.f=0.2))
polygon(f_cho, border="black", col=adjustcolor("darkblue",alpha.f=0.2))
polygon(macs_def, border="black", col=adjustcolor("chartreuse3",alpha.f=0.2))
polygon(macs_cho, border="black", col=adjustcolor("darkgreen",alpha.f=0.2))

legend(1100, 0.0055, c("FSeq default parameters", "FSeq chosen parameters", "Macs2 default parameters", "Macs2 chosen parameters"), lty=c(1,1), lwd=c(4.5,4.5),col=c("deepskyblue", "darkblue", "chartreuse3", "darkgreen"), cex=2)


plot(f_cho, main="Peak width distribution", cex.axis=2, cex.lab=2, xlab="", ylab="", cex.main=3, xlim=c(0, 500))

polygon(f_def, border="black", col=adjustcolor("deepskyblue",alpha.f=0.2))
# abline(v=median(as.numeric(fseq_fedault[,3]-fseq_fedault[,2])), col="deepskyblue", lwd=3)

polygon(f_cho, border="black", col=adjustcolor("darkblue",alpha.f=0.2))
# abline(v=median(as.numeric(fseq_chosen[,3]-fseq_chosen[,2])), col="darkblue", lwd=3)

polygon(macs_def, border="black", col=adjustcolor("chartreuse3",alpha.f=0.2))
# abline(v=median(as.numeric(macs_default[,3]-macs_default[,2])), col="chartreuse3", lwd=3)

polygon(macs_cho, border="black", col=adjustcolor("darkgreen",alpha.f=0.2))
# abline(v=median(as.numeric(macs_chosen[,3]-macs_chosen[,2])), col="darkgreen", lwd=3)

legend(300, 0.0055, c("FSeq default parameters", "FSeq chosen parameters", "Macs2 default parameters", "Macs2 chosen parameters"), lty=c(1,1), lwd=c(4.5,4.5),col=c("deepskyblue", "darkblue", "chartreuse3", "darkgreen"), cex=2)

