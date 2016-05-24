f_def_input <- read.table ("fseq_chosen.counts.txt", header=FALSE)
f_cho_input <- read.table ("fseq_default.counts.txt", header=FALSE)
m_def_input <- read.table ("macs.default.counts.txt", header=FALSE)
m_cho_input <- read.table ("macs.chosen.counts.txt", header=FALSE)

f_def <- density(as.numeric(f_def_input[,3]))
f_cho <- density(as.numeric(f_cho_input[,3]))
m_def <- density(as.numeric(m_def_input[,3]))
m_cho <- density(as.numeric(m_cho_input[,3]))


par(mfrow=c(1,2))
plot(f_cho, main="Peak intensity distribution", cex.axis=2, cex.lab=2, xlab="Normalized peak intensity", ylab="", cex.main=3, xlim=c(0,100))

polygon(f_def, border="black", col=adjustcolor("deepskyblue",alpha.f=0.2))
polygon(f_cho, border="black", col=adjustcolor("darkblue",alpha.f=0.2))
polygon(m_def, border="black", col=adjustcolor("chartreuse3",alpha.f=0.2))
polygon(m_cho, border="black", col=adjustcolor("darkgreen",alpha.f=0.2))

legend(70, 0.25, c("FSeq default", "FSeq chosen", "Macs2 default", "Macs2 chosen"), lty=c(1,1), lwd=c(3.5,3.5),col=c("deepskyblue", "darkblue", "chartreuse3", "darkgreen"), cex=1.5)


plot(f_cho, main="Peak intensity distribution", cex.axis=2, cex.lab=2, xlab="Normalized peak intensity", ylab="", cex.main=3, xlim=c(0,30))

polygon(f_def, border="black", col=adjustcolor("deepskyblue",alpha.f=0.2))
polygon(f_cho, border="black", col=adjustcolor("darkblue",alpha.f=0.2))
polygon(m_def, border="black", col=adjustcolor("chartreuse3",alpha.f=0.2))
polygon(m_cho, border="black", col=adjustcolor("darkgreen",alpha.f=0.2))

legend(21, 0.25, c("FSeq default", "FSeq chosen", "Macs2 default", "Macs2 chosen"), lty=c(1,1), lwd=c(3.5,3.5),col=c("deepskyblue", "darkblue", "chartreuse3", "darkgreen"), cex=1.5)

