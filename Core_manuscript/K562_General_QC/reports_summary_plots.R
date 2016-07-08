data <- read.table("reports_summary.txt", header=TRUE, sep='\t')

# ============ BOXPLOT ============

pdf("K562.three_repl.barplot.QC.pdf")
barplot(as.integer(data[1,]), main="Number of mapped reads after filtering", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[2,]), main="Average fragment width", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[3,]), main="Number of called peaks", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[4,]), main="Number of peaks per chromosome", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[5,]), main="Number of peaks normalized by total mapped reads", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[6,]), main="Number of peaks normalized by reads mapped per chr", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[7,]), main="Reads in peaks normalized by total mapped reads", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[8,]), main="Signal to noise ratio", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
barplot(as.integer(data[9,]), main="Average peak width", names.arg=c("fresh", "", "", "frozen", "", "", "fix_3d", "", "", "fix_7d", "", ""), col=c("lightsteelblue1", "lightsteelblue1", "lightsteelblue1", "mediumpurple", "mediumpurple", "mediumpurple", "darkslateblue", "darkslateblue", "darkslateblue", "royalblue4", "royalblue4", "royalblue4"))
dev.off()

# ============ BOXPLOT ============

pdf("K562.three_repl.boxplot.QC.pdf")
boxplot(as.integer(data[1,seq(1,3)]), as.integer(data[1,seq(4,6)]), as.integer(data[1,seq(7,9)]), as.integer(data[1,seq(10,12)]), main="Number of mapped reads after filtering", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[2,seq(1,3)]), as.integer(data[2,seq(4,6)]), as.integer(data[2,seq(7,9)]), as.integer(data[2,seq(10,12)]), main="Average fragment width", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[3,seq(1,3)]), as.integer(data[3,seq(4,6)]), as.integer(data[3,seq(7,9)]), as.integer(data[3,seq(10,12)]), main="Number of called peaks", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[4,seq(1,3)]), as.integer(data[4,seq(4,6)]), as.integer(data[4,seq(7,9)]), as.integer(data[4,seq(10,12)]), main="Number of peaks per chromosome", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[5,seq(1,3)]), as.integer(data[5,seq(4,6)]), as.integer(data[5,seq(7,9)]), as.integer(data[5,seq(10,12)]), main="Number of peaks normalized by total mapped reads", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[6,seq(1,3)]), as.integer(data[6,seq(4,6)]), as.integer(data[6,seq(7,9)]), as.integer(data[6,seq(10,12)]), main="Number of peaks normalized by reads mapped per chr", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[7,seq(1,3)]), as.integer(data[7,seq(4,6)]), as.integer(data[7,seq(7,9)]), as.integer(data[7,seq(10,12)]), main="Reads in peaks normalized by total mapped reads", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[8,seq(1,3)]), as.integer(data[8,seq(4,6)]), as.integer(data[8,seq(7,9)]), as.integer(data[8,seq(10,12)]), main="Signal to noise ratio", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
boxplot(as.integer(data[9,seq(1,3)]), as.integer(data[9,seq(4,6)]), as.integer(data[9,seq(7,9)]), as.integer(data[9,seq(10,12)]), main="Average peak width", names=c("fresh", "frozen", "fix_3d", "fix_7d"), col=c("lightsteelblue1", "mediumpurple", "darkslateblue", "royalblue4"))
dev.off()














