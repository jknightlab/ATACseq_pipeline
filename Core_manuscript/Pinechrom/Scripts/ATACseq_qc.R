#Pipeline for qc
library(ggplot2)

args<-commandArgs(TRUE)

peaks_length_hist <- args[1]
peaks_width_hist <- args[2]
stats <- args[3]
fragment_insertSize_hist <- args[4]
output_filename <- args[5]

## Reading in the data
combined.table <- read.delim(stats, header=TRUE)
chromosome.names <- combined.table[,1]

chromosome.names <- factor (chromosome.names, levels = chromosome.names)
peak.per.chr <- combined.table[,2]
reads.mapping.peaks.per.chr <- combined.table[,3]
reads.per.chr <- combined.table[,4]
reads.mapping.offpeaks.same.length.per.chr <- combined.table[,5]

peak.width.table <- read.delim(peaks_width_hist, header=FALSE)

peak.length.table <- read.delim(peaks_length_hist, header=FALSE)

insert.size.table <- read.delim(fragment_insertSize_hist, header=TRUE)
# Converting read counts to "normalized" read counts
insert.size.table$Normalised.read.density.1000 <- 1000*insert.size.table$All_Reads.fr_count/sum(insert.size.table$All_Reads.fr_count)


## ===================================================
## 
##                Creating the QC plots

## 
## ===================================================

# Periodicity plot

plot1 <- ggplot(insert.size.table,
    aes(x=insert_size, y=Normalised.read.density.1000)) + theme_bw() +
    geom_line(stat="identity", colour = "red", alpha=1) +
    scale_x_continuous(name="Fragment length (bp)",
    breaks=seq(0,800,by=200))+
    scale_y_continuous(name="Normalised insert size density (x10^-3)") + 
    ggtitle("Normalised insert size density")+
    theme(plot.title = element_text(size=20,face="bold"),
        legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15),
        axis.title.x = element_text(face="bold",size=18),
        axis.text.y = element_text(hjust = 1, size=15),
        axis.title.y = element_text(face="bold",size=18))

# Periodicity plot on normalized data

plot2 <-ggplot(insert.size.table, aes(x=insert_size, y=All_Reads.fr_count)) + theme_bw() +
    geom_line(stat="identity", colour = "red") + 
    scale_y_log10(name="Normalised read density") + 
    scale_x_continuous(name="Fragment length (bp)", breaks=seq(0,1000,by=200))+
    ggtitle("Normalised insert size density\n  (y axis in log10 scale)")+
    theme(plot.title = element_text(size=20,face="bold"),
        legend.title=element_blank(), axis.text.x = element_text(hjust = 1, size=15),
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15),
        axis.title.y = element_text(face="bold",size=18))

# Peak width distribution

x_axis_breaks <- length(seq(from = 200, 
    to = max(peak.width.table$V1), 
    by = 100))
plot3<-hist (peak.width.table$V1, 
    breaks=x_axis_breaks,col="blue", 
    main="Peak width distribution",
    xlab="length (bp)", ylab="Number of peaks")

# Peak length distribution

# missing plot4


# Number of peaks per chromosome

plot5 <- ggplot(data=combined.table, 
    aes(x=combined.table$Chromosome,y=combined.table$Called_peaks,
    fill=combined.table$Chromosome)) + geom_bar(stat="identity") +
    xlab("Chromosome") + ylab("Number of peaks") + 
    ggtitle("Number of peaks per chromosome")+
    theme(plot.title = element_text(size=20,face="bold"),
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=15),
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15),
        axis.title.y = element_text(face="bold",size=18))


# Correlation between the number of peaks and the chromosome length

# Do not include?


# Q1 -- number of peaks per chromosome divided by number of reads per chromosome

table.num_peaks_per_mapped_reads_per_chrom <- data.frame(
    chromosome=chromosome.names,
    peaks.per.reads=((peak.per.chr/reads.per.chr)*100))

plot6 <- ggplot(data=table.num_peaks_per_mapped_reads_per_chrom,
    aes(x=table.num_peaks_per_mapped_reads_per_chrom$chromosome,
    y=table.num_peaks_per_mapped_reads_per_chrom$peaks.per.reads,
    fill=table.num_peaks_per_mapped_reads_per_chrom$chromosome))+
    geom_bar(stat="identity")+
    xlab("Chromosome") + ylab("Number of peaks per number of PE reads (%)") +
    ggtitle("Number of peaks per number of\n PE reads (%) in each chromosome\n")+
    theme(plot.title = element_text(size=16,face="bold"),
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=15),
        axis.title.x = element_text(face="bold",size=15),
        axis.text.y = element_text(hjust = 1, size=15),
        axis.title.y = element_text(face="bold",size=14))


# Q2 -- number of reads mapping to peaks dvided by number of reads mapped to the entire chromosome

table.num_reads_in_peaks_per_total_reads_in_chr <- data.frame(
    chromosome=chromosome.names,
    reads.in.peaks.per.total.reads=((
    reads.mapping.peaks.per.chr/reads.per.chr)*100))

plot7 <- ggplot(data=table.num_reads_in_peaks_per_total_reads_in_chr,
    aes(x=table.num_reads_in_peaks_per_total_reads_in_chr$chromosome,
    y=table.num_reads_in_peaks_per_total_reads_in_chr$reads.in.peaks.per.total.reads,
    fill=table.num_reads_in_peaks_per_total_reads_in_chr$chromosome)) + geom_bar(stat="identity") +
    xlab("Chromosome") + ylab("Number of PE reads in peaks\n per total number fragments (%)") +
    ggtitle("Number of PE reads in peaks per number of\n PE reads in each chromosome")+
    theme(plot.title = element_text(size=20,face="bold"),
        legend.title=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size=15),
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15),
        axis.title.y = element_text(face="bold",size=18))


# Q4 -- number of reads mapping to peaks divided by number of reads mapping off-peaks to the regions of peak length

table.reads.in.peaks.per.reads.off.peak.same.length <- data.frame (
    chromosome=chromosome.names,
    reads.in.peaks.per.reads.off.peaks=(
    reads.mapping.peaks.per.chr/reads.mapping.offpeaks.same.length.per.chr))

plot8 <- ggplot(data=table.reads.in.peaks.per.reads.off.peak.same.length,
    aes(x=table.reads.in.peaks.per.reads.off.peak.same.length$chromosome,
    y=table.reads.in.peaks.per.reads.off.peak.same.length$reads.in.peaks.per.reads.off.peaks,
    fill=table.reads.in.peaks.per.reads.off.peak.same.length$chromosome)) +
    geom_bar(stat="identity") +
    xlab("Chromosome") + 
    ylab("Fragments fold enrichment over off peak regions") + 
    ggtitle("Fold enrichment of fragments in peaks\n versus off-peak regions") +
    theme(plot.title = element_text(size=20,face="bold"),
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=15),
        axis.title.x = element_text(face="bold",size=20),
        axis.text.y = element_text(hjust = 1, size=15),
        axis.title.y = element_text(face="bold",size=18))



# Reporting the plots

pdf (output_filename)
plot1
plot2

x_axis_breaks <- length(seq(from = 200,
    to = max(peak.width.table$V1),
    by = 100))
plot3 <- hist (peak.width.table$V1,
    breaks=x_axis_breaks,col="midnightblue",
    main="Peak width distribution",
    xlab="width (bp)", ylab="Number of peaks",
    xlim=c(0,5000))

x_axis_breaks <- length(seq(from = 200,
    to = max(peak.length.table$V1),
    by = 100))
plot4 <- hist (peak.length.table$V1,
    breaks=x_axis_breaks,col="darkslateblue",
    main="Peak height distribution",
    xlab="height (bp)", ylab="Number of peaks",
    xlim=c(0,10000))

plot5
plot6
plot7
plot8
dev.off()


