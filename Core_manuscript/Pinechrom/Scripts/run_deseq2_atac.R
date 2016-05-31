library(DESeq2)

args<-commandArgs(TRUE)

input_counts <- args[1]
input_colnames <- args[2]
output_header <- args[3]

# Generating filenames
filename_count_table <- paste(output_header, ".normalized_counts.txt", sep="")
filename_de_plots <- paste(output_header, ".DE_plots.pdf", sep="")
filename_de_results <- paste(output_header, ".DE_results.txt", sep="")


## -------------------------------
## 
##     Reading the input data          
## 
## -------------------------------

# Input -- a file containing count data and info about conditions
countData <- read.table(input_counts, sep='\t', header=TRUE)
colData <- read.table(input_colnames, sep='\t', header=TRUE)

print ("Data was successfully loaded.")

# Converting input data into a deseq object
dds <- DESeqDataSetFromMatrix(countData = countData,
    colData = colData,
    design = ~ condition)

# Filtering to remove low counts
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Tell deseq which samples to compaire
dds$condition <- factor(dds$condition, 
    levels=c("untreated","treated"))

print ("DESeq2 object was created.")

## -------------------------------
## 
##       Running DE analysis
## 
## -------------------------------

# Running DE analysis
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=T)
write.table(as.data.frame(normalized_counts),
    file=filename_count_table, quote=FALSE, sep='\t')

print ("Differential expression analysis was finished.")

# Extracting the results matrix
# The results are shown for treated vs untreated, fold change log2 (treated/untreated)
res <- results(dds)

# alpha is the FDR cutoff, we can change it
res05 <- results(dds, alpha=0.05)


## -------------------------------
## 
##           Extra plots          
## 
## -------------------------------
library("ggplot2")
library("pheatmap")
library("RColorBrewer")

print ("Started to generate plots.")
pdf(filename_de_plots)

# plot MA
plotMA(res, main="DESeq2", ylim=c(-2,2))
resMLE <- results(dds, addMLE=TRUE)
plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))

# plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# plot counts in ggplot
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",
    returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) +
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_log10(breaks=c(25,100,400))

# plot heatmap of counts
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
    cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
    cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
    cluster_cols=FALSE, annotation_col=df)

# plot sample to sample distance
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors)

# plot PCA
plotPCA(rld, intgroup=c("condition", "type"))
data <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

print ("Finished generating plots.")

## -------------------------------
## 
##        Exporting results          
## 
## -------------------------------

print ("Preparing to print the list of differentially expressed regions.")

# Sort results on adjusted pvalue
resOrdered <- res[order(res$padj),]

# Write output in a csv file
write.table(as.data.frame(resOrdered), file=filename_de_results, 
    quote=FALSE, sep='\t')

print ("Finished printing DESeq2 output.")

