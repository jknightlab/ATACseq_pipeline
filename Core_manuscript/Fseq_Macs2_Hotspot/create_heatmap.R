library(reshape2)
library(gplots)

args <- commandArgs(T)
input_parameter <- args[1]
output_name <- args[2]
title <- args[3]

nba <- read.table(input_parameter, header=TRUE, sep='\t')

nba_matrix <- data.matrix(nba)

my_palette <- colorRampPalette(c("royalblue", "deeppink4"))(n = 500)

pdf(output_name)
heatmap.2(nba_matrix, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', density.info="none", main=title, tracecol=NA, col=my_palette, labRow=rownames(nba_matrix), margins=c(10, 10), cexRow=2, cexCol=2, keysize = 1)
dev.off()

