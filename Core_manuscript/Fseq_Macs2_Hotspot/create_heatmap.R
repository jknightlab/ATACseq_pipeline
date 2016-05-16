library(reshape2)
library(gplots)

args <- commandArgs(T)
input_parameter <- args[1]
output_name <- args[2]
title <- args[3]

nba <- read.table(input_parameter, header=TRUE, sep='\t')

nba_matrix <- data.matrix(nba)

my_palette <- colorRampPalette(c("royalblue", "deeppink4"))(n = 50)

pdf(output_name)
heatmap.2(nba_matrix, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', density.info="none", main=title, tracecol=NA, col=my_palette, labRow=rownames(nba_matrix), margins=c(10, 10), cexRow=2, cexCol=2, keysize = 1)
dev.off()



heatmap.2(nba_matrix, Rowv=FALSE,Colv=FALSE, 
          col = my_palette,
          scale="none",
          margins=c(10, 4),
          trace='none', 
          symkey=FALSE, 
          symbreaks=FALSE, 
          dendrogram='none',
          cexRow=2,
          cexCol=2,
          density.info='none', 
          denscol=tracecol,
          key=2,
          key.xlab = "",
          key.title = "",
          key.xtickfun = NULL,
          key.par=list(mar=c(3.5,0,3,0)),
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1.5, 5), lwid=c(1, 10, 1))


