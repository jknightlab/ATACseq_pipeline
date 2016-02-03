##### Convert reps-intrsct.all Into Intersect, Union, and Median Peak Boundries
## JHendry 2016/01/21
## Idea here is to take the reps-intrsct.all file, containing:
##	- intersecting regions boundries (col1-3)
##  - replicate 1 peak in intersection (col4-7)
##  - replicate 2 peak in intersection (col8-11)
##  - replicate 3 peak in intersection (col12-15)
## and generate three .bed files:
##	- sample.intrsct.bed, region of rep1/rep2/rep3 peaks intersection
##	- sample.union.bed, region of rep1/rep2/rep3 peaks union
##	- sample.median.bed, region defined by median start and end boundry of rep1/rep2/rep3


### Enable command line arguments
args <- commandArgs(T)

### Load reps-intrsct.all
# It's columns are defined below
colsAll <- c("intrsctChr", "intrsctStart", "intrsctEnd", "rep1", "rep1Chr", "rep1Start", "rep1End", "rep2", "rep2Chr", "rep2Start", "rep2End", "rep3", "rep3Chr", "rep3Start", "rep3End")
repsAll <- read.table(file = paste(args[2], "/", args[1], sep = ""), header = F, col.names = colsAll)

### Compute min, max, and median 5' and 3' boundries as needed
# Grab 5' (Start) boundry and 3' (End) boundry for all replicates
repStarts <- as.matrix(repsAll[ , c(6, 10, 14)])
repEnds <- as.matrix(repsAll[ , c(7, 11, 15)])

# Compute
repsAll <- within(repsAll, {
	minStart <- apply(repStarts, 1, min)
	maxEnd <- apply(repEnds, 1, max)
	medStart <- apply(repStarts, 1, median)
	medEnd <- apply(repEnds, 1, median)
})

### Generate and export intresect, union and median boundries as seperate .bed files
repsIntrsct <- repsAll[ , c("intrsctChr", "intrsctStart", "intrsctEnd")]
repsUnion <- repsAll[ , c("intrsctChr", "minStart", "maxEnd")]
repsMedian <- repsAll[ , c("intrsctChr", "medStart", "medEnd")]
repsAll <- within(repsAll, rm(minStart, maxEnd, medStart, medEnd))

write.table(repsAll, paste(args[2], "/", args[2], ".all.bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(repsIntrsct, paste(args[2], "/", args[2], ".intrsct.bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(repsUnion, paste(args[2], "/", args[2], ".union.bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(repsMedian, paste(args[2], "/", args[2], ".median.bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)

