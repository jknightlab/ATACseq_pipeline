#!/bin/bash

#$ -N intrNPFs
#$ -P jknight.prjc -q short.qc
#$ -o stdout_intrsctNPFs -e sterr_intrsctNPFs -j y
#$ -cwd -V

### Generate Consensus Peaks from Three NPF Files
## JHendry, 2016/12/20

## Idea is take folder containing three .npf files corresponding to
## peaks either from MACS2, Fseq, Hotspot, representing different
## replicates from one sample type (cell/preparation method) and
## determine the peaks called in all three replicates.
## Standard for overlap is 1bp.
## Afterwards, report the 5' and 3' boundries of the intersect, 
## union, and median across the 3 replicates for the overlapping peaks.
##
## Schematic Eg.
##  Rep1            !===========:
##  Rep2                 *============!
##  Rep3               :======*
##
##  intersect            *====*  
##  median             :========:  
##  union           !=================!
##
## Run:
##  qsub run_intrDir-npf.sh <dir-with-replicate-npfs>
##  Requires Rscipt generateIntrsctUnnMdn.R
## Output includes:
##  - .all.bed
##  - .intrsct.bed
##  - .union.bed
##  - .median.bed


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"

### Define input directory, trim .npf to first three columns
# Below is modified such that all temporary files (.temp.)
# are produced within 'repDir'. This allows for run_intrDir-npf.sh
# to be run simultaneously on multiple folders without
# risking mixing of temporary files in the home directory.

repDir=$1
arrReps=($(ls $repDir)) # /* not included, and so only file names (NOT path)

cd $repDir

cat ${arrReps[0]} | awk '{print $1 "\t" $2 "\t" $3}' > rep1.temp.npf
cat ${arrReps[1]} | awk '{print $1 "\t" $2 "\t" $3}' > rep2.temp.npf
cat ${arrReps[2]} | awk '{print $1 "\t" $2 "\t" $3}' > rep3.temp.npf

echo "Input Directory: "$repDir
echo "Replicate 1 File: "${arrReps[0]}
echo "  No. Peaks: "`cat rep1.temp.npf | wc -l`
echo "Replicate 2 File: "${arrReps[1]}
echo "  No. Peaks: "`cat rep2.temp.npf | wc -l`
echo "Replicate 3 File: "${arrReps[2]}
echo "  No. Peaks: "`cat rep3.temp.npf | wc -l`
echo ""

### Serially intersect replicates 1, 2 and 3
bedtools intersect \
-a rep1.temp.npf \
-b rep2.temp.npf \
> rep1-rep2.temp.intrsct

echo "Replicate 1 v. 2 interescting peaks:"`cat rep1-rep2.temp.intrsct | wc -l`

bedtools intersect \
-a rep1-rep2.temp.intrsct \
-b rep3.temp.npf \
> rep1-rep2-rep3.temp.intrsct

echo "Three-way intersecting peaks:" `cat rep1-rep2-rep3.temp.intrsct | wc -l`
echo ""

bedtools intersect -wb \
-a rep1-rep2-rep3.temp.intrsct \
-b rep1.temp.npf rep2.temp.npf rep3.temp.npf \
-names rep1 rep2 rep3 \
> rep-all.temp.intrsct

### Grab overlapping peaks for each relicate
cat rep-all.temp.intrsct | awk '$4 == "rep1"' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7}' > rep1-all.temp.intrsct
cat rep-all.temp.intrsct | awk '$4 == "rep2"' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7}' > rep2-all.temp.intrsct
cat rep-all.temp.intrsct | awk '$4 == "rep3"' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7}' > rep3-all.temp.intrsct

### Merge intersection and overlapping rep1, rep2 and rep3 peaks into one file by row
paste rep1-rep2-rep3.temp.intrsct rep1-all.temp.intrsct rep2-all.temp.intrsct rep3-all.temp.intrsct > reps-intrsct.temp.all

cd .. # return to home directory for pathing in R script

### Generate intersect, median, and union consensus boundries
module load R/3.1.3
Rscript generate-IntrsctUnnMdn.R reps-intrsct.temp.all $repDir

intrsctFile=$(echo $repDir"/"$repDir".intrsct.bed")
medFile=$(echo $repDir"/"$repDir".median.bed")
unnFile=$(echo $repDir"/"$repDir".union.bed")


echo "No. Consensus Peaks, Pre-merge: "
echo "  Intersect: "`cat $intrsctFile | wc -l`
echo "  Median: "`cat $medFile | wc -l`
echo "  Union:  "`cat $unnFile | wc -l`

### Before merging, must ensure .bed files are sorted by chr then start
intrsctFileSorted=$(echo $intrsctFile".sorted")
medFileSorted=$(echo $medFile".sorted")
unnFileSorted=$(echo $unnFile".sorted")

sort -k1.4,1.5V -k2,2n $intrsctFile > $intrsctFileSorted
sort -k1.4,1.5V -k2,2n $medFile > $medFileSorted
sort -k1.4,1.5V -k2,2n $unnFile > $unnFileSorted 

bedtools merge -i $intrsctFileSorted > $intrsctFile
bedtools merge -i $medFileSorted > $medFile
bedtools merge -i $unnFileSorted > $unnFile

echo "No. Consensus Peaks, Post-merge: "
echo "  Intersect: "`cat $intrsctFile | wc -l`
echo "  Median: "`cat $medFile | wc -l`
echo "  Union:  "`cat $unnFile | wc -l`

### Remove temporary files
rm $repDir/*.sorted
rm $repDir/*.temp.*


echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""

