#!/bin/bash

#$ -N compSamps
#$ -P jknight.prjc -q short.qc
#$ -o stdout_compSamps -e sterr_compSamps.log -j y
#$ -cwd -V

### Perform Basic Comparison Between Two Samples 
## JHendry, 2016/01/14

## Idea is take consensus peak files from two samples and automatedly
## perform a series of basic analyses. Before run_compareSamps.sh
## can be run, the two samples in question must have already been
## intersected with run_intrsctSamps.sh, and the intersect file
## must exist in ./comparisons.
## The analyses performed include:
##  From compareSamples.R:
##  - Comparing global overlap in peak calls between two samples
##    by determining number of intersecting peaks that satisfy
##    a given threshold (#bps, reciprocal percent, both).
##  From analyzePeaks:R
##  - Side-by-side of total number of peaks, peaks/chromosome, length
##  and density distribution of peaks for both samples.

## Run:
##  qsub run_compareSamps-v2.sh <dir-samp1> <dir-samp2> <intrsct/median/union>
## Requires:
##  compareSamples.R
##  analyzePeaks.R
## Output:  	
##  (A) All files begin with the prefix <dir-samp1>-<dir-samp2> and are deposited
##  in a folder ./comparisons/<dir-samp1>-<dir-samp2>-<intrsct/median/union>.
##  (B) Between the prefix and suffix, the type of threshold used for
##  looking at overlap is specified, either:
##   - <prefix>-OverlapByBp-<suffix> : by absolute # bps in overlap
##   - <prefix>-OverlapByPer-<suffix> : by reciprocal percent of overlap
##   - <prefix>-OverlapByHq-<suffix> : by >20bp and reciprocal percent
##  (C) Suffix of files include:
##     -absTable.txt : absolute numbers for overlaps
##     -perTable.txt : percent values for overlaps (calc. from -absTable.txt)
##     .pdf : stacked barplot to visualize -absTable.txt
##  Also generated is comparisons of total number, length and density
##  of peaks, with suffixes:
##  _numPeaks.pdf : total number of peaks, number peaks/chromosome.
##  _widthPeaks.pdf : length (~width) distribution of peaks.
##  _densityPeaks.pdf : density distribution of peaks.

echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"

samp1Dir=$1
samp2Dir=$2
intrsctType=$3

samp1=$(echo $samp1Dir"/"$samp1Dir"."$intrsctType".bed")
samp2=$(echo $samp2Dir"/"$samp2Dir"."$intrsctType".bed")
intersectFile=$(echo "comparisons/"$samp1Dir"-"$samp2Dir"."$intrsctType".intrsct")


echo "Sample 1:" $samp1Dir
echo "Sample 2:" $samp2Dir
echo "Intersection Type:" $intrsctType
echo "Intersect File:" $intersectFile
echo ""
echo "Peaks in Sample 1:"`cat $samp1 | wc -l`
echo "Peaks in Sample 2:"`cat $samp2 | wc -l`
echo "Peaks in Intersect File:"`cat $intersectFile | wc -l`
echo ""
echo ""


### Load R module and run R scripts
module load R/3.1.3
Rscript sampCompare.R $samp1 $samp2 $intersectFile
Rscript analyzePeaks.R $samp1 $samp2

### Generate comparisons directory and move all outputs there
mkdir comparisons/$samp1Dir-$samp2Dir-$intrsctType
mv $samp1Dir-$samp2Dir* comparisons/$samp1Dir-$samp2Dir-$intrsctType


echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""

