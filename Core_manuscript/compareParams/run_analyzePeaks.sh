#/bin/bash

#$ -N aPeaks
#$ -P jknight.prjc -q short.qc
#$ -o stdout_aPeaks.log -e sterr_aPeaks.log -j y

#$ -cwd -V

### Compute Basic Statistics On a Directory of Bed Files
## JHendry, 2016/02/03
##
##
## This script runs run_analyzePeaks-n.R on a user-
## specified directory. The output includes statistics
## on the total number of peaks, number of basepairs,
## peak width and peak densities, and is deposited in
## a folder within the comparisons directory.


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"

### Define Input Directory
bedDir=$1

echo "Input File Directory:"$bedDir
echo "Running run_analyzePeaks-n.R"
module load R/3.1.3
Rscript run_analyzePeaks-n.R $bedDir
echo "run_analyzePeaks-n.R complete."

echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
