#/bin/bash

#$ -N wrapFseq
#$ -P jknight.prjc -q short.qc
#$ -o stdout_wrapFseq -e sterr_wrapFseq.log -j y
#$ -cwd -V

### Run F-seq Across a Range of Length and Threshold Parameter Values
## JHendry, 2016/12/01
##
##
## This script wraps run_Fseq-parallel.sh, allowing for F-seq to
## be run on a given sample type (fresh, frozen, fixed) over a 
## range of user-specified length (-l) and threshold (-t)
## parameters. To define a range of values, the following 
## notation is used:
##
##                 minVal-maxVal:stepSize
## e.g.            100-2000:100 or 2-16:2
##
## When ranges are specified for both the length and threshold 
## parameters, F-seq is run across all possible combinations.
## This script also supports fixed value, e.g. if one is
## interested in running a range threshold values against
## a single length value.
##
## A series of directories is created reflecting the range
## of parameter values across which F-seq was run.
## See run_Fseq-parallel.sh for details.

## Inputs:
## $1 ---> Sample Input (fresh, frozen, fixed3d, fixed7d)
##    * Note that bed files for all replicates belonging to the
##      specified sample are processed
## $2 ---> Length Range (or fixed value)
##    * e.g. 200-2000:200, or 400
## $3 ---> Threshold Range (or fixed value)
##    * e.g. 2-16:2, or 8


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"

sampInput=$1
lengthInput=$2
threshInput=$3

### Evolve inputs into appropriate arrays
# Format for length and threshold inputs:
# start-end:step

echo "Sample Input:"$sampInput
sampDirs=$(ls -d bedFiles/$sampInput*)
numSamps=$(echo $sampDirs | wc -w)
echo "Sample Directories:"
echo $sampDirs
echo ""

echo "Length Input:"$lengthInput
lengthStep=${lengthInput#*:}
lengthRange=${lengthInput%:*}
lengthStart=${lengthRange%-*}
lengthEnd=${lengthRange#*-}
echo "Length Start:"$lengthStart
echo "Length End:"$lengthEnd
echo "Length Step:"$lengthStep
lengthValues=$(seq -w $lengthStart $lengthStep $lengthEnd)
numLength=$(echo $lengthValues | wc -w)
echo "Therefore, generating "$numLength "length values:"
echo $lengthValues
echo ""

echo "Threshold Input:"$threshInput
threshStep=${threshInput#*:}
threshRange=${threshInput%:*}
threshStart=${threshRange%-*}
threshEnd=${threshRange#*-}
echo "Threshold Start:"$threshStart
echo "Threshold End:"$threshEnd
echo "Threshold Step:"$threshStep
threshValues=$(seq -w $threshStart $threshStep $threshEnd)
numThresh=$(echo $threshValues | wc -w)
echo "Therefore, generating" $numThresh "threshold values:"
echo $threshValues
echo ""
echo "No. Jobs to Submit:" `echo $[$numSamps*$numLength*$numThresh]`
echo "----------------------------------------------------------------------"
echo ""

for samp in $sampDirs;
do
  for lengthVal in $lengthValues;
  do
    for threshVal in $threshValues;
    do
      echo "Submitting: qsub run_Fseq-parallel.sh" $samp $lengthVal $threshVal
      qsub run_Fseq-parallel.sh $samp $lengthVal $threshVal
    done
  done
done
echo ""
echo "----------------------------------------------------------------------"
echo "No. of jobs successfully submitted:" `qstat | wc -l`
echo ""
echo ""

echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
