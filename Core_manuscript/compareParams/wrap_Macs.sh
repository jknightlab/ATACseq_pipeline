#/bin/bash

#$ -N wrapFseq
#$ -P jknight.prjc -q short.qc
#$ -o stdout_wrapFseq -e sterr_wrapFseq.log -j y
#$ -cwd -V

### Run MACS2 Across a Range of User-specified  Parameter Values
## JHendry, 2016/12/01
##
## This script is a wrap for run_Macs.sh, which is called
## across or a range of user-specified parameter values.
## Follows the same scheme as wrap_Fseq.sh, but for MACS2 
## parameters, --q-value and --extsize. I.e., ranges are
## specified as follows:
##
##                minVal-maxVal:stepSize
## e.g.                  50-500:50
## 
## One important feature has been added to the MACS2
## wrapper. Since --q-value specifies a FDR rate, it is
## impractical to have a fixed step size (as interesting
## FDRs are often seperated by orders of magnitude). 
## This has been addressed by allowing the user to
## input "st" as the step size for the q-value, which
## produces q-values as follows at .1 and .5 between
## the specified minVal-maxVal:
##
## e.g.             0.001-0.1:st
## produces  0.001 0.005 0.01 0.05 0.1
   
## Inputs:
## $1 ---> Sample (fresh, frozen, fixed3d, fixed7d)
## $2 ---> q-value Range (or fixed value)
## $3 ---> extsize Range (or fixed value)


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"

sampInput=$1
qInput=$2
extInput=$3

### Evolve inputs into appropriate arrays
# Format for length and threshold inputs:
# start-end:step

echo "Sample Input:"$sampInput

echo "Q-value Input:"$qInput
qStep=${qInput#*:}
qRange=${qInput%:*}
qStart=${qRange%-*}
qEnd=${qRange#*-}
echo "Q-value Start:"$qStart
echo "Q-value End:"$qEnd
echo "Q-value Step:"$qStep

if [ "$qStep" != "st" ]; then
  qValues=$(seq -w $qStart $qStep $qEnd)
elif [ "$qStep" == "st" ]; then
  ## Somewhat complicated operation to generate nice q-value thresholds
  ## 0.001 0.005 0.01 0.05 ... etc.
  echo "Special Step Option Detected."
  qAll=($(seq $qStart $qStart $qEnd) $qEnd)
  qNumOrdMag=$[${#qStart} - ${#qEnd} + 1]
  echo "Q-values span" $qNumOrdMag "orders of magnitude."
  i=0
  for ord in `seq $qNumOrdMag`; do
    ordMag=$[10**($ord - 1)]
    ##echo ${qAll[$ordMag - 1]} 
    qValues[i]=`echo ${qAll[$ordMag - 1]}`
    (( i++ ))
    ##echo ${qAll[$[5*$ordMag]-1]}
    qValues[i]=`echo ${qAll[$[5*$ordMag]-1]}`
    (( i++ ))
  done
  qValues=${qValues[@]}
fi
numQs=$(echo $qValues | wc -w)
echo "Therefore, generating "$numQs "q-values:"
echo $qValues
echo ""

echo "Extsize Input:"$extInput
extStep=${extInput#*:}
extRange=${extInput%:*}
extStart=${extRange%-*}
extEnd=${extRange#*-}
echo "Extsize Start:"$extStart
echo "Extsize End:"$extEnd
echo "Extsize Step:"$extStep
extValues=$(seq -w $extStart $extStep $extEnd)
numExt=$(echo $extValues | wc -w)
echo "Therefore, generating" $numExt "extsize values:"
echo $extValues
echo ""
echo "No. Jobs to Submit:" `echo $[$numQs*$numExt]`
echo "----------------------------------------------------------------------"
echo ""
i=1
for qVal in $qValues; do
  for extVal in $extValues; do
    echo $i "Submitting: qsub run_Fseq-parallel.sh" $sampInput $qVal $extVal
    qsub run_Macs.sh $sampInput $qVal $extVal
    (( i++ ))
  done
done
echo ""
echo "----------------------------------------------------------------------"
echo "No. of jobs successfully submitted:" $i

echo ""
echo ""
echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
