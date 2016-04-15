#/bin/bash

#$ -N wrapNPFs
#$ -P jknight.prjc -q short.qc
#$ -o stdout_wrapNPFs -e sterr_wrapNPFs.log -j y
#$ -cwd -V

### Intersect NPF Files of Peak Calls From F-seq or MACS2
### Producing Consensus Peak Calls Across a Range of
### Algorithm Parameter Values
## JHendry, 2016/12/01
##
##
## This is a wrap for run_intrsctNPFs.sh -- it allows the 
## computation of consensus peak calls for NPF files
## across a range of either F-seq or MACS2 parameter
## values.  In essence, it was built to be run after 
## wrap_Fseq.sh or wrap_Macs.sh, in order to produce
## consensus peak calls for all of the NPF files
## produced by those scripts.
## The input is very similar to those scripts,
## except that the first command line argument is 
## what algorithm was run.
##  An example workflow:
##
## . wrap_Macs.sh fresh 0.001-0.1:st 50-500:50
## . wrap_intrsctNPFs.sh macs fresh 0.001-0.1:st 50-500:50
##
## Would produce consensus peak calls for the entire
## matrix of MACS2 parameter values defined above.

## Inputs:
## $1 ---> algorithm name (either "fseq" or "macs")
## $2 ---> sample type (fresh, frozen...)
## $3 ---> length or q-value input (dep. on specified algorithm)
## $4 ---> threshold or extsize input
 

echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"


algoInput=$1
if [ "$algoInput" == "fseq" ]; then
  echo "Algorithm Input:" $algoInput
  shift
  sampInput=$1
  lengthInput=$2
  threshInput=$3

  ### Evolve inputs into appropriate arrays
  # Format for length and threshold inputs:
  # start-end:step

  echo "Sample Input:" $sampInput
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
  echo "No. Jobs to Submit:" `echo $[$numLength*$numThresh]`
  echo "----------------------------------------------------------------------"
  echo ""
  i=1
  for lengthVal in $lengthValues;
  do
    for threshVal in $threshValues;
    do
      sampDirName=$(echo "fseql"$lengthVal"t"$threshVal"_"$sampInput)
      echo $i "Submitting: qsub run_intrsctNPFs.sh" $sampDirName
      qsub run_intrsctNPFs.sh $sampDirName
      (( i++ ))
    done
  done

  echo ""
  echo "----------------------------------------------------------------------"
  echo ""
  echo ""
elif [ "$algoInput" == "macs" ]; then
echo "Algorithm Input:" $algoInput
  shift
  sampInput=$1
  qInput=$2
  extInput=$3

  ### Evolve inputs into appropriate arrays
  # Format for length and threshold inputs:
  # start-end:step

  echo "Sample Input:" $sampInput

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
  i=0
  for qVal in $qValues;
  do
    for extVal in $extValues;
    do
      qValRmvZeroDot=${qVal##*.}
      sampDirName=$(echo "macsq"$qValRmvZeroDot"e"$extVal"_"$sampInput)
      echo $i "Submitting: qsub run_intrsctNPFs.sh" $sampDirName
      qsub run_intrsctNPFs.sh $sampDirName
      (( i++ ))
    done
  done

  echo ""
  echo "----------------------------------------------------------------------"
  echo ""
  echo ""
else
  echo "Algorithm input not recognized."
fi


echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
