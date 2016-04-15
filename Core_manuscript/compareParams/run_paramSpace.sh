#/bin/bash

#$ -N paramSp
#$ -P jknight.prjc -q short.qc
#$ -o stdout_paramSpace -e sterr_paramSpace.log -j y
#$ -cwd -V

### Analyze Basic Peak Statistics Across a Range of Algorithm
## Parameters, for MACS2 or F-seq
## JHendry, 2016/12/01

echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"


algoInput=$1
sampInput=$2
pOneInput=$3
pTwoInput=$4
consensusType=$5

### Want to slowly build a list of file names that satisfy the input conditions
if [ "$algoInput" == "fseq" ]; then
  echo "Input Algorithm:" $algoInput
  allSamps=$(ls -d fseql*t*_$sampInput)
elif [ "$algoInput" == "macs" ]; then 
  allSamps=$(ls -d macsq*e*_$sampInput)
else
  echo "Input Algorithm:" $algoInput
  echo "Algorithm type not recognized."
fi

### Extract Paramter 1 and 2 Min & Max
# pOne -- parameter 1
# pTwo -- parameter 2
pOneMin=${pOneInput%-*}
pOneMax=${pOneInput#*-}

pTwoMin=${pTwoInput%-*}
pTwoMax=${pTwoInput#*-}


### Extract All Samples Inside Paramter 1 and Parameter 2 Range
i=1
for curSamp in $allSamps;
do
  if [ "$algoInput" == "fseq" ]; then
    pOneVal=$(echo ${curSamp#*l} | cut -c1-4)
    pTwoVal=$(echo ${curSamp#*t} | cut -c1-2)
  elif [ "$algoInput" == "macs" ]; then
    pOneVal=$(echo ${curSamp#*q} | cut -c1-4)
    pTwoVal=$(echo ${curSamp#*e} | cut -c1-3)
  else
    echo "Algorithm type not recognized."
  fi
  if [ $pOneVal -ge $pOneMin -a $pOneVal -le $pOneMax ] && [ $pTwoVal -ge $pTwoMin -a $pTwoVal -le $pTwoMax ]
  then
    allSampsInParams[i]=$(echo $curSamp"/"$curSamp"."$consensusType".bed")
    ((i++))
  fi
done

### Input fseqInParams to Rscript
module load R/3.1.3
Rscript run_paramSpace-n.R `echo ${allSampsInParams[@]}`

### Generate comparisons directory and move all outputs there

if [ "$algoInput" == "fseq" ]; then
  outputDir=$(echo "comparisons/fseql"$pOneInput"t"$pTwoInput"_"$sampInput"-"$consensusType)
elif [ "$algoInput" == "macs" ]; then
  outputDir=$(echo "comparisons/macsq"$pOneInput"e"$pTwoInput"_"$sampInput"-"$consensusType)
else
  echo "Algorithm type not recognized."
fi

if [ ! -d "$outputDir" ]
then
  mkdir $outputDir
fi
mv paramSpace* $outputDir


echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
