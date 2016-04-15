#/bin/bash

#$ -N compAnnot
#$ -P jknight.prjc -q short.qc
#$ -o stdout_compAnnot -e sterr_compAnnot.log -j y
#$ -cwd -V

### Compare a Set of Sample Files to an Annotation Directory
## JHendry, 2016/02/09

# Eventually want several input options for sample files
# for now, just fseq /w range


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"

### Part 1. Determine Format Of Input
## A. Range
##  -r fseq fresh 200-1200 2-16 intrsct
## B. Directory
##  -d sampleSets/sampleDir
##
## ... and annotation directory:
##   -a annotation/ENCODEsegmentation

echo "Dectected $# input parameters."
while [ "$#" != "0" ]; do
  if [ "$1" == "-r" ]; then
    echo "Processing samples as parameter range (-r option)."
    runStyle="-r"
    algoInput=$2
    sampInput=$3
    pOneInput=$4
    pTwoInput=$5
    consensusType=$6
    echo "Inputs:"
    echo "Algorithm Input:" $algoInput
    echo "Sample Input:" $sampInput
      if [ "$algoInput" == "macs" ]; then
        echo "FDR Input (--qvalue):" $pOneInput
        echo "Extsize Input (--extsize):" $pTwoInput
      elif [ "$algoInput" == "fseq" ]; then
        echo "Length Input (-l):" $pOneInput
        echo "Threshold Input (-t):" $pTwoInput
      fi
    echo "Consensus Peak Type:" $consensusType
    echo ""
    shift 6
  elif [ "$1" == "-d" ]; then
    echo "Processing samples from named directory (-d option)."
    runStyle="-d"
    sampDir=$2
    echo "Directory:" $sampDir
    echo ""
    shift 2
  elif [ "$1" == "-a" ]; then
    echo "Found annotation directory (-a option)."
    annotDir=$2
    echo "Annotation Directory:" $annotDir 
    echo ""
    shift 2
  else
    echo "Unrecognized input parameter."
    echo "Skipping."
    echo ""
    shift
  fi
done

### Generate, For Either Range (-r) or Directory (-d) Input Method
## (1) sampFiles -- A List of Sample File Names
## (2) outputDir -- Directory To Dump bedtools intersect Outputs
if [ "$runStyle" == "-r" ]; then
  annotDirBase=${annotDir##*/}
  ### Want to slowly build a list of file names that satisfy the input conditions
  if [ "$algoInput" == "fseq" ]; then
    allSamps=$(ls -d fseql*t*_$sampInput)
    outputDir=$(echo "comparisons/fseql"$pOneInput"t"$pTwoInput"_"$sampInput"-"$consensusType"-vs-"$annotDirBase)  
  elif [ "$algoInput" == "macs" ]; then 
    allSamps=$(ls -d macsq*e*_$sampInput)
    outputDir=$(echo "comparisons/macsq"$pOneInput"e"$pTwoInput"_"$sampInput"-"$consensusType"-vs-"$annotDirBase)  
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


  # (1) Assign discovered files to sampFiles
  sampFiles=${allSampsInParams[@]}
  
  # (2) Generate ouputDir based on range (-r) inputs
  if [ ! -d "$outputDir" ]
  then
    echo "Generating output directory:" $outputDir
    echo ""
    mkdir $outputDir
    mkdir $outputDir"/sampleBeds"
    mkdir $outputDir"/intersectBeds"
  fi
elif [ "$runStyle" == "-d" ]; then
  ### Extract Samples from Sample Directory &
  # (1) Assign discovered files to sampFiles
  sampFiles=$(ls $sampDir/*)  
  
  # (2) Generate ouputDir based on directory name (-d)
  annotDirBase=${annotDir##*/}
  sampDirBase=${sampDir##*/}
  outputDir=$(echo "comparisons/"$sampDirBase"-vs-"$annotDirBase)
  if [ ! -d "$outputDir" ]
  then
    echo "Generating output directory:" $outputDir
    echo ""
    mkdir $outputDir
    mkdir $outputDir"/sampleBeds"
    mkdir $outputDir"/intersectBeds"
  fi
else
  echo "Unrecognized run style."
fi

### Find Annotation Files & Annotation Key
# Same method for both (-r) and (-d) options
annotFiles=$(ls $annotDir/*.bed)
annotKey=$(echo $annotDir"/annotKey.txt")


### Move Sample Files into $outputDir/sampleBeds
# Necessary for later analysis with compareToAnnot.R
cp $sampFiles $outputDir"/sampleBeds"

echo "Input Processing Complete."
numSamps=$(echo $sampFiles | wc -w)
echo "Discovered" $numSamps "sample files."
echo "...vs..."
numAnnots=$(echo $annotFiles | wc -w)
echo "Discovered" $numAnnots "annotation files:"
printf ${annotFiles[@]}
echo ""
echo ""


### Iteratively Intersect sampFiles vs. annotFiles
## 1.Compute Number of Intersections Expected Before
## 2. What type of intersection files should I build?
##  - all? intersect?
echo "Performing" $[$numSamps*$numAnnots] "intersections."
echo "--------------------------------------------------------------------"
i=1
for annotFile in $annotFiles;
do
  for sampFile in $sampFiles;
  do
    sampFileNameExt=${sampFile##*/}
    sampFileName=${sampFileNameExt%%.*}
    annotFileNameExt=${annotFile##*/}
    annotFileNameNoExt=${annotFileNameExt%%.*}
    # Grab Condensed Annotation File Name from annotKey
    annotFileName=$(cat $annotKey | awk -v val=$annotFileNameNoExt '$1 == val' | awk '{print $2}')
    intrsctFileName=$(echo $outputDir"/intersectBeds/"$sampFileName"-"$annotFileName".bed")
    
    onlyIntrsctFileName=$(echo $intrsctFileName".only.temp")
    allIntrsctFileName=$(echo $intrsctFileName".all.temp")
    annotIntrsctFileName=$(echo $intrsctFileName".annot.temp")
    sampIntrsctFileName=$(echo $intrsctFileName".samp.temp")
    
    echo $i":" 
    echo "Intersecting" $sampFileNameExt "and" $annotFileNameExt "..."

    ### Primary Intersect
    # Here the output is a three column .bed
    # chr intrsctStart intrsctEnd
    echo "Finding intersections..."
    bedtools intersect \
    -a $annotFile \
    -b $sampFile \
    > $onlyIntrsctFileName
    echo "Getting original peaks..."
    # Secondary Intersect
    # Here, I intersect the intersect bed with the original files
    # and use -wb option, to get original peak/annot
    bedtools intersect \
    -a $onlyIntrsctFileName \
    -b $sampFile $annotFile \
    -names $sampFileName $annotFileName \
    -wb \
    > $allIntrsctFileName
    echo "Extracting Sample & Annotation peaks..."
    awk -v annot=$annotFileName '{if($4==annot){print $4"\t"$5"\t"$6"\t"$7}}' $allIntrsctFileName > $annotIntrsctFileName
    awk -v samp=$sampFileName '{if($4==samp){print $4"\t"$5"\t"$6"\t"$7}}' $allIntrsctFileName > $sampIntrsctFileName
    paste $onlyIntrsctFileName $sampIntrsctFileName $annotIntrsctFileName > $intrsctFileName
    echo " Output File:" ${intrsctFileName##*/}
    echo ""
    ((i++)) 
   done
done
echo "--------------------------------------------------------------------"
echo "Completed" $[$i - 1] "intersections."
echo ""
echo ""

### At This Point, should have:
# In outputDir
# 1. sampleBeds -- directory of all samples
# 2. intersectBeds -- directory of all intersects
# As well as the annotDir
# These are all fed to compareToAnnot.R, along with the runType
# e.g.
# Rscript compareToAnnot.R <outputDir> <annotDir> <runType>
# & "intersectBeds" & "sampleBeds" are assumed.

### Load R module and run R scripts
echo "Running run_compareToAnnot-n.R"
module load R/3.1.3
Rscript run_compareToAnnot-n.R $outputDir $annotDir $runStyle
echo "run_compareToAnnot-n.R complete."
echo ""

echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
