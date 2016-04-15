#/bin/bash

#$ -N fseq
#$ -P jknight.prjc -q short.qc
#$ -o stdout_Fseq -e sterr_Fseq.log -j y
#$ -cwd -V

### Run F-seq on a BED File With User-Specified Parameters
## JHendry, 2016/12/01

## Idea here is to take directory containing .bed file and generate
## an .npf file representing peak calls from F-seq algorithm
## run under certain parameters.
## Default for F-seq is to produce a .npf file for each chromosome:
## e.g. 1.npf 2.npf ... X.npf
## These files are concatenated into one .npf file, and
## the resultant file named to match the original .bed file.
## Finally, this .npf file is moved into a new directory
## which is named as follows:
##
## fseq<l-setting><t-setting>_<sample-name>
##
## This folder is created IF it does not already exist.
## The ultimate result is the original folder containing the .bed
## file is returned to its original state (all intermediary 
## NPF files are removed, F-seq can be run again) and an .npf 
## file representing peaks is moved into an informatively 
## named directory for subsequent analysis.

## To achieve this, inputs are:
## $1 ---> .bed file containing directory
## $2 ---> value for F-seq -l parameter
## $3 ---> value for F-seq -t parameter

## Note: Chromosomal NPF files (1.npf, 2.npf etc.) produced by
## the F-seq algorithm are deposited in a temporary subdirectory 
## created within the .bed file containing directory. 
## This allows for run_Fseq.sh to be run multiple times simultaneously 
## on the same .bed file containing directory without the chromosomal NPF 
## files (1.npf.. etc.) being confused.


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"


### Define Input Directory, Parameter Values,  Get Input File Name
bedDir=$1
bedFile=$(ls $bedDir/*.bed)
bedFile=${bedFile##*/}

lengthVal=$2
threshVal=$3


echo "Input File Directory:"$bedDir
echo "File Name:"$bedFile
echo "Feature Length Value (-l): "$lengthVal
echo "Threshold Value (-t): "$threshVal
echo ""

### Temporary Subdirectory to Output 1.npf, 2.npf, ... X.npf files
npfDirPrefix=$(echo "fseql"$lengthVal"t"$threshVal)
tempNpfDir=$(echo $bedDir"/"$npfDirPrefix)

if [ ! -d "$tempNpfDir" ]; then
  mkdir $tempNpfDir
  echo "Making Temporary Output Directory:"$tempNpfDir
  echo ""
fi

### Run F-seq
fseq -d $bedDir -o $tempNpfDir \
-f 0 \
-l $lengthVal \
-t $threshVal \
-of npf \
-v

echo ""
echo "F-seq Complete"
echo "Catenating Output..."
echo ""

### Make Destination Directory for NPF file
npfDirSuffixClrPath=${bedDir##*/}
npfDirSuffix=${npfDirSuffixClrPath%_*}
npfDir=$(echo $npfDirPrefix"_"$npfDirSuffix)

if [ ! -d "$npfDir" ]; then
  mkdir $npfDir
  echo "Making Destination Directory:"$npfDir
  echo ""
fi

npfFileName=$(echo $npfDir"/"$bedFile | sed 's/.bed$/.npf/')
cat $tempNpfDir/*.npf > $npfFileName
rm -r $tempNpfDir

echo "Output File:" $npfFileName
echo "Output Directory:" $npfDir
echo ""

echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
