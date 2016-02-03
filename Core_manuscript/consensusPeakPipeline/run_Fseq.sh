#/bin/bash

#$ -N fseq
#$ -P jknight.prjc -q short.qc
#$ -o stdout_Fseq -e sterr_Fseq.log -j y
#$ -cwd -V

### Run F-seq with Default Parameters on Explicitly Named Directory
## JHendry, 2016/12/01

## Idea here is to take directory containing .bed file and generate
## an .npf file representing peak calls from F-seq algorithm
## run under certain parameters.
## Default for F-seq is to produce a .npf file for each chromosome:
## e.g. 1.npf 2.npf ... X.npf
## These files must all be concatenated into one .npf file, and
## the resultant file named to match the original .bed file.
## Finally, this .npf file must be moved into a new directory
## which will be named as follows:
##
## fseq<l-setting><t-setting>_<sample-name>
##
## This folder is created IF it does not already exist.
## The ultimate result is the original folder containing the .bed
## file is returned to its original state (such that F-seq can be
## run again) and an .npf file representing peaks is moved
## into an informatively named directory for subsequent analysis.

## To achieve this, inputs are:
## $1 ---> .bed file containing directory
## $2 ---> value for F-seq -l parameter
## $3 ---> value for F-seq -t parameter

## NB: YOU CAN ONLY RUN THIS SCRIPT ON A DIRECTORY **ONCE**
## AT A TIME--OTHERWISE THE DEFAULT OUTPUT FILES FROM
## FSEQ (1.npf, 2.npf, etc.) WILL BE CONFUSED AND CATENATED
## INTO ONE OUTPUT FILE.

echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"


### Define Input Directory, Get Input File Name
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

### Run F-seq
fseq -d $bedDir -o $bedDir \
-f 0 \
-l $lengthVal \
-t $threshVal \
-of npf \
-v

echo ""
echo "F-seq Complete"
echo "Catenating Output..."
echo ""

### If necessary, make .npf directory
npfDirPrefix=$(echo "fseql"$lengthVal"t"$threshVal)
npfDirSuffixClrPath=${bedDir##*/}
npfDirSuffix=${npfDirSuffixClrPath%_*}
npfDir=$(echo $npfDirPrefix"_"$npfDirSuffix)

if [ ! -d "$npfDir" ]; then
  mkdir $npfDir
fi

npfFileName=$(echo $npfDir"/"$bedFile | sed 's/.bed$/.npf/')
cat $bedDir/*.npf > $npfFileName
rm $bedDir/*.npf

echo "Output File:" $npfFileName
echo "Output Directory:" $npfDir
echo ""

echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
