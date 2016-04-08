#/bin/bash

#$ -N fseq
#$ -P jknight.prjc -q short.qc
#$ -e sterr_Fseq.log -j y
#$ -cwd -V


bedFile=$1
lengthVal=$2
threshVal=$3

bedDir=`echo $bedFile | sed 's/\.bed$/\_output/g'`
mkdir $bedDir
mv $bedFile $bedDir
tempDir="$bedDir/TEMP.fseq.length_$lengthVal.tresh_$threshVal"
mkdir $tempDir

fseq -d $bedDir -o $tempDir \
-f 0 \
-l $lengthVal \
-t $threshVal \
-of npf \
-v

npfFile=`echo $bedFile | sed -e "s/\.bed$/\.fseq.length_\$lengthVal.tresh_\$threshVal.narrowPeak/g"`
rm -rf $bedDir/$npfFile

for i in `echo {1..22} X Y`
do
    cat $tempDir/chr$i.npf >> $bedDir/$npfFile
    echo done for chr$i
done

rm -rf $tempDir

