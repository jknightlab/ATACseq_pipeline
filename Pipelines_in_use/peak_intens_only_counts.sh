#!/bin/bash

#$ -N peak_intens
#$ -P jknight.prjc -q long.qc
#$ -e stderr_L.log -j y
#$ -cwd -V

echo "START : `date`"

PeakFile=$1 # bed file containing peaks
BamFile=$2 # bam file needs to be indexed
FileName=$3 

for i in `cat $PeakFile | awk '{print $1 ":" $2 "-" $3}'`; 
do
    reads=`samtools view $BamFile $i -F 4 -c` #count total number of reads which are not unmapped in that location
    echo -e "$i\t$reads" | awk '{print $1 "\t" $2}' >> $FileName.counts.txt
done


echo "END : `date`"
