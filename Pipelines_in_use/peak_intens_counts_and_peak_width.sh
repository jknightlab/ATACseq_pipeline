#!/bin/bash

#$ -N peak_intens
#$ -P jknight.prjc -q long.qc
#$ -e stderr_L.log -j y
#$ -cwd -V

echo "START : `date`"

PeakFile=$1 # bed file containing peaks
BamFile=$2 # bam file needs to be indexed

for i  in `cat $PeakFile | awk '{print $1 ":" $2 "-" $3}'`;
do
    reads=`samtools view $BamFile $i -F 4 -c` #count total number of reads which are not unmapped in that location
    bp_in_peak=`echo $i | sed 's/\:/\t/g' | sed 's/\-/\t/g' | awk '{print $3-$2}'`  #$3-$2 is calculating number of bp in the peak
    echo -e "$i\t$reads\t$bp_in_peak" | awk '{print $1 "\t" $2 "\t" $3}' >> $PeakFile.counts.txt
done


echo "END : `date`"
