#/bin/bash

#$ -N macs
#$ -P jknight.prjc -q short.qc
#$ -e sterr_Macs.log -j y
#$ -cwd -V

MACS2="/apps/well/python/2.7.8/bin/macs2"

Input_File=$1
Prefix=$2
length=$3
thresh=$4

$MACS2 callpeak \
    --nomodel \
    -t $Input_File \
    --name $Prefix.macs.length_$length.thresh_$thresh.narrowPeak \
    --outdir $Prefix.macs2.length_$length.thresh_$thresh \
    --extsize $length \
    --qvalue $thresh \
    --nolambda \
    --shift 5 \
    --keep-dup all \
    --slocal 10000 \
    --SPMR \
    --bdg

cat $Prefix.macs2.length_$length.thresh_$thresh/$Prefix.macs.length_$length.thresh_$thresh.narrowPeak_peaks.narrowPeak | \
    awk '{print "chr" $0}' > $Prefix.macs.length_$length.thresh_$thresh.narrowPeak

rm -rf $Prefix.macs2.length_$length.thresh_$thresh

