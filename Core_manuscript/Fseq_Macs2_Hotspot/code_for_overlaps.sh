## This code should be sufficient to generate intersections of peaks and their overlaps with the annotation

for l in `echo 100 200 400 600 800 1000 2000`
do
    for t in `echo 2 4 6 8 10 12 14 16`
    do
        fresh1="dnase_rep1.filtered.bedpe.fragments.fseq.length_$l.tresh_$t.narrowPeak"
        fresh2="dnase_rep2.filtered.bedpe.fragments.fseq.length_$l.tresh_$t.narrowPeak"
        fresh3="dnase_rep3.filtered.bedpe.fragments.fseq.length_$l.tresh_$t.narrowPeak"
        bedtools intersect -f 0.3 -r -a $fresh1 -b $fresh2 |  bedtools intersect -f 0.3 -r -a - -b $fresh3 | bedtools sort -i - | bedtools merge -i - > dnase.fseq.length_$l.tresh_$t.intrsct.bed;
        echo Done for $l $t
    done
done



for l in `echo 10 100 200 400 600 800 1000 2000`
do
    for t in `echo 0.001 0.01 0.1 0.99`
    do
        fresh1="dnase_rep1.macs.length_$l.thresh_$t.narrowPeak"
        fresh2="dnase_rep2.macs.length_$l.thresh_$t.narrowPeak"
        fresh3="dnase_rep3.macs.length_$l.thresh_$t.narrowPeak"
        bedtools intersect -f 0.3 -r -a $fresh1 -b $fresh2 | bedtools intersect -f 0.3 -r -a - -b $fresh3 | bedtools sort -i - | bedtools merge -i - > dnase.macs.length_$l.thresh_$t.intrsct.bed;
        echo Done for $l $t
    done
    echo
done



annotation="seg_on.bed"
for i in `ls dnase.fseq.*.intrsct.bed`
do
    name=`echo $i | sed s/100.t/0100.t/g | sed s/200.t/0200.t/g | sed s/400.t/0400.t/g | sed s/600.t/0600.t/g | sed s/800.t/0800.t/g | sed s/h_2.intrsct.bed/h_02.intrsct.bed/g | sed s/h_4.intrsct.bed/h_04.intrsct.bed/g | sed s/h_6.intrsct.bed/h_06.intrsct.bed/g | sed s/h_8.intrsct.bed/h_08.intrsct.bed/g`
    all_bases=`bedtools merge -i $i | awk '{sum += $3-$2} END {print sum}'`
    all_TSS=`bedtools sort -i $annotation | bedtools merge -i - | awk '{sum += $3-$2} END {print sum}'`
    called_TSS=`bedtools intersect -f 0.99 -a $i -b $annotation | bedtools sort -i - | bedtools merge -i - | awk '{sum += $3-$2} END {print sum}'`
    called_nonTSS=`echo $all_bases $called_TSS | awk '{print $1-$2}'`
    echo -e "$name\tall_bases=$all_bases\tall_TSS=$all_TSS\tcalled_TSS=$called_TSS\tcalled_nonTSS=$called_nonTSS" | sed 's/all_bases=//g' | sed 's/all_TSS=//g' | sed 's/called_TSS=//g' | sed 's/called_nonTSS=//g' | sed s/dnase.fseq.length/l/g | sed s/tresh/t/g | sed s/.intrsct.bed//g
done


for i in `ls dnase.macs.*.intrsct.bed`
do
    name=`echo $i | sed s/10.t/0010.t/g | sed s/100.t/0100.t/g | sed s/200.t/0200.t/g | sed s/400.t/0400.t/g | sed s/600.t/0600.t/g | sed s/800.t/0800.t/g | sed s/h_2.intrsct.bed/h_02.intrsct.bed/g | sed s/h_4.intrsct.bed/h_04.intrsct.bed/g | sed s/h_6.intrsct.bed/h_06.intrsct.bed/g | sed s/h_8.intrsct.bed/h_08.intrsct.bed/g`
    all_bases=`bedtools merge -i $i | awk '{sum += $3-$2} END {print sum}'`
    all_TSS=`bedtools sort -i $annotation | bedtools merge -i - | awk '{sum += $3-$2} END {print sum}'`
    called_TSS=`bedtools intersect -f 0.99 -a $i -b $annotation | bedtools sort -i - | bedtools merge -i - | awk '{sum += $3-$2} END {print sum}'`
    called_nonTSS=`echo $all_bases $called_TSS | awk '{print $1-$2}'`
    echo -e "$name\tall_bases=$all_bases\tall_TSS=$all_TSS\tcalled_TSS=$called_TSS\tcalled_nonTSS=$called_nonTSS" | sed 's/all_bases=//g' | sed 's/all_TSS=//g' | sed 's/called_TSS=//g' | sed 's/called_nonTSS=//g' | sed s/dnase.fseq.length/l/g | sed s/tresh/t/g | sed s/.intrsct.bed//g
done | sort -n -k 1


