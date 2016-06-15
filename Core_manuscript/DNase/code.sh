## Generating intersections

fresh3="./dnase.rep3.pinechrom_general/dnase.rep3.fseq.chr.narrowPeak"
fresh2="./dnase.rep2.pinechrom_general/dnase.rep2.fseq.chr.narrowPeak"
fresh1="./dnase.rep1.pinechrom_general/dnase.rep1.fseq.chr.narrowPeak"

bedtools intersect \
    -f 0.3 -r \
    -a $fresh1 \
    -b $fresh2 | \
    bedtools intersect \
    -f 0.3 -r \
    -a - \
    -b $fresh3 | \
    bedtools sort \
    -i - | \
    bedtools merge \
    -i - > \
    dnase.fseq.intrsct.bed



## How annotation was created:

we="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-WE.bed"
e="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-E.bed"
tss="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-TSS.bed"
ctcf="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-CTCF.bed"

cat $we $e $tss $ctcf | bedtools sort -i - | bedtools merge -i - > segmentation.on_target.bed

pf="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-PF.bed"
r="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-R.bed"
t="/home/irina/Work/ATAC-Seq/Git/ATACseq_pipeline/Core_manuscript/Peaks/Annotation/wgEncodeAwgSegmentationCombinedK562-T.bed"

cat $pf $r $t | bedtools sort -i - | bedtools merge -i - > segmentation.off_target.bed



## Calculating overlap with the annotation:

annotation="segmentation.off_target.bed"

for sample in `ls *intrsct.bed`
do
    bases_in_sample=`cat $sample | awk '{sum += $3-$2} END {print sum}'`
    bases_in_overlap=`bedtools intersect -f 0.3 -r -a $sample -b $annotation | bedtools sort -i - | bedtools merge -i - | awk '{sum += $3-$2} END {print sum}'`
    percent_overlap=`echo $bases_in_sample $bases_in_overlap | awk '{print ($2*100)/$1}'`
    echo -e "$sample\t$percent_overlap"
done


## Distribution of peaks across categories of the annotation

for i in `ls Annotation/wgEncodeAwgSegmentationCombinedK562-*bed`
do
    for j in `ls *unique.bed`
    do
        annotation=`echo $i | sed s/.*K562-//g`
        peaks=`bedtools intersect -f 0.3 -r -a $j -b $i | bedtools intersect -u -a - -b $j | bedtools sort -i - | bedtools merge -i - | wc -l`
        bases=`bedtools intersect -f 0.3 -r -a $j -b $i | bedtools intersect -u -a - -b $j | bedtools sort -i - | bedtools merge -i - | awk '{sum += $3-$2} END {print sum}'`
        echo -e "$annotation\t$j\t$peaks\t$bases"
    done
    echo
done

