#!/bin/bash

#$ -N iri_pipeline
#$ -P jknight.prjc -q long.qc
#$ -o stdout_pipeline.log -e stderr_L.log -j y
#$ -cwd -V

echo "START : `date`"

module load R/3.1.3
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/

# Pre-installed tools

ZCAT="/bin/zcat"
FASTQC="/well/jknight/Irina/Programs/FastQC/fastqc"
CUTADAPT="/apps/well/python/2.7.8/bin/cutadapt"
BWA="/apps/well/bwa/0.7.10/bwa"
BWA_INDEX="/well/jknight/reference/mapping/bwa_g/Gencode19.GRCh37.genome.PGF.decoy.ERCC"
BOWTIE="/apps/well/bowtie2/2.2.5/bowtie2"
BOWTIE_INDEX="/well/jknight/reference/GRCh37/bowtie2/hs37d5c2_cox"
SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"
PICARD="/apps/well/picard-tools/1.111/picard-Xmx3g"
BEDTOOLS="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools"
BED_TO_WIG="/well/jknight/software/rescomp/bin/bedGraphToBigWig"
CHROM_SIZE="/well/jknight/reference/hs37_cox.chrom.sizes"
MACS2="/apps/well/python/2.7.8/bin/macs2"
BAMTOOLS="/apps/well/bamtools/2.3.0/bin/bamtools"
ATAC_QC="/well/jknight/Scripts/ATACseq_qc.R"
DIST_BETWEEN_PEAKS="/well/jknight/Scripts/distance_between_regions.pl"
EXTRACT_FRAGMENTS="/well/jknight/Scripts/bed_extract_fragments.pl"
BLACKLIST="/well/jknight/ATACseq/ATACseq_001/Analysis/Peak_calling/blacklist_filtering/wgEncodeDacMapabilityConsensusExcludable.bed"


# Input data

INPUT_READ1=$1
INPUT_READ2=$2


# Names of the files that will be generated while the script is running

FILENAME=$3
READ1=`echo $FILENAME.read1`
READ2=`echo $FILENAME.read2`


# ======================================================================
# =                                                                    =
# =                              ANALYSIS                              =
# =                                                                    =
# ======================================================================

mkdir $FILENAME.atac_analysis

# Extracting raw data from a zip archive

echo -n "`date`: unzipping the raw fastq files... "
$ZCAT $INPUT_READ1 > $FILENAME.atac_analysis/$READ1.fastq
$ZCAT $INPUT_READ2 > $FILENAME.atac_analysis/$READ2.fastq
echo "done."

cd $FILENAME.atac_analysis

# Quality control on raw data

mkdir $READ1.fastqc_output
mkdir $READ2.fastqc_output

echo -n "`date`: running FASTQC on raw fastq files... "
$FASTQC \
          $READ1.fastq \
          -f fastq \
          -o $READ1.fastqc_output

$FASTQC \
          $READ2.fastq \
          -f fastq \
          -o $READ2.fastqc_output
echo "done."

# Adapter trimming, removing low quality bases

echo -n "`date`: removing adapters... "
$CUTADAPT \
          -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
          -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
          --overlap 10 \
          --minimum-length=20 \
          -o $READ1.trimmed.fastq \
          -p $READ2.trimmed.fastq \
          $READ1.fastq $READ2.fastq
echo "done."

# Quality control on processed (cleaned) data

mkdir $READ1.trimmed.fastqc_output
mkdir $READ2.trimmed.fastqc_output

echo -n "`date`: running FASTQC on trimmed files... "
$FASTQC \
          $READ1.trimmed.fastq \
          -f fastq \
          -o $READ1.trimmed.fastqc_output

$FASTQC \
          $READ2.trimmed.fastq \
          -f fastq \
          -o $READ2.trimmed.fastqc_output
echo "done."

# Aligning processed data

echo -n "`date`: running the alignment... "
## $BOWTIE \
##           -p 2 \
##           -X 2000 -I 0 \
##           --no-discordant --no-mixed \
##           -x $BOWTIE_INDEX \
##           -1 $READ1.trimmed.fastq \
##           -2 $READ2.trimmed.fastq \
##           -S $FILENAME.sam \
##           > $FILENAME.bowtie.log \
##           2> $FILENAME.bowtie.err


$BWA mem \
          -t 3 \
          $BWA_INDEX \
          $READ1.trimmed.fastq \
          $READ2.trimmed.fastq > \
          $FILENAME.sam
echo "done."

# Manipulations on sam file: creating bam, sorting, indexing

echo -n "`date`: converting sam to bam... "
$SAMTOOLS view -bSh $FILENAME.sam > $FILENAME.bam
echo "done."

echo -n "`date`: sorting non-filtered bam $FILENAME.bam ... "
$SAMTOOLS sort \
          $FILENAME.bam \
          $FILENAME.sorted
echo "done."

echo -n "`date`: indexing non-filtered sorted bam $FILENAME.sorted ... "
$SAMTOOLS index \
          $FILENAME.sorted.bam
echo "done."

# Removing intermediate files
rm $READ1.trimmed.fastq
rm $READ2.trimmed.fastq
rm $READ1.fastq
rm $READ2.fastq
rm $FILENAME.sam
rm $FILENAME.bam

# Post-alignment quality control, statistical report

echo -n "`date`: generating insertSizeMetric report with picard..."
$PICARD CollectInsertSizeMetrics.jar \
          VALIDATION_STRINGENCY=LENIENT \
          ASSUME_SORTED=true \
          HISTOGRAM_FILE=$FILENAME.picard_histogram \
          INPUT=$FILENAME.sorted.bam \
          OUTPUT=$FILENAME.picard_insertSizeMetric_report
echo "done."

# Post-alignment data filtering/cleaning

# 1. Removing duplicates

echo -n "`date`: filtering step 1 -- removing duplicates..."
$PICARD MarkDuplicates.jar \
          REMOVE_DUPLICATES=true \
          INPUT=$FILENAME.sorted.bam \
          OUTPUT=$FILENAME.nodup.bam \
          METRICS_FILE=$FILENAME.picard_metrics \
          ASSUME_SORTED=true
echo "done."


# 2. Removing reads with low mapping quality (MAPQ < 30)

echo -n "`date`: filtering step 2 -- selecting reads with high mapping quality > 30..."
$BAMTOOLS filter \
          -mapQuality ">30" \
          -in $FILENAME.nodup.bam \
          -out $FILENAME.nodup.MAPQ30.bam
echo "done."


# 3. Removing non-uniquely mapped reads

echo -n "`date`: filtering step 3 -- removing non-uniquely mapped reads..."
$SAMTOOLS view \
          -h $FILENAME.nodup.MAPQ30.bam | \
          grep -v "XA" | \
          $SAMTOOLS view -bS - > \
          $FILENAME.nodup.MAPQ30.uniq.bam
echo "done."


# 4. Removing reads that were not properly paired (wrong orientation, based on the bitscore flag 0x62)

echo -n "`date`: filtering step 4 -- removing not properly paired reads..."

## AC
# $SAMTOOLS view \
#           -h -b \
#           -f 0x62 \
#           $FILENAME.nodup.MAPQ30.uniq.bam > \
#           $FILENAME.nodup.MAPQ30.uniq.proper_paired_read1.bam

# $SAMTOOLS view \
#           -h -b \
#           -f 0x92 \
#           $FILENAME.nodup.MAPQ30.uniq.bam > \
#           $FILENAME.nodup.MAPQ30.uniq.proper_paired_read2.bam

# $PICARD MergeSamFiles.jar \
#           INPUT=$FILENAME.nodup.MAPQ30.uniq.proper_paired_read1.bam \
#           INPUT=$FILENAME.nodup.MAPQ30.uniq.proper_paired_read2.bam \
#           SORT_ORDER=coordinate \
#           OUTPUT=$FILENAME.nodup.MAPQ30.uniq.proper_paired.bam

## AC
$SAMTOOLS view \
	  -h \
	  -b \
	  -f 3 \
	  -F 12 \
	  $FILENAME.nodup.MAPQ30.uniq.bam > \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.bam


echo "done."

## AC
# rm $FILENAME.nodup.MAPQ30.uniq.proper_paired_read1.bam
# rm $FILENAME.nodup.MAPQ30.uniq.proper_paired_read2.bam

# 5. Removing reads mapped to mitochondrial DNA

echo -n "`date`: filtering step 5 -- removing reads mapped to non-conventional chromosomes and mitochondrial DNA... "
$SAMTOOLS view \
          -h $FILENAME.nodup.MAPQ30.uniq.proper_paired.bam | \
          grep -v -i "chrM\|chrMT\|GL\|NC\|hs" | \
          $SAMTOOLS view -bS - > \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.bam

$SAMTOOLS index \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.bam
echo "done."


# Post-alignment quality control, statistical report on filtered file

echo -n "`date`: generating insertSizeMetric report on filtered data..."
$PICARD CollectInsertSizeMetrics.jar \
          VALIDATION_STRINGENCY=LENIENT \
          ASSUME_SORTED=true \
          HISTOGRAM_FILE=$FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.picard_histogram \
          INPUT=$FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.bam \
          OUTPUT=$FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.picard_insertSizeMetric_report
echo "done."


# Creating a bigWig file from a filtered bam file

echo -n "`date`: generating the bigwig file for the filtered bam file... "
$BEDTOOLS genomecov \
          -bg \
          -ibam $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.bam > \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.coverage

## AC
LC_COLLATE=C
sort -k1,1 -k2,2n \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.coverage > \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.sorted.coverage

$BED_TO_WIG \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.sorted.coverage \
          $CHROM_SIZE \
          $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM.bw
echo "done."

# Changing names of semi-final files

echo -n "`date`: renaming filtered files... "
for i in `ls $FILENAME.nodup.MAPQ30.uniq.proper_paired.noM*`
do
    j=`echo $i | sed s/nodup.MAPQ30.uniq.proper_paired.noM/filtered/g`
    mv $i $j
done
echo "done."


# Calling peaks

echo -n "`date`: calling peaks on $FILENAME.filtered.bam ... "
$MACS2 callpeak \
          --nomodel \
          -t $FILENAME.filtered.bam \
          --name $FILENAME.macs2 \
          --outdir $FILENAME.macs2_results \
          --nolambda \
          --shift 5 \
          --keep-dup all \
          --slocal 10000 \
          --SPMR \
          --bdg
echo "done."


# ======================================================================
# =                                                                    =
# =                         QUALITY CONTROL                            =
# =                                                                    =
# ======================================================================

# Making sure the peak file contains "chrN" as chromosome name
cat $FILENAME.macs2_results/$FILENAME.macs2_peaks.narrowPeak | \
    sed s/chr//g | \
    awk '{print "chr" $0}' > \
    $FILENAME.macs2_results/$FILENAME.macs2_peaks.chr.narrowPeak

# Filtering peaks -- removing regions from "black list"
$BEDTOOLS intersect \
    -v \
    -a $FILENAME.macs2_results/$FILENAME.macs2_peaks.chr.narrowPeak \
    -b $BLACKLIST \
    > $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak

# Distribution of the peak length  -- for the periodicity plot
 
$DIST_BETWEEN_PEAKS \
    $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak > \
    $FILENAME.peaks_length.hist


# Distribution of the peak width -- for the histogram of peak widths

cat $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
    awk '{print $3-$2}' > \
    $FILENAME.peaks_width.hist

# QUALITY CONTROL 1 -- Called_peaks / Fragments_per_chromosome

# Number of called peaks per chromosome
echo -e "Chromosome\tCalled_peaks" > $FILENAME.temp_peaks_per_chrom
for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> $FILENAME.temp_peaks_per_chrom
    cat $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
    grep -cP "^chr$i\t" >> \
    $FILENAME.temp_peaks_per_chrom
done

# Number of reads mapped per chromosome -- creating a BEDPE file

# creating bam file sorted on names
$SAMTOOLS sort \
    -n \
    -T temp.bam \
    $FILENAME.filtered.bam \
    -o $FILENAME.filtered.name_sort.bam

# converting bam to bedpe

## AC: output file gets overwritten
# $BEDTOOLS bamtobed \
#     -bedpe \
#     -i $FILENAME.filtered.name_sort.bam > \
#     $FILENAME.filtered.bedpe.bed

# $BEDTOOLS bamtobed \
#     -bedpe \
#     -i $FILENAME.filtered.name_sort.bam | \
#     grep -v 'GL' | \
#     grep -v 'NC' | \
#     grep -v 'hs' | \
#     grep -P '\t\+\t\-' > \
#     $FILENAME.filtered.bedpe.bed

# $BEDTOOLS sort \
#     -i $FILENAME.filtered.bedpe.bed > \
#     $FILENAME.filtered.bedpe.name_sorted.bed

# # Selecting beginning and end of each fragment
# cat $FILENAME.filtered.bedpe.name_sorted.bed | \
#     sed s/chr//g | awk '{print "chr" $1 "\t" $2 "\t" $6 "\t" $7 }' > \
#     $FILENAME.filtered.bedpe.fragments.bed

getFragmentBed_pe_version.pl \
    $FILENAME.filtered.name_sort.bam \
    $FILENAME.filtered.bedpe.bed

$BEDTOOLS \
    sort \
    -i $FILENAME.filtered.bedpe.bed \
    $FILENAME.filtered.bedpe.fragments.bed 


# Number of reads mapped per chromosome
echo -e "Chromosome\tFragments_per_chrom" > $FILENAME.temp_fragm_per_chrom

for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> \
        $FILENAME.temp_fragm_per_chrom
    cat $FILENAME.filtered.bedpe.fragments.bed | \
        grep -cP "^chr$i\t" >> \
        $FILENAME.temp_fragm_per_chrom
done

# QUALITY CONTROL 2 -- Fragments_in_peaks / Fragments_per_chromosome

# Overlapping fragments with the peak file
$BEDTOOLS intersect \
    -f 0.10 \
    -wa \
    -a $FILENAME.filtered.bedpe.fragments.bed \
    -b $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
    uniq > \
    $FILENAME.filtered.bedpe.fragments_to_peaks.bed

echo -e "Chromosome\tFragments_in_peaks" > $FILENAME.temp_fragm_per_peaks

for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> \
        $FILENAME.temp_fragm_per_peaks
    cat $FILENAME.filtered.bedpe.fragments_to_peaks.bed | \
        grep -cP "^chr$i\t" >> \
        $FILENAME.temp_fragm_per_peaks
done

# QUALITY CONTROL 4 -- Fragments_in_peaks / Fragments_off_peaks_same_size

# Generating off peak genome-wide file
$BEDTOOLS subtract \
    -a $CHROM_SIZE \
    -b $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak > \
    $FILENAME.no_peak.bed

# Generating regions outside peaks of peak length
echo > $FILENAME.off_peaks.peak_size.bed
for i in `echo {1..22} X Y`
do
    cat $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
        grep -P "^chr$i\t" > \
        peaks.$i.bed
    cat $FILENAME.no_peak.bed | \
        grep -P "^chr$i\t" > \
        no_peaks.$i.bed
    $EXTRACT_FRAGMENTS \
        peaks.$i.bed \
        no_peaks.$i.bed >> \
        $FILENAME.off_peaks.peak_size.bed
    rm peaks.$i.bed no_peaks.$i.bed
done

# Getting reads mapped to regions outside peaks of peak length
$BEDTOOLS intersect \
    -f 0.10 \
    -wa \
    -a $FILENAME.filtered.bedpe.fragments.bed \
    -b $FILENAME.off_peaks.peak_size.bed | \
    uniq > \
    $FILENAME.fragments_to_peaks.peak_length.bed

echo -e "Chromosome\tFragments_per_off_peaks_peak_length" > \
    $FILENAME.temp_fragm_per_off_peaks_peak_length
for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> \
        $FILENAME.temp_fragm_per_off_peaks_peak_length
    cat $FILENAME.fragments_to_peaks.peak_length.bed | \
        grep -cP "^chr$i\t" >> \
        $FILENAME.temp_fragm_per_off_peaks_peak_length
done

paste \
          <(awk '{print $0}' $FILENAME.temp_peaks_per_chrom) \
          <(awk '{print $2}' $FILENAME.temp_fragm_per_peaks) \
          <(awk '{print $2}' $FILENAME.temp_fragm_per_chrom) \
          <(awk '{print $2}' $FILENAME.temp_fragm_per_off_peaks_peak_length) > $FILENAME.stat.txt

rm $FILENAME.temp_*
cat $FILENAME.filtered.picard_insertSizeMetric_report | \
    grep -v -P '.*\t.*\t.*\t' | \
    grep -P '^\d|^insert' > \
    $FILENAME.filtered.insertSize_hist

echo -n "`date`: creating plots with quality metrics... "
Rscript $ATAC_QC \
    $FILENAME.peaks_length.hist \
    $FILENAME.peaks_width.hist \
    $FILENAME.stat.txt \
    $FILENAME.filtered.insertSize_hist \
    $FILENAME.QC.pdf
echo "done."

echo "END : `date`"

