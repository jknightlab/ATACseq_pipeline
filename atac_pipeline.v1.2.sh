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
SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"
PICARD="/apps/well/picard-tools/1.111/picard-Xmx3g"
BEDTOOLS="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools"
BEDGRAPH_TO_BIGWIG="/well/jknight/software/rescomp/bin/bedGraphToBigWig"
GENOME_COVERAGE_BED="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/genomeCoverageBed"
CHROM_SIZE="/well/jknight/reference/hs37_cox.chrom.sizes"
MACS2="/apps/well/python/2.7.8/bin/macs2"
BAMTOOLS="/apps/well/bamtools/2.3.0/bin/bamtools"
ATAC_QC="/well/jknight/Scripts/ATACseq_qc.R"
DIST_BETWEEN_PEAKS="/well/jknight/Scripts/distance_between_regions.pl"
EXTRACT_FRAGMENTS="/well/jknight/Scripts/bed_extract_fragments.pl"
BLACKLIST="/well/jknight/ATACseq/ATACseq_001/Analysis/Peak_calling/blacklist_filtering/wgEncodeDacMapabilityConsensusExcludable.bed"
ATAC_QC="/well/jknight/Scripts/ATACseq_qc.R"


# Input data, fastq files

INPUT_READ1=$1
INPUT_READ2=$2


# Names of the files that will be generated while the script is running

FILENAME=$3
READ1=`echo $FILENAME.read1`
READ2=`echo $FILENAME.read2`


# ======================================================================
# =                                                                    =
# =                        PRIMARY ANALYSIS                            =
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
    -o $READ1.fastqc_output \
    2> $FILENAME.fastqc_stderr.txt

$FASTQC \
    $READ2.fastq \
    -f fastq \
    -o $READ2.fastqc_output \
    2> $FILENAME.fastqc_stderr.txt
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
    $READ1.fastq $READ2.fastq \
    2> $FILENAME.cutadapt_stderr.txt
echo "done."

# Quality control on processed (cleaned) data

mkdir $READ1.trimmed.fastqc_output
mkdir $READ2.trimmed.fastqc_output

echo -n "`date`: running FASTQC on trimmed files... "
$FASTQC \
    $READ1.trimmed.fastq \
    -f fastq \
    -o $READ1.trimmed.fastqc_output \
    2> $FILENAME.fastqc_stderr.txt

$FASTQC \
    $READ2.trimmed.fastq \
    -f fastq \
    -o $READ2.trimmed.fastqc_output \
    2> $FILENAME.fastqc_stderr.txt
echo "done."

# Aligning processed data
echo -n "`date`: running the alignment... "
$BWA mem \
    -t 3 \
    $BWA_INDEX \
    $READ1.trimmed.fastq \
    $READ2.trimmed.fastq > \
    $FILENAME.sam \
    2> $FILENAME.bwa_stderr.txt
echo "done."

# Manipulations on sam file: creating bam, sorting, indexing

echo -n "`date`: converting sam to bam... "
$SAMTOOLS view \
    -bSh $FILENAME.sam > \
    $FILENAME.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: sorting alignments in orer to remove duplicates... "
$SAMTOOLS sort \
    $FILENAME.bam \
    $FILENAME.sorted \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: filtering step 1 -- removing duplicates using SAMTOOLS 0.1.19... "
/apps/well/samtools/0.1.19/bin/samtools rmdup \
    $FILENAME.sorted.bam \
    $FILENAME.nodup.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: filtering step 2 -- selecting reads with high mapping quality > 30..."
$BAMTOOLS filter \
    -mapQuality ">30" \
    -in $FILENAME.nodup.bam \
    -out $FILENAME.nodup.MAPQ30.bam \
    2> $FILENAME.bamtools_stderr.txt
echo "done."

echo -n "`date`: filtering step 3 -- removing non-uniquely mapped reads..."
$SAMTOOLS view \
    -h $FILENAME.nodup.MAPQ30.bam | \
    grep -v "XA" | \
    $SAMTOOLS view -bS - > \
    $FILENAME.nodup.MAPQ30.uniq.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: filtering step 4 -- removing not properly paired reads... "
$SAMTOOLS view \
    -h \
    -b \
    -f 3 \
    -F 12 \
    $FILENAME.nodup.MAPQ30.uniq.bam > \
    $FILENAME.nodup.MAPQ30.uniq.proper_paired.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: Generating sam with reads mapped to conventional chromosomes... "
$SAMTOOLS view \
    -H $FILENAME.nodup.MAPQ30.uniq.proper_paired.bam | \
    grep chr | \
    grep -v chrM > \
    $FILENAME.filtered.non_sorted.sam

$SAMTOOLS view \
    $FILENAME.nodup.MAPQ30.uniq.proper_paired.bam | \
    awk '$3 ~ /chr/ {print $0}' | \
    grep -v chrM >> \
    $FILENAME.filtered.non_sorted.sam

$SAMTOOLS view \
    -bS $FILENAME.filtered.non_sorted.sam > \
    $FILENAME.filtered.non_sorted.bam

$SAMTOOLS sort \
    $FILENAME.filtered.non_sorted.bam \
    $FILENAME.filtered \
    2> $FILENAME.samtools_stderr.txt
rm $FILENAME.filtered.non_sorted.sam
rm $FILENAME.filtered.non_sorted.bam
echo "done."

# Creating WIGGLE file
head -24 $CHROM_SIZE | sed s/chr//g | awk '{print "chr" $1 "\t" $2}' > $FILENAME.chrom_sizes.txt

echo -n "`date`: Generating bedgraph to create a bigwig file... "
$GENOME_COVERAGE_BED \
    -bg \
    -ibam $FILENAME.filtered.bam \
    -split \
    -g $FILENAME.chrom_sizes.txt > \
    $FILENAME.bedgraph \
    2> $FILENAME.bedtools_stderr.txt
echo "done."

echo -n "`date`: Sorting bedgraph to create a bigwig file... "
LC_COLLATE=C
sort -k1,1 -k2,2n $FILENAME.bedgraph > $FILENAME.sorted.bedgraph
echo "done."

echo -n "`date`: Generating bigwig from bedgraph... "
$BEDGRAPH_TO_BIGWIG \
    $FILENAME.sorted.bedgraph \
    $FILENAME.chrom_sizes.txt $FILENAME.bw \
    2> $FILENAME.bedtools_stderr.txt
echo "done".

echo -n "`date`: sorting filtered bam file by name... "
$SAMTOOLS sort \
    -n \
    -T $FILENAME.temp.bam \
    $FILENAME.filtered.bam \
    -o $FILENAME.filtered.name_sort.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: Generating bedpe from bam... "
$BEDTOOLS bamtobed \
    -bedpe \
    -i $FILENAME.filtered.name_sort.bam > \
    $FILENAME.filtered.bedpe.bed \
    2> $FILENAME.bedpe_stderr.txt
echo "done."


echo -n "`date`: Fixing the flags... "
$SAMTOOLS fixmate \
    $FILENAME.filtered.name_sort.bam \
    $FILENAME.fixed.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: Sorting fixed file by name... "
$SAMTOOLS sort \
    -n \
    -T $FILENAME.temp.bam \
    $FILENAME.fixed.bam \
    -o $FILENAME.fixed.name_sorted.bam \
    2> $FILENAME.samtools_stderr.txt
echo "done."

echo -n "`date`: BEDPE Creating bedpe file... "
$BEDTOOLS bamtobed \
    -bedpe \
    -i $FILENAME.fixed.name_sorted.bam  > \
    $FILENAME.filtered.bedpe.bed \
    2> $FILENAME.bedpe_stderr.txt
echo "done."

echo -n "`date`: Sorting bedpe file... "
$BEDTOOLS \
    sort \
    -i $FILENAME.filtered.bedpe.bed > \
    $FILENAME.filtered.bedpe.fragments.bed \
    2> $FILENAME.bedpe_stderr.txt
echo "done."

echo -n "`date`: calling peaks on $FILENAME.filtered.bam ... "
$MACS2 callpeak \
    --nomodel \
    -t $FILENAME.fixed.bam \
    --name $FILENAME.macs2 \
    --outdir $FILENAME.macs2_results \
    --nolambda \
    --shift 5 \
    --keep-dup all \
    --slocal 10000 \
    --SPMR \
    --bdg 2> \
    $FILENAME.macs2_stderr.txt
echo "done."

echo -n "`date`: Making sure the peak file contains "chrN" as chromosome name... "
cat $FILENAME.macs2_results/$FILENAME.macs2_peaks.narrowPeak | \
    sed s/chr//g | \
    awk '{print "chr" $0}' > \
    $FILENAME.macs2_results/$FILENAME.macs2_peaks.chr.narrowPeak
echo "done."

echo -n "`date`: Filtering peaks -- removing regions from black list... "
$BEDTOOLS intersect \
    -v \
    -a $FILENAME.macs2_results/$FILENAME.macs2_peaks.chr.narrowPeak \
    -b $BLACKLIST \
    > $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak \
    2> $FILENAME.bedtools_stderr.txt
echo "done."

echo -n "`date`: creating a list of peak widths -- for the histogram of peak width... "
cat $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
    awk '{print $3-$2}' > \
    $FILENAME.peaks_width.hist
echo "done."

echo "QUALITY CONTROL 1 -- Called_peaks / Fragments_per_chromosome"

echo -n "`date`: QC1 -- counting number of called peaks per chromosome... "
echo -e "Chromosome\tCalled_peaks" > $FILENAME.temp_peaks_per_chrom
for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> $FILENAME.temp_peaks_per_chrom
    cat $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
    grep -cP "^chr$i\t" >> \
    $FILENAME.temp_peaks_per_chrom
done
echo "done."

echo -n "`date`: QC1 -- counting number of reads mapped per chromosome... "
echo -e "Chromosome\tFragments_per_chrom" > $FILENAME.temp_fragm_per_chrom
for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> \
        $FILENAME.temp_fragm_per_chrom
    cat $FILENAME.filtered.bedpe.fragments.bed | \
        grep -cP "^chr$i\t" >> \
        $FILENAME.temp_fragm_per_chrom
done
echo "done."

echo "QUALITY CONTROL 2 -- Fragments_in_peaks / Fragments_per_chromosome"
echo -n "`date`: QC2 -- overlapping fragments with the peak file... "
$BEDTOOLS intersect \
    -f 0.10 \
    -wa \
    -a $FILENAME.filtered.bedpe.fragments.bed \
    -b $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak | \
    uniq > \
    $FILENAME.filtered.bedpe.fragments_to_peaks.bed \
    2> $FILENAME.bedtools_stderr.txt
echo "done."

echo -n "`date`: QC2 -- counting number of fragments in peaks per chromosome... "
echo -e "Chromosome\tFragments_in_peaks" > $FILENAME.temp_fragm_per_peaks
for i in `echo {1..22} X Y`
do
    echo -n -e chr$i"\t" >> \
        $FILENAME.temp_fragm_per_peaks
    cat $FILENAME.filtered.bedpe.fragments_to_peaks.bed | \
        grep -cP "^chr$i\t" >> \
        $FILENAME.temp_fragm_per_peaks
done
echo "done."

echo "QUALITY CONTROL 4 -- Fragments_in_peaks / Fragments_off_peaks_same_size"

echo -n "`date`: Generating off peak genome-wide file... "
cat $CHROM_SIZE | sed s/chr//g | awk '{print "chr" $1 "\t1\t" $2}' > $FILENAME.temp_chrom
$BEDTOOLS subtract \
    -a $FILENAME.temp_chrom \
    -b $FILENAME.macs2_results/$FILENAME.filtered.narrowPeak > \
    $FILENAME.no_peak.bed \
    2> $FILENAME.bedtools_stderr.txt
echo "done."

echo -n "`date`: Generating regions outside peaks of peak length... "
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
echo "done."

echo -n "`date`: Getting reads mapped to regions outside peaks of peak length... "
$BEDTOOLS intersect \
    -f 0.10 \
    -wa \
    -a $FILENAME.filtered.bedpe.fragments.bed \
    -b $FILENAME.off_peaks.peak_size.bed | \
    uniq > \
    $FILENAME.fragments_to_peaks.peak_length.bed \
    2> $FILENAME.bedtools_stderr.txt
echo "done."

echo -n "`date`: Counting the fragments mapped off peaks... "
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
echo "done."

echo -n "`date`: gathering all intermediate files into $FILENAME.stat.txt... "
paste \
    <(awk '{print $0}' $FILENAME.temp_peaks_per_chrom) \
    <(awk '{print $2}' $FILENAME.temp_fragm_per_peaks) \
    <(awk '{print $2}' $FILENAME.temp_fragm_per_chrom) \
    <(awk '{print $2}' $FILENAME.temp_fragm_per_off_peaks_peak_length) > \
    $FILENAME.stat.txt
rm $FILENAME.temp*
echo "done."

echo -n "`date`: generating insertSizeMetric report on filtered data..."
$PICARD CollectInsertSizeMetrics.jar \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true \
    HISTOGRAM_FILE=$FILENAME.filtered.picard_histogram \
    INPUT=$FILENAME.filtered.bam \
    OUTPUT=$FILENAME.filtered.picard_insertSizeMetric_report \
    2> $FILENAME.picard_stderr.txt
echo "done."

echo -n "`date`: parsing picard output... "
cat $FILENAME.filtered.picard_insertSizeMetric_report | \
    grep -v -P '.*\t.*\t.*\t' | \
    grep -P '^\d|^insert' > \
    $FILENAME.filtered.insertSize_hist
echo "done."

echo -n "`date`: creating plots with quality metrics... "
Rscript $ATAC_QC \
    $FILENAME.peaks_width.hist \
    $FILENAME.stat.txt \
    $FILENAME.filtered.insertSize_hist \
    $FILENAME.QC.pdf \
    2> $FILENAME.R_stderr.txt
echo "done."

echo "END : `date`"


