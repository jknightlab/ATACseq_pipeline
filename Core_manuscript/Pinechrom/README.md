PineChrom -- ATAC-Seq analysis pipeline
--------------------------------------

**Pi**peli**ne** for **chrom**atin analysis using ATAC-Seq, a pipeline to
analyze ATAC-Seq data. It consists of multipe independant modules that can be
run separately. If you want to run just one module, make sure that you follow
the instructions on how your input files should be called. **Pinechrom**
consists of three modules:

- `pinechrom_general` is run per sample, it performs read alignment, filtering
   of mapped reads, quality checks and generates statistical overview of the
   quality, as well as some informative plots;
- `pinechrom_compare` is run for a pair of conditions, each condition can be
   represented by one or multiple samples. This module creates consensus peak
   lists for each condition, extracts peaks unique per condition, runs
   differential expression analysis on the common peaks and checks for the
   existences linear correlation between the intensities common peaks;
- `pinechrom_genes`

#### Prerequisites

In order for script to run without errors and recognize all
necessary tools, please specify the paths to the following
tools (within the script):

- zcat
- FastQC
- cutadapt
- bwa version 0.7.10
- bwa index files of human reference hg19
- samtools version 1.2
- picard version 1.1 and higher
- bedtools version 2.2 and higher
- bedGraphToBigWig
- file containing chromosome names (as specified in the bwa index file) and
  chromosome lengths
- macs2
- bamtools version 2.3.0
- blacklist region -- regions that should be removed from consideration due to
  a number of reasons

You also need to specify path to the following in-house developed scripts:

- ATACseq\_qc.R
- bed\_extract\_fragments.pl



#### Module "pinechrom_general"

This module of pinechrom, pinechrom_general, will run per-sample analysis on
your ATAC-Seq samples. Starting from fastq or bam files, it will perform
quality of mapped reads, statistical assesement of the quality of the sample,
peak calling and statistical assesement of called peaks.

pinechrom_general is flexible and can run different set of analysis steps
depending on your input. You can run it in three modes:
```
    pinechrom_general --prefix PREFIX --fastq /path/to/fq/read1.fastq.gz,/path/to/fq/read2.fastq.gz
```
In case you want to start your analysis from raw fastq files, perform quality
checks of fastq files, run the alignment, filter good-quality reads and call
peaks, you should use this command. **--prefix** will be the name of the folder
that pinechrom_general will create and run all the analysis in, as well as the
prefix of all output file names. Please make sure that your fastq files:
  - are gzipped;
  - are separated by comma without spaces (`file1,file2` NOT `file1, file2`);
  - you specify the *full* path to your fastq files, even if they are located
    in the same directory as the executable of pinechrom.
```
    pinechrom_general --prefix PREFIX --bam /path/to/bam/sample.bam
```
In case you already generated bam files with mapped reads and want to proceed
with filtering reads based on their quality, calling peaks and performing
statistical analysis, you should use this command. You should specify the *full*
path to your bam file using **--bam**.
```
    pinechrom_general --prefix PREFIX --bam_filt /path/to/bam/sample.filt.bam
```
In case you performed filtering mapped reads yourself and use pinechrom to
generate statistical reports, call peaks and assess the quality of the peaks,
you should use this command. The pipeline will start with creating a bigwig
file for your input filtered bam file, perform statistical analysis, call peaks
and perform statistical analysis on them. You should specify the *full* path to
your filtered bam file using **--bam_filt**.

The following procedures are carried out:

- removing adapter sequences from the raw fastq files (`--fastq`);
- aligning trimmed fastq files to the human reference hg19 with `BWA`
  (`--fastq`);
- filtering the alignment file: removing duplicates, alignments with mapping
  quality < 30, non-uniquely mapped reads, not properly paired reads, reads
  mapping to non-conventional chromosomes or mitochondrial DNA (`--fastq` /
  `--bam`);
- generating BigWig file with peaks of ATAC signal (`--fastq`, `--bam`,
  `--bam_filt`);
- calling peaks on filtered alignments using macs2 (`--fastq`, `--bam`,
  `--bam_filt`);
- running quality controls, such as number of mapped reads at every stage of
  filtering the alignments, generating plots with number of reads mapped to
  peaks, signal to noise ratio, peak width (`--fastq`, `--bam`, `--bam_filt`);


#### Module "pinechrom_compare"

This module of pinechrom, pinechrom_compare, will compare peaks called under
each of the two conditions. Starting from bam files and bed files containing
called peaks, it will (if applicable) generate a consensus peak file for each
condition, compare two consensus peak files and generate a list of peaks
called only under the first condition, list a list of peaks called only under
the second condition, a list of peaks common for both conditions and
differentially expressed and a list of peaks common for both conditions and
expressed at the same level.

pinechrom_compare requires bam files containing mapped reads (per replicate)
and bed files containing called peaks (also per replicate).

|     |     |
| --- | --- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/Pinechrom/pinechrom_general_schema.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/Pinechrom/pinechrom_compare_schema.png) |


----------------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
