PineChrom -- ATAC-Seq analysis pipeline
--------------------------------------

**Pi**peli**ne** for **chrom**atin analysis using ATAC-Seq,
a pipeline to analyze ATAC-Seq data.

### Usage

qsub ./run_atac_pipeline.sh fastq1.gz fastq2.gz output_prefix

A folder `output_prefix.atac_analysis` will be generated
and output files with common prefix `output_prefix` will
be created inside it.


### Requirements

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
- file containing chromosome names (as specified in the bwa index file) and chromosome lengths
- macs2
- bamtools version 2.3.0
- blacklist region -- regions that should be removed from consideration due to a number of reasons

You also need to specify path to the following in-house developed scripts:

- ATACseq_qc.R
- bed_extract_fragments.pl


### Overview of the pipeline

The following procedures are carried out:

- removing adapter sequences from the raw fastq files;
- aligning trimmed fastq files to the human reference hg19 with `BWA`
- filtering the alignment file (removing duplicates, alignments with 
  mapping quality < 30, non-uniquely mapped reads, not properly
  paired reads, reads mapping to non-conventional chromosomes or
  mitochondrial DNA)
- generating BigWig file with peaks of ATAC signal
- calling peaks on filtered alignments using macs2
- running quality controls, such as number of mapped reads at every
  stage of filtering the alignments, generating plots with number of
  reads mapped to peaks, signal to noise ratio, peak width


### Future plans

In the nearest future we are planning to implement:

- generation of heatmap showing enrichment of different fragment lengths
  over different types of regulatory elements
- annotating fragments and calculating enrichment of each category with
  `cGAT`


---
designed by I. Pulyakhina
