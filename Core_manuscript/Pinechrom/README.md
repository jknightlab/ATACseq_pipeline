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
   represented by one or multiple samples. This module creates consensus peak lists for each condition, extracts peaks unique per condition, runs differential expression analysis on the common peaks and checks for the existences linear correlation between the intensities common peaks;
- `pinechrom_genes`


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
