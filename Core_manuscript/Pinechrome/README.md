PineChrom -- ATAC-Seq analysis pipeline
--------------------------------------

**Pi**peli**ne** for **chrom**atin analysis using ATAC-Seq,
a pipeline to analyze ATAC-Seq data.

This module of **PineChrom** is starting the analysis with
the filtered *bam* files. It generates *bigwig* files, calls
peaks, performs quality control per sample and generates some
statistical reports as tables and figures.

These are the instructions on how to run the pipeline:

```
pinechrom <BAM_FILENAME> <BAM_LOCATION> <FOLDER>
```

where:
```
<BAM_FILENAME>  - name of the bam file without ".bam" (e.g., if
                  the bam file is called "sample1.bam", specify
                  "sample" as BAM_FILENAME)
<BAM_LOCATION>  - full path to the bam file
<FOLDER>        - name of the output folder; the pipeline will
                  create all intermediate and final results and
                  files in that folder
```


----------------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk

