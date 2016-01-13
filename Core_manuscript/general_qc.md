Primary quality control of K562 samples
---------------------------------------

Hereby we explored 12 ATAC samples generated for different conditions of
chromatin for K562 cell line:

- fresh cells (3 replicates)
- frozen cells (3 replicates)
- fixed cells, day 3 after fixation (3 replicates)
- fixed cells, day 7 after fixaton (3 replicates)

The following measures were calculated by running `pinechrom` pipeline:

- Number of mapped reads after filtering
- Average fragment width
- Number of called peaks
- Number of peaks per chromosome
- Number of peaks normalized by total mapped reads
- Number of peaks normalized by reads mapped per chr
- Reads in peaks normalized by total mapped reads
- Signal to noise ratio
- Average peak width

The plots were generated with
[these commands](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/reports_summary_plots.R)
from
[this file](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/reports_summary.txt).

Single PDF file containing all plots in high resolution can be found
[here](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/K562.three_repl.boxplot.QC.pdf).



