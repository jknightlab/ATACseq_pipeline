Comparing ATAC and DNase samples
---------------------------------------

#### General statistics and parameters

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_va_atac.png)

#### Correlation between common peaks

|  macs2-macs2  | macs2-fseq       |
| ------- | ------ |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/macs2_dnase_atac_common_peaks.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_common_peaks.log.png) |
|  macs2-fseq  | macs2-fseq       |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_TSS_common_peaks.log.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_nonTSS_common_peaks.log.png) |


#### Examples

|         |        |
| ------- | ------ |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_example1.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_example2.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_example3.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_example4.png) |

#### Downsampling DNase bam files

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/atac_dnase_downsampled.png)

#### Overlap with the annotation

| ATAC | DNase |
| ---- | ----- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/ATAC_k562_annotation_all_peaks_pie.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_k562_annotation_all_peaks_pie.png) |
|      |       |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/ATAC_annotated_peaks_classification_pie.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/DNase_annotated_peaks_classification_pie.png) |

**Generalized annotation**

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/ATAC_dnase_generalized_annotated_peaks_classification_pie.png)

**Common peaks between ATAC and downsampled DNase**
![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/downsampled.png)

**Peaks unique for DNase or ATAC**

|     |
| --- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_unique_peak_annotation.png) |
|     |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/unique_peak_width_distr.png) |

Code to generate the last figure with overlaped densities:
```
> plot (atac_dens, xlim=c(0,2000), cex.main=2, main="Unique peak width", cex.axis=2, xlab="Peak width", cex.lab=2, ylab="")
> polygon(atac_dens, col=adjustcolor("mediumpurple",alpha.f=0.3), border="mediumpurple")
> lines (dnase_dens)
> polygon(dnase_dens, col=adjustcolor("aquamarine4",alpha.f=0.3), border="aquamarine4")
> legend(1500, 0.0035, c("ATAC", "DNase"), lty=c(1,1), col=c("mediumpurple", "aquamarine4", lwd=c(5, 5)))
```

#### Fresh ATAC against downsampled ATAC

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/fresh_vs_downsampled_frozen.png)

------------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
