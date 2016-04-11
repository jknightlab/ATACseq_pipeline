Comparing the performance of peak callers
------------------------------------------

#### Calling peaks

For F-Seq, we called peaks varying parameters **length** and **threshold** in the following range:

- Length = 100 or 200 or 400 or 600 or 800 or 1000 or 2000
- Threshold = 2 or 4 or 6 or 8 or 10 or 12 or 14 or 16

[Code](https://raw.githubusercontent.com/jknightlab/ATACseq_pipeline/master/Core_manuscript/Fseq_Macs2_Hotspot/run_fseq.sh).

For MACS2

For Hotspot


#### Generating consensus

As we had three replicate for each condition (fresh, frozen, etc.), we decided to generate a consensus peak file for each parameter set of each peakcaller. We require two regions to have at least 25% of their length reciprocally overlapping:

```
bedtools intersect -f 0.25 -r -a First_replicate -b Second_replicate | \
bedtools intersect -f 0.25 -r -a - -b Third_replicate > \
consensus_peak_file.bed
```

#### Generating QC metrics

We used the following measurements as our QC metrics to the performance of each peak caller to each other:

- Number_of_peaks
- Average_peak_width
- Number of peaks_below_200
- Number of peaks_above_200
- Per cent of peaks_below_200
- Per cent of peaks_above_200
- True Positive Rate, TPR
- False Discovery Rate, FDR
- Specificity
- Sensitivity
- A combination of sensitivity and specificity, F-Score

[Code](https://raw.githubusercontent.com/jknightlab/ATACseq_pipeline/master/Core_manuscript/Fseq_Macs2_Hotspot/QC_metrics_compare_peakcallers.sh).

More info about some QC metrics:

- True positive rate is calculated as the number of correctly predicted TSS divided by all TSS. Formula: `TPR=TP/(TP+FN)`
- False discovery rate is calculated as the fraction of falsely identified peaks, non-TSS peaks. Formula: `FDR=FP/(TP+FP)`
- Specificity is calculated as a fraction of identified TSS out of all annotated TSS
- Sensitivity is calculated as a fraction of peaks identified as TSS out of all called peaks
- F-Score is calculated as a combination of Sensitivity and Specificity:
```
	beta = 1
	FScore=(1 + beta*beta) * (Sens*Spec) / (beta*beta*Sens + Spec)
```

#### Generating heatmaps for each parameter

For each peak caller, we are testing two parameters. This means that the visualization needs to be 2-dimentional, e.g., through a heatmap.

All the input files -- bed files containing consensus peaks called with each paramater set -- are all called in a similar way:
`condition.aligner.length_Param1Value.tresh_Param2Value.bed.qc.txt`. This makes it easier to collect information for all parameter sets:

For F-Seq:
```
for l in 100 200 400 600 800 1000 2000
do
	echo -n -e "parameters\t"
	ls QC_metrics/fresh.fseq.length_$l.* | \
	    sed 's/.*tresh\_/tresh\=/g' | \
	    sed 's/\.bed.*//g' | \
	    tr '\n' '\t' | \
	    sed 's/\t$//g'
	echo
	echo -n -e "Len=$l\t"
	cat QC_metrics/fresh.fseq.length_$l.* | \
	    grep FDR | \
	    awk '{print $2}' | \
	    sed s/\%//g | \
	    tr '\n' '\t' | \
	    sed 's/\t$//g'
	echo
	echo
done | \
    sort -r | \
    uniq | \
    awk '{print $1 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > \
    fresh.fseq.fdr.unsrt

./sort_file.sh fresh.fseq.fdr.unsrt fresh.fseq.fdr.txt
```

`fresh.fseq.fdr.txt` can now be fed to the heatmap script and a heatmap will be generated in a separate *pdf* file.

[Code](https://raw.githubusercontent.com/jknightlab/ATACseq_pipeline/master/Core_manuscript/Fseq_Macs2_Hotspot/create_heatmap.R).











-------------------------------------------------
developed by Irina Pulyakhina irina@well.ox.ac.uk
