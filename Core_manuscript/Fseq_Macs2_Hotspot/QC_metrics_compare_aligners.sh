#!/bin/bash

# Folder "Consensus" should contain bed files wit hpeak calling results
for i in `ls Consensus/`
do
	Filename=`echo Consensus/$i`
	Report=`echo QC_metrics/$i.qc.txt`
	TSS_file="/well/bsgjknight/Jason_analysis/atac_k562/annotations/ENCODEsegmentation/wgEncodeAwgSegmentationCombinedK562-TSS.bed"

	echo -n -e "Number_of_peaks\t" > $Report
	cat $Filename | wc -l >> $Report

	All_Peaks=`cat $Filename | wc -l`
	echo -n -e "Average_peak_width\t" >> $Report
	cat $Filename | awk -v var="$All_Peaks" '{sum += $3-$2} END {print sum/var}' >> $Report

	echo -n -e "#peaks_below_200\t" >> $Report
	cat $Filename | awk '$3-$2 < 200 {print $0}' | awk -v var="$All_Peaks" '{sum += $3-$2} END {print sum/var}' >> $Report

	echo -n -e "#peaks_above_200\t" >> $Report
	cat $Filename | awk '$3-$2 >= 200 {print $0}' | awk -v var="$All_Peaks" '{sum += $3-$2} END {print sum/var}' >> $Report

	narrowPeaks=`cat $Filename | awk '$3-$2 < 200 {print $0}' | wc -l`
	widePeaks=`cat $Filename   | awk '$3-$2 >= 200 {print $0}' | wc -l`

	echo -n -e "%peaks_below_200\t" >> $Report
	echo "$narrowPeaks $widePeaks" | awk '{print 100*($1/($1+$2)) "%"}' >> $Report

	echo -n -e "%peaks_above_200\t" >> $Report
	echo "$narrowPeaks $widePeaks" | awk '{print 100*($2/($1+$2)) "%"}' >> $Report

	Peaks_as_TSS=`bedtools intersect -f 0.25 -wa -u -a $Filename -b $TSS_file | wc -l`
	TSSs_recover=`bedtools intersect -f 0.25 -wa -u -b $Filename -a $TSS_file | wc -l`
	TSSs_all=`cat $TSS_file | wc -l`

	## True positive rate -- correctly predicted TSS divided by all TSS
	## TPR=TP/(TP+FN)
	TPR=`echo "$TSSs_recover $TSSs_all" | awk '{print 100*($1/($1+$2)) "%"}'`
	echo -e "TPR\t$TPR" >> $Report

	## False discovery rate -- falsely identified peaks, non-TSS peaks
	## FDR=FP/(TP+FP)
	FDR=`echo $Peaks_as_TSS $All_Peaks | awk '{print 100*(($2-$1)/$2) "%"}'`
	echo -e "FDR\t$FDR" >> $Report

	## Specificity
	## % TSSs out of the Annotation
	Spec=`echo "$TSSs_recover $TSSs_all" | awk '{print 100*($1/($1+$2)) "%"}'`
	echo -e "Specificity\t$Spec" >> $Report

	## Sencitivity
	## % TSSs out of all Peaks
	Sens=`echo "$Peaks_as_TSS $All_Peaks" | awk '{print 100*($1/($1+$2)) "%"}'`
	echo -e "Sensitivity\t$Sens" >> $Report

	## F-Score
	b=1
	FScore=`echo $b $Sens $Spec | awk '{print (1+$1*$1)*($2*$3)/($1*$1*$2+$3)}'`
	echo -e "F-Score\t$FScore" >> $Report

	echo Done for file $i
done

