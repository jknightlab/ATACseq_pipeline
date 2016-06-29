Comparing called peaks in fresh, frozen and fixed ATAC
----------------------------------

**Legend for the names of the samples**

| Sample type                            | Short name |
| -------------------------------------- | ---------- |
| Fresh                                  | F          |
| Frozen                                 | Fr         |
| Fixed day 1 inactivated only           | Fxa1       |
| Fixed day 1 inactivated de-crosslinked | Fx1        |
| Fixed day 3 inactivated de-crosslinked | Fx3        |
| Fixed day 7 inactivated de-crosslinked | Fx7        |


#### General statistics

Here we analyze some basic metrics and statistics to compare samples and get
idea of how close each material is to another and whether the data generated
from any material type is meaningful at all.

|  parameter                   | F          | Fr         | Fxa1       | Fx1        | Fx3        | Fx7        |
| ---------------------------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| Number of filtered reads     | 11,860,451 | 15,979,159 | 33,108,409 | 26,462,788 | 36,800,607 | 49,502,051 |
| Average fragment width       | 494        | 517        | 599        | 640        | 376        | 381        |
| Number of peaks              | 17,655     | 25,464     | 23,919     | 22,844     | 21,822     | 22,675     |
| Normalized number of peaks   | 148.86     | 159.36     | 72.24      | 2.42       | 1.02       | 0.78       |
| Number of bases in peaks     | 20,401,332 | 30,954,626 | 37,835,051 | 35,027,062 | 24,281,923 | 24,272,305 |
| [Signal to noise ratio](https://raw.githubusercontent.com/jknightlab/ATACseq_pipeline/master/Core_manuscript/FSeq_compare_samples/s2n_calculations.txt)  | 5.724 | 4.692 | 3.535 | 4.722 | 0.724 | 0.44 |
| Overlap between replicates   | 61.22%     | 63.14%     | 37.17%     | 46.97%     | 25.95%     | 19.49%     |
| Bases on target, %           | 81.09%     | 72.88%     | 60.4%      | 64.24%     | 22.82%     | 2.95%      |
| Bases off target. %          | 1.33%      | 2.69%      | 9.45%      | 8.49%      | 51.7%      | 79.78%     |
| Bases on TSS, %              | 60.36%     | 48.32%     | 44.13%     | 48.04%     | 15.84%     | 1.21%      |
| Bases mapped on Enhancers, % | 15%        | 16.22%     | 8.75%      | 8.75%      | 4.41%      | 0.57%      |
| Bases mapped on CTCF, %      | 2.04%      | 3.5%       | 3.17%      | 3%         | 1.3%       | 0.95%      |
| Bases mapped on Weak Enh., % | 0.77%      | 1.48%      | 0.4%       | 0.37%      | 0.44%      | 0.09%

# OLD

**Overlap between different material types**

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/overlap_between_types.png)


**Overlap with the annotation**

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/overlap_of_types_with_annotation.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/on-off_target.png)


**Screenshots from UCSC**

|             |            |
| ----------- | ---------- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/hgt_genome_euro_5481_d95ac0.pdf.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/hgt_genome_euro_5956_d95f90.pdf.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/hgt_genome_euro_5ec0_d96600.pdf.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/hgt_genome_euro_681e_d96f20.pdf.png) |



**Observations**

- Number of peaks called in fresh samples and normalized to the number of
  mapped and filtered reads is slightly higher than for frozen samples, almost
  two times higher than in samples fixed for 1 day (both inact only and inact
  decross), more than 20 times higher than in samples fixed for 3 days or 7
  days.

- Average peak width of 550-570 is very similar for fresh and frozen samples,
  much wider (780-800 bp) for both fixed_day1 samples, wider (670-680 bp) for
  fixed_day3 and fixed_day7.

- Both fresh and frozen samples have over 65% bases overlap with DNaseI peaks
  (slightly higher overlap for fresh samples), only 40% for both fixed_day1
  samples. Surprisingly the overlap between fixed_day3 with DNaseI is much
  higher than for fixed_day7 (47% against 8%) and even higher than for
  fixed_day1 (47% against 40%).

- Overlap between replicates is similarly high for fresh and frozen samples
  (76-77%), slightly lower for both fixed_day1 samples (70-72%) and much lower
  for fixed_day3 and fixed_day7 (below 50%).

- Number of bases mapped on target is higher for frozen samples (27%) than for
  fresh (17%). This percentage is comparably low for all fixed samples (7-10%).

- Number of bases mapped off target is a bit lower for fresh than for frozen
  samples (0.02% against 0.03%). Fixed_day1 samples have 10 times higher peak
  calling of target, fixed_day3 and fixed_day7 -- 100 times higher.

- Signal-to-noise ratio is the highest for fresh samples, similarly high for
  frozen and fixed_day1 inact+decross, much lower for fixed_day1_inact only
  (2.06 against 1.39).

- Overlap between samples of different types of material: the highest between
  fixed_day1 inact+decross and fixed_day1_inact (78.7%); fresh and frozen
  (72.1%); fixed_day1 inact+decross and fresh/frozen (60.7% / 64%).

- Looking at the overlap with annotation, fixed_day3 still looks a bit like
  fresh samples, but fixed_day7 looks very different.

- Looking at the overlap with annotation, fresh samples have higher percentage
  of TSS and frozen samples have higher percentage of enhancers, weak enhancers
  and CTCFs.

- Looking at the overlap with annotation, both fixed_day1 samples look almost
  identical.

- Looking at the images from UCSC, it looks like the signal is the same in all
  samples, sometimes stronger in fresh, sometimes stronger in frozen samples.
  And it gets gradually diluted from fixed_day1 to fixed_day3 and fixed_day7.

- Fresh and frozen samples look the best. Some metrics look better for fresh,
  some look better for frozen.

- Fixed_day1 inact+decross looks slightly better than fixed_day1 inact only.

- Fixed_day3 and day7 contain very noisy, diluted signal.


#### Unique peaks

This is what happens with the distribution of peak categories in fresh ATAC
sample when we subtract peaks overlapping between fresh and frozen samples,
fresh and fixed day1, etc. A big difference happens when we subtract frozen
peaks from fresh peaks, after that, when we remove fixed peaks, not much
change. Also, the main change happens in the number of TSS, they are over 50%
of peaks that we loose when we subtract frozen peaks.

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/fresh_minus_frozen-fixed.png)


This figure shows the number and percentage of peaks with different annotations
which are unique in pairwise comparisons. E.g., fresh vs frozen means that
these peaks are remaining in a fresh sample when we subtract frozen, but we
have no information whether they are also present in any over peak list. We can
see that subtracting fresh/frozen/fixed_day1 peaks from fixed_day3/fixed_day7
leaves us with almost no peaks. We can also see that fresh and frozen samples
look more similar to each other than to fixed_day1 samples. Both fresh and
frozen samples have a lot more enhancers than fixed_day1.

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/unique_pairwise_peaks.png)


On this plot we can see the annotation of peaks unique for each material type
compared to all other types of material. Firstly, we can see that regardless of
a very similar number of called peaks in fresh (~11,000), frozen (~13,000) and
both fixed_day1 samples (~12,000), frozen samples contain many more unique
peaks than fresh or fixed samples. We can also see than fixed day1 samples
contain more unique peaks than fresh; and that fixed day1 inact_only contains
more unique peaks than fixed_day1 inact+decrosslink. We can see that unique
frozen peaks are enriched for enhancers and weak enhancers, while unique
fixed_day1 peaks are enriched for TSS (and less for CTCF).


![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/unique_peak_annotation.png)


#### Peaks overlapping between different types of material

This bar plot contains the annotation-based distribution of peaks common
between each two types of materials. We can once again appreciated that frozen
samples are enriched with enhancers and the fraction of overlapping enhancers
is the highest for frozen samples. It is higher when we compare fresh and
frozen samples than when we compare fixed and frozen samples. We can also
appreciate that the number of overlapping peaks is very similar when we compare
overlaps of fresh or frozen samples with fresh/frozen/fixed_day1. However, when
we look at fixed_day1 overlaps with the other fixed_Day1 sample and
fresh/frozen, two fixed_day1 samples are more similar than a fixed_day1 and a
fresh/frozen sample.  From the Y axis we can also appreciate that both
fixed_day3 and fixed_day7 have a significantly lower overlap with other types
of mterial. 


![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/overlap_between_material_types.png)


**Peak intensity per group**

This figure shows the distribution of peak intensities for each of the main
categories of RE -- CTCF, Enhancers, Weak Enhancers and TSS. We can appreciate
that both fresh and frozen samples contain many more peaks with higher coverage
than both fixed day1 samples. We can also appreciate that for each category
fixed_day1 inact+decross has more peaks with higher intensity than fixed_day1
inact_only.

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_intensity_on-target.png)

On this density plot we see the distribution of peak intensity for all peaks
annotated as Enhancers called in frozen samples and peaks annotated as
enhancers found uniquely in frozen samples. We can see almost no
difference/shift in the two distributions, just a slight distribution to higher
intensity values for unique Enhancers.

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_intensity_frozen_enh.png)



#### Peak width for each of the main annotation categories

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/width_distr_of_RE.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_width_average.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_width_CTCF.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples//peak_width_Enh.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_width_WE.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_width_TSS.png)







#### ToDo:

- + visual analysis of samples
- + number of called peaks
- + number of called peaks normalized by library size
- + overlap between replicates
- + average peak intensity
- + normalized signal-to-noise ratio
- + peak width distribution
- + % overlap with DNaseI peaks
- + % overlap with annotation "on-target"
- + % overlap with annotation "off-target"
- + % overlap with each functional group of the annotation
- + average peak intensity per group
- + average peak width per group
- + number of peaks overlapping between two material types in a pairwise manner
- + pairwise overlap between different types of material and its description (height, width, overlap with functional categories)
- + similar analysis of peaks specific for each type of material (height, width, overlap with functional categories)
- sequence analysis of common peaks
- sequence analysis of unique peeaks
- Downsample different types of material and see how the distribution of peaks across categories changes

-------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
