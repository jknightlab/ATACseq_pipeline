Comparing called peaks in fresh, frozen and fixed ATAC
----------------------------------



|    | fresh | frozen | fixed_day1_inact_only | fixed_day1_inact_decrosslink | fixed_day3_inact_decrosslink | fixed_day7_inact_decrosslink |
| -- | ----- | ------ | --------------------- | ---------------------------- | ---------------------------- | ---------------------------- |
| Number of filtered reads | 11,860,451 | 15,979,159 | 33,108,409 | 26,462,788 | 36,800,607 | 49,502,051 |
| Average fragment width | 494 | 517 | 599 | 640 | 376 | 381 |
| Number of peaks | 10211 | 13602 | 12,987 | 12,949 | 2,043 | 945 |
| Normalized number of peaks | 86.09 | 85.12 | 39.23 | 48.93 | 5.55 | 1.91 |
| Average peak width | 554 | 577 | 786 | 799 | 669 | 672 |
| Number of bases in peaks | 5,660,092 | 7,854,369 | 10,202,624 | 10,349,877 | 1,366,587 | 635,152 |
| Overlap with DNaseI, peaks | 6,640 | 8,571 | 5,410 | 5,275 | 974 | 69 |
| Overlap with DNaseI, bases | 924,498 | 1,198,802 | 764,260 | 748,556 | 130,493 | 9,271 |
| Overlap with DNaseI, % peaks | 16.33 | 15.26 | 7.49 | 7.23 | 9.55 | 1.46 |
| Overlap with DNaseI, % bases | 65.03 | 63.01 | 41.66 | 40.74 | 47.67 | 7.3 |
| Overlap between replicates | 76.02 | 76.65 | 71.12 | 71.76 | 48.09 | 35.05 |


![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_width_distr.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_intensity.png)


ToDo:

- + number of called peaks
- + number of called peaks normalized by library size
- + overlap between replicates
- + average peak intensity
- normalized signal-to-noise ratio
- + peak width distribution
- + % overlap with DNaseI peaks
- % overlap with annotation "on-target"
- % overlap with annotation "off-target"
- % overlap with each functional group of the annotation
- average peak intensity per group
- average peak width per group
- number of peaks overlapping between two material types in a pairwise manner
- pairwise overlap between different types of material and its description (height, width, overlap with functional categories)
- similar analysis of peaks specific for each type of material (height, width, overlap with functional categories)
- sequence analysis of common peaks
- sequence analysis of unique peeaks

-------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
