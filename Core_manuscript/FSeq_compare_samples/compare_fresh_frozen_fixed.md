Comparing called peaks in fresh, frozen and fixed ATAC
----------------------------------

#### General statistics

Here we analyze some basic metrics and statistics to compare samples and get
idea of how close each material is to another and whether the data generated
from any material type is meaningful at all.

|  param  | fresh | frozen | fixed day1 inact only | fixed day1 inact decrosslink | fixed day3 inact decrosslink | fixed day7 inact decrosslink |
| ---------------------------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| Number of filtered reads     | 11,860,451 | 15,979,159 | 33,108,409 | 26,462,788 | 36,800,607 | 49,502,051 |
| Average fragment width       | 494        | 517        | 599        | 640        | 376        | 381        |
| Number of peaks              | 10,211     | 13,602     | 12,987     | 12,949     | 2,043      | 945        |
| Normalized number of peaks   | 86.09      | 85.12      | 39.23      | 48.93      | 5.55       | 1.91       |
| Average peak width           | 554        | 577        | 786        | 799        | 669        | 672        |
| Number of bases in peaks     | 5,660,092  | 7,854,369  | 10,202,624 | 10,349,877 | 1,366,587  | 635,152    |
| Overlap with DNaseI, peaks   | 6,640      | 8,571      | 5,410      | 5,275      | 974        | 69         |
| Overlap with DNaseI, bases   | 924,498    | 1,198,802  | 764,260    | 748,556    | 130,493    | 9,271      |
| Overlap with DNaseI, % peaks | 16.33      | 15.26      | 7.49       | 7.23       | 9.55       | 1.46       |
| Overlap with DNaseI, % bases | 65.03      | 63.01      | 41.66      | 40.74      | 47.67      | 7.3        |
| Overlap between replicates   | 76.02      | 76.65      | 71.12      | 71.76      | 48.09      | 35.05      |
| On target, peaks             | 1,425      | 2,969      | 1,567      | 1,527      | 350        | 420        |
| On target, bases             | 960,648    | 2,160,639  | 762,077    | 808,954    | 99,956     | 66,124     |
| Off target, peaks            | 9          | 12         | 47         | 57         | 331        | 465        |
| Off target, bases            | 1,320      | 2,135      | 39,000     | 52,737     | 182,924    | 222,286    |
| On target, peaks, %          | 13.96%     | 21.83%     | 12.07%     | 11.79%     | 17.13%     | 44.44%     |
| On target, bases, %          | 16.97%     | 27.51%     | 7.47%      | 7.82%      | 7.31%      | 10.41%     |
| Off target, peaks, %         | 0.09%      | 0.09%      | 0.36%      | 0.44%      | 16.2%      | 49.21%     |
| Off target, bases, %         | 0.02%      | 0.03%      | 0.38%      | 0.51%      | 13.39%     | 35%        |
| [Signal to noise ratio](https://raw.githubusercontent.com/jknightlab/ATACseq_pipeline/master/Core_manuscript/FSeq_compare_samples/s2n_calculations.txt)  | 3.49  | 2.9   | 2.7   | 3.11  | 4.74   | 2.95|


|    |        |
| -- | ------ |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_width_distr.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/peak_intensity.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/overlap_between_types.png) |   |


![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/overlap_of_types_with_annotation.png)

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/FSeq_compare_samples/on-off_target.png)





ToDo:

- visual analysis of samples
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
- average peak intensity per group
- average peak width per group
- number of peaks overlapping between two material types in a pairwise manner
- pairwise overlap between different types of material and its description (height, width, overlap with functional categories)
- similar analysis of peaks specific for each type of material (height, width, overlap with functional categories)
- sequence analysis of common peaks
- sequence analysis of unique peeaks

-------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
