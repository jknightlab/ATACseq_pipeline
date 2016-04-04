Primary quality control of K562 samples
---------------------------------------


![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/P160113/general_metrics.png)

----------------------------------------

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/P160113/S2N_example.png)

----------------------------------------

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/P160113/S2N_example2.png)

----------------------------------------

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/P160113/overlap_within_group.png)

----------------------------------------

![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/P160113/overlap_between_groups.png)


#### Comparison to DNase data

This table contains the information about the similarity between
ATAC fixed day1 (inact. decrosslink. or inact. only) and ENCODE
DNase of K562. The table shows that, even though the total number
of called peaks in fixed day 1 inact decrosslink was lower (22,692)
than in inact only (24,697), higher fraction of peaks could be
found back in the DNase data with any minimal required overlap tested.

```
I   -- Required number of overlapping base, % of region length
II  -- number of overlapping peaks between DNase and Fixed day 1 INACT DECROSSLINK
III -- number of overlapping peaks normalized by the total number of peaks in Fixed day 1 INACT DECROSSLINK
IV  -- number of overlapping peaks between DNase and Fixed day 1 INACT ONLY
V   -- number of overlapping peaks normalized by the total number of peaks in Fixed day 1 INACT ONLY
VI  -- Ratio between number of overlapping peaks in Fixed day 1 INACT DECROSSLINK divided by INACT ONLY
```

| I    | II    | III     | IV    | V       | VI    |
| ---- | ----- | ------- | ----- | ------- | ----- |
| 0.1  | 21115 | 9305.04 | 22243 | 9006.36 | 1.033 |
| 0.15 | 20721 | 9131.41 | 21853 | 8848.44 | 1.032 |
| 0.2  | 20189 | 8896.97 | 21302 | 8625.34 | 1.031 |
| 0.25 | 19475 | 8582.32 | 20567 | 8327.73 | 1.031 |
| 0.3  | 18584 | 8189.67 | 19669 | 7964.13 | 1.028 |
| 0.35 | 17558 | 7737.53 | 18577 | 7521.97 | 1.029 |
| 0.4  | 16463 | 7254.98 | 17430 | 7057.54 | 1.028 |
| 0.45 | 15305 | 6744.67 | 16226 | 6570.03 | 1.027 |
| 0.5  | 14117 | 6221.14 | 14910 | 6037.17 | 1.03  |

### Conclusion about fixed day 1 inact+decrosslinked OR inact only

**Observations**

Overall inact+decross and inact_only data for fixed_day1 looks very
comparable, however, there are slight differences:

1. Inact+decross data has a lightly higher signal-to-noise ratio
(average across three replicates **5.3**) than inact_only
(average **4.3**).
2. Inact+decross data has a higher normalized average peak intensity
(average **58**) than inact_only (average **53**).
3. Less peaks (normalized by the library size) was called in
inact+decross data compared to inact_only data.
4. Three replicates of inact+decross have higher number of overlapping
peaks between three replicates (above 50% of called peaks) than for
inact_only (below 50%).
5. When compared to fresh or frozen peaks (which we consider gold standard
at the moment), inact+decross have higher overlap (above 60% of called peaks) 
than inact_only (below 60%).
6. When compared to DNase data analyzed with 
[Pinechrom](https://github.com/jknightlab/ATACseq_pipeline/tree/master/Core_manuscript/Pinechrom)
(jknight pipeline to analyze ATAC data), inact+decross shows higher overlap
with DNase peaks than inact_only.

**Conclusion**

Fixed day 1 inact and decrosslinked data looks slightly more similar to
the current gold standards (fresh/frozen ATAC, DNase) and generates cleaner
signal (higher signal-to-noise ratio). Therefore, when one has to choose
between inact_decross or inact_only, we suggest to choose
**inactivated and decrosslinked**.


-------------------
#### Developed by Irina Pulyakhina irina@well.ox.ac.uk
