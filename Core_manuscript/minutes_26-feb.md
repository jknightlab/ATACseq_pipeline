Performance of peak callers on ATAC
----------------------------------------------

**Intermediate results**

- Fseq
- fresh K562
- intersect
- ChromM annotation
- Active promoters and Strong enhancers

1) number of Peaks

- length value of Fseq doesn't seem to effect Enh as much as it affects Prom.
- when we are looking at the number of peaks for the same t value and different
  l values, we loose more peaks (50%) for Prom and 25% for Enh.

2) Looking at SP: bp level; when we increase t parameter we get enrichment for
promoters and do not have any enrichment in enhancers.

3) we want to push for more specificity -- to make sure that our peaks, almost
all of them, are meaningful. Beta 0.5 is good for now (to select the parameter
of Fseq).

4) we Do enrich for DNase peaks as well.

5) F score 0.5 for DNase: first you call TSSs and when you lower the T value
you enrich for artificial peaks and possibly DNase peaks which are not in TSSs
(were not annotated as TSSs)



ToDo
-- signal-to-noise ratio (Irina)
-- sens and spec on peak level -- take 25% overlap (Jason)
-- make sure that sens and spec on peak level is working (Jason)
-- jacard (Jason)
-- frozen
-- macs
-- coverage around TSSs

-- plot

