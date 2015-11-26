ATAC-Seq technical paper
--------------------------------------

Hereby the analysis plan suggested for the paper.

### Overview of the pipeline

The main steps of the primary analysis should contain:

- pre-alignment cleaning of the data (trimming adapters,
removing low quality bases);
- alignment;
- post-alignment cleaning of the data;
- calling peaks on cleaned data using different peak callers
and different parameters;
- annotating peaks (using cell-line specific functional
annotation).

Suggested filtering steps to perform on aligned data:

- remove duplicates;
- remove non-uniquely mapped reads;
- remove reads that were not paired properly;
- remove reads mapped to mitochondrial DNA.


### Peak callers and changing parameters

Current papers on ATAC-Seq used `macs2` to call peak, however,
other peak callers have not been explored (except for ZINBA,
which did not generate convincing results). Therefore, it is
suggested to compare the performance of different peak callers
to be able to draw a conclusion on the best peak caller for
ATAC-Seq data.

Suggested peak callers:

- `macs2`
- `F-Seq`
- `HotSpot`

Also, no data on the effect of different `macs2` parameters on
the results is currently available. It is therefore suggested to
run peak calling with different sets of parameters to be able
to decide what is the best (it is possible that different parameters
will be most suitable for different biological questions/types of
ATAC-Seq data).

Suggested parameters to explore (at the moment -- only for
`macs2`):

- `--call-summits`
- `--bw`
- `--extsize`
- `--slocal`
- `--llocal`
- `--broad`


### Comparison of the samples

It is suggested to compare some technical aspects
of the samples generated per sample (number of peaks,
average length of peaks, average coverage, etc) and to
assess peaks which were differentially called in each
pair of samples (using `macs2 diffpeak`, please scroll
down for more detail).

#### Technical comparison of quality

Metrics to compare:

** before functional annotation **

- mapped reads after filtering
- number of peaks normalized by the number of mapped reads
- number of reads in peaks normalized by the number of mapped reads
- signal-to-noise ratio (number of reads in peaks divided by number
of reads mapped outside peaks to the regions of equal length)
- average peak width
- average peak height (normalized by the number of reads used for
peak calling)

** after functional annotation **

- number of peaks in each functional group (normalized by the
number of reads used for peak calling)
- average peak height in each functional group (normalized by the
number of reads used for peak calling)
- average peak width in each functional group
- enrichment score in each functional group (using `cGAT`)


#### Biological comparison -- differential peak calling

First, look at the influence of different parameters of
`macs2 diffpeak` on the results:

- `--peak-min-len`
- `--diff-min-len`
- `--ignore-duplicate-peaks`

Second, 














