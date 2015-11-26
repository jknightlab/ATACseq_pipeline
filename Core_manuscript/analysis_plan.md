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

macs2 diffpeak















