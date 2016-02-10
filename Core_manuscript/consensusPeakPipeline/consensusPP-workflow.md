## Workflow for Consensus Peak Pipeline
---------------------------------------
JHendry, 2016/02/03
---------------------------------------

When we have multiple technical or biological replicates generated for one
sample, we want to get the consensus of peaks called in each replicate, e.g., a
peak present in all replicates. Later on we will compare such consensus peak
lists between conditions/samples.


#### 1. Create consensus peak calls from narrow peak files.
```
Script: run_intrsctNPFs.sh
Requirements: generate-IntrsctUnnMdn.R
```
  
  Three narrow peak files (.npf), contained within a single directory are
  combined into consensus peak files (.bed). Three types of consensus
  peak calls are generated:
 ```
    .intrsct.bed (3'/'5 boundaries by intersect of three .npfs)
    .median.bed (3'/5' boundaries by median of three .npfs)
    .union.bed (3'/5' boundaries by union of three .npfs)
```

Intersection, median and union of peaks are defined as follows:

```
## Schematic Eg.
##  Rep1            !===========:
##  Rep2                 *============!
##  Rep3               :======*
##
##  intersect            *====*  
##  median             :========:  
##  union           !=================!
```
  See header of `run_intrsct_NPFs.sh` for more detail.
  
  Run as:
```
  qsub run_intrsctNPFs.sh <dir-with-npfs>
```

#### 2. Interesect consensus peak calls from two samples.

After we generate consensus peak list for each of the two samples we wan to
compare, we assess the number of peaks shared/common for the two samples.

``` 
  Script: run_intrsctSamps.sh
  Requirements: none
```

  Two consensus peak files (.bed) are interesected and the output
  is deposited in ./comparison. The output file is named to
  reflect the samples interesected, and the type of consensus
  boundary used. See header of `run_intrsctSamps.sh` for more detail.
  
  Run as:
```
  qsub run_intrsctSamps.sh <dir-samp1> <dir-samp2> <type boundary:intrsct/median/union>
```
  
#### 3. Perform a basic statistical and graphical analysis of two samples.

After we identify peaks shared between conditions, we perform some statistical
analysis and generate plots to assess the quality of the data.

``` 
  Script: run_compareSamps.sh
  Requirements: analyzePeaks.R, sampCompare.R
```

  Two consensus peak files are compared and several graphical and
  statistical analysis performed. This script only works on samples
  that have *already* been intersected (e.g., `run_intrsctSamps.sh`).
  Output includes information regarding extent of global overlap
  between peak calls, total number of peaks called, width and
  density of peak calls.
  
`compareSamples.R` will compare global overlap in peak calls between two
samples by determining number of intersecting peaks that satisfy a given
threshold (#bps, reciprocal percent, both).

`analyzePeaks.R` will plot side-by-side of total number of peaks,
peaks/chromosome, length and density distribution of peaks for both samples.

See header of `run_compareSamps.sh` for more detail.
  
  Run as:
```
  qsub run_compareSamps.sh <dir-samp1> <dir-samp2> <type boundary:intrsct/median/union>
```  
