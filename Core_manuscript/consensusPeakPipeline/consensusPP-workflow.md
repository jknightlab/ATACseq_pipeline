## Workflow for Consensus Peak Pipeline
JHendry, 2016/02/03

#### 1. Create consensus peak calls from narrow peak files.
```
Script: run_intrsctNPFs.sh
Requirements: generate-IntrsctUnnMed.R
```
  
  Three narrow peak files (.npf), contained within a single directory are
  combined into consensus peak files (.bed). Three types of consensus
  peak calls are generated:
 ```
    .intrsct.bed (3'/'5 boundaries by intersect of three .npfs)
    .median.bed (3'/5' boundaries by median of three .npfs)
    .union.bed (3'/5' boundaries by union of three .npfs)
```
  See header of run_intrsct_NPFs.sh for more detail.
  
  Run as:
```
  qsub run_intrsctNPFs.sh <dir-with-npfs>
```

#### 2. Interesect consensus peak calls from two samples.
``` 
  Script: run_intrsctSamps.sh
  Requirements: none
```

  Two consensus peak files (.bed) are interesected and the output
  is deposited in ./comparison. The output file is named to
  reflect the samples interesected, and the type of consensus
  boundary used. See header of run_intrsctSamps.sh for more detail.
  
  Run as:
```
qsub run_intrsctSamps.sh <dir-samp1> <dir-samp2> <type boundary:intrsct/median/union>
```
  
#### 3. Perform a basic statistical and graphical analysis of two samples.
``` 
  Script: run_compareSamps.sh
  Requirements: analyzePeaks.R, sampCompare.R
```

  Two consensus peak files are compared and several graphical and
  statistical analysis performed. This script only works on samples
  that have *already* been intersected (e.g., run_intrsctSamps.sh).
  Output includes information regarding extent of global overlap
  between peak calls, total number of peaks called, width and
  density of peak calls. See header of run_compareSamps.sh
  for more detail.
  
  Run as:
```
  qsub run_compareSamps.sh <dir-samp1> <dir-samp2> <type boundary:intrsct/median/union>
```  
