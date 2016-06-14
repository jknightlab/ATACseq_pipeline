Comparing ATAC, DNaseI and FAIRE samples
------------------------------------------

#### Workflow

1. We downloaded DNaseI data from
[here](https://www.encodeproject.org/experiments/ENCSR000EKN/).
Three replicates for K562 were available. Description of samples from the
ENCODE website: K562 (Homo sapiens, adult 53 year female), immortalized cell
line, isogenic replicates, DNase-seq on human K562 treated with SAHA at 1uM for
72 hour. Treaed with 1 μM vorinostat (CHEBI:45716) for 72 hours.

We downloaded FAIRE data from
[here](https://www.encodeproject.org/experiments/ENCSR000DCK/).
Only two replicates for K562 data were available.  Description of samples from
the ENCODE website: K562 (Homo sapiens, adult 53 year female), immortalized
cell line, isogenic replicates, FAIRE-seq on human K562, treated with 100 uM
hydroxyurea for 72 hours.
```
ENCFF000TLN.bam -> faire.rep1.bam
ENCFF000TLP.bam -> faire.rep2.bam
```

2. We downsampled DNaseI and FAIRE data to the amount of reads in fresh
samples.

|                    | Replicate 1 | Replicate 2 | Replicate 3 |
| ------------------ | ----------- | ----------- | ----------- |
| Fresh ATAC         |  11,785,502 |  12,841,034 |  10,954,818 |
| Frozen ATAC        |  15,377,840 |  15,690,118 |  16,869,520 |
| DNaseI             | 166,566,404 | 160,212,700 | 176,522,007 |
| DNaseI downsampled |  11,782,306 |  12,830,604 |  10,967,282 |
| FAIRE              |  64,997,252 |  73,682,365 |      -      |
| FAIRE downsampled  |  11,782,011 |  12,845,269 |      -      |

Downsampling:
```
samtools view -b -s 0.1813 faire.rep1.bam > faire.rep1.fseq_downs.bam
samtools view -b -s 0.1743 faire.rep2.bam > faire.rep2.fseq_downs.bam
```

3. We analyzed DNaseI and FAIRE data with
[pinechrom](https://github.com/jknightlab/ATACseq_pipeline/tree/master/Core_manuscript/Pinechrom)
in the same way we analyzed ATAC data (using F-Seq as a peak caller).


#### Comparison of general metrics and stats

Code:
```
fresh3="./dnase.rep3.pinechrom_general/dnase.rep3.fseq.chr.narrowPeak"
fresh2="./dnase.rep2.pinechrom_general/dnase.rep2.fseq.chr.narrowPeak"
fresh1="./dnase.rep1.pinechrom_general/dnase.rep1.fseq.chr.narrowPeak"

bedtools intersect \
    -f 0.3 -r \
    -a $fresh1 \
    -b $fresh2 | \
    bedtools intersect \
    -f 0.3 -r \
    -a - \
    -b $fresh3 | \
    bedtools sort \
    -i - | \
    bedtools merge \
    -i - > \
    dnase.fseq.intrsct.bed
```

| Parameter                         | Fresh                 | Frozen                | DNaseI                | FAIRE              |
| --------------------------------- | --------------------- | --------------------- | --------------------- | ------------------ |
| Number of filtered reads          |            11,860,451 |            15,979,159 |            11,860,064 |         12,313,640 |
| Number of peaks                   |                10,211 |                13,602 |                27,493 |              1,015 |
| Normalized number of peaks        |                 86.09 |                 85.12 |                231.81 |               8.24 |
| Number of bases in peaks          |             5,660,092 |             7,854,369 |            13,750,508 |            596,555 |
| Average peak width                |                   554 |                   577 |                   500 |                588 |
| Overlap with fresh ATAC, % bases  |                    NA | 76.84% (5,188,735 bp) | 43.47% (4,118,590 bp) |   0.41% (9,628 bp) |
| Overlap with frozen ATAC, % bases | 76.84% (5,188,735 bp) |                    NA | 48.63% (5,092,382 bp) |   0.29% (9,012 bp) |
| Overlap with DNaseI, % bases      | 43.47% (4,118,590 bp) | 48.63% (5,092,382 bp) |                    NA | 1.55% (156,001 bp) |
| Overlap with FAIRE, % bases       |      0.41% (9,628 bp) |      0.29% (9,012 bp) |    1.55% (156,001 bp) |                 NA |
| Overlap between replicates, %     |                76.02% |                76.65% |                79.25% |             52.74% |
| On target, bases, %               |                16.97% |                27.51% |                50.87% |              0.10% |
| Off target, bases, %              |                 0.02% |                 0.03% |                   XXX |                XXX |
| Signal to noise ratio             |                  2.43 |                  2.04 |                  3.22 |                XXX |













OLD

#### Correlation between common peaks

|  macs2-macs2  | macs2-fseq       |
| ------- | ------ |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/macs2_dnase_atac_common_peaks.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_common_peaks.log.png) |
|  macs2-fseq  | macs2-fseq       |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_TSS_common_peaks.log.png) | ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/DNase/dnase_atac_nonTSS_common_peaks.log.png) |



------------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
