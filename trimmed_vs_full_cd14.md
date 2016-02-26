Trimmed vs full-length reads
----------------------------------

In order to see how much read length influences the outcome
of an ATAC-Seq experiment, we run `pinechrom` on two samples:
a full length CD14 monocytes sample sequenced at 80 bp length
and the same sample with shorter, trimmed reads of 55 bp each.
3' end of each end of a read pair was trimmed.

#### Results: general characteristics

As can be seen from this table, trimming the fastq files leads
to fewer mapped reads, hence the reduced number of called peaks.
Also, as expected, when peak calling is done on shorter reads,
average peak width decreases (267 bp vs 261 bp).

|                                   | 80 bp      | 55 bp | Difference|
| --------------------------------- | ---------- | ---------- | ---- |
| mapped reads after filt           | 29,867,693 | 29,197,636 | 2.27 |
| average fragm width               | 545.545    | 534.34     | 2.08 |
| total called peaks                | 46,015     | 45,455     | 1.22 |
| average peaks per chrom           | 2000.65    | 1976.3     | 1.22 |
| peaks norm mapped reads           | 66.98      | 67.6871    | 1.04 |
| peaks norm mapped reads per chrom | 3254.81    | 3340.02    | 2.58 |
| reads in peaks norm mapped reads  | 3945.15    | 3962.36    | 0.44 |
| signal to noise                   | 5.99       | 6.07       | 1.28 |
| average peak width                | 267.06     | 261.21     | 2.22 |

#### Results: number of overlappping peaks

When we look at the number of overlapping peaks, we get over
96% overlap which is very concinsing considering the fact that
we trimmed 30% of read length and had 2.5% less mapped reads.

| bedtools par   | #peaks | % in 80bp | % in 55bp |
|:--------------:|:------:|:---------:|:---------:|
| at least 1 bp  | 45,288 | 98.41%    | 99.63%    |
| 10% reciprocal | 45,267 | 98.37%    | 99.58%    |
| 20% reciprocal | 45,111 | 98.03%    | 99.23%    |
| 30% reciprocal | 44,866 | 97.5%     | 98.7%     |
| 40% reciprocal | 44,615 | 96.95%    | 98.14%    |
| 50% reciprocal | 44,286 | 96.23%    | 97.42%    |

Names of columns: `bedtools par` -- bedtools parameter;
`#peaks` -- number of overlapping peaks; `% in 80bp` --
% of all peaks called in 80bp files; `% in 55bp`-- % of
all peaks called in 55bp files.

#### Results: reads which make a difference

The difference between the number of mapped reads before and
after trimming is around 700,000 reads. One may think that
after trimming some reads become too short and are impossible
to map uniquely. However, when we compare the names of mapped
and unmapped reads, we get *1,203,435* mapped only when they
are full length, and *533,378* mapped only when they are
trimmed. One explanation can be that some reads had bad quality
bases on their 3' end and after trimming those bases we manage
to map them. Another explanation can be that some reads are of
poor/medium quality, but after we remove 25 bases we are able to
map the remaining "bad" reads and produce artificial alignments.
The second explanation is not likely, as we filtered mappings
based on quality and only selected high quality alignments.


-------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
