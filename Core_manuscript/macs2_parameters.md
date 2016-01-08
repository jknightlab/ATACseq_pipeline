Peak calling with MACS2
-------------------------------

Hereby we investigate peak calling with `macs2` and the influence of its
parameters on the results.

#### MACS2 parameters

There are several parameters which might influence the results directly and
drastically:
- `--call-summits`: if set, MACS will use a more sophisticated signal
  processing approach to find subpeak summits in each enriched peak region.
- `--slocal`: the small nearby region in basepairs to calculate dynamic lambda.
  This is used to capture the bias near the peak summit region. Invalid if
  there is no control data. If you set this to 0, MACS will skip slocal lambda
  calculation. *Note* that MACS will always perform a d-size local lambda
  calculation. The final local bias should be the maximum of the lambda value
  from d, slocal, and llocal size windows. DEFAULT: 1000
- `--extsize`: the arbitrary extension size in bp. When nomodel is true, MACS
  will use this value as fragment size to extend each read towards 3' end, then
  pile them up.  DEFAULT: 200.
- `--bw`: band width for picking regions to compute fragment size. This value
  is only used while building the shifting model. DEFAULT: 300


After running peak calling with a range of values for each of the parameters
(or "TRUE"/"FALSE" for `--call-summits`), we discovered that only `extsize` and
`call-summits` directly influences the number, width and height of called
peaks.


|     |     |
| --- | --- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_macs2_number_of_peaks.png) |   ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_macs2_peak_width.png) |


#### What is a shifting model and `extsize`

A parameter `shiftsize`, later replaced by `extsize`, was initially designed
for single-end ChIP-Seq data.

The simplest approach for calling enriched regions in ChIP-seq data is to take
a direct census of mapped tag sites along the genome and allow every contiguous
set of base pairs with more than a threshold number of tags covering them to
define an enriched sequence region. Such approach would only work effectively
for very precise experiments with extremely strong ChIP enrichment. However,
overall it is not satisfactory due to inherent complexities of the signals as
well as experimental noise and/or artifacts. Additional information present in
the data is now used to help discriminate true positive signals from various
artifacts. Because immunoprecipitated DNA fragments were initially typically
sequenced as single-ended reads, the tags are expected to come on average
equally frequently from each strand, thus giving rise to 2 related
distributions of stranded reads. The corresponding individual strand
distributions will occur upstream and downstream, shifted from the source point
("summit") by half-the average sequenced fragment length, which is typically
referred to as the "shift".

The assumption behind shifting is not that the binding occurs in the middle of
the fragments. It is that the binding occurs within the fragment. The most
likely binding location is determined not from a single fragment, but from the
superposition of multiple fragments.

The goal of the shifting is to find the midpoint of each sequenced fragment.
The distribution of these midpoints should be centered on the binding site of
the protein.

MACS empirically models the shift size of ChIP-Seq tags, and uses it to improve
the spatial resolution of predicted binding sites.

#### Shifting for paired-end data

For SE reads, MACS shifts and extends the reads to build the whole genome
profile. The shift and extend distance is a fixed value estimated from the
double peak pattern or provided by command line parameters. This might be
inaccurate if the wrong shift/extend distance were used, and might cause the
peak appear as doublets. So the peak heights are sensitive to the shift/extend
values. However, paired-end sequencing provides fragment size information,
therefore build whole genome profile is straightforward, no shift or extend
involved.

However, we still need to check the influence of `extsize` on the results.
Suggestions on how to use `extsize` from `macs2` developers:

- to find enriched cutting sites such as some DNAse-Seq datasets, use
  '--nomodel --shift -100 --extsize 200';
- to find nucleosome positioning, use '--nomodel --shift 37 --extsize 73'.


#### Examples of identified peaks

A UCSC session with coverage of one K562 fresh sample
(20150818_Bulk_ATAC_K562_1.dedup.q30f3F12) and peaks identified using different
extsize values can be found
[here](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=pulyakhina&hgS_otherUserSessionName=Core_narrowPeak_extsizes).
From the session, we can see that none of the values works optimally, possibly
due to relatively low coverage. Sometimes smaller extsize helps to identify to
separate peaks while larger extsize values merge multiple peaks into one big
peak. Sometimes smaller extsize subdivides one peak into multiple small ones,
while larger extsize identifies one big peak (as it is supposed to be
identified).

Examples:

|     |     |
| --- | --- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_ASCC3.png) |   ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_intergenic.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_RPE_prom.png) |   ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_triple_peak_low_coverage.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_double_peak.png) |   ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_mult_peaks.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_extsize_SPOPL_double_peak.png) |   |



#### Filtering peaks

**Filtering using macs2 outcome**
As we could see from the UCSC session mentioned above, the major points of
confusion are regions with low/insufficient coverage. Perhaps if we select
significant peaks, peak calling results will be more optimal.

Looking at the narrowPeak file
[description](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) (output of
`macs2`), multiple columns can be used for filtering:

- Column 5 -- score (?)
- Column 7 -- signalValue (measurement of overall average enrichment for the
  region)
- Column 9 -- -log10 qValue
- Column 10 -- peak (point-source called for this peak)

Correlation between columns 5/7 or 5/9 or 7/9 is 0.99. Correlation between
column 5/7/9 and column 10 is below 0.18. So columns 5, 7 or 9 should have
equal contribution to the results and filtering. However, trying to filter
based on either of these columns didn't work out: there is no distribution of
the values (they are rather descrete, so it is hard to choose a sensible
cutoff).

**Filtering using coverage**

Coverage seems to be one of the major issues when calling peaks; regions with
low coverage tend to produce lots of small falty peaks. So, to avoid that, we
created a bedgraph from a bigwig file and selected regions with coverage above
20. When we look at called peaks in regions with coverage > 20, peaks tend to
be less ambiguous (but still not optimal). In that case, `extsize 100` seems
to produce the best results.


Examples:

|     |     |
| --- | --- |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_sign_coverage20_example2.png) |   ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_sign_coverage20_example3.png) |
| ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_sign_coverage20_example4.png) |   ![alt text](https://github.com/jknightlab/ATACseq_pipeline/blob/master/Core_manuscript/fresh1_sign_coverage20_examples1.png) |

#### Conclusion on extsize and filtering peaks

- Select reigons with coverage > 20 (has to be specific per samples).
- Select extsize around 100.





