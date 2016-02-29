Overview of histone marks and their analysis
----------------------------------------------

A histone modification is a covalent post-translational modification (PTM) to
histone proteins which includes methylation, phosphorylation, acetylation,
ubiquitylation, and sumoylation. The PTMs made to histones can impact gene
expression by altering chromatin structure or recruiting histone modifiers.


#### Types of histone marks

Two histone modifications are particularly **associated with active
transcription**:

- `Trimethylation of H3 lysine 4 (H3K4Me3)`. It is an excellent mark of active
  promoters and the level of this histone modification at a gene’s promoter is
  broadly correlated with transcriptional activity of the gene.

- `Trimethylation of H3 lysine 36 (H3K36Me3)`. This trimethylation occurs in
  the body of active genes and is deposited by the methyltransferase Set2.
  This protein associates with elongating RNA polymerase II, and H3K36Me3 is
  indicative of actively transcribed genes. `H3K36Me3` is recognised by the
  Rpd3 histone deacetylase complex, which removes acetyl modifications from
  surrounding histones, increasing chromatin compaction and repressing spurious
  transcription. Increased chromatin compaction prevents transcription factors
  from accessing DNA, and reduces the likelihood of new transcription events
  being initiated within the body of the gene. This process therefore helps
  ensure that transcription is not interrupted.

Three histone modifications are particularly **associated with repressed
genes**:

- `Trimethylation of H3 lysine 27 (H3K27Me3)`. This histone modification is
  depositied by the polycomb complex PRC2. It is a clear marker of gene
  repression, and is likely bound by other proteins to exert a repressive
  function.

- `Di and tri-methylation of H3 lysine 9 (H3K9Me2/3)`. `H3K9Me2/3` is a
  well-characterised marker for heterochromatin, and is therefore strongly
  associated with gene repression. The formation of heterochromatin in yeast
  *Schizosaccharomyces pombe* is initiated by recruitment of the RNA-induced
  transcriptional silencing complex to double stranded RNAs produced from
  centromeric repeats. RITS recruits the Clr4 histone methyltransferase which
  deposits `H3K9Me2/3` (histone methylation). `H3K9Me2/3` serves as a binding
  site for the recruitment of Swi6 (heterochromatin protein 1 or HP1, another
  classic heterochromatin marker) which in turn recruits further repressive
  activities including histone modifiers such as histone deacetylases and
  histone methyltransferases.

- `Trimethylation of H4 lysine 20 (H4K20Me3)`. This modification is tightly
  associated with heterochromatin, although its functional importance remains
  unclear. This mark is placed by the Suv4-20h methyltransferase, which is at
  least in part recruited by heterochromatin protein 1.

Activating acetylation marks (H3K27ac and H3K9ac) are roughly as informative as
activating methylation marks (H3K4me3 and H3K4me2).  H3K79me2 and H3K36me3
marks, both of which mark gene bodies, likely reflecting recruitment of
modification enzymes by polymerase isoforms H3K79me2 occurs preferentially at
the 5' ends of gene bodies and H3K36me3 occurs more 3'.

#### ENCODE histon modifications for K562

| Mark     | Effect         |
| -------- | -------------- |
| H3k4me1  | activation     |
| H3k4me2  | activation     |
| H3k4me3  | activation     |
| H3k9ac   | activation     |
| H3k9me1  | activation     |
| H3k9me3  | **repression** |
| H3k27ac  | activation     |
| H3k27me3 | **repression** |
| H3k36me3 | activation     |
| H3k79me2 | activation     |
| H4k20me1 | activation     |

#### Studying histone marks using NGS

Histone modification was first detected on a genome wide level through the
coupling of chromatin immunoprecipitation (ChIP) technology with DNA
microarrays, termed ChIP-Chip (Barski et al. 2007). However instead of
isolating a DNA-binding transcription factor or enhancer protein through
chromatin immunoprecipitation, the proteins of interest are the modified
histones themselves. First, histones are cross-linked to DNA in vivo through
light chemical treatment (e.g., formaldehyde). The cells are next lysed,
allowing for the chromatin to be extracted and fragmented, either by sonication
or treatment with a non-specific restriction enzyme (e.g., micrococcal
nuclease). Modification-specific antibodies in turn, are used to
immunoprecipitate the DNA-histone complexes (Kouzarides 2007). Following
immunoprecipitation, the DNA is purified from the histones, amplified via PCR
and labeled with a fluorescent tag (e.g., Cy5, Cy3). The final step involves
hybridization of labeled DNA, both immunoprecipitated DNA and
non-immunoprecipitated onto a microarray containing immobilized gDNA. Analysis
of the relative signal intensity allows the sites of histone modification to be
determined (Gibson 2009 229-230; Russell 2010 p. 532).

In order to study histone modifications on a truly genome level, other
high-throughput methods were coupled with the chromatin immunoprecipitation,
namely: SAGE: serial analysis of gene expression (ChIP-SAGE), PET: paired end
ditag sequencing (ChIP-PET) and more recently, next-generation sequencing
(ChIP-Seq). ChIP-seq follows the same protocol for chromatin
immunoprecipitation but instead of amplification of purified DNA and
hybridization to a microarray, the DNA fragments are directly sequenced using
next generation parallel re-sequencing. It has proven to be an effective method
for analyzing the global histone modification patterns and protein target
sites, providing higher resolution than previous methods (Barski et al. 2007;
Gibson 2009 229-232).









----------------------------------
Designed by Irina Pulyakhina, irina@well.ox.ac.uk
