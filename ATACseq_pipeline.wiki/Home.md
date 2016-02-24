# ATAC-Seq

## Overview

 ATACseq is a technique used to identify open versus closed
chromatin regions (depleted on nucleosomes and high chromatin
structure) in a context specific way (in particular cell types
from controls and patients). The output of ATACseq is similar
to the information provided by the well-established DNase I
hypersensitivity but the former one has the advantage of reduced
number of cells to be implemented (50,000 cells). This enables
to apply it in clinical samples and allows comparison of the
chromatin landscape between healthy and diseased individuals.
High-density sequencing using Hiseq output will also allow to
develop footprints for transcription binding sites. Overlapping
of ATACseq information in cells from patients with other
annotation information such as K27Ac, K4me1,K4me3, RNAseq, DHS,
eQTLs...will help to prioritise those SNPs within loci identified
to be associated to disease in different GWAS which may have a
functional role in regulation of gene expression. These candidate
SNPs will be suitable for later confirmation using the
high-throughput chromosome conformation capture technique known as
Capture-C.

| **Project ID**               |  ATAC-seq  |
|:-----------------------------|:-----------|
| **Status**                   | Ongoing    |
| **Dates**                    | 2014-2017  |
| **Publication**              | Unpublished|
| **Ethics**                   | Psoriasis and ankylosing spondylitis patients recruitment |

| Sample information           | &nbsp;|
|:-----------------------------|:------|
| **Organism**                 | Homo sapience  |
| **Organism part**            | Blood and skin |
| **Cell type**                | T cells (CD4+ and CD8+), monocytes, keratinocytes and others; K562 (collaboration with the Core) |
| **Experimental conditions**  | healthy controls and patients |
| **Phenotype**                | Psoriasis and ankylosing spondylitis patients |
| **Sample source**            | ATACseq_001 folder contains Miseq NGS data from blood-isolated human monocytes using CD14+ Miltenyi magnetic beads and cultured in RPMI 20% FBS for 24 hours  |


## Available data
| Sequencing project ID | Data type | Platform | Repository name *1 |Repository ID *2| Data release status | Analysis *3|Script Repository|
|:---------------------:|:---------:|:--------:|:------------------:|:--------------:|:-------------------:|:------------------------------:|:---------------:|
| XXX      | ATAC sequencing | Illumina Miseq | rescomp1: /well/jknight/Irina/ATAC | ???  | NA | analysis with PineChrom, ATAC-Seq pipeline | https://github.com/jknightlab/ATACseq_pipeline |
| P150493  | ATAC sequencing | Illumina Hiseq2500 | rescomp1: /well/bsgjknight | ???  | NA | analysis with PineChrom and Core's ATAC pipeline | https://github.com/jknightlab/ATACseq_pipeline/tree/master/Core_manuscript |


- *1 - Name of public data repository to which data has been submitted
- *2 - ID assigned by repository
- *3 - Brief description of analysis. Include links to scripts in Git repository and pages describing results where possible.


| Contact information          | &nbsp;|
|:-----------------------------|:------|
| **Contact1, name and email** | Irina Pulyakhina, irina@well.ox.ac.uk      |
| **Contact2, name and email** | Alicia Lledo Lara, alicia@well.ox.ac.uk    |
| **Contact3, name and email** | Anna Sanniti, asanniti@well.ox.ac.uk       |
| **Contact4, name and email** | Julian Knight, julian.knight@well.ox.ac.uk‎ |
