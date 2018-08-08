There are at least [three hard problems](https://martinfowler.com/bliki/TwoHardThings.html) when working with genome-wide sequencing data: [naming things](https://www.quora.com/Why-is-naming-things-hard-in-computer-science-and-how-can-it-can-be-made-easier), [rerunning modified code](http://lab.loman.net/uncategorized/2013/02/15/lomans-law-of-bioinformatics/) (and off-by-one errors). The pipelines in this repository try to address this by using snakemake, and organising data files using the following schema.

Every data file has four parts: *dataset*, *step*, *suffix* and *prefix*. These are used in file names as follows:
`~/relmapping/{prefix}/{steps}/{dataset}.{steps}{suffix}`

Example of an SPMR-normalised coverage track an ATAC-seq sample:
- *dataset*: `HS232_L1_JA6_N2_atac_S2`
- *steps*: `tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.macs2_pe_lt200`
- *suffix*: `_treat_pileup.bw`
- *prefix*: `atac`

*Dataset* is a chunk of the original fastq file name ("as it came off from the sequencer") that is sufficient to identify the raw data.

*Steps* is a dot-separated string of snakemake rules that were sequentially applied to the raw data obtain the file. Typically, each step runs one main "tool", preceded/followed by tailored pre-/post-processing. For example, [`macs2_pe_lt200`](https://github.com/jurgjn/relmapping/blob/master/workflows/macs2.snakefile#L58-L82)  takes a .bam-file with a paired-end alignment, size-selects for reads <200bp, uses MACS2 to call peaks+calculate normalised coverage, and converts the resulting signal to BigWig using kentUtils.

*Suffix* typically just has the file extension. In cases where the step produces multiple output files (e.g. MACS2 peak calls and coverage tracks) it specifies something to distinguish these different outputs.

*Prefix* specifies a top-level directory to e.g. separate processing for each assay: `atac`, `lcap`, `scap`, ...
