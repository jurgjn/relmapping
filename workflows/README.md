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

Once datasets/processing for a particular purpose have settled down, they'll also be assigned shorter names that make semantic sense in the context of the project/paper. They are also cleanly re-processed (or linked) under a separate prefix. In practice, as the data may still further change, *prefix* and *dataset* will also include a three-digit number indicating the year and the week of processing (e.g. `atac814` would refer to year 2018, week 14). The relationships between "raw" and "semantic" dataset identifiers are kept in [workflows/config.yaml](/workflows/config.yaml#L240-L301).  For example, processed coverage data for wild-type L1 stage,  biological replicate #2, as defined on week 14 of year 2018 would be stored as:

- *dataset*: `HS232_L1_JA6_N2_atac_S2` => `atac814_wt_l1_rep2`
- *steps*: `tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.macs2_pe_lt200` => `atac_ce10_spmr_pe`
- *suffix*: `_treat_pileup.bw` => `.bw`
- *prefix*: `atac` => `atac814`

So the final file name would be: `atac814/atac_ce10_spmr_pe/atac814_wt_l1_rep2.atac_ce10_spmr_pe.bw`. 

(This approach makes it possible to use "semantic names" in the downstream analysis code, and re-run the downstream analyses on updated data with minor changes to the input file names.)
