# Long cap alignment and processing

Current default steps: `trim20.bwa_pe.rm_rRNA.p10_keep.filled_fwd` (`.filled_rev`)

1. Trim both reads to 20bp; follow (Chen et al. 2013) and align using `bwa` in paired-end mode (`trim20.bwa_pe`).
2. Remove rRNA (`rm_rRNA`); filter for mapq10, properly paired, remove mitochondrial, modENCODE blacklisted regions (`p10_keep`)
    - Ribosomal genes have a lot of copies; hence ribosomal RNA contamination is estimated by counting all, including multi-mapping, reads.
3. To monitor transcription elongation as in (Chen et al. 2013), calculate paired-end read coverage with regions between read pairs filled-in (`filled_fwd`, `filled_rev`).
    - With paired-end long cap, the 5' end of the RNA fragment corresponds to the second sequencing read. As `bedtools` determines the strand of a paired-end fragment from the first read, the strandedness argument is set to the opposite of the desired strand.
    (Chen et al. 2013)

(4. TODO Normalisation using DESeq `sizeFactors` -- specific to the particular set of samples...)

## QC metrics
`c_r1`, `c_r2`, `trim20.bwa_pe.c_r1` Number of reads in the read1/read2 fastq file, trimmed-and-aligned .bam file. These should all equal each other as bwa includes all (mapped and unmapped) reads in the output.

`trim20.bwa_pe.c_aln_r1` Number of reads that aligned (anywhere) on the worm genome.

`trim20.bwa_pe.rm_rRNA.c_r1` Number of reads left after masking for ribosomal regions

`trim20.bwa_pe.rm_rRNA.p10_keep.c_r1` Number of reads where (both) reads aligned, are properly paired for paired-end, mapq 10, do not map to mitochondrial DNA or [modENCODE blacklists](https://www.encodeproject.org/comparative/regulation/#Wormset5).

`trim20_1M.bwa_pe_ecoli.c_aln_r1` Number of reads (from a sample of 1M reads randomly chosen from the raw fastq file) that aligned to the *e. coli* genome.
