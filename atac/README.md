# ATAC-seq processing and peak calling
Brackets refer to corresponding processing rules in [workflows/](/workflows/).
1. Trim adapters/bad-quality reads using `trim_galore` (`tg_pe` for paired-end, `tg_se` for single-end)
2. Align with `bwa sampe` (BWA-backtrack) (`bwa_pe` for paired-end, `tg_se` for single-end)
3. Keep alignments where (both) reads aligned, are properly paired for paired-end, mapq 10, do not map to mitochondrial DNA or [modENCODE blacklists](https://www.encodeproject.org/comparative/regulation/#Wormset5) (`p10_keep`)
    * For developmental/aging maps, pool technical replicates (e.g. `atac_dm_emb_rep1/2`, `atac_am_d3_rep1/2`)
    * Radomly sample pooled data for pseudoreplicates for IDR peak calling (e.g. `atac_dm_emb_pep1/2`, `atac_am_d3_pep1/2`)
4. Calculate coverage tracks using `macs2` by treating the alignments as:
    * Paired-end, size-selecting for reads <300, <200, <150, <100 (`macs2_pe_ltXXX`)
    * Single-end, extending tags by a fixed amount from the 5' end (`macs2_se_extsizeXXX`)
    * Single-end, focusing on [where the cut sites are](https://github.com/taoliu/MACS/issues/145) by extending and shifting the tags, thereby [effectively smoothing the cut site coverage signal](https://groups.google.com/forum/#!topic/macs-announcement/4OCE59gkpKYs) (`macs2_se_extsizeXXX_shiftmYYY`)
5. Peak calling with [jurgjn/yapc](https://github.com/jurgjn/yapc)
