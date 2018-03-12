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
5. Peak calling (`concave_idr`)
    1. Identify representative peak locations (refered to as oracle peaks in IDR2) via concave regions on mean coverage calculated across all stages and replicates
    2. Determine statistical significance by running IDR on stage-specific concave regions, but using the representative peak locations from the previous step as the oracle peak list
        * Peaks are ranked by the mean second derivative of the coverage
        * As advised in the IDR documentation, this done using replicates and pseudoreplicates
    3. Aggregate peak calls from the individual stages, use an IDR cutoff around 0.01/0.03/0.05 ([on globalIDR, not localIDR](https://groups.google.com/forum/#!topic/idr-discuss/FY2K5VKx8AQ)) in any stage as a threshold for statistical significance

## Additional notes/references

### IDR
- Broad overview of the IDR workflow -- replicates, pseudoreplicates, self-pseudoreplicates, etc -- in [(Landt et al. 2012)](https://dx.doi.org/10.1101%2Fgr.136184.111) around [Fig 7](http://genome.cshlp.org/content/22/9/1813/F7.expansion.html).  
- There are two IDR implementations: [an initial R version](https://sites.google.com/site/anshulkundaje/projects/idr), and a [v2 python rewrite](https://github.com/nboley/idr/) ([available in bioconda](https://anaconda.org/bioconda/idr)). The latter [has several improvements](https://groups.google.com/forum/#!topic/idr-discuss/A7PaMnzoFwg) and is being used by the current version of this pipeline.
- `idr` script in IDR v2 is the equivalent of `batch-consistency-analysis.r` in the initial R version; [it should be used in a similar fashion](https://groups.google.com/forum/#!topic/idr-discuss/_a_GKfw7kwM) (replicates, pseudoreplicates, self-pseudoreplicates)
- An explanation of [globalIDR vs localIDR](https://groups.google.com/forum/#!topic/idr-discuss/FY2K5VKx8AQ)
    > Use the global IDR for thresholding. Essentially the local IDR is akin to the posterior prob. of a peak belonging to the irreproducible noise component. The global IDR is analogous to a multiple hypothesis correction on a p-value to compute an FDR.
- Replicate-specific peaks/peak scores can be replaced with [a metric calculated from a replicate-specific coverage track at regions defined by the oracle peaks](https://groups.google.com/forum/#!topic/idr-discuss/HuS9RyLQsi4)
- Interpreting IDR v2 [diagnostic plots](https://groups.google.com/forum/#!topic/idr-discuss/KCUPDALSBmI)
- Weak but highly correlated peaks [can be problematic with IDR](https://sites.google.com/site/anshulkundaje/projects/idr#TOC-CALL-PEAKS-ON-INDIVIDUAL-REPLICATES)
    > truncate the number of peaks to the top 100k-125k. Using more than this simply increases the running time of the IDR analysis with no advantage. Infact using more peaks with MACS2 can cause problems with the IDR model because MACS2 seems to produce strange highly correlated peak scores for very weak and noisy detections. This can confuse the IDR model.
