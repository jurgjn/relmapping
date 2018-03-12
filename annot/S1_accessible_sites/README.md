# S1 Identification of accessible sites

[`S1a_accessible_sites.tsv`](S1a_accessible_sites.tsv) locations and signal metrics of accessible sites
- `atac_mode` location of peak accessibility within an accessible site; useful for aligning a group of sites for downstream analyses such as heatmaps or aggregate plots
- `atac_source` source data used to define the location of the accessible site; one of `atac_wt_pe`, `atac_wt_se` or `atac_glp1_se`
- `atac_%stage_height` stage-specific peak accessibility within an accessible site; used as input for cluster analyses
- `atac_%stage_%rep_count` stage/replicate-specific read counts; used as input for statistical tests of differential accessibility

[`S1b_differential_accessibility_wt.tsv`](S1b_differential_accessibility_wt.tsv) differential accessibility tests
- two-sided tests for differential read counts from `%stage1` to `%stage2` (all pairwise combinations of wild-type stages, `n=15`), as reported by [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html):
  - `atac_%stage1_to_%stage2_log2FoldChange` log2 fold change
  - `atac_%stage1_to_%stage2_padj` FDR-corrected p-value
- `atac_%stage1_to_%stage2_incr` statistically significant increase in accessibility from `%stage1` to `%stage2` with at least two-fold change (`log2FoldChange > 1`) at `padj < 0.05`
- `atac_%stage1_to_%stage2_decr` statistically significant decrease in accessibility from `%stage1` to `%stage2` with at least two-fold change (`log2FoldChange < -1`) at `padj < 0.05`

[`S1c_accessibility_patterns_wt.tsv`](S1c_accessibility_patterns_wt.tsv) summary metrics characterising developmental patterns
- `atac_is_dynamic` either set to:
  - `dynamic` if the peak tested positive for any change (increase or decrease) between any two wild-type stages
  - `constitutive` if no change between any stages
- `atac_peak_accessibility` peak accessibility, obtained by coding the stage as an integer (emb=1, L1=2, ..., ya=6), and calculating the weighted average with stage-specific peak heights as weights; useful for a coarse relative ranking of dynamic accessibility patterns with a lower value corresponding to "earlier" peak accessibility, and a higher value corresponding to "later" peak accessibility
- `atac_accessibility_cluster_id` cluster id of clustering peak heights of developmentally dynamic peaks using _k_-means:
  - constitutive peaks (`atac_is_dynamic` set to `constitutive`) are set to `0`
  - clusters id-s of clusters containing dynamic peaks are ranked by mean peak accessibility of all sites within the cluster (as defined by `atac_peak_accessibility`)

[`S1c_accessibility_patterns_wt.bed`](`S1c_accessibility_patterns_wt.bed`) peaks coloured by whether they are identified as `dynamic` (black) or `constitutive` (orange)
