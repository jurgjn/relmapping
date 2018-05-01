# Regulatory annotation of accessible sites
[`Fig2D2_regulatory_annotation_Apr27_ce10.bed`](Fig2D2_regulatory_annotation_Apr27_ce10.bed) contains the summary annotation in an IGV-friendly .bed-file (track "Regulatory elements" on the screen shot above)
- colour shows the regulatory type:
  - `coding_promoter`: red
  - `pseudogene_promoter`: orange
  - `non-coding_RNA`: black
  - `unknown_promoter`: yellow
  - `putative_enhancer`: green
  - `other_element`: blue
- `name` ("labels below peaks") shows relevant gene(s), depending on the regulatory type of the accessible site:
  - `promoter_locus_id_fwd / promoter_locus_id_rev` for bidirectional promoters (either coding or non-coding)
  - `promoter_locus_id_fwd` or `promoter_locus_id_rev` for unidirectional promoters (either coding or non-coding)
  - `(associated_locus_id)` for all other sites
- `annot_%strand` strand-specific transcription pattern (shown when hovering above a specific peak with the mouse)

[`Fig2D2_regulatory_annotation_Apr27_ce11.bed`](Fig2D2_regulatory_annotation_Apr27_ce11.bed) contains the summary annotation as previously, except that the peak coordinates are given in `ce11` as detailed below

[`Fig2D2_regulatory_annotation.tsv`](Fig2D2_regulatory_annotation.tsv) contains an expanded set of metrics for all sites
- `chrom`, `start`, `end` coordinates of the accessible site in `ce10`
  - `chrom_ce11`, `start_ce11`, `end_ce11` coordinates of the accessible site in `ce11`, based on lifting over the position with peak accessibility
- `annot` summarises strand-specific patterns, already aggregated across stages, into a final regulatory element type:
  - `coding_promoter` on either strand => `coding_promoter`
  - `pseudogene_promoter` on either strand => `pseudogene_promoter`
  - `non-coding_RNA` on either strand => `non-coding_RNA`
  - `unknown_promoter` on either strand => `unknown_promoter`
  - `transcription_initiation` on either strand => `putative_enhancer`
  - all remaining sites => `other_element`
- `annot_%strand` strand/stage-specific transcription patterns across all stages by finding the 'highest' annotation across all stages, using the ranking below:
  1. `coding_promoter`
  2. `pseudogene_promoter`
  3. `non-coding_RNA`
  4. `unknown_promoter`
  5. `transcription_initiation`
  6. `no_transcription`

  (For example, a site called as `coding_promoter` in `wt_l2` will be called a `coding_promoter` in the aggregated annotation, regardless of what it was annotated as in other stages.)
- `promoter_gene_id_%strand`, `promoter_locus_id_%strand`, `promoter_gene_biotype_%strand` WormBase gene id, locus id, and biotype (`protein_coding`, `pseudogene`) of the associated promoter assignment
- `associated_gene_id`, `associated_locus_id` list of genes that overlap the site. The overlap regions include the outron region. If a site overlaps multiple genes, all overlaps are reported, separated by commas. This is used to tentatively assign enhancers to their putative target genes, based on the assumption that most enhancers target the closest gene (e.g. [Andersson et al. 2014](https://doi.org/10.1038/nature12787)).
- `tss_%strand` position of the transcription initiation mode within `-125:+125` of peak accessibility; for sites without reproducible short cap, defaults to 60bp downstream of peak accessibility (exprapolated based on the mean transcription initiation pattern at all sites)
- `lcap_%stage_%rep_%strand_ucount`, `_dcount` number reads with 5' ends within -250:-75 upstream (+75:+250 downstream) of the accessible site
- `lcap_%stage_%strand_baseMean`, `_log2FoldChange`, `_lfcSE`, `_stat`, `_pvalue`, `_padj` [DESeq2 output](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis) of a one-side test for an increase in downstream counts against upstream counts ("jump test")
- `lcap_%stage_%strand_passed` TRUE or FALSE based on whether the site has a statistically significant increase in long cap (padj<0.1 and log2FoldChange>1.5)
- `lcap_%stage_%strand_passed_incr` TRUE or FALSE based on whether the site passed the "incr" test
- `lcap_%stage_%strand_passed` TRUE or FALSE based on whether the site passed at least one of the jump/incr tests
