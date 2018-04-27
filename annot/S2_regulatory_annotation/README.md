# S2 Regulatory annotation of accessible sites
![Fig_ztf-11.png](Fig_ztf-11.png)

[`S2_regulatory_annotation.bed`](S2_regulatory_annotation.bed) contains the summary annotation in an IGV-friendly .bed-file (track "Regulatory elements" on the screen shot above)
- colour shows the regulatory type:
  - `coding_promoter`: red
  - `pseudogene_promoter`: orange
  - `non-coding_RNA`: black
  - `unknown_promoter`: yellow
  - `putative_enhancer`: green
  - `other_element`: blue
- `name` ("labels below peaks") shows relevant gene(s), depending on the regulatory type of the accessible site:
  - `promoter_gene_id_fwd / promoter_gene_id_rev` for bidirectional promoters (either coding or non-coding)
  - `promoter_gene_id_fwd` or `promoter_gene_id_rev` for unidirectional promoters (either coding or non-coding)
  - `(associated_gene_id)` for all other sites
- `strand` ("arrows") shows the "transcriptional strandedness" around the site (see below)
- `annot_fwd` and `annot_rev` show the strand-specific transcription pattern (shown when hovering above a specific peak with the mouse)

[`S2_regulatory_annotation.tsv`](S2_regulatory_annotation.tsv) contains a more extensive set of metrics for all sites, described further below

## Outline of the methodology
We aim to annotate accessible sites with regulatory type based on patterns of initiating and elongating transcription at the site, and the positioning of the site relative to transcript annotations.
1. Calculate a series of metrics/tests on strand/stage specific coverage patterns
2. Using these metrics, annotate the strand/stage specific pattern: `coding_promoter`, `pseudogene_promoter`, `non-coding_RNA`, `unknown_promoter`, `transcription_initiation` or `no_transcription`
3. Aggregate strand/stage specific patterns across stages, and then across the two strands, to determine the final type: `coding_promoter`, `pseudogene_promoter`, `non-coding_RNA`, `unknown_promoter`, `putative_enhancer` or `other_element`

## Individual tests
All tests/metrics are performed/calculated separately for every site and on both strands. Tests that use short/long cap data are performed on both strands, and all stages.
- Tests for a **local increase in elongation**:
  - estimate transcription elongation upstream, and downstream of a site by counting 5' ends of long cap reads within `-250:-75` and `+75:+250` of peak accessibility
  -  `jump` method: test for an increase in downstream signal against upstream signal using DESeq2 (one-sided test, significance thresholds set to `log2FoldChange > 1.5`, and `padj < 0.1`)
  - `incr` method:
    - upstream counts `0` (zero) in both replicates
    - downstream counts `>=1` in both replicates individually
    - downstream counts `>=3` when pooled across replicates
- Test for **reproducible transcription initiation**: check for a non-singleton short cap stack within 125bp of peak accessibility (using short cap signal pooled across all replicates/stages, and filtered for reproducibility)
- Identification of the **transcription initiation mode**:
  - sites with reproducible short cap are assigned the position with maximum short cap signal
  - sites without reproducible short cap are assigned an extrapolated, "best-guess" position of 60bp downstream of peak accessibility (motivated by the global genome-wide short cap distribution at accessible sites)
- Identification of **proximal exons**:
  - **closest first exon**: closest first exon of a `coding_promoter` or `pseudogene` within 250bp upstream to anywhere downstream of peak accessibility
  - **closest other exon**: closest non-first exon of a `coding_promoter` or `pseudogene`

## Stage/strand-specific annotation
- `coding_promoter`, `pseudogene_promoter`:
  - **joint criteria**, regardless of confidence levels:
    - site is either "close to a first exon" or "away from any other exon"
      - "close to a first exon" = 5' end of the closest first exon is within 250bp of peak accessibility
      - "away from any other exon" = 5' end of the closest other exon is further than 250bp away from peak accessibility
    - region between peak accessibility and the 5' end of the closest exon 1 does not contain the 5' end of a non-first exon
    - transcription initiation mode is located "properly relative to the TSS"
      - "properly relative to the TSS" = transcription initiation mode is either upstream of the annotated transcript start site, or, in the presence of a UTR, up to 250bp downstream of the annotated transcript start site, within the UTR
    - "if distal, site has continuous transcription": for sites further than 250bp upstream of the closest first exon, check for continuous long cap coverage from 250bp downstream of peak accessibility to the 5' end of the closest first exon
  - **non-low confidence** `coding_promoter` or `pseudogene_promoter` if any of the following:
    - has transcription initiation, and passes the jump test
    - has transcription initiation, and passes the incr test
    - passes the jump and incr tests
  - **low-confidence** `coding_promoter` or `pseudogene_promoter` if:
    - target gene does not get any promoters using the other methods
    - site is either intergenic, or within the first 250bp of a gene
    - site has a long cap log2FoldChange > 1
    - site has the highest log2FoldChange out of all sites that fulfil all other criteria
- `non-coding_RNA`: closest downstream first exon of a `tRNA`, `snoRNA`, `miRNA`, `snRNA` or `rRNA` is within 250bp of peak accessibility
- `unknown_promoter`:
  - site is at least 250bp away from the 5' end of the any (first or other) exon
  - site is intergenic: peak accessibility is not within a (sense-strand) gene body
  - has transcription initiation
  - passes the jump test
- `transcription_initiation`: has transcription initiation
- `no_transcription`: was not assigned to any previous category

## Aggregating the annotation across stages and strands + additional metrics
- `annot_%strand` strand/stage-specific transcription patterns across all stages by finding the 'highest' annotation across all stages, using the ranking below:

  1. `coding_promoter`
  2. `pseudogene_promoter`
  3. `non-coding_RNA`
  4. `unknown_promoter`
  5. `transcription_initiation`
  6. `no_transcription`

  (For example, a site called as `coding_promoter` in `wt_l2` will be called a `coding_promoter` in the aggregated annotation, regardless of what it was annotated as in other stages.)

- `annot` summarises strand-specific patterns, already been aggregated across stages, into a final regulatory element type:
  - `coding_promoter` on either strand => `coding_promoter`
  - `pseudogene_promoter` on either strand => `pseudogene_promoter`
  - `non-coding_RNA` on either strand => `non-coding_RNA`
  - `unknown_promoter` on either strand => `unknown_promoter`
  - `transcription_initiation` on either strand => `putative_enhancer`
  - all remaining sites => `other_element`

- `strand` ("arrows in the final .bed-file"):
  - if annotated as `coding_promoter`, `pseudogene_promoter` or `non-coding_RNA` on at least one strand, annotate as bidirectional, forward or reverse based on those two annotations. For example, `coding_promoter` on fwd and `non-coding_RNA` on rev is assigned as bidirectional.
  - otherwise (=not annotated as `coding_promoter`, `pseudogene_promoter` or `non-coding_RNA` on either strand), the strandedness is defined using long cap, adapting the approach from ([Chen et al. 2013](https://doi.org/10.1101/gr.153668.112)). Specifically, compare long cap read counts at the site, pooled across stages, and assign the strand with higher counts. Sites with equal long cap counts are annotated as unstranded (`.`).

- `promoter_gene_id_%strand` for `coding_promoter`, `pseudogene_promoter`, and `non-coding_RNA`, set to the WormBase `gene_id` of downstream gene linked by transcription. Set to `.` for other types of elements

- `associated_gene_id` WormBase `gene_id` of the closest downstream or overlapping gene that is located on the strand with higher levels of transcriptional elongation at the site (adapted from [Chen et al. 2013](https://doi.org/10.1101/gr.153668.112))). The aim of this is to identify the likeliest target gene for an enhancer, based on the assumption that enhancers are most likely to interact with the closest gene (e.g. [Andersson et al. 2014](https://doi.org/10.1038/nature12787)). For a site that overlaps multiple genes, all overlaps are reported, separated by commas.
