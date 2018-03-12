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
We aim to annotate accessible sites with regulatory type based on patterns of initiating and elongating transcription at the site.
1. Calculate a series of metrics/tests on strand/stage specific coverage patterns
2. Using these metrics, annotate the strand/stage specific pattern: `coding_promoter`, `pseudogene_promoter`, `non-coding_RNA`, `unknown_promoter`, `transcription_initiation` or `no_transcription`
3. Aggregate strand/stage specific patterns across stages, and then across the two strands, to determine the final type: `coding_promoter`, `pseudogene_promoter`, `non-coding_RNA`, `unknown_promoter`, `putative_enhancer` or `other_element`

## Individual tests
All tests are performed on strand-specific data from one stage.
- Tests for a **local increase in elongating transcription**:
  - estimate transcription elongation upstream, and downstream of a accessible site using two methods (5' counts, and median signal levels); use results with larger fold change in subsequent tests:
    - count 5' ends of long cap reads within `-300:-50` and `+50:+300` of `atac_mode`
    - measure median long cap levels within `-300:-50` and `+50:+300` of `atac_mode`
  - `jump` method:
    - test for an increase in downstream signal against upstream signal using DESeq2 (one-sided test)
    - significance thresholds set to `log2FoldChange > 1.5`, and `padj < 0.1`
  - `incr` method:
    - upstream counts `0` (zero) in both replicates
    - downstream counts `>=1` in both replicates individually
    - downstream counts `>=3` when pooled across replicates
- Identification of the closest first exon of the **putative (downstream) target gene**:
  - find closest first exon with the 5' end downstream of `atac_mode - 200`
  - uses first exons with `transcript_biotype` set to `protein_coding`, `pseudogene`, `tRNA`, `snoRNA`, `miRNA`, `snRNA` or `rRNA` ([`/WS260_ce10/WS260_ce10.transcripts.annot.gtf.gz`](/WS260_ce10/WS260_ce10.transcripts.annot.gtf.gz))
- Test for **continuous transcription elongation** between the accessible site, and the closest downstream gene:
  - calculate the maximum gap (zero-coverage region) in long cap signal, pooled across replicates, within `atac_mode + 300` to 5' end of the first exon of the putative target gene (`maxgap`).
  - test passes if `maxgap=0`
- Test for **proximity to 5' end of non-first exons**:
  - (*Motivation: long cap RNA-seq consists of a mixture of spliced and unspliced RNA such that transcription elongation at active gene bodies is higher on exons than on introns. Therefore, sites near 5' ends of active exons are very likely to test positive in the `jump`/`incr` tests. However, it is not possible to determine whether this increase is caused by transcription elongation that originates from the accessible site, or mature transcripts.*)
  - find closest up- or downstream non-first exon (derived from the same set as first exons)
  - calculate the distance between `atac_mode`, and the 5' end of the exon (*=location of long cap increase caused by mature mRNA*)
  - long cap jump is considered unreliable if the distance is under 300bp (up- or downstream)
- Short cap test for **reproducible transcription initiation**:
  - calculate maximum short cap signal (pooled across stages) within 100bp of `atac_mode`
  - sites with at non-singleton short cap stack (*minimum criteria to define a TIC in ([Chen et al. 2013](https://doi.org/10.1101/gr.153668.112))*) are considered to have reproducible transcription initiation

## Stage/strand-specific annotation
- Sites at least 300bp away from exons:
  - elongation (`jump` or `incr`) and continuous transcription (`maxgap=0`) to a `protein_coding` or `pseudogene` gene => `coding_promoter` or `pseudogene_promoter`
  - elongation (`jump` or `incr`) but transcription not continuous (`maxgap>0`) => `unknown_promoter`
  - has transcription initiation => `transcription_initiation`
  - remaining sites => `no_transcription`
- 5' end of first exon within 300bp up- or downstream of `atac_mode`, and gene is `tRNA`, `snoRNA`, `miRNA`, `snRNA` or `rRNA` (=other than `protein_coding`) => `non-coding_RNA`
- 5' end of first exon within 300bp downstream of `atac_mode`, and gene is `protein_coding` or `pseudogene`:
  - elongation via `incr` => `coding_promoter` or `pseudogene_promoter`
  - elongation via `jump` and reproducible transcription initiation => `coding_promoter` or `pseudogene_promoter`
  - has transcription initiation => `transcription_initiation`
  - remaining sites => `no_transcription`
- 5' end of first exon is within 300bp upstream of `atac_mode`:
  - elongation via `incr` or `jump`; has a transcription initiation mode within 50bp of annotated TSS => `coding_promoter` or `pseudogene_promoter`
  - has transcription initiation => `transcription_initiation`
  - remaining sites => `no_transcription`
- Remaining sites (sites within 300bp of an other exon):
  - (*No attempt at any type of promoter annotation as elongation tests are unreliable.*)
  - has transcription initiation => `transcription_initiation`
  - remaining sites => `no_transcription`

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
