# S3 Genomic regions

[`S3_genomic_regions.bed`](S3_genomic_regions.bed) has the entire `ce10` genome partitioned into `outronic` (red), `exonic` (orange), `intronic` (yellow), `intergenic` (light blue) or `mixed` (black) based on the regulatory annotation, and WS260_ce10 genomic annotations.

[`S3_genomic_regions.tsv`](S3_genomic_regions.tsv) has the same contents as a tab-separated table.

Outline of the partitioning steps:
1. `exonic`: union of all annotated exons sharing a `gene_id`
2. `intronic`: sequence between `exonic` regions
3. `gene_end`: sequence flanking 100bp at the 3' end of the last exon of a gene
4. `outronic`: sequence spanning from a `coding_promoter` to the closest downstream exon
5. `intergenic`: all remaining regions
6. `mixed`: all regions that have previously been assigned more than one type, e.g. `exonic` and `intronic`
