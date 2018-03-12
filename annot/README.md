# Figures
## Fig 1 Overview of the project

(Suppl) Biological reproducibility of long cap data. To broadly assess biological reproducibility, we calculated gene-level read counts, and performed hierarchial clustering of top 10,000 most highly expressed genes. We used logged pseudocounts as input, and correlation distance as the metric. Biological replicates prepared from the same stage show high correlation,  and are identified as closest to each other. ([Figure](FigA_clustering/reproducibility_lcap.pdf))

## Fig 2 Annotation of regulatory elements
(C) Number of sites identified as unidirectional forward, unidirectional reverse, and bidirectional promoters of protein-coding genes
([Figure](FigA_mapping/coding_promoter_by_type.pdf), [Source data](FigA_mapping/coding_promoter_by_type.tsv))

(D) Distribution of promoters per gene
([Figure](FigA_mapping/promoters_per_gene.pdf), [Source data](FigA_mapping/promoters_per_gene.tsv))

(E) Distribution of enhancers per gene
([Figure](FigA_mapping/enhancers_per_gene.pdf), [Source data](FigA_mapping/enhancers_per_gene.tsv))

(F) Number of enhancers assigned to a gene, split by the number of promoters assigned to the gene
([Figure](FigA_mapping/npromoter_vs_nenhancer.pdf))

(Figure 2 -- Supplement 2) [Examples of unknown promoters](Fig2S2/)
1. Probable protein-coding promoters not captured due to weak long cap signal
2. Antisense promoters at protein-coding genes

(Figure 2 -- Supplement 3a) Comparison of our ATAC-seq, (Daugherty et al., 2017), and (Ho et al., 2017) to modENCODE/modERN TFBS clusters [Figure](Fig2S3_overlaps/Fig2S3a_compareto_TFBS.pdf)
- modENCODE/modERN clusters were defined by extending all peaks by 200bp on each side of the summit, and clustering using a single-linkage approach

(Figure 2 -- Supplement 3b) Comparison of our ATAC-seq to existing sets of accessible sites.
- (Daugherty et al., 2017) ATAC-seq peaks from Supplemental Table S3 ([Figure](Fig2S3_overlaps/Fig2S3b_Daugherty2017.pdf)
- (Ho et al., 2017) DNase-seq from GEO accession GSE97425 ([Figure](Fig2S3_overlaps/Fig2S3b_Ho2017.pdf)
- Subplots from left-to-right:
    1. Venn diagram showing the overlap between accessible sites from this study, and previously defined sets of accessible sites
    2. modENCODE/modERN filtered peak call coverage at four sets of regions, as defined by the overlap
    3. Transcription initiation at four sets of regions, as defined by the overlap
    4. Exon coverage at four sets of regions, as defined by the overlap

(Figure 2 -- Supplement 4) Screen shots of stage-specific genome-wide accessibility profiles from Daugherty et al., 2017 (ATAC-seq), Ho et al., 2017 (DNase-seq), and this study (ATAC-seq)
- Regions taken from (Daugherty et al., 2017), Fig 3b, 3d, 3f
- ([Daugherty et al., Fig 3b chrI:994,900-1,011,200 C54G6.3](Fig2S4_atac_screen_shots/3b_C54G6.2.png))
- ([Daugherty et al., Fig 3d chrX:13,007,000-13,015,000 nhr-25](Fig2S4_atac_screen_shots/3d_nhr-25.png))
- ([Daugherty et al., Fig 3f chrX:2,111,000-2,120,000 swip-10](Fig2S4_atac_screen_shots/3f_swip-10.png))

## Fig 3 Properties of promoters/enhancers
(A) Heatmap of histone marks at `coding_promoter`, and `putative_enhancer`, sorted by H3K4me3 levels
([Figure](FigB_chromatin/hmods_prom_enh_sort_H3K4me3.png))

(B) Inr, TATA, CpG profiles on `coding_promoter`, and `putative_enhancer`, binned into three based on H3K4me3 levels as used in (A)
([Figure](FigB_motifs/Inr_TATA_CpG_sortbyH3K4me3.png))

(Figure 3 -- Supplement 1) Same setup as main figure, except sorted by CV 

(A) Heatmap of histone marks ([Figure](FigB_chromatin/hmods_prom_enh_sort_CV.png))

(B) Inr, TATA, CpG profiles ([Figure](FigB_motifs/Inr_TATA_CpG_sortbyCV.png))

(Figure 3 -- Supplement 2) Comparison of the transcription-based annotation (rows) to chromatin states (columns), as determined by the chromatin state that overlaps the mode of the accessible site. Colours indicate the fraction of annotations that overlap a particular state, out of all annotations of that type.
- [Ho2014_EE](FigB_chromatin_states/annot_vs_Ho2014_EE.png)
- [Ho2014_L3](FigB_chromatin_states/annot_vs_Ho2014_L3.png)
- [Evans2016_EE](FigB_chromatin_states/annot_vs_Evans2016_EE.png)
- [Evans2016_L3](FigB_chromatin_states/annot_vs_Evans2016_L3.png)

(Figure 3 -- Supplement 2) [Screen shots of examples](FigB_screen_shots/) -- broad H3K4me3 domain with multiple elements; promoter in enhancer state; promoter in regulated chromatin...
