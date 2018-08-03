# Batch processing (for jadb)
1. Place a three column text file with the datasets to-be-processed to `batch/datasets.txt`, e.g.:
```
Dataset              Assay    Organism
HS086_Caro-A         atac     ce
HS086_Caro-E         atac     ce
HS476_MS51           atac     ce
HS476_MS52           atac     ce
HS645_ATAC_TG_CB2    atac     cb
```
* `Dataset` is used to locate the sequenced reads, e.g. for the first dataset, the gzipped fastq file is assumed to be at `samples/HS086_Caro-A.r1.fq.gz` (add `.r2.fq.gz` for paired-end data). The pipeline will automatically detect whether a dataset is single-end or paired-end based on the presence/absence of the second read file.
* `Assay` should be set to `atac` (other assays -- `scap`, `lcap`, ... will be added after the initial pipeline is working).
* `Organism` should be set to `ce` for _C. elegans_, or `cb` for _C. briggsae_.
2. Run the following command (in `~/relmapping`):
```
snakemake --cluster sbatch --jobs 50 batch
```
3. Once finished, a table of processed files can be found at  `batch/processed_files.txt`:
```
Dataset           Processing                   Scale Genome Resolution Filetype File
HS086_Caro-A      atac_ce10_spmr_se            SPMR  ce10   10bp       bw       HS086_Caro-A.atac_ce10_spmr_se.bw
HS086_Caro-A      atac_ce10_spmr_pe            SPMR  ce10   10bp       bw       HS086_Caro-A.atac_ce10_spmr_pe.bw
HS086_Caro-E      atac_ce10_spmr_se            SPMR  ce10   10bp       bw       HS086_Caro-E.atac_ce10_spmr_se.bw
HS086_Caro-E      atac_ce10_spmr_pe            SPMR  ce10   10bp       bw       HS086_Caro-E.atac_ce10_spmr_pe.bw
HS476_MS51        atac_ce10_spmr_se            SPMR  ce10   10bp       bw       HS476_MS51.atac_ce10_spmr_se.bw
HS476_MS52        atac_ce10_spmr_se            SPMR  ce10   10bp       bw       HS476_MS52.atac_ce10_spmr_se.bw
HS645_ATAC_TG_CB2 atac_CB4_spmr_se_preliminary SPMR  CB4    10bp       bw       HS645_ATAC_TG_CB2.atac_CB4_spmr_se_preliminary.bw
HS645_ATAC_TG_CB2 atac_CB4_spmr_pe_preliminary SPMR  CB4    10bp       bw       HS645_ATAC_TG_CB2.atac_CB4_spmr_pe_preliminary.bw
```
