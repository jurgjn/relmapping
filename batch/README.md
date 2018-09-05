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

## Description for default processing

Brackets refer to corresponding processing rules in [workflows/](/workflows/).

### ATAC-seq (atac_ce10_spmr_se)

1. Reads were trimmed using trim_galore (`tg_se`), and aligned using bwa in single-end mode (`bwa_se`).
2. Low-quality (q<10, `rm_unmapped/rm_q10`), mitochondrial (`rm_chrM`), and [modENCODE blacklisted](https://www.encodeproject.org/comparative/regulation/#Wormset5) (`rm_blacklist`) reads were discarded.
3. Normalised coverage was calculated using MACS2 by extending/shifting the tags, and keeping duplicates (`macs2_se_extsize150_shiftm75_keepdup_all`).
    - This leads to high signal [where the cut sites are](https://github.com/taoliu/MACS/issues/145) by extending and shifting the tags, thereby [effectively smoothing the cut site coverage signal](https://groups.google.com/forum/#!topic/macs-announcement/4OCE59gkpKYs)
4. Signal was smoothed by binning at 10bp resolution (`bin10`).

Paired-end ATAC-seq data is also processed as follows (`atac_ce10_spmr_pe`).

1. Reads were trimmed using trim_galore (`tg_pe`), and aligned using bwa in single-end mode (`bwa_pe`).
2. Low-quality (q<10, `rm_unmapped_pe/rm_q10`), mitochondrial (`rm_chrM`), and [modENCODE blacklisted](https://www.encodeproject.org/comparative/regulation/#Wormset5) (`rm_blacklist`) reads were discarded.
3. Normalised coverage was calculated using MACS2 in paired-end mode, size-selecting for <200bp fragments (`macs2_pe_lt200`).
4. Signal was smoothed by binning at 10bp resolution (`bin10`).
