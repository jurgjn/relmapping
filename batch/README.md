# Batch processing (for jadb)
1. Write a table with the datasets to-be-processed to `batch/samples.txt`, e.g.:
```
Dataset              Assay    Organism
HS086_Caro-A         atac     ce
HS086_Caro-E         atac     ce
```
* `Dataset` is used to locate the sequenced reads, e.g. for the first dataset, the gzipped fastq file is assumed to be at `samples/HS086_Caro-A.r1.fq.gz` (add `.r2.fq.gz` for paired-end data). The pipeline will automatically detect whether a dataset is single-end or paired-end based on the presence/absence of the second read file.
* `Assay`  should be set to `atac` (other data types will be added after the initial pipeline is working).
* `Organism` should be set to `ce` (other organisms (briggsae) will be added after the initial pipeline is working).
2. Run the following command:
```
smj 100 jadb -n
```
3. Once finished, a table of processed files can be found at  `batch/processed_files.txt`:
```
Dataset       Processing         Scale  Genome  Resolution  Filetype  File 
HS086_Caro-A  atac_ce10_spmr_se  SPMR   ce10    10bp        bw        HS086_Caro-A.atac_ce10_spmr_se.bw
HS086_Caro-E  atac_ce10_spmr_se  SPMR   ce10    10bp        bw        HS086_Caro-A.atac_ce10_spmr_se.bw
```
