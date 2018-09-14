# Briggsae annotation (do not run -- draft)

1. Specify conditions (stages), and files of raw reads for two (biological) replicates in the `annot_cb` section of [workflows/config.yaml](workflows/config.yaml), e.g.:
```
annot_cb:
  atac_samples:
    wt_ya:
      rep1:
        - samples_raw/ATAC-JA60-YA_S1_L001_R1_001.fastq.gz
        - samples_raw/ATAC-JA60-YA_S1_L002_R1_001.fastq.gz
      rep2:
        - samples_raw/ATAC-TG-CB2_S5_L001_R1_001.fastq.gz
        - samples_raw/ATAC-TG-CB2_S5_L002_R1_001.fastq.gz
    #glp1_ya: # Uncomment & add files as above, once data becomes available...
  lcap_samples:
    wt_ya:
      # Mock setup where lanes are treated as replicate
      # Once second replicate becomes available:
      # 1) add both lanes to rep1
      # 2) add new data as rep2
      rep1_read1:
        - samples_raw/LCRNA-cb-YA_S14_L001_R1_001.fastq.gz
      rep1_read2:
        - samples_raw/LCRNA-cb-YA_S14_L001_R2_001.fastq.gz
      rep2_read1:
        - samples_raw/LCRNA-cb-YA_S14_L002_R1_001.fastq.gz
      rep2_read2:
        - samples_raw/LCRNA-cb-YA_S14_L002_R2_001.fastq.gz
    #glp1_ya: # Uncomment & add files as above, once data becomes available...
```

2. Run the target (first to test with -n; then remove -n to start the run).
```
$ smc 20 annot_cb -n
```
