rule atac_cb3: # alignments & coverage tracks
    input:
        expand(pf('{sample}', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}', 'tg_se.bwa_se_cb3.rm_unmapped.rm_chrM.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac'), sample=config['atac_cb3']),

rule atac_cb3_yapc_spr:
    """
    Self-pseudoreplicate peak calls:
    pre: snakemake --use-conda --cores 12 -s sample_prp.snakefile atac_cb3_spr -n
    run: smc 10 atac_cb3_yapc_spr --keep-going --restart-times 5 -n
    """
    input:
        expand(pf('{sample}', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10.macs2_pe_lt200.yapc_spr_pe', '.tsv', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}', 'tg_se.bwa_se_cb3.rm_unmapped.rm_chrM.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.yapc_spr_se', '.tsv', 'atac'), sample=config['atac_cb3']),
