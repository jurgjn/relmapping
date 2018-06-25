rule atac_cb3:
    input:
        expand(pf('{sample}', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}', 'tg_se.bwa_se_cb3.rm_unmapped.rm_chrM.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'), sample=config['atac_cb3']),
