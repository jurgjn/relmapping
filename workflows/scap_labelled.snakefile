l_lid = [config['scap'][bid]['library_series_id'] for bid in config['scap']]

rule scap_lsid_bw:
    input:
        lambda wildcards: pf(scap_bid_from_lsid(wildcards.lsid), wildcards.step, '.bw', 'scap'),
    output:
        pf('{lsid}', '{step}.lsid', '.bw', 'scap'),
    shell:
        'cp {input} {output}'

rule txn_init_fwd: #https://bitbucket.org/snakemake/snakemake/issues/397/unable-to-set-utime-on-symlink-your-python
    input:
        lambda wildcards: pf(scap_bid_from_lsid(wildcards.lid), 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{lid}', 'txn_init_fwd', '.bw', 'scap_labelled'),
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule txn_init_rev: #https://bitbucket.org/snakemake/snakemake/issues/397/unable-to-set-utime-on-symlink-your-python
    input:
        lambda wildcards: pf(scap_bid_from_lsid(wildcards.lid), 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{lid}', 'txn_init_rev', '.bw', 'scap_labelled'),
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap_labelled:
    input:
        expand(pf('{lid}', 'txn_init_fwd', '.bw', 'scap_labelled'), lid=l_lid),
        expand(pf('{lid}', 'txn_init_rev', '.bw', 'scap_labelled'), lid=l_lid),
