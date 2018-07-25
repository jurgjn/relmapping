
rule processed_tracks:
    # smj 100 processed_tracks --keep-going --restart-times 3 -n
    # smc 4 processed_tracks -n
    input:
        #'processed_tracks/atac_ce10_stats.tsv',
        #expand(pf('{pid}', 'atac_ce10_spmr_se', '.bw', 'processed_tracks'), pid=[* df_atac().query('is_se')['pid'] ]),
        #expand(pf('{pid}', 'atac_ce10_spmr_pe', '.bw', 'processed_tracks'), pid=[* df_atac().query('is_pe')['pid'] ]),
        #'processed_tracks/scap_ce10_stats.tsv',
        #expand(pf('{pid}', 'scap_ce10_init_fwd', '.bw', 'processed_tracks'), pid=[* df_scap()['pid'] ]),
        #expand(pf('{pid}', 'scap_ce10_init_rev', '.bw', 'processed_tracks'), pid=[* df_scap()['pid'] ]),
        #'processed_tracks/lcap_ce10_stats.tsv',
        #expand(pf('_{bid}', 'lcap_ce10_linear_fwd', '.bw', 'processed_tracks'), bid=config['lcap_raw']),
        #expand(pf('_{bid}', 'lcap_ce10_linear_rev', '.bw', 'processed_tracks'), bid=config['lcap_raw']),
        'processed_tracks/dmnase_ce10_stats.tsv',
        #expand(pf('{label}', 'dmnase_ce10_spmr', '.bw', 'processed_tracks'), label=d_dnase_mnase_label_bid.keys()),
