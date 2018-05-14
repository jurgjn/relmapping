
rule atac_cfp1:
    input:
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.fsizes', '.txt', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt150', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_keepdup_all', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        #expand(pf('{sample}', 'tg_pe.bwa_pe_ce11.rm_unmapped_pe.rm_chrM', '.bam', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_pe.bwa_pe_ce11.rm_unmapped_pe.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),
        expand(pf('{sample}', 'tg_pe.bwa_pe_ce11.rm_unmapped_pe.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac_cfp1'), sample=config['atac_cfp1_805']),

rule atac_cfp1_qc:
    input:
        expand(pf('{sample}', '{step}', '.txt', 'atac_cfp1'), sample=config['atac_cfp1_805'],
            step=[
                'c_r1', # Total reads
                'tg_pe.bwa_pe.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr_rmdup',
                ]),
    output:
        'atac_cfp1/atac_qc_counts.tsv',
        'atac_cfp1/atac_qc_passed.tsv',
    run:
        df = pd.DataFrame()
        df.index.name = 'sample'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            df.ix[bid, step] = read_int(pf(bid, step, suffix, prefix))
        df.to_csv(output[0], sep='\t')

        def pct_(a, b): return list(map(lambda a_i, b_i: '%.01f%%' % (100.0 * a_i / b_i,), a, b))
        def loss_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * (a_i - b_i) / a_i,), df[col_a], df[col_b]))
        def keep_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * b_i / a_i,), df[col_a], df[col_b]))

        df_ = pd.DataFrame()
        df_['raw_reads'] = df['c_r1'].astype(int)
        df_['mapped'] = keep_pct_(
                'tg_pe.bwa_pe.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.c_r1',
        )
        df_['not mitochondrial'] = keep_pct_(
                'tg_pe.bwa_pe.rm_unmapped_pe.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
        )
        df_['not blacklisted'] = keep_pct_(
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.c_r1',
        )
        df_['mapq10'] = keep_pct_(
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_r1',
        )
        df_['useful_reads'] = df['tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_r1'].astype(int)
        
        df_['complexity'] = (df['tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr_rmdup'].astype(int) \
                           / df['tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr'].astype(int)).map('{:,.02f}'.format)

        df_.sort_index().to_csv(output[1], sep='\t')

"""
Merging across different samples...

# N2 live:	HS086-Caro-A  + HS086-Caro-E
scripts/bigWiggleTools.ipy write_bg atac_cfp1/merged_tracks/atac_N2_live.bw mean \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS086_Caro-A.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS086_Caro-E.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw

# set2 live:  HS086-Caro-I  +  HS086-Caro-K
scripts/bigWiggleTools.ipy write_bg atac_cfp1/merged_tracks/atac_set2_live.bw mean \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS097_Caro-I.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS097_Caro-K.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw

# N2 frozen:    HS322-JA33 + HS298_JA26
scripts/bigWiggleTools.ipy write_bg atac_cfp1/merged_tracks/atac_N2_frozen.bw mean \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS322_crb55_EMB-JA33.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS298_JA26_N2_atac_S1.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw

# cfp1 frozen:  HS333-RC22  + HS333-JA30
scripts/bigWiggleTools.ipy write_bg atac_cfp1/merged_tracks/atac_cfp1_frozen.bw mean \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS333_crb55_RC22.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw \
atac_cfp1/tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300/HS333_crb55_JA30.tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300_treat_pileup.bw
"""