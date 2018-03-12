l_HS491 = ['HS491_metset_ATAC_rep1', 'HS491_metset_ATAC_rep2', 'HS491_tm548_ATAC_rep1', 'HS491_tm548_ATAC_rep2']
rule atac_HS491:
    input:
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac_HS491'), bid=l_HS491),

rule atac_HS491_qc:
    input:
        expand(pf('{bid}', '{step}', '.txt', 'atac_HS491'), bid=l_HS491,
            step=[
                'c_r1', # Total reads
                'tg_se.bwa_se.c',
                'tg_se.bwa_se.rm_unmapped.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c',
                ]),
    output:
        'atac_HS491/atac_qc_counts.tsv',
        'atac_HS491/atac_qc_passed.tsv',
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
                'tg_se.bwa_se.c',
                'tg_se.bwa_se.rm_unmapped.c',
        )
        df_['not mitochondrial'] = keep_pct_(
                'tg_se.bwa_se.rm_unmapped.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c',
        )
        df_['not blacklisted'] = keep_pct_(
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c',
        )
        df_['mapq10'] = keep_pct_(
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c',
        )
        df_['useful_reads'] = df['tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c'].astype(int)
        df_.sort_index().to_csv(output[1], sep='\t')
