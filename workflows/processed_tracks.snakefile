df_atac_ = None
def df_atac():
    global df_atac_

    def processed_id(bid, lid):
        if lid != lid or lid == '':
            return '_%(bid)s' % locals()
        return '%(lid)s_%(bid)s' % locals()
    def is_se(bid): return os.path.isfile('samples/%(bid)s.r1.fq.gz' % locals())
    def is_pe(bid): return os.path.isfile('samples/%(bid)s.r2.fq.gz' % locals())

    if df_atac_ is None:
        fp_ = 'processed_tracks/metadata/Worm Regulatory Mapping Data Sets - Libraries (ATAC-, DNase-, MNase-seq).tsv'
        df_ = pd.read_csv(fp_, sep='\t').query('(Enzyme == "Tn5") & (Genome == "ce10")')[['Bioinformatics ID(s)', 'Library series ID']].rename(columns={'Bioinformatics ID(s)': 'bid', 'Library series ID': 'lid'}).reset_index(drop=True)
        df_['pid'] = [ *map(processed_id, df_['bid'], df_['lid']) ]
        df_['is_pe'] = [ *map(is_pe, df_['bid']) ]
        df_['is_se'] = [ *map(is_se, df_['bid']) ]
        #df_ = df_.query('bid == "HS298_JA26_N2_atac_S1" | bid == "HS491_metset_ATAC_rep1" | bid == "HS491_tm548_ATAC_rep1"')
        df_atac_ = df_

    return df_atac_

def atac_ce10_spmr_se_input_(wildcards):
    pid_ = wildcards.pid
    df_ = df_atac().query('pid == @pid_')
    assert(len(df_) == 1)
    return pf(df_.iloc[0]['bid'], 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac')

def atac_ce10_spmr_pe_input_(wildcards):
    pid_ = wildcards.pid
    df_ = df_atac().query('pid == @pid_')
    assert(len(df_) == 1)
    return pf(df_.iloc[0]['bid'], 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac')

rule atac_ce10_spmr_se:
    input:
        atac_ce10_spmr_se_input_
    output:
        pf('{pid}', 'atac_ce10_spmr_se', '.bw', 'processed_tracks')
    shell:
        'scripts/bigWiggleTools.ipy write {output} scale 0.1 bin 10 {input}'

rule atac_ce10_spmr_pe:
    input:
        atac_ce10_spmr_pe_input_
    output:
        pf('{pid}', 'atac_ce10_spmr_pe', '.bw', 'processed_tracks')
    shell:
        'scripts/bigWiggleTools.ipy write {output} scale 0.1 bin 10 {input}'

rule atac_processed_stats:
    input:
        expand(pf('{bid}', '{step}', '.txt', 'atac'), bid=df_atac().query('is_se')['bid'].tolist(),
            step=[
                'c_r1', # Total reads
                'tg_se.bwa_se.c', # After adapter trimming
                'tg_se.bwa_se.rm_unmapped.c', # Mapped reads
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c', # Non-mitochondrial
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c', # Non-blacklisted
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c', # Q<10-aligned
                ]),
    output:
        'atac/atac_stats_raw_counts.tsv', # raw read counts at each step
        'atac/atac_stats_frac_passed.tsv',
        'processed_tracks/atac_ce10_stats.tsv', # percentages that passed each step

    run:
        """
        df = make_dnase_mnase_samples(expanded=True)[['Bioinformatics ID(s)', 'Library series ID', 'geo_id', 'Digestion conditions (homogenized)']].set_index('Bioinformatics ID(s)') #pd.DataFrame()
        df.index.name = 'dataset_id'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            df.ix[bid, step] = read_int(pf(bid, step, suffix, prefix))
        df.to_csv(output[0], sep='\t')
        """
        df = df_atac().query('is_se')[['bid', 'lid']].set_index('bid')
        df.index.name = 'dataset_id'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            df.ix[bid, step] = read_int(pf(bid, step, suffix, prefix))
        df.to_csv(output[0], sep='\t')

        def pct_(a, b): return list(map(lambda a_i, b_i: '%.01f%%' % (100.0 * a_i / b_i,), a, b))
        def loss_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * (a_i - b_i) / a_i,), df[col_a], df[col_b]))
        def keep_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * b_i / a_i,), df[col_a], df[col_b]))

        df_ = pd.DataFrame()
        df_['raw_reads'] = df['c_r1'].astype(int).map(yp.f_uk)
        df_['mapped'] = keep_pct_(
                'tg_se.bwa_se.c',
                'tg_se.bwa_se.rm_unmapped.c',
        )
        df_['not_mitochondrial'] = keep_pct_(
                'tg_se.bwa_se.rm_unmapped.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c',
        )
        df_['not_blacklisted'] = keep_pct_(
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c',
        )
        df_['mapq10'] = keep_pct_(
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c',
        )
        df_['useful_reads'] = df['tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c'].astype(int).map(yp.f_uk)
        df_.sort_index().to_csv(output[1], sep='\t')
        
        df_.insert(loc=0, column='library_series_id', value=df['lid'])
        d_ = dict(zip(config['atac814'].values(), config['atac814'].keys()))
        df_.insert(loc=1, column='geo_id', value=[ *map(lambda bid: d_.get(bid, ''), df_.index) ])

        df_.sort_index().to_csv(output[2], sep='\t')

rule processed_tracks:
    # smj 100 processed_tracks --keep-going --restart-times 3 -n --quiet
    # smc 4 processed_tracks -n
    input:
        'processed_tracks/atac_ce10_stats.tsv',
        #expand(pf('{pid}', 'atac_ce10_spmr_se', '.bw', 'processed_tracks'), pid=[* df_atac().query('is_se')['pid'] ]),
        #expand(pf('{pid}', 'atac_ce10_spmr_pe', '.bw', 'processed_tracks'), pid=[* df_atac().query('is_pe')['pid'] ]),
        #'processed_tracks/scap_ce10_stats.tsv',
        #expand(pf('{pid}', 'scap_ce10_init_fwd', '.bw', 'processed_tracks'), pid=[* df_scap()['pid'] ]),
        #expand(pf('{pid}', 'scap_ce10_init_rev', '.bw', 'processed_tracks'), pid=[* df_scap()['pid'] ]),
        #'processed_tracks/lcap_ce10_stats.tsv',
        #expand(pf('_{bid}', 'lcap_ce10_linear_fwd', '.bw', 'processed_tracks'), bid=config['lcap_raw']),
        #expand(pf('_{bid}', 'lcap_ce10_linear_rev', '.bw', 'processed_tracks'), bid=config['lcap_raw']),
        #'processed_tracks/dmnase_ce10_stats.tsv',
        #expand(pf('{label}', 'dmnase_ce10_spmr', '.bw', 'processed_tracks'), label=d_dnase_mnase_label_bid.keys()),
