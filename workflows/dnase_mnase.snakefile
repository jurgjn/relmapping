def make_dnase_mnase_samples(expanded=False):
    fp_inp = 'processed_tracks/metadata/Worm Regulatory Mapping Data Sets - Libraries (ATAC-, DNase-, MNase-seq).tsv'
    fp_out = 'dnase_mnase/dnase_mnase_samples.tsv'
    df_ = pd.read_csv(fp_inp, sep='\t').query('(Strain == "N2") & (Enzyme == "dnase" | Enzyme == "mnase")').sort_values(['Enzyme', 'Stage'])
    df_ = df_.loc[df_['Bioinformatics ID(s)'] == df_['Bioinformatics ID(s)']] # All sequenced libraries
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA0-N2-dnase-test')] # Early test data, partial digestion course
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA11-N2-dnase-S2.1xtraHS082')] # Pippen/beads tests
    df_ = df_.loc[~(df_['Library series ID'] == 'L2-JAX-N2-dnase-SX')] # Pippen/beads tests
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA12-N2-mnase-S2')] # Djem's time course
    df_.reset_index(drop=True, inplace=True)
    df_.to_csv(fp_out, sep='\t', index=None)

    # Expand DNase/MNase samples table to have one row for each dataset (instead of each digestion series)
    if expanded:
        # Check that Bioinformatics ID(s) and annotations are consistent
        for i, r in df_.iterrows():
            lid = r['Library series ID']
            l_bid = str(r['Bioinformatics ID(s)']).split(';')
            l_cnd = str(r['Digestion conditions (homogenized)']).split(';')
            assert len(l_bid) == len(l_cnd), 'bid/annot lists do not match for library %s' % (lid,) 

        # Expand dataset_id, annotation labels
        # Adapted from: https://stackoverflow.com/questions/17116814/pandas-how-do-i-split-text-in-a-column-into-multiple-rows/17116976#17116976
        s1 = df_['Bioinformatics ID(s)'].str.split(';').apply(pd.Series, 1).stack()
        s1.index = s1.index.droplevel(-1)
        s1.name = 'Bioinformatics ID(s)'
        s2 = df_['Digestion conditions (homogenized)'].str.split(';').apply(pd.Series, 1).stack()
        s2.index = s2.index.droplevel(-1)
        s2.name = 'Digestion conditions (homogenized)'
        assert (s1.index == s2.index).all()
        del df_['Bioinformatics ID(s)']
        del df_['Digestion conditions (homogenized)']
        df_ = df_.join(pd.concat([s1, s2], axis=1)).reset_index(drop=True)

        # Attach GEO ID (from config.yaml)
        df_geo_ = pd.read_csv('dnase_mnase819_geo/dnase_mnase_geo1_samples.tsv', sep='\t')
        d_ = {k: v for k, v in zip(df_geo_['bid'], df_geo_['title'])}
        df_.insert(loc=2, column='geo_id', value=[ *map(lambda bid: d_.get(bid, ''), df_['Bioinformatics ID(s)']) ])

    return df_

def l_dnase_mnase_samples():
    #if not(os.path.isfile('dnase_mnase/dnase_mnase_samples.tsv')):
    #    make_dnase_mnase_samples()
    #l_bid_dnase_mnase = []
    #for l_bid_ in pd.read_csv('dnase_mnase/dnase_mnase_samples.tsv', sep='\t')['Bioinformatics ID(s)']:
    #    for bid in l_bid_.split(';'):
    #        l_bid_dnase_mnase.append(bid)
    #return l_bid_dnase_mnase
    return make_dnase_mnase_samples(expanded=True)['Bioinformatics ID(s)'].tolist()

def dnase_mnase_label_bid_make():
    dnase_mnase_label_bid = {}
    df_ = make_dnase_mnase_samples()#.head(5)
    for i, r in df_.iterrows():
        libid = r['Library series ID'].rstrip()
        l_bid = r['Bioinformatics ID(s)'].rstrip().split(';')
        l_annot = r['Digestion conditions (homogenized)'].rstrip().split(';')
        for bid, annot in zip(l_bid, l_annot):
            label = '%s/%s_%s_%s' % (libid, libid, annot, bid)
            #print(label, bid, libid, annot)
            #print(pf(label, 'spmr_lt300', '.bw', 'dnase_mnase_labelled'))
            dnase_mnase_label_bid[label] = bid
    return dnase_mnase_label_bid

d_dnase_mnase_label_bid = dnase_mnase_label_bid_make()

rule dmnase_stats:
    input:
        expand(pf('{bid}', '{step}', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples(),
            step=[
                'c_r1',
                'tg_pe.bwa_pe.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.c_r1',
                'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_r1',
                #'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr',
                #'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr_rmdup',
                ]),

    output:
        'dnase_mnase/dmnase_stats_raw_counts.tsv', # raw read counts at each step
        'dnase_mnase/dmnase_stats_frac_passed.tsv',
        'processed_tracks/dmnase_ce10_stats.tsv', # percentages that passed each step
    run:
        df = make_dnase_mnase_samples(expanded=True)[['Bioinformatics ID(s)', 'Library series ID', 'geo_id', 'Digestion conditions (homogenized)']].set_index('Bioinformatics ID(s)') #pd.DataFrame()
        df.index.name = 'dataset_id'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            df.ix[bid, step] = read_int(pf(bid, step, suffix, prefix))
        df.to_csv(output[0], sep='\t')

        def pct_(a, b): return list(map(lambda a_i, b_i: '%.01f%%' % (100.0 * a_i / b_i,), a, b))
        def loss_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * (a_i - b_i) / a_i,), df[col_a], df[col_b]))
        def keep_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * b_i / a_i,), df[col_a], df[col_b]))

        df_ = pd.DataFrame()
        df_['raw_reads'] = df['c_r1'].astype(int).astype(int).map(yp.f_uk)
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
        df_['useful_reads'] = df['tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_r1'].astype(int).astype(int).map(yp.f_uk)
        #df_['complexity'] = (df['tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr_rmdup'].astype(int) \
        #                   / df['tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.c_fr'].astype(int)).map('{:,.02f}'.format)
        df_.sort_index().to_csv(output[1], sep='\t')

        df_.insert(loc=0, column='library_series_id', value=df['Library series ID'])
        df_.insert(loc=1, column='digest_conditions', value=df['Digestion conditions (homogenized)'])
        df_.insert(loc=2, column='geo_id', value=df['geo_id'])
        df_.sort_index().to_csv(output[2], sep='\t')

rule dnase_mnase:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'c_r2', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'tg_pe.bwa_pe', '.bam', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        # GEO -- md5 checksums of fastq files
        expand(pf('{bid}', 'md5sum_r1', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'md5sum_r2', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        # GEO -- read lengths
        expand(pf('{bid}', 'readlen_r1', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'readlen_r2', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        # GEO -- fragment sizes
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.fsizes', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),

rule dmnase_ce10_spmr:
    input:
        lambda wildcards: pf('%s', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'dnase_mnase') % (d_dnase_mnase_label_bid[wildcards.label],),
    output:
        pf('{label}', 'dmnase_ce10_spmr', '.bw', 'processed_tracks'),
    shell:
        'scripts/bigWiggleTools.ipy write {output} scale 0.1 bin 10 {input}'

def dnase_mnase819_geo_read1_(wildcards):
    fn_ = wildcards.raw_file_1
    df_ = pd.read_csv('dnase_mnase819_geo/dnase_mnase_geo1_samples.tsv', sep='\t').query('raw_file_1 == @fn_')
    assert(len(df_) == 1)
    return 'samples/%s.r1.fq.gz' % (df_.iloc[0]['bid'],)

def dnase_mnase819_geo_read2_(wildcards):
    fn_ = wildcards.raw_file_2
    df_ = pd.read_csv('dnase_mnase819_geo/dnase_mnase_geo1_samples.tsv', sep='\t').query('raw_file_2 == @fn_')
    assert(len(df_) == 1)
    return 'samples/%s.r2.fq.gz' % (df_.iloc[0]['bid'],)

rule dnase_mnase819_geo_read1:
    input:
        dnase_mnase819_geo_read1_,
    output:
        'dnase_mnase819_geo/reads/{raw_file_1}'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule dnase_mnase819_geo_read2:
    input:
        dnase_mnase819_geo_read2_,
    output:
        'dnase_mnase819_geo/reads/{raw_file_2}'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

def dnase_mnase819_geo_processed_(wildcards):
    fn_ = wildcards.raw_file_2
    df_ = pd.read_csv('dnase_mnase819_geo/dnase_mnase_geo1_samples.tsv', sep='\t').query('raw_file_2 == @fn_')
    assert(len(df_) == 1)
    return 'samples/%s.r2.fq.gz' % (df_.iloc[0]['bid'],)

rule dnase_mnase819_geo_processed:
    input:
        dnase_mnase819_geo_processed_,
    output:
        'dnase_mnase819_geo/reads/{raw_file_1}'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule dnase_mnase819_geo:
    input:
        expand('dnase_mnase819_geo/reads/{raw_file_1}', raw_file_1=itertools.islice(pd.read_csv('dnase_mnase819_geo/dnase_mnase_geo1_samples.tsv', sep='\t')['raw_file_1'].tolist(), None)),
        expand('dnase_mnase819_geo/reads/{raw_file_2}', raw_file_2=itertools.islice(pd.read_csv('dnase_mnase819_geo/dnase_mnase_geo1_samples.tsv', sep='\t')['raw_file_2'].tolist(), None)),

