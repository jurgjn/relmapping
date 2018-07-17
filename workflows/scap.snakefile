# Short cap-specific processing
def scap_samples_unpooled():
    return [k for k, v in config['scap'].items() if not('pooled_from' in v.keys())]

def scap_samples_unpooled_rm_qcfail():
    return [k for k, v in config['scap'].items() if not('qc_fail' in v.keys()) and not('pooled_from' in v.keys())]

def scap_samples_unpooled_rm_qcfail_lsid():
    return [v['library_series_id'] for k, v in config['scap'].items() if not('qc_fail' in v.keys()) and not('pooled_from' in v.keys())]

def scap_samples_pooled_only():
    return [k for k, v in config['scap'].items() if 'pooled_from' in v.keys()]

def scap_samples_all():
    return config['scap'].keys()

def scap_bid_from_lsid(lsid):
    d_lsid_bid = {config['scap'][bid]['library_series_id']: bid for bid in scap_samples_unpooled()}
    return d_lsid_bid[lsid]

rule firstbp_fwd_scap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output:
        temp(pf('{bid}', '{step}.firstbp_fwd', '.tmp', '{prefix}')),
        pf('{bid}', '{step}.firstbp_fwd', '.bw', '{prefix}'),
    params: ce10_chroms='shared/ce10.chroms'
    shell:
        '''
        samtools view -b {input} | bedtools genomecov -ibam stdin -bga -5 -strand + > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} {params.ce10_chroms} {output[1]}
        '''

rule firstbp_rev_scap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output:
        temp(pf('{bid}', '{step}.firstbp_rev', '.tmp', '{prefix}')),
        pf('{bid}', '{step}.firstbp_rev', '.bw', '{prefix}'),
    params: ce10_chroms='shared/ce10.chroms'
    shell:
        '''
        samtools view -b {input} | bedtools genomecov -ibam stdin -bga -5 -strand - > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} {params.ce10_chroms} {output[1]}
        '''

rule firstbp_fwd_ce11_scap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output:
        temp(pf('{bid}', '{step}.firstbp_fwd_ce11', '.tmp', '{prefix}')),
        pf('{bid}', '{step}.firstbp_fwd_ce11', '.bw', '{prefix}'),
    params: ce11_chroms='shared/ce11.chroms'
    shell:
        '''
        samtools view -b {input} | bedtools genomecov -ibam stdin -bga -5 -strand + > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} {params.ce11_chroms} {output[1]}
        '''

rule firstbp_rev_ce11_scap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output:
        temp(pf('{bid}', '{step}.firstbp_rev_ce11', '.tmp', '{prefix}')),
        pf('{bid}', '{step}.firstbp_rev_ce11', '.bw', '{prefix}'),
    params: ce11_chroms='shared/ce11.chroms'
    shell:
        '''
        samtools view -b {input} | bedtools genomecov -ibam stdin -bga -5 -strand - > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} {params.ce11_chroms} {output[1]}
        '''

rule bw_neg:
    input:
        pf('{bid}', '{step}', '.bw', '{prefix}'),
    output:
        pf('{bid}', '{step}.neg', '.bw', '{prefix}'),
    shell:
        'scripts/bigWiggleTools.ipy write_bg {output} "scale -1" {input}'

rule bw_neg_ce11:
    input:
        pf('{bid}', '{step}', '.bw', '{prefix}'),
    output:
        pf('{bid}', '{step}.neg_ce11', '.bw', '{prefix}'),
    shell:
        'scripts/bigWiggleTools_ce11.ipy write_bg {output} "scale -1" {input}'

rule keep:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.keep', '.bam', '{prefix}'),
        pf('{bid}', '{step}.drop', '.bam', '{prefix}'),
    shell:
        # 4 = read unmapped
        'samtools view -b -F 4 -L shared/ce10_keep.bed -U {output[1]} {input} > {output[0]}'

rule scap_1M_stats:
    input:
        expand(pf('{bid}_1M', '{step}', '.txt', 'scap'), bid=scap_samples_unpooled(),
            step=[
                'c_r1', # Total reads
                'tg_se_q0.c_r1', # Adapter-trimmed reads
                'tg_se.c_r1', # Adapter-trimmed + quality-trimmed reads
                'tg_se.bwa_se.c', # Total reads
                #'tg_se.bwa_se_ecoli.c', # All reads in the fastq file -- sanity check
                'tg_se.bwa_se_ecoli.c_aln',
                'tg_se.bwa_se_ecoli.c_q10',
                'tg_se.bwa_se.c_aln', # Aligned reads
                'tg_se.bwa_se.c_chrM', # Aligned chrM
                'tg_se.bwa_se.c_blacklist', # Aligned chrM
                'tg_se.bwa_se.keep.c', # Aligned excl chrM, blacklist
                'tg_se.bwa_se.q10_keep.c',
                'tg_se.bwa_se.q10_keep.c_atac',
                'trim20.bwa_se.q10_keep.c_atac', # Alternative brute-force trim
                #'tg_se_q10.bwa_se.q10_keep.c',
                #'tg_se_q15.bwa_se.q10_keep.c',
                #'tg_se_q20.bwa_se.q10_keep.c',
                #'tg_se_len16.bwa_se.q10_keep.c'
                ]),
    output:
        'scap/scap_1M_stats.tsv'
    run:
        df = pd.DataFrame()
        df.index.name = 'sample'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            df.ix[bid, step] = '%d' % (read_int(pf(bid, step, suffix, prefix)),)
            #df.ix[bid, step] = read_int(pf(bid, step, suffix, prefix))
        #df.sort_index(axis=1, inplace=True)
        df.sort_index(axis=0, inplace=True)
        df.to_csv(output[0], sep='\t')

rule tics:
    # naive TIC caller that merges anything within 50bp (including singletons -- unlike Chen 2013)
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.tics', '.bed', 'scap'),
    shell: '''
        bigWigToBedGraph {input[0]} stdout | awk '(int($4)) > 0' | bedtools merge -i stdin -d 50 -o sum -c 4 | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,"tic",$4,"+"}}' > {output}
        bigWigToBedGraph {input[1]} stdout | awk '(int($4)) > 0' | bedtools merge -i stdin -d 50 -o sum -c 4 | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,"tic",$4,"-"}}' >> {output}
        sort -k1,1 -k2,2n -k3,3n -k6,6 {output} -o {output}
    '''

rule tics_top5:
    input:
        pf('{bid}', '{step}.tics', '.bed', 'scap'),
    output:
        pf('{bid}', '{step}.tics_top5', '.bed', 'scap'),
        pf('{bid}', '{step}.c_rm_tics_top5', '.txt', 'scap'),
    run:
        df = pd.read_csv(input[0], names=['chrom', 'start', 'end', 'name', 'score', 'strand'], sep='\t')
        df['frac'] = (100.0 * df['score'] / df['score'].sum()).map('{:,.01f}%'.format)
        df.sort_values('score', inplace=True, ascending=False)
        df.head(5).to_csv(output[0], sep='\t', index=False)
        with open(output[1], 'w') as fh:
            print(df['score'].sum() - df.head(5)['score'].sum(), file=fh)

rule tics_top20:
    input:
        pf('{bid}', '{step}.tics', '.bed', 'scap'),
    output:
        pf('{bid}', '{step}.tics_top20', '.bed', 'scap'),
        pf('{bid}', '{step}.c_rm_tics_top20', '.txt', 'scap'),
    run:
        df = pd.read_csv(input[0], names=['chrom', 'start', 'end', 'name', 'score', 'strand'], sep='\t')
        df['frac'] = (100.0 * df['score'] / df['score'].sum()).map('{:,.01f}%'.format)
        df.sort_values('score', inplace=True, ascending=False)
        df.head(20).to_csv(output[0], sep='\t', index=False)
        with open(output[1], 'w') as fh:
            print(df['score'].sum() - df.head(20)['score'].sum(), file=fh)

rule tics_top100:
    input:
        pf('{bid}', '{step}.tics', '.bed', 'scap'),
    output:
        pf('{bid}', '{step}.tics_top100', '.bed', 'scap'),
        pf('{bid}', '{step}.c_rm_tics_top100', '.txt', 'scap'),
    run:
        df = pd.read_csv(input[0], names=['chrom', 'start', 'end', 'name', 'score', 'strand'], sep='\t')
        df['frac'] = (100.0 * df['score'] / df['score'].sum()).map('{:,.01f}%'.format)
        df.sort_values('score', inplace=True, ascending=False)
        df.head(100).to_csv(output[0], sep='\t', index=False)
        with open(output[1], 'w') as fh:
            print(df['score'].sum() - df.head(100)['score'].sum(), file=fh)

rule c_capstack:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.c_capstack', '.txt', 'scap'),
    run:
        dm_step = 'dm160615'
        df_sites = pd.read_csv('../HTSProcessing/%(dm_step)s/hs_sites.atac.tsv' % locals(), sep='\t')
        (fp_fwd, fp_rev) = (input[0], input[1])
        df_sites['_bid'] = list(map(max, zip(
            read_regions(fp_fwd, df_sites['chrom'].tolist(), (df_sites['start'] - 50).tolist(), (df_sites['end'] + 50).tolist(), f=np.nanmax),
            read_regions(fp_rev, df_sites['chrom'].tolist(), (df_sites['start'] - 50).tolist(), (df_sites['end'] + 50).tolist(), f=np.nanmax)
        )))
        n_ = sum(map(lambda x: 1 if x >= 2 else 0, df_sites['_bid']))
        with open(output[0], 'w') as fh:
            print(n_, file=fh)

rule c_chen_kruesi_consensus:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.c_chen_kruesi_consensus', '.txt', 'scap'),
    run:
        flank_len = 50
        (fp_fwd, fp_rev) = (input[0], input[1])
        df_ = pd.read_csv('../HTSProcessing/_chen_kruesi_l3_consensus_tss.tsv', sep='\t')
        def nanargmax_(a):
            try:
                return np.nanargmax(a)
            except ValueError:
                return float('NaN')

        df_['scap_mode_fwd'] = df_['start'] - flank_len + list(read_regions(fp_fwd, df_['chrom'].tolist(), (df_['start'] - flank_len).tolist(), (df_['end'] + flank_len).tolist(), f=nanargmax_))
        df_['scap_mode_rev'] = df_['start'] - flank_len + list(read_regions(fp_rev, df_['chrom'].tolist(), (df_['start'] - flank_len).tolist(), (df_['end'] + flank_len).tolist(), f=nanargmax_))
        df_fwd = df_.query('strand == "+"')
        df_rev = df_.query('strand == "-"')
        #print(df_fwd[['chrom', 'start', 'end', 'strand', 'scap_mode_fwd', 'scap_mode_rev']].head(20))
        #print(df_rev[['chrom', 'start', 'end', 'strand', 'scap_mode_fwd', 'scap_mode_rev']].head(20))
        frac_ = (len(df_fwd.query('start == scap_mode_fwd')) + len(df_rev.query('start == scap_mode_rev'))) / float(len(df_))
        with open(output[0], 'w') as fh:
            print('%d' % (100*frac_,), file=fh)

rule c_kruesi_l3:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.c_kruesi_l3', '.txt', 'scap'),
    run:
        flank_len = 50
        (fp_fwd, fp_rev) = (input[0], input[1])
        df_ = pd.read_csv('../HTSProcessing/_kruesi_l3_tss.tsv', sep='\t')
        def nanargmax_(a):
            try:
                return np.nanargmax(a)
            except ValueError:
                return float('NaN')

        df_['scap_mode_fwd'] = df_['start'] - flank_len + list(read_regions(fp_fwd, df_['chrom'].tolist(), (df_['start'] - flank_len).tolist(), (df_['end'] + flank_len).tolist(), f=nanargmax_))
        df_['scap_mode_rev'] = df_['start'] - flank_len + list(read_regions(fp_rev, df_['chrom'].tolist(), (df_['start'] - flank_len).tolist(), (df_['end'] + flank_len).tolist(), f=nanargmax_))
        df_fwd = df_.query('strand == "+"')
        df_rev = df_.query('strand == "-"')
        #print(df_fwd[['chrom', 'start', 'end', 'strand', 'scap_mode_fwd', 'scap_mode_rev']].head(20))
        #print(df_rev[['chrom', 'start', 'end', 'strand', 'scap_mode_fwd', 'scap_mode_rev']].head(20))
        frac_ = (len(df_fwd.query('start == scap_mode_fwd')) + len(df_rev.query('start == scap_mode_rev'))) / float(len(df_))
        with open(output[0], 'w') as fh:
            print('%d' % (100*frac_,), file=fh)

rule c_TCA:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.c_TCA', '.txt', 'scap'),
    run:
        fp_ = 'annot/S1_accessible_sites/S1a_accessible_sites.tsv'
        df_sites = pd.read_csv(fp_, sep='\t').query('atac_source != "atac_glp1_se"')

        fp_fwd = str(input[0])
        fp_rev = str(input[1])

        flank_len = 50
        df_sites['scap_mode_raw_fwd'] = df_sites['start'] - flank_len + list(map(yp.nanargmax_median, yp.read_regions(fp_fwd, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len)))
        df_sites['scap_mode_raw_rev'] = df_sites['start'] - flank_len + list(map(yp.nanargmax_median, yp.read_regions(fp_rev, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len)))
        df_sites['scap_mode_count_fwd'] = list(map(lambda c: int(np.nanmax(c)) if not(np.isnan(np.nanmax(c))) else 0, yp.read_regions(fp_fwd, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len)))
        df_sites['scap_mode_count_rev'] = list(map(lambda c: int(np.nanmax(c)) if not(np.isnan(np.nanmax(c))) else 0, yp.read_regions(fp_rev, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len)))

        df_sites['scap_mode_seq_fwd'] = list(map(lambda s: s.upper(), yp.read_seq(df_sites['chrom'], df_sites['scap_mode_raw_fwd'] - 2, df_sites['scap_mode_raw_fwd'] + 1)))
        df_sites.loc[df_sites['scap_mode_count_fwd'] == 0, 'scap_mode_seq_fwd'] = 'NNN'
        n_TCA_fwd = sum(df_sites['scap_mode_seq_fwd'] == 'TCA')
        #n_exp = len(df_sites) / (4**3)
        #r_enr = n_TCA / n_exp
        #print(n_TCA, n_exp, r_enr)

        df_sites['scap_mode_seq_rev'] = list(map(lambda s: s.upper(), yp.read_seq(df_sites['chrom'], df_sites['scap_mode_raw_rev'], df_sites['scap_mode_raw_rev'] + 3)))
        df_sites.loc[df_sites['scap_mode_count_rev'] == 0, 'scap_mode_seq_rev'] = 'NNN'
        n_TCA_rev = sum(df_sites['scap_mode_seq_rev'] == 'TGA') # reverse complement of TCA
        #n_exp = len(df_sites) / (4**3)
        #r_enr = n_TCA / n_exp
        #print(n_TCA, n_exp, r_enr)

        with open(output[0], 'w') as fh:
            print(n_TCA_fwd + n_TCA_rev, file=fh)

rule scap_qc:
    input:
        expand(pf('{bid}', '{step}', '.txt', 'scap'), bid=config['scap'].keys(),
            step=[
                'c_r1', # Total reads
                'tg_se.c_r1', # Adapter-trimmed + quality-trimmed reads (#1)
                'tg_se.bwa_se.c', # Adapter-trimmed + quality-trimmed reads (#2)
                'tg_se.bwa_se.rm_unmapped.c', # Aligned reads
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c', # Aligned chrM
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c', # Blacklisted
                ####
                #'tg_se_q0.c_r1', # Adapter-trimmed reads
                #'tg_se.bwa_se.c', # Total reads
                #'tg_se.bwa_se_ecoli.c', # All reads in the fastq file -- sanity check
                #'tg_se.bwa_se_ecoli.c_aln',
                #'tg_se.bwa_se_ecoli.c_q10',
                #'tg_se.bwa_se.keep.c', # Aligned excl chrM, blacklist
                #'tg_se.bwa_se.q10_keep.c',
                #'tg_se.bwa_se.q10_keep.c_atac',
                #'tg_se_q10.bwa_se.q10_keep.c',
                #'tg_se_q15.bwa_se.q10_keep.c',
                #'tg_se_q20.bwa_se.q10_keep.c',
                #'tg_se_len16.bwa_se.q10_keep.c'
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.c', # Highly abundant, non-coding
                #'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_rm_tics_top5',
                #'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_rm_tics_top20',
                #'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_rm_tics_top100',
                #'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_capstack',
                #'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_chen_kruesi_consensus',
                #'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_kruesi_l3',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_TCA',                
                ]),
    output:
        'scap/scap_qc_counts.tsv',
        'scap/scap_qc_leaks.tsv',
    run:
        df = pd.DataFrame()
        df.index.name = 'bid'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            #df.ix[bid, step] = '%d' % (read_int(pf(bid, step, suffix, prefix)),)
            df.ix[bid, step] = read_int(pf(bid, step, suffix, prefix))
        #df.sort_index(axis=1, inplace=True)
        #df.sort_index(axis=0, inplace=True)
        df.to_csv(output[0], sep='\t')

        """

        l_raw = df['c_r1'].astype(int) # Raw reads, as they came off the sequencer
        l_tgr = df['c_r1'] - df['tg_se.c_r1'] # Reads removed by trim_galore
        l_unm = df['tg_se.c_r1'] - df['tg_se.bwa_se.c_aln'] # Reads not mapped to the worm genome
        l_chm = df['tg_se.bwa_se.c_chrM'] # Reads mapped to chrM
        l_blk = df['tg_se.bwa_se.c_blacklist'] # Reads mapped to chrM

        #rm_tg -- either adapter dimer, or low-quality
        ##pct_(l_tgr, l_raw)
        """
        def pct_(a, b): return list(map(lambda a_i, b_i: '%.01f%%' % (100.0 * a_i / b_i,), a, b))
        def loss_pct_(col_a, col_b): return list(map(lambda a_i, b_i: '%.02f%%' % (100.0 * (a_i - b_i) / a_i,), df[col_a], df[col_b]))

        df_ = pd.DataFrame()
        df_['raw_reads'] = df['c_r1'].astype(int).map(yp.f_uk)
        df_['adapter_or_lowqual'] = loss_pct_('c_r1', 'tg_se.c_r1')
        df_['unmapped'] = loss_pct_('tg_se.bwa_se.c', 'tg_se.bwa_se.rm_unmapped.c')
        df_['mitochondrial'] = loss_pct_('tg_se.bwa_se.rm_unmapped.c','tg_se.bwa_se.rm_unmapped.rm_chrM.c')
        df_['blacklist'] = loss_pct_('tg_se.bwa_se.rm_unmapped.rm_chrM.c', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c')
        df_['mapq_lt10'] = loss_pct_('tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c')
        df_['non_coding'] = loss_pct_('tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.c')        
        df_['useful_reads'] = df['tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.c'].astype(int).map(yp.f_uk)
        #df_['tics_top5'] = loss_pct_('tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_rm_tics_top5')
        #df_['tics_top20'] = loss_pct_('tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_rm_tics_top20')
        df_['inr_sites'] = df['tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c_TCA'].astype(int).map(yp.f_uk)

        """
        df_['rm_tg_frac'] = l_tgr
        df_['not_mapped'] = l_unm
        df_['mitochondrial'] = l_chm
        df_['blacklist'] = l_blk

        #samtools mapped:
        #chrM -- samtools mapped chrM
        #mapq0 -- samtools mapped not chrM mapq zero
        #mapq_low -- samtools mapped not chrM
        #highly_abundant
        #useful
        """
        df_.to_csv(output[1], sep='\t')

def scap_wt_qc_pass():
    df_ = pd.DataFrame(config['scap']).transpose().query('(strain == "N2") & (qc_fail != qc_fail)')
    df_.index.name = 'bid'
    lsid_counts = collections.Counter(df_['library_series_id'])
    for k, v in lsid_counts.items():
        assert v == 1, 'mutliple occurrences of lsid=%s' % (k,)
    return df_

rule scap_samples:
    input:
        'scap/scap_qc_leaks.tsv',
    output:
        'scap/scap_samples.tsv',
    run:
        l_col = ['library_series_id', 'collection_id', 'strain', 'stage', 'enzyme', 'rna_type', 'rna_purification', 'qc_fail']#, 'pooled_from', 'scap649', 'scap740',
        # 'original_name',
        # 'manual_demultiplex_of'
        df = scap_wt_qc_pass()
        df_stats = pd.read_csv(str(input), sep='\t', index_col='bid')
        df.ix[config['scap'].keys()].merge(df_stats, left_index=True, right_index=True, how='left').to_csv(str(output), sep='\t')

rule plot_scap:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
        pf('chen13_scap_rep1a', '{step}.firstbp_fwd', '.bw', 'scap'),
        pf('chen13_scap_rep1a', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.plot_scap', '.pdf', 'scap'),
    run:
        def mean_log2_fwd(x, a):
            x_ = np.mean(np.log2(x + 1), a)
            return x_ / x_.sum()

        def mean_log2_rev(x, a):
            x_ = np.mean(np.log2(x + 1), a)
            return -x_ / x_.sum()
            #return -np.mean(np.log2(x + 1), a)

        fp_ = 'annot/Fig1D1_accessible_sites/Fig1D1_accessible_sites.tsv'
        df_sites = pd.read_csv(fp_, sep='\t')
        bid = wildcards.bid
        #bid = 'HS252_Yan_S2'
        #bid_ref = 'chen13_scap_rep1a'
        #step = 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10'
        #step_fwd = step + '.firstbp_fwd'
        #step_rev = step + '.firstbp_rev'
        #fp_fwd = pf(bid, step_fwd, '.bw', 'scap')
        #fp_rev = pf(bid, step_rev, '.bw', 'scap')
        #fp_ref_fwd = pf(bid_ref, step_fwd, '.bw', 'scap')
        #fp_ref_rev = pf(bid_ref, step_rev, '.bw', 'scap')
        fp_fwd = input[0]
        fp_rev = input[1]
        fp_ref_fwd = input[2]
        fp_ref_rev = input[3]

        gdf = yp.GenomicDataFrame(df_sites)
        gdf.add_track(bid + '_fwd', fp_fwd , bin_size=1, flank_len=400, memoized=False)
        gdf.add_track(bid + '_rev', fp_rev , bin_size=1, flank_len=400, memoized=False)
        gdf.add_track('chen13_scap_rep1a_fwd', fp_ref_fwd , bin_size=1, flank_len=400, memoized=False)
        gdf.add_track('chen13_scap_rep1a_rev', fp_ref_rev , bin_size=1, flank_len=400, memoized=False)

        bid_ref = 'chen13_scap_rep1a'
        plt.figure(figsize=(8, 8))
        plt.title(bid)
        gdf.t[bid + '_fwd'].plot(f=mean_log2_fwd, label='fwd: %s' % (bid,), color=yp.RED)
        gdf.t[bid + '_rev'].plot(f=mean_log2_rev, label='rev: %s' % (bid,), color=yp.BLUE)
        gdf.t[bid_ref + '_fwd'].plot(f=mean_log2_fwd, label='fwd: %s' % (bid_ref,), color=yp.RED, alpha=0.5)
        gdf.t[bid_ref + '_rev'].plot(f=mean_log2_rev, label='rev: %s' % (bid_ref,), color=yp.BLUE, alpha=0.5)
        plt.gca().legend(loc='upper left')
        plt.axhline(0, color='k')
        plt.savefig(str(output), bbox_inches='tight')

rule scap:
    input:
        'scap/scap_samples.tsv',
        expand(pf('{bid}', 'c_r1', '.txt', 'scap'), bid=config['scap'].keys()),
        expand(pf('{bid}', 'tg_se.bwa_se', '.bam', 'scap'), bid=config['scap'].keys()),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist', '.bam.bai', 'scap'), bid=config['scap'].keys()),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10', '.bam.bai', 'scap'), bid=config['scap'].keys()),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd', '.bw', 'scap'), bid=config['scap'].keys()),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev', '.bw', 'scap'), bid=config['scap'].keys()),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.plot_scap', '.pdf', 'scap'), bid=config['scap'].keys()),
        #expand(pf('{lsid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.lsid', '.bw', 'scap'), lsid=scap_samples_unpooled_rm_qcfail_lsid()),
        #expand(pf('{lsid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.lsid', '.bw', 'scap'), lsid=scap_samples_unpooled_rm_qcfail_lsid()),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_fwd.rm_exonic', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_rev.rm_exonic', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10_reversed.firstbp_fwd', '.bw', 'scap'), bid=scap_samples_unpooled_rm_qcfail()),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10_reversed.firstbp_rev', '.bw', 'scap'), bid=scap_samples_unpooled_rm_qcfail()),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.tics', '.bed', 'scap'), bid=['scap541_emb_l3_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.tics_top100', '.bed', 'scap'), bid=['scap541_emb_l3_ya']),

rule scap815_sum_by_stage:
    input:
        pf('scap815_{stage}_rep1', '{step}', '.bw', 'scap815'),
        pf('scap815_{stage}_rep2', '{step}', '.bw', 'scap815'),
    output:
        pf('scap815_{stage}', '{step}.sum_by_stage', '.bw', 'scap815'),
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]}
    '''

rule scap815_sum_by_stage_ce11:
    input:
        pf('scap815_{stage}_rep1', '{step}', '.bw', 'scap815'),
        pf('scap815_{stage}_rep2', '{step}', '.bw', 'scap815'),
    output:
        pf('scap815_{stage}', '{step}.sum_by_stage_ce11', '.bw', 'scap815'),
    shell: '''
        scripts/bigWiggleTools_ce11.ipy write_bg {output[0]} sum {input[0]} {input[1]}
    '''

rule scap815_gt0x2_all:
    input:
        expand(pf('scap815_{stage}', '{{step}}', '.bw', 'scap815'), stage=config['stages_wt_rep_scap815']),
    output:
        pf('scap815_wt_all', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_all', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_all', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} {input[12]} {input[13]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]} gt 0 {input[2]} gt 0 {input[3]} gt 0 {input[4]} gt 0 {input[5]} gt 0 {input[6]} gt 0 {input[7]} gt 0 {input[8]} gt 0 {input[9]} gt 0 {input[10]} gt 0 {input[11]} gt 0 {input[12]} gt 0 {input[13]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_all_ce11:
    input:
        expand(pf('scap815_{stage}', '{{step}}', '.bw', 'scap815'), stage=config['stages_wt_rep_scap815']),
    output:
        pf('scap815_wt_all', '{step}.gt0x2_ce11', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_all', '{step}.gt0x2_ce11', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_all', '{step}.gt0x2_ce11', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools_ce11.ipy write_bg {output[0]} sum {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} {input[12]} {input[13]}
        scripts/bigWiggleTools_ce11.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]} gt 0 {input[2]} gt 0 {input[3]} gt 0 {input[4]} gt 0 {input[5]} gt 0 {input[6]} gt 0 {input[7]} gt 0 {input[8]} gt 0 {input[9]} gt 0 {input[10]} gt 0 {input[11]} gt 0 {input[12]} gt 0 {input[13]}
        scripts/bigWiggleTools_ce11.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_emb:
    input:
        expand(pf('scap815_wt_emb_{rep}', '{{step}}', '.bw', 'scap815'), rep=['rep1', 'rep2', 'rep3', 'rep4']),
    output:
        pf('scap815_wt_emb', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_emb', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_emb', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]} {input[2]} {input[3]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]} gt 0 {input[2]} gt 0 {input[3]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_l1:
    input:
        expand(pf('scap815_wt_l1_{rep}', '{{step}}', '.bw', 'scap815'), rep=['rep1', 'rep2']),
    output:
        pf('scap815_wt_l1', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_l1', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_l1', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_l2:
    input:
        expand(pf('scap815_wt_l2_{rep}', '{{step}}', '.bw', 'scap815'), rep=['rep1', 'rep2']),
    output:
        pf('scap815_wt_l2', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_l2', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_l2', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_l3:
    input:
        expand(pf('scap815_wt_l3_{rep}', '{{step}}', '.bw', 'scap815'), rep=['rep1', 'rep2']),
    output:
        pf('scap815_wt_l3', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_l3', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_l3', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_l4:
    input:
        expand(pf('scap815_wt_l4_{rep}', '{{step}}', '.bw', 'scap815'), rep=['rep1', 'rep2']),
    output:
        pf('scap815_wt_l4', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_l4', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_l4', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_gt0x2_ya:
    input:
        expand(pf('scap815_wt_ya_{rep}', '{{step}}', '.bw', 'scap815'), rep=['rep1', 'rep2']),
    output:
        pf('scap815_wt_ya', '{step}.gt0x2', '_sum.bw', 'scap815'), # sum of signal across all replicates
        pf('scap815_wt_ya', '{step}.gt0x2', '_gt0_count.bw', 'scap815'), # number of replicates with non-zero signal
        pf('scap815_wt_ya', '{step}.gt0x2', '.bw', 'scap815'), # summed signal from base pairs with non-zero coverage in at least two replicates
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} sum {input[0]} {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[1]} sum gt 0 {input[0]} gt 0 {input[1]}
        scripts/bigWiggleTools.ipy write_bg {output[2]} mult {output[0]} gt 1 {output[1]}
    '''

rule scap815_read1:
    input:
        'samples/scap815_{replicate}.r1.fq.gz'
    output:
        'scap815_geo/reads/scap_{replicate}.fastq.gz'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_fwd:
    input:
        pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_fwd.sum_by_stage', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_fwd/scap_{stage}_fwd.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_rev:
    input:
        pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.sum_by_stage.neg', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_rev/scap_{stage}_rev.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_fwd_ce11:
    input:
        pf('scap815_{stage}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_fwd_ce11.sum_by_stage_ce11', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_ce11_fwd/scap_{stage}_ce11_fwd.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_rev_ce11:
    input:
        pf('scap815_{stage}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11.sum_by_stage_ce11.neg_ce11', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_ce11_rev/scap_{stage}_ce11_rev.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_all_fwd:
    input:
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_fwd.gt0x2', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_fwd/scap_wt_all_fwd.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_all_rev:
    input:
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.gt0x2.neg', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_rev/scap_wt_all_rev.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_all_ce11_fwd:
    input:
        pf('scap815_wt_all', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_fwd_ce11.gt0x2_ce11', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_ce11_fwd/scap_wt_all_ce11_fwd.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815_all_ce11_rev:
    input:
        pf('scap815_wt_all', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11.gt0x2_ce11.neg_ce11', '.bw', 'scap815'),
    output:
        'scap815_geo/tracks_ce11_rev/scap_wt_all_ce11_rev.bw'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap815:
    input:
        expand(pf('scap815_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_fwd', '.bw', 'scap815'), sample=techreps_collapse(config['scap815'].keys(), include_raw=True)),
        expand(pf('scap815_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev', '.bw', 'scap815'), sample=techreps_collapse(config['scap815'].keys(), include_raw=True)),
        expand(pf('scap815_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.c', '.txt', 'scap815'), sample=techreps_collapse(config['scap815'].keys(), include_raw=True)),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_fwd.sum_by_stage', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.sum_by_stage', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.sum_by_stage.neg', '.bw', 'scap815'), stage=config['stages_wt']),
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_fwd.gt0x2', '.bw', 'scap815'),
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.gt0x2', '.bw', 'scap815'),
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.gt0x2.neg', '.bw', 'scap815'),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_fwd.gt0x2', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.gt0x2', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.rm_q10.firstbp_rev.gt0x2.neg', '.bw', 'scap815'), stage=config['stages_wt']),
        # GEO submission -- reads
        expand('scap815_geo/reads/scap_{replicate}.fastq.gz', replicate=config['scap815'].keys()),
        # GEO submission -- coverage tracks
        expand('scap815_geo/tracks_fwd/scap_{stage}_fwd.bw', stage=config['stages_wt'] + ['wt_all']),
        expand('scap815_geo/tracks_rev/scap_{stage}_rev.bw', stage=config['stages_wt'] + ['wt_all']),

rule scap815_ce11:
    input:
        expand(pf('scap815_{sample}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_fwd_ce11', '.bw', 'scap815'), sample=techreps_collapse(config['scap815'].keys(), include_raw=True)),
        expand(pf('scap815_{sample}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11', '.bw', 'scap815'), sample=techreps_collapse(config['scap815'].keys(), include_raw=True)),
        pf('scap815_wt_all', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_fwd_ce11.gt0x2_ce11', '.bw', 'scap815'),
        pf('scap815_wt_all', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11.gt0x2_ce11', '.bw', 'scap815'),
        pf('scap815_wt_all', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11.gt0x2_ce11.neg_ce11', '.bw', 'scap815'),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_fwd_ce11.sum_by_stage_ce11', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11.sum_by_stage_ce11', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_non_coding_ce11.rm_q10.firstbp_rev_ce11.sum_by_stage_ce11.neg_ce11', '.bw', 'scap815'), stage=config['stages_wt']),
        # GEO submission -- coverage tracks
        expand('scap815_geo/tracks_ce11_fwd/scap_{stage}_ce11_fwd.bw', stage=config['stages_wt'] + ['wt_all']),
        expand('scap815_geo/tracks_ce11_rev/scap_{stage}_ce11_rev.bw', stage=config['stages_wt'] + ['wt_all']),

rule scap815_mapq0:
    input:
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.firstbp_fwd.gt0x2', '.bw', 'scap815'),
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.firstbp_rev.gt0x2', '.bw', 'scap815'),
        pf('scap815_wt_all', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.firstbp_rev.gt0x2.neg', '.bw', 'scap815'),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.firstbp_fwd.gt0x2', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.firstbp_rev.gt0x2', '.bw', 'scap815'), stage=config['stages_wt']),
        expand(pf('scap815_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_non_coding.firstbp_rev.gt0x2.neg', '.bw', 'scap815'), stage=config['stages_wt']),

def df_scap():
    def processed_id(bid, lid):
        if lid != lid or lid == '':
            return '_%(bid)s' % locals()
        return '%(lid)s_%(bid)s' % locals()
    df_ = pd.DataFrame(config['scap']).transpose()
    df_.index.name = 'bid'
    df_.rename(columns={'library_series_id': 'lid'}, inplace=True)
    df_ = df_.reset_index()
    df_['pid'] = [ *map(processed_id, df_['bid'], df_['lid']) ]
    return df_#[['bid', 'pid', 'lid']]

def scap_ce10_init_fwd_input_(wildcards):
    pid_ = wildcards.pid
    df_ = df_scap().query('pid == @pid_')
    assert(len(df_) == 1)
    return pf(df_.iloc[0]['bid'], 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_fwd', '.bw', 'scap')

def scap_ce10_init_rev_input_(wildcards):
    pid_ = wildcards.pid
    df_ = df_scap().query('pid == @pid_')
    assert(len(df_) == 1)
    return pf(df_.iloc[0]['bid'], 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_rev', '.bw', 'scap')

rule scap_ce10_init_fwd:
    input:
        scap_ce10_init_fwd_input_,
    output:
        pf('{pid}', 'scap_ce10_init_fwd', '.bw', 'processed_tracks'),
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap_ce10_init_rev:
    input:
        scap_ce10_init_rev_input_,
    output:
        pf('{pid}', 'scap_ce10_init_rev', '.bw', 'processed_tracks'),
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule scap_processed:
    input:
        expand(pf('{pid}', 'scap_ce10_init_fwd', '.bw', 'processed_tracks'), pid=[* df_scap()['pid'] ]),
        expand(pf('{pid}', 'scap_ce10_init_rev', '.bw', 'processed_tracks'), pid=[* df_scap()['pid'] ]),
