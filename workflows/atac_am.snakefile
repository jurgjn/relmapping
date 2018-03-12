l_stage_am = ['ya', 'd3', 'd7', 'd10', 'd14']

# Processing steps specific to ATAC-seq
"""
rule atac_am1_pep1:
    input:
        pf('atac_am1_{stage}_rep1', '{step}', '.bam', 'atac_am'),
        pf('atac_am1_{stage}_rep2', '{step}', '.bam', 'atac_am'),
    output:
        pf('atac_am1_{stage}_pep1', '{step}', '.bam', 'atac_am'),
    threads: 4
    shell: '''
        samtools merge -u --threads {threads} - {input[0]} {input[1]} | samtools view -h -S --threads {threads} - \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (rand() < .5)' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule atac_am1_pep2:
    input:
        pf('atac_am1_{stage}_rep1', '{step}', '.bam', 'atac_am'),
        pf('atac_am1_{stage}_rep2', '{step}', '.bam', 'atac_am'),
    output:
        pf('atac_am1_{stage}_pep2', '{step}', '.bam', 'atac_am'),
    threads: 4
    shell: '''
        samtools merge -u --threads {threads} - {input[0]} {input[1]} | samtools view -h -S --threads {threads} - \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (.5 <= rand())' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''
"""

rule atac_am1_am_concave150_idr_raw_oracle: # Generate oracle peak list -- all stages/replicate -- mean coverage + identify concave regions
    input:
        expand(pf('atac_am1_{stage}_rep1', '{{step}}', '_treat_pileup.bw', 'atac_am'), stage=l_stage_am),
        expand(pf('atac_am1_{stage}_rep2', '{{step}}', '_treat_pileup.bw', 'atac_am'), stage=l_stage_am),
    output:
        temp(pf('atac_am1', '{step}.concave150_idr_raw', '.bdg', 'atac_am')),
        pf('atac_am1', '{step}.concave150_idr_raw', '.bw', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr_raw', '.bed', 'atac_am'),
    shell: '''
        wiggletools write_bg {output[0]} mean {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        bigWigToBedGraph {output[1]} stdout | scripts/concave_regions 150 > {output[2]}
    '''

rule concave150_idr_raw_D2:
    input:
        pf('{bid}', '{step}', '_treat_pileup.bw', 'atac_am'),
    output:
        pf('{bid}', '{step}.concave150_idr_raw', '_D2.bw', 'atac_am'),
    shell: '''
        scripts/bigWigConvolve {input} D2_s3x150_1E6 {output}
    '''

rule concave150_idr_raw_oracle_D2: # Annotate oracle peak list by D2 from a specific stage
    input:
        pf('atac_am1', '{step}.concave150_idr_raw', '.bed', 'atac_am'),
        pf('atac_am1_{stage}_{rep}', '{step}.concave150_idr_raw', '_D2.bw', 'atac_am'),
    output:
        pf('atac_am1', '{step}.concave150_idr_raw', '_{stage}_{rep}_D2.bed', 'atac_am'),
    run:
        df = pd.read_csv(input[0], sep='\t', names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        df['score'] = list(map(lambda x: -x, read_regions(input[1], df.chrom, map(int, df.start), map(int, df.end), np.mean)))
        #df['score'] = list(map(lambda x: -x, read_regions(input[1], df.chrom, map(int, df.start), map(int, df.end), np.min)))
        df.to_csv(output[0], header=False, index=False, sep='\t')

rule concave200_idr_raw_oracle_D2: # Annotate oracle peak list by D2 from a specific stage
    input:
        pf('atac_am1', '{step}.concave200_idr_raw', '.bed', 'atac_am'),
        pf('atac_am1_{stage}_{rep}', '{step}.concave200_idr_raw', '_D2.bw', 'atac_am'),
    output:
        pf('atac_am1', '{step}.concave200_idr_raw', '_{stage}_{rep}_D2.bed', 'atac_am'),
    run:
        df = pd.read_csv(input[0], sep='\t', names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        df['score'] = list(map(lambda x: -x, read_regions(input[1], df.chrom, map(int, df.start), map(int, df.end), np.mean)))
        #df['score'] = list(map(lambda x: -x, read_regions(input[1], df.chrom, map(int, df.start), map(int, df.end), np.min)))
        df.to_csv(output[0], header=False, index=False, sep='\t')

rule concave150_idr_raw: # IDR analysis for a single stage+rep/pep; using pooled peaks as the oracle peak list
    input:
        pf('atac_am1', '{step}.concave150_idr_raw', '_{stage}_{rep}1_D2.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr_raw', '_{stage}_{rep}2_D2.bed', 'atac_am'),
    output:
        pf('atac_am1', '{step}.concave150_idr_raw', '_{stage}_{rep}_idr.bed', 'atac_am'),
    params:
        idr_args = '--input-file-type bed --rank score --peak-merge-method max'
    shell: '''
        idr {params.idr_args} \
            --samples \
                <(cat {input[0]} | sort -k5nr,5nr | head -n 100000 | sort -k1,1 -k2,2n) \
                <(cat {input[1]} | sort -k5nr,5nr | head -n 100000 | sort -k1,1 -k2,2n) \
            --output-file {output}
        sort -k1,1 -k2,2n -o {output} {output}
    '''

rule concave150_idr_merge: # Merge stage+rep/pep reproducibility calls
    input:
        pf('atac_am1', '{step}.concave150_idr_raw', '.bed', 'atac_am'), # Generate oracle peak list
        expand(pf('atac_am1', '{{step}}.concave150_idr_raw', '_{stage}_rep_idr.bed', 'atac_am'), stage=l_stage_am),
        expand(pf('atac_am1', '{{step}}.concave150_idr_raw', '_{stage}_pep_idr.bed', 'atac_am'), stage=l_stage_am),
    output:
        pf('atac_am1', '{step}.concave150_idr', '.tsv', 'atac_am'),
    run:
        df = pd.read_csv(input[0], sep='\t', names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
        for fp_ in input[1:1+12]:
            col_prefix = 'atac_' + parse_pf(fp_)[2][1:-8] # TODO danger zone
            print(col_prefix)
            names = [
                'chrom', 'start', 'end', 'name',
                '%(col_prefix)s_scaledIDR' % locals(), 'strand',
                '%(col_prefix)s_localIDR' % locals(), '%(col_prefix)s_globalIDR' % locals(),
                '%(col_prefix)s1_start' % locals(), '%(col_prefix)s1_end' % locals(), '%(col_prefix)s1_score' % locals(),
                '%(col_prefix)s2_start' % locals(), '%(col_prefix)s2_end' % locals(), '%(col_prefix)s2_score' % locals(),
            ]
            df_ = pd.read_csv(fp_, sep='\t', names=names)
            df = df.merge(df_, how='left', on=['chrom', 'start', 'end', 'name', 'strand'])
        df.to_csv(output[0], header=True, index=False, sep='\t')

def concave_idr_th_(input, output, th=0.05):
    df_inp = pd.read_csv(input[0], sep='\t')
    th_globalIDR = -math.log(th, 10)

    df_out = df_inp[['chrom', 'start', 'end']].copy()
    df_attr = pd.DataFrame()
    for stage in l_stage_am:
        n_rep = len(df_inp.query('atac_%(stage)s_rep_globalIDR > %(th_globalIDR)s' % locals()))
        n_pep = len(df_inp.query('atac_%(stage)s_pep_globalIDR > %(th_globalIDR)s' % locals()))
        if n_pep > n_rep: # pseudoreplicate-based peak calls
            df_attr['atac_%(stage)s_idr' % locals()] = list(map(lambda globalIDR: 10**(-globalIDR), df_inp['atac_%(stage)s_pep_globalIDR' % locals()]))
        else: # replicate-based peak calls
            df_attr['atac_%(stage)s_idr' % locals()] = list(map(lambda globalIDR: 10**(-globalIDR), df_inp['atac_%(stage)s_rep_globalIDR' % locals()]))
        n_out = len(df_attr.query('atac_%(stage)s_idr < %(th)s' % locals()))
        print(stage, n_rep, n_pep, n_out)

    def pack_gfftags(df):
        def pack_row(r): return (";".join([("%s=%s" % (k, v)).replace(" ", "%20") for k, v in zip(df.columns, r)]))
        return df.apply(pack_row, axis=1, reduce=False, raw=True)
    df_out['name'] = pack_gfftags(df_attr[['atac_%(stage)s_idr' % locals() for stage in l_stage_am]])
    df_out['score'] = df_attr[['atac_%(stage)s_idr' % locals() for stage in l_stage_am]].min(axis=1)
    #df_inp[['atac_%(stage)s_rep_scaledIDR' % locals() for stage in l_stage_am] + ['atac_%(stage)s_pep_scaledIDR' % locals() for stage in l_stage_am]].max(axis=1)
    df_out['strand'] = '.'
    #df_out['thickStart'] = list(map(int, df_inp[['start', 'end']].mean(axis=1)))
    #df_out['thickEnd'] = df_out['thickStart'] + 1

    def idr_score_threshold(idr_score): return(int(-125*math.log(idr_score, 2)))

    trackline='#track gffTags=on'# useScore=1'
    with open(output[0], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.001').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    with open(output[1], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.003').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    with open(output[2], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.005').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    with open(output[3], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.01').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    with open(output[4], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.05').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    with open(output[5], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.3').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    with open(output[6], 'w') as fh:
        print(trackline, file=fh)
        df_out.query('score < 0.8').to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

rule concave150_idr_th:
    input:
        pf('atac_am1', '{step}.concave150_idr', '.tsv', 'atac_am'),
    output:
        pf('atac_am1', '{step}.concave150_idr', '_0.001.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr', '_0.003.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr', '_0.005.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr', '_0.01.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr', '_0.05.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr', '_0.3.bed', 'atac_am'),
        pf('atac_am1', '{step}.concave150_idr', '_0.8.bed', 'atac_am'),
    run:
        concave_idr_th_(input, output)

rule concave_idr_tracks: # Stage-specific mean coverage
    input:
        pf('atac_am1_{stage}_rep1', '{step}', '_treat_pileup.bw', 'atac_am'),
        pf('atac_am1_{stage}_rep2', '{step}', '_treat_pileup.bw', 'atac_am'),
    output:
        temp(pf('atac_am1_{stage}', '{step}.concave_idr', '.bdg', 'atac_am')),
        pf('atac_am1_{stage}', '{step}.concave_idr', '.bw', 'atac_am'),
    shell: '''
        wiggletools write_bg {output[0]} mean {input[0]} {input[1]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
    '''

#def atac_am_samples():
#    #return list(sorted(config['atac_am'].keys()))
#    return ['HS446_MS1', 'HS446_MS2', 'HS446_MS3', 'HS446_MS4', 'HS446_MS5',
#            'HS453_MS6', 'HS453_MS7', 'HS456_MS8', 'HS456_MS9', 'HS456_MS10',
#            'HS456_MS11', 'HS456_MS12', 'HS456_MS13', 'HS456_MS14', 'HS456_MS15',
#            'HS456_MS16', 'HS456_MS17', 'HS456_MS18', 'HS456_MS19',
#            'HS457_MS19', 'HS457_MS20', 'HS457_MS21', 'HS457_MS22', 'HS457_MS23', 'HS457_MS24', 'HS457_MS25',
#            'HS457_MS26', 'HS457_MS27', 'HS457_MS28', 'HS457_MS29', 'HS457_MS30',
#            'HS462_MS37', 'HS462_MS38', 'HS462_MS39', 'HS462_MS40', 'HS462_MS41', 'HS462_MS42', 'HS462_MS43',
#            ]

#l_macs2_step = [
#    'macs2_se_extsize150_shiftm75', 'macs2_se_extsize150_shiftm75_keepdup_all', 'macs2_se_extsize150_shiftm75_keepdup_auto',
#    'macs2_se_extsize200_shiftm100', 'macs2_se_extsize200_shiftm100_keepdup_all', 'macs2_se_extsize200_shiftm100_keepdup_auto',
#]
l_macs2_step = ['macs2_se_extsize150_shiftm75_keepdup_all']

def atac_am_samples():
    #fp = 'atac_am/atac_am_sample_sheet.tsv'
    #df = pd.read_csv(fp, sep='\s+', comment='#')
    #return set(df['sample_name'].tolist()) | set(df['label'].tolist())
    return []

rule atac_am:
    input:
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep', '.bam', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep', '.bam.bai', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep.{macs2_step}', '_treat_pileup.bw', 'atac_am'), bid=atac_am_samples(), macs2_step=l_macs2_step),

rule atac_am_concave_idr:
    input:
        # Generate oracle peak list
        expand(pf('atac_am1', 'tg_se.bwa_se.q10_keep.{macs2_step}.concave150_idr_raw', '.bed', 'atac_am'), macs2_step=l_macs2_step),
        # Annotate oracle peaks by second derivative from specific stage/rep
        expand(pf('atac_am1', 'tg_se.bwa_se.q10_keep.{macs2_step}.concave150_idr_raw', '_{stage}_{rep}_D2.bed', 'atac_am'), macs2_step=l_macs2_step, stage=l_stage_am, rep=l_rep),
        # Stage/rep-specific IDR analysis
        expand(pf('atac_am1', 'tg_se.bwa_se.q10_keep.{macs2_step}.concave150_idr_raw', '_{stage}_{rep}_idr.bed', 'atac_am'), macs2_step=l_macs2_step, stage=l_stage_am, rep=['rep', 'pep']),
        # Pool stage-rep-specific analyses
        expand(pf('atac_am1', 'tg_se.bwa_se.q10_keep.{macs2_step}.concave150_idr', '.tsv', 'atac_am'), macs2_step=l_macs2_step),
        # Select fixed cutoffs
        expand(pf('atac_am1', 'tg_se.bwa_se.q10_keep.{macs2_step}.concave150_idr', '_0.001.bed', 'atac_am'), macs2_step=l_macs2_step),
        # Coverage tracks, pooled across stages (for browsing)
        expand(pf('atac_am1_{stage}', 'tg_se.bwa_se.q10_keep.{macs2_step}.concave_idr', '.bw', 'atac_am'), macs2_step=l_macs2_step, stage=l_stage_am),

rule atac_am_qc:
    input:
        # Full samples
        expand(pf('{bid}', 'c_r1', '.txt', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se.c_r1', '.txt', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se.bwa_se.c_chrM', '.txt', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep.c', '.txt', 'atac_am'), bid=atac_am_samples()),
        # 1M subsampled for quick QC
        expand(pf('{bid}', 'tg_se_1M.c_r1', '.txt', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se_1M.bwa_se_ecoli.c', '.txt', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se_1M.bwa_se_ecoli.c_q10', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_1M.bwa_se.c_chrM', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_1M.bwa_se.q10_keep.c', '.txt', 'atac_am'), bid=atac_am_samples()),
        # trim_galore with lower threshold -- neglible increase in reads aligned (<1%)
        #expand(pf('{bid}', 'tg_se_q10.c_r1', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_q10.bwa_se_ecoli.c_q10', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_q10.bwa_se.c_chrM', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_q10.bwa_se.q10_keep.c', '.txt', 'atac_am'), bid=atac_am_samples()),
    output:
        'atac_am/atac_am_qc.tsv'
    run:
        df = pd.DataFrame()
        df.index.name = 'sample'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            #df.ix[bid, 'label'] = config['atac_am'][bid]['label']
            df.ix[bid, step] = '%d' % (read_int(pf(bid, step, suffix, prefix)),)
        df.sort_index(axis=0, inplace=True)
        df.to_csv(output[0], sep='\t')

"""
rule atac_am_summary_pooled:
    input:
        # Full samples
        expand(pf('{bid}', 'tg_se.bwa_se.c', '.txt', 'atac_am'), bid=atac_am_samples_pooled()),
        #expand(pf('{bid}', 'tg_se.bwa_se_ecoli.c_q10', '.txt', 'atac_am'), bid=atac_am_samples_pooled()),
        expand(pf('{bid}', 'tg_se.bwa_se.c_chrM', '.txt', 'atac_am'), bid=atac_am_samples_pooled()),
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep.c', '.txt', 'atac_am'), bid=atac_am_samples_pooled()),
        # 1M subsampled for quick QC
        #expand(pf('{bid}', 'tg_se_1M.c_r1', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_1M.bwa_se_ecoli.c_q10', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_1M.bwa_se.c_chrM', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_1M.bwa_se.q10_keep.c', '.txt', 'atac_am'), bid=atac_am_samples()),
        # trim_galore with lower threshold -- neglible increase in reads aligned (<1%)
        #expand(pf('{bid}', 'tg_se_q10.c_r1', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_q10.bwa_se_ecoli.c_q10', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_q10.bwa_se.c_chrM', '.txt', 'atac_am'), bid=atac_am_samples()),
        #expand(pf('{bid}', 'tg_se_q10.bwa_se.q10_keep.c', '.txt', 'atac_am'), bid=atac_am_samples()),
    output:
        'atac_am/atac_am_summary_pooled.tsv'
    run:
        df = pd.DataFrame()
        df.index.name = 'sample'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            #df.ix[bid, 'label'] = config['atac_am'][bid]['label']
            df.ix[bid, step] = '%d' % (read_int(pf(bid, step, suffix, prefix)),)
        df.sort_index(axis=0, inplace=True)
        df.to_csv(output[0], sep='\t')
"""

rule atac_am_728:
    input:
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep.macs2_se_extsize150_shiftm75_keepdup_auto', '_treat_pileup.bw', 'atac_am'), bid=atac_am_samples()),
        expand(pf('{bid}', 'tg_se.bwa_se.q10_keep.macs2_se_extsize150_shiftm75_keepdup_1', '_treat_pileup.bw', 'atac_am'), bid=atac_am_samples()),
