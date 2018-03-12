
def lcap_am_samples():
    #fp = 'lcap_am/lcap_am_sample_sheet.tsv'
    #df = pd.read_csv(fp, sep='\t', comment='#')
    #return df['sample_name'].tolist()
    return []

rule lcap_am_filled_pooled:
    input:
        pf('lcap_am1_{stage}_rep1', '{step}', '.bw', 'lcap_am'),
        pf('lcap_am1_{stage}_rep2', '{step}', '.bw', 'lcap_am'),
    output:
        temp(pf('lcap_am1_{stage}', '{step}.pooled', '.bdg', 'lcap_am')),
        pf('lcap_am1_{stage}', '{step}.pooled', '.bw', 'lcap_am'),
    shell: '''
        wiggletools write_bg {output[0]} sum {input[0]} {input[1]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
    '''

rule lcap_am_rafts_pooled:
    # Raw rafts of continuous read coverage; pooled across biological replicates; annotated with raw read counts
    input:
        pf('lcap_am1_{stage}_rep1', '{step}', '.bam', 'lcap_am'),
        pf('lcap_am1_{stage}_rep2', '{step}', '.bam', 'lcap_am'),
    output:
        temp(pf('lcap_am1_{stage}', '{step}.rafts_pooled', '_namesorted.bam', 'lcap_am')),
        pf('lcap_am1_{stage}', '{step}.rafts_pooled', '_raw.bed', 'lcap_am'),
    threads: 8
    shell: '''
        samtools sort -n --threads {threads} <(samtools cat {input[0]} {input[1]}) > {output[0]}
        bedtools bamtobed -bedpe -mate1 -i {output[0]} \
        | awk -F'\\t' -v OFS='\\t' '{{print $1,$2<$5?$2:$5,$3<$6?$6:$3,$7,$8,$10}}' | sort -k 1,1 -k 2,2n -k 3,3n -k6,6 \
        | bedtools merge -s -c 1 -o count -i stdin > {output[1]}
    '''

rule lcap_am_rafts_pooled_filtered:
    # Raw rafts of continuous read coverage; pooled across biological replicates; annotated with raw read counts
    input:
        pf('{bid}', '{step}.rafts_pooled', '_raw.bed', 'lcap_am'),
    output:
        pf('{bid}', '{step}.rafts_pooled', '_fwd.bed', 'lcap_am'),
        pf('{bid}', '{step}.rafts_pooled', '_rev.bed', 'lcap_am'),
    shell: '''
        awk -F'\\t' -v OFS='\\t' '$4=="+"' {input} \
        | bedtools merge -d 100 -c 5 -o sum \
        | awk -F'\\t' -v OFS='\\t' 'int($4)>1 {{print $1,$2,$3,"raft_" $1 "_" $2 "_" $3,$4,"+"}}' \
        > {output[0]}
        awk -F'\\t' -v OFS='\\t' '$4=="-"' {input} \
        | bedtools merge -d 100 -c 5 -o sum \
        | awk -F'\\t' -v OFS='\\t' 'int($4)>1 {{print $1,$2,$3,"raft_" $1 "_" $2 "_" $3,$4,"-"}}' \
        > {output[1]}
    '''

rule lcap_am:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'c_r2', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe', '.bam.bai', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep', '.bam.bai', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.filled_fwd', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.filled_rev', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.unfilled_fwd', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.unfilled_rev', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.lt200_filled_fwd', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.lt200_filled_rev', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.htseq_exon_counts', '.tsv', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.startbp_fwd', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.startbp_rev', '.bw', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('lcap_am1_{stage}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.filled_fwd.pooled', '.bw', 'lcap_am'), stage=l_stage_am),
        expand(pf('lcap_am1_{stage}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.filled_rev.pooled', '.bw', 'lcap_am'), stage=l_stage_am),
        expand(pf('lcap_am1_{stage}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.rafts_pooled', '_raw.bed', 'lcap_am'), stage=l_stage_am),
        expand(pf('lcap_am1_{stage}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.rafts_pooled', '_fwd.bed', 'lcap_am'), stage=l_stage_am),
        expand(pf('lcap_am1_{stage}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.rafts_pooled', '_rev.bed', 'lcap_am'), stage=l_stage_am),

rule lcap_am_qc:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'c_r2', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.c_r1', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.c_r1', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_r1', '.txt', 'lcap_am'), bid=lcap_am_samples()),
        # Full samplesc_aln_r1
        #expand(pf('{bid}', 'tg_pe_1M', '.r1.fq.gz', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'tg_pe_1M', '.r2.fq.gz', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'tg_pe', '.r1.fq.gz', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'tg_pe', '.r2.fq.gz', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'tg_pe_1M', '.r1.fq.gz', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.c_aln_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.c_app_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA.c_aln_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.c_aln_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_fr', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_fr_rmdup', '.txt', 'lcap'), bid=lcap_samples()),
        # Complexity at outron regions
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.outrons.c_fr', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.outrons.c_fr_rmdup', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20_1M.bwa_pe_ecoli.c_aln_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.c_r1_rRNA', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.c_r1_rRNA_q10', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20_1M.bwa_se_ecoli.c_q10', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20_1M.bwa_pe_ecoli.c_aln_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20_1M.bwa_pe_ecoli.c_p10', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.c_aln_r1', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.p10_keep.c', '.txt', 'lcap'), bid=lcap_samples()),
        #expand(pf('{bid}', 'trim20.bwa_pe.c_chrM', '.txt', 'lcap'), bid=lcap_samples()),
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
        'lcap_am/lcap_am_qc.tsv'
    run:
        df = pd.DataFrame()
        df.index.name = 'sample'
        for (bid, step, suffix, prefix) in map(parse_pf, input):
            try:
                df.ix[bid, 'label'] = config['samples_lcap'][bid]['dm_label']
            except KeyError:
                df.ix[bid, 'label'] = bid
            except TypeError:
                df.ix[bid, 'label'] = bid
            df.ix[bid, step] = '%d' % (read_int(pf(bid, step, suffix, prefix)),)
        #df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.complexity'] = (df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_fr_rmdup'].astype(int) / df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_fr'].astype(int)).map('{:,.02f}'.format)
        #df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.outrons.complexity'] = (df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.outrons.c_fr_rmdup'].astype(int) / df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.outrons.c_fr'].astype(int)).map('{:,.02f}'.format)
        df.sort_index(axis=0, inplace=True)
        df.sort_values('label', inplace=True)
        df['raw_reads'] = df['c_r1']
        df['useable_reads'] = df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_r1']
        #df['useable_fraction'] = (df['trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_r1'].astype(int) / df['c_r1'].astype(int)).map('{:,.02f}'.format)
        #df['needed_useable'] = 20000000 - df['useable_reads'].astype(int)
        #df['needed_raw'] = (df['needed_useable'].astype(float) / df['useable_fraction'].astype(float)).astype(int)
        df.to_csv(output[0], sep='\t')
