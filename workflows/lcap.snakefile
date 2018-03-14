# Long cap-specific processing

rule htseq_counts_lcap:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.htseq_counts', '_namesorted.bam', '{prefix}')),
        pf('{bid}', '{step}.htseq_counts', '.tsv', '{prefix}'),
    params:
        gff_file = '<(gunzip -c WS260_ce10/WS260_ce10.transcripts.gtf.gz)'
    shell:
        '''
        samtools sort -n {input} > {output[0]}
        htseq-count -f bam -r name -a 10 -s reverse -m intersection-strict --type=exon {output[0]} {params.gff_file} > {output[1]}
        '''
        #'htseq-count -f bam -r pos -a 10 -s reverse -m intersection-strict {input} {params.gff_file} > {output}' # htseq-count alignment buffer overflows with some long cap libraries when coordinate-sorted

rule filled_fwd_lcap:
    '''
    Proper samtools args for mapped paired-end alignments:
    -f 3 => read paired, mapped in proper pair
    -F 12 => NOT read unmapped, NOT mate unmapped
    '''
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.filled_fwd', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.filled_fwd', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -b -u {input} | bedtools genomecov -ibam stdin -bg -pc -strand - > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule filled_rev_lcap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.filled_rev', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.filled_rev', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -b -u {input} | bedtools genomecov -ibam stdin -bg -pc -strand + > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule unfilled_fwd_lcap:
    '''
    Proper samtools args for mapped paired-end alignments:
    -f 3 => read paired, mapped in proper pair
    -f 128 => second in sequencing
    -F 12 => NOT read unmapped, NOT mate unmapped
    '''
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.unfilled_fwd', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.unfilled_fwd', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -b -u {input} | bedtools genomecov -ibam stdin -bg -strand + > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule unfilled_rev_lcap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.unfilled_rev', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.unfilled_rev', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -b -u {input} | bedtools genomecov -ibam stdin -bg -strand - > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule lt200_filled_fwd_lcap:
    '''
    Proper samtools args for mapped paired-end alignments:
    -f 3 => read paired, mapped in proper pair
    -F 12 => NOT read unmapped, NOT mate unmapped
    '''
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.lt200_filled_fwd', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.lt200_filled_fwd', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -f 3 -F 12 -h {input} \
        | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 200*200)' \
        | samtools view -b -u - \
        | bedtools genomecov -ibam stdin -bg -pc -strand - > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule lt200_filled_rev_lcap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.lt200_filled_rev', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.lt200_filled_rev', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -f 3 -F 12 -h {input} \
        | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 200*200)' \
        | samtools view -b -u - \
        | bedtools genomecov -ibam stdin -bg -pc -strand + > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule startbp_fwd_lcap:
    '''
    Proper samtools args for mapped paired-end alignments (+select for second reads to capture 5' ends of the fragments):
    -f 131 => read paired, mapped in proper pair, second in sequencing
    -F 12 => NOT read unmapped, NOT mate unmapped
    Strand of the second-read is the 'correct one' (=does not need to be flipped)
    '''
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.startbp_fwd', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.startbp_fwd', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -b -f 128 -u {input} | bedtools genomecov -ibam stdin -bg -5 -strand + > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule startbp_rev_lcap:
    input: pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.startbp_rev', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.startbp_rev', '.bw', '{prefix}'),
    shell:
        '''
        samtools view -b -f 128 -u {input} | bedtools genomecov -ibam stdin -bg -5 -strand - > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule rafts:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.rafts', '_namesorted.bam', '{prefix}')),
        pf('{bid}', '{step}.rafts', '_raw.bed', '{prefix}'),
    threads: 8
    shell: '''
        samtools sort -n --threads {threads} {input} > {output[0]}
        bedtools bamtobed -bedpe -mate1 -i {output[0]} \
        | awk -F'\\t' -v OFS='\\t' '{{print $1,$2<$5?$2:$5,$3<$6?$6:$3,$7,$8,$10}}' | sort -k 1,1 -k 2,2n -k 3,3n -k6,6 \
        | bedtools merge -s -c 1 -o count -i stdin > {output[1]}
    '''

def lcap_samples():
    #return ['chen13_emb_lcap1_1M', 'chen13_emb_lcap2_1M']
    #return ['HS434_K19_LibA', 'HS434_K19_LibB', 'HS434_K19_LibC']
    #return ['HS482_JA8_lcRNA_S1r', 'HS482_L06_N2_R1_lcRNA', 'HS482_L06_N2_R2_lcRNA', 'HS482_L07_DM_R1_lcRNA', 'HS482_L07_DM_R2_lcRNA', 'HS482_P49_YA_lcRNA', 'HS482_P50_D3_lcRNA', 'HS482_P62_D3_lcRNA', 'HS482_P62_YA_lcRNA',]
        #list(sorted(config['samples_lcap'].keys())) +\
        #['HS327_AM132_1', 'HS327_AM132_2', 'HS327_AM132_3', 'HS327_AM132_4'] +\
        #['HS460_crb55_MSR01_YA_rep2', 'HS460_crb55_MSR02_YA_rep2ca', 'HS460_crb55_MSR03_D10_rep2', 'HS460_crb55_MSR04_D7_rep1', 'HS460_crb55_MSR05_YA_rep1', 'HS460_crb55_MSR06_D7_rep2'] +\
    return ['HS485_JA22_lcRNA_S1r2', 'HS485_JA26_lcRNA_S1r', 'HS485_JA2_lcRNA_S1r2', 'HS485_JA3_lcRNA_S1r2', 'HS485_JA4_lcRNA_S1r2',  'HS485_JA5_lcRNA_S1r2','HS485_JA6_lcRNA_S1r2', 'HS485_JA7_lcRNA_S1r2', 'HS485_JA9_lcRNA_S1r2',]

#def lcap_process_bid_annot():
#    return [config['samples'][bid]['annot'] for bid in lcap_process_bid()]

"""
c_rRNA can be 2x larger than c_rRNA_q10, as rRNA is expected to have a significant
amount of repetitive genomic sequence due to multiple copies...
"""
rule c_r1_RNA:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_r1_rRNA', '.txt', '{prefix}')
    shell: 'samtools view -c -f 64 -L shared/WS253_ce10.rRNA.bed {input} > {output}'

rule c_r1_rRNA_q10:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_r1_rRNA_q10', '.txt', '{prefix}')
    shell: 'samtools view -c -f 64 -q 10 -L shared/WS253_ce10.rRNA.bed {input} > {output}'

rule c_r1_RNA_broad:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_r1_RNA_broad', '.txt', '{prefix}')
    shell: 'samtools view -c -f 64 -L shared/WS253_ce10.rRNA_broad.bed {input} > {output}'

rule c_app_r1: # "aligned in a proper pair" -f 67 = paired, mapped in proper pair, first in pair -F 12 = (not) unmapped (read+mate)
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_app_r1', '.txt', '{prefix}')
    shell: 'samtools view -c -f 67 -F 12 {input} > {output}'

rule outrons:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.outrons', '.bam', '{prefix}'),
    shell:
        'samtools view -b -L shared/dm160615_outrons.bed {input} > {output}'

rule lcap:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'lcap'), bid=config['lcap_raw']),
        expand(pf('{bid}', 'c_r2', '.txt', 'lcap'), bid=config['lcap_raw']),
        #expand(htp('.bid/{bid}_1M.r1.fq.gz'), bid=config['lcap_raw']),
        #expand(htp('.bid/{bid}_1M.r2.fq.gz'), bid=config['lcap_raw']),
        #expand(pf('{bid}', 'trim14.bwa_pe', '.bam', 'lcap'), bid=config['lcap_raw']),
        #expand(pf('{bid}', 'trim20.bwa_pe', '.bam', 'lcap'), bid=config['lcap_raw']),
        #expand(pf('{bid}', 'trim20.bwa_pe', '.bam.bai', 'lcap'), bid=config['lcap_raw']),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap'), bid=config['lcap_raw']),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev', '.bw', 'lcap'), bid=config['lcap_raw']),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.htseq_counts', '.tsv', 'lcap'), bid=config['lcap_HS569']),

rule lcap_nreads_vs_outron_coverage_sample:
    input:
        expand('samples/{bid}_{sub}.r1.fq.gz', bid=['HS352_JA1_lcRNA_S1r', 'HS352_JA7_lcRNA_S1r'], sub=['1M', '5M', '10M', '20M']),
        expand('samples/{bid}_{sub}.r2.fq.gz', bid=['HS352_JA1_lcRNA_S1r', 'HS352_JA7_lcRNA_S1r'], sub=['1M', '5M', '10M', '20M']),

rule lcap_qc:
    input:
        expand(pf('{bid}', '{step}', '.txt', 'lcap'), bid=config['lcap_raw'],
            step=[
                'c_r1', # Total reads
                'trim20.bwa_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1',
                ]),
    output:
        'lcap/lcap_qc_counts.tsv',
        'lcap/lcap_qc_passed.tsv',
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
        df_['raw_reads'] = df['c_r1'].astype(int).map(yp.f_uk)
        df_['mapped'] = keep_pct_(
                'trim20.bwa_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.c_r1',
        )
        df_['not_mitochondrial'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
        )
        df_['not_rRNA'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1',
        )
        df_['not_blacklisted'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1',
        )
        df_['mapq10'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1',
        )
        df_['useful_reads'] = df['trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1'].astype(int).map(yp.f_uk)
        df_.sort_index().to_csv(output[1], sep='\t')

rule lcap_nreads_vs_outron_coverage:
    input:
        'lcap/lcap_qc.tsv'
    output:
        'lcap/lcap_nreads_vs_outron_coverage.pdf'
    run:
        dm_step = 'dm160615'
        # Load unstranded promoter annotation
        fp_prom = '../HTSProcessing/%(dm_step)s/hs_sites.promoters_Dec14_unstranded.tsv' % locals()
        df_prom = pd.read_csv(fp_prom, sep='\t')
        #print df_prom['is_promoter_fwd'].value_counts()
        #print df_prom['is_promoter_rev'].value_counts()
        print(len(df_prom.query('is_promoter_fwd & is_promoter_rev')),\
            len(df_prom.query('is_promoter_fwd & ~is_promoter_rev')),\
            len(df_prom.query('~is_promoter_fwd & is_promoter_rev')), 'bi/fwd/rev promoters')

        def frac_fwd(bid, dm_step = 'dm160615'):
            print('frac_fwd', bid)
            fp_ = pf(bid, 'trim20.bwa_pe.rm_rRNA.p10_keep.filled_fwd', '.bw', 'lcap')
            fp_out = '../HTSProcessing/%(dm_step)s/hs_sites.promoters_Dec14' % locals()
            q_ = 'prom_summary == "high_confidence_promoter" | prom_summary == "low_confidence_promoter"'
            df_prom_fwd = pd.read_csv(fp_out + '_fwd.tsv', sep='\t').query(q_).query('exon1_start - end > 200')
            #df_prom_rev = pd.read_csv(fp_out + '_rev.tsv', sep='\t').query(q_).query('exon1_start - end > 200')
            #print len(df_prom_fwd), len(df_prom_rev)
            chroms = df_prom_fwd['chrom'].tolist()
            starts = df_prom_fwd['end'].tolist()
            ends = df_prom_fwd['exon1_start'].tolist()
            df_prom_fwd['outron_min'] = list(
                read_regions(fp_, chroms, starts, ends, np.min))
            return sum(df_prom_fwd['outron_min'] >= 1) / float(len(df_prom_fwd))

        df_qc = pd.read_csv(input[0], sep='\t')#.head(4)
        #print(df_qc)
        df_qc['frac_fwd'] = list(map(frac_fwd, df_qc['sample'].tolist()))
        df_qc_ = df_qc.copy()#.query('colour == "r" | colour == "g"')
        y_col = 'trim20.bwa_pe.rm_rRNA_broad.p10_keep.c_r1'
        plt.figure(figsize=(12,8))
        plt.scatter(df_qc_[y_col], df_qc_['frac_fwd'])#, color=df_qc_['colour'])
        for i, r in df_qc_.iterrows():
            plt.gca().annotate(r['label'], (r[y_col], r['frac_fwd']), size=4)
        #plt.gca().set_ylim(0.0,1)
        plt.gca().set_xlabel('Useful reads')
        plt.gca().set_ylabel('Fwd promoters with continuous signal')
        plt.savefig(output[0], bbox_inches='tight')

"""
rule lcap728:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'c_r2', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_fwd', '.bw', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_fwd_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_rev', '.bw', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_rev_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10', '.bam.bai', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.htseq_counts', '.tsv', 'lcap728'), bid=lcap728_samples),
        pf('lcap728', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.sizefactors_lcap728', '_counts.tsv', 'lcap728'),
        pf('lcap728', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.sizefactors_lcap728', '_sizefactors.tsv', 'lcap728'),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev', '.bw', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_fwd', '.bw', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_rev', '.bw', 'lcap728'), bid=lcap728_samples),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.rafts', '_raw.bed', 'lcap728'), bid=lcap728_samples),

rule lcap728_qc:
    input:
        expand(pf('{bid}', '{step}', '.txt', 'lcap728'), bid=l_lcap_pooled,
            step=[
                'c_r1', # Total reads
                'trim20.bwa_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1',
                ]),
    output:
        'lcap728/lcap_qc_counts.tsv',
        'lcap728/lcap_qc_passed.tsv',
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
                'trim20.bwa_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.c_r1',
        )
        df_['not_mitochondrial'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
        )
        df_['not_rRNA'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1',
        )
        df_['not_blacklisted'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1',
        )
        df_['mapq10'] = keep_pct_(
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1',
                'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1',
        )
        df_['useful_reads'] = df['trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1'].astype(int)
        df_.sort_index().to_csv(output[1], sep='\t')
"""

rule lcap805:
    input:
        expand(pf('lcap805_{bid}', 'c_r1', '.txt', 'lcap805'), bid=config['stages_rep']),
        expand(pf('lcap805_{bid}', 'c_r2', '.txt', 'lcap805'), bid=config['stages_rep']),
        # Coverage tracks, q10
        expand(pf('lcap805_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap805'), bid=config['stages_rep']),
        expand(pf('lcap805_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev', '.bw', 'lcap805'), bid=config['stages_rep']),
        # WS260_ce10 exon counts
        expand(pf('lcap805_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.htseq_counts', '.tsv', 'lcap805'), bid=config['stages_rep']),
        # Normalised coverage, q10
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        ######
        #expand(pf('{bid}', 'trim20.bwa_pe.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_fwd', '.bw', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_fwd_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_rev', '.bw', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_rev_sfnorm', '.bw', 'lcap728'), bid=l_lcap_pooled),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.c_r1', '.txt', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10', '.bam.bai', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.htseq_counts', '.tsv', 'lcap728'), bid=lcap728_samples),
        #pf('lcap728', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.sizefactors_lcap728', '_counts.tsv', 'lcap728'),
        #pf('lcap728', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.sizefactors_lcap728', '_sizefactors.tsv', 'lcap728'),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_fwd', '.bw', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_rev', '.bw', 'lcap728'), bid=lcap728_samples),
        #expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.rafts', '_raw.bed', 'lcap728'), bid=lcap728_samples),

rule lcap808_read1: #https://bitbucket.org/snakemake/snakemake/issues/397/unable-to-set-utime-on-symlink-your-python
    input:
        'samples/lcap808_{sample}.r1.fq.gz'
    output:
        'lcap808_geo/reads/lcap_{sample}.read1.fastq.gz'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule lcap808_read2: #https://bitbucket.org/snakemake/snakemake/issues/397/unable-to-set-utime-on-symlink-your-python
    input:
        'samples/lcap808_{sample}.r2.fq.gz'
    output:
        'lcap808_geo/reads/lcap_{sample}.read2.fastq.gz'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule sizefactors_lcap808_counts:
    input:
        expand(pf('lcap808_{sample}', '{{step}}.htseq_counts', '.tsv', 'lcap808'), sample=config['stages_rep']),
    output:
        pf('lcap808', '{step}.sizefactors_lcap808', '_counts.tsv', 'lcap808'),
    run:
        def counts_(fp): return pd.read_csv(fp, sep='\t', names=('gene_id', 'counts'))['counts'].tolist()
        df = pd.DataFrame(collections.OrderedDict([(parse_pf(input_)[0], counts_(input_)) for input_ in input]),
            index=pd.read_csv(input[0], sep='\t', names=('gene_id', 'counts'))['gene_id'].tolist())\
            .drop(['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique'])\
            .to_csv(output[0], sep='\t')

rule sizefactors_lcap808_deseq:
    input:
        pf('lcap808', '{step}.sizefactors_lcap808', '_counts.tsv', 'lcap808')
    output:
        pf('lcap808', '{step}.sizefactors_lcap808', '_sizefactors.tsv', 'lcap808')
    script:
        'sizefactors_lcap808_deseq.R'

rule filled_fwd_lcap808_sfnorm:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('lcap808', '{step}.sizefactors_lcap808', '_sizefactors.tsv', 'lcap808'),
    output:
        temp(pf('{bid}', '{step}.filled_fwd_sfnorm', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.filled_fwd_sfnorm', '.bw', '{prefix}'),
    shell:
        '''
        sample=$(basename "{input[0]}" | cut -f1 -d".")
        sf=$(cat {input[1]} | awk -F'\\t' -v OFS='\\t' -v sample=$sample '($1 == sample) {{print 1/$2}}')
        samtools view -b -u {input[0]} | bedtools genomecov -ibam stdin -bg -pc -strand - -scale $sf > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule filled_rev_lcap808_sfnorm:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('lcap808', '{step}.sizefactors_lcap808', '_sizefactors.tsv', 'lcap808'),
    output:
        temp(pf('{bid}', '{step}.filled_rev_sfnorm', '.bedGraph', '{prefix}')),
        pf('{bid}', '{step}.filled_rev_sfnorm', '.bw', '{prefix}'),
    shell:
        '''
        sample=$(basename "{input[0]}" | cut -f1 -d".")
        sf=$(cat {input[1]} | awk -F'\\t' -v OFS='\\t' -v sample=$sample '($1 == sample) {{print 1/$2}}')
        samtools view -b -u {input[0]} | bedtools genomecov -ibam stdin -bg -pc -strand + -scale $sf > {output[0]}
        sort -k1,1 -k2,2n {output[0]} -o {output[0]}
        bedGraphToBigWig {output[0]} <(samtools view -H {input} | awk 'substr($2,1,2)=="SN" {{print substr($2, 4),substr($3, 4)}}') {output[1]}
        '''

rule mean_by_stage:
    input:
        pf('{bid}_rep1', '{step}', '.bw', '{prefix}'),
        pf('{bid}_rep2', '{step}', '.bw', '{prefix}'),
    output:
        pf('{bid}', '{step}.mean_by_stage', '.bw', '{prefix}'),
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} mean {input[0]} {input[1]}
    '''

rule lcap808:
    input:
        # raw read counts
        expand(pf('lcap808_{bid}', 'c_r1', '.txt', 'lcap808'), bid=config['stages_rep']),
        expand(pf('lcap808_{bid}', 'c_r2', '.txt', 'lcap808'), bid=config['stages_rep']),
        # Coverage tracks, q10
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap808'), bid=config['stages_rep']),
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev', '.bw', 'lcap808'), bid=config['stages_rep']),
        # Coverage tracks, q10, mean by stage
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd.mean_by_stage', '.bw', 'lcap808'), bid=config['stages']),
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev.mean_by_stage', '.bw', 'lcap808'), bid=config['stages']),
        # WS260_ce10 exon counts
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.htseq_counts', '.tsv', 'lcap808'), bid=config['stages_rep']),
        # Startbp tracks for jump/incr tests
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_fwd', '.bw', 'lcap808'), bid=config['stages_rep']),
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_rev', '.bw', 'lcap808'), bid=config['stages_rep']),
        # raw reads & final tracks; final "geo" names
        expand('lcap808_geo/reads/lcap_{sample}.read1.fastq.gz', sample=list(config['lcap808'].keys())),
        expand('lcap808_geo/reads/lcap_{sample}.read2.fastq.gz', sample=list(config['lcap808'].keys())),
        # calculate sizeFactors from gene-level read counts using DESeq2
        pf('lcap808', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.sizefactors_lcap808', '_counts.tsv', 'lcap808'),
        pf('lcap808', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.sizefactors_lcap808', '_sizefactors.tsv', 'lcap808'),
        # normalise tracks by sizeFactors
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd_sfnorm', '.bw', 'lcap808'), bid=config['stages_rep']),
        expand(pf('lcap808_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap808'), bid=config['stages_rep']),

