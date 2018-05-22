
rule alignment_raw_reads_by_rep:
    input:
        pf('atac728_{sample}', 'tg_se.bwa_se', '.bam', 'atac728'),
    output:
        pf('atac_{sample}', 'alignment_raw_reads_by_rep', '.bam', 'atac_geo'),
    shell:
        'ln -s `pwd`/{input} `pwd`/{output}'

rule alignment_q10_reads_by_rep:
    input:
        pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10', '.bam', 'atac728'),
    output:
        pf('atac_{sample}', 'alignment_q10_reads_by_rep', '.bam', 'atac_geo'),
    shell:
        'ln -s `pwd`/{input} `pwd`/{output}'

rule coverage_spmr_all_reads_by_rep:
    input:
        pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac728'),
    output:
        pf('atac_{sample}', 'coverage_spmr_all_reads_by_rep', '.bw', 'atac_geo'),
    shell:
        'ln -s `pwd`/{input} `pwd`/{output}'

rule coverage_spmr_q10_reads_by_rep:
    input:
        pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac728'),
    output:
        pf('atac_{sample}', 'coverage_spmr_q10_reads_by_rep', '.bw', 'atac_geo'),
    shell:
        'ln -s `pwd`/{input} `pwd`/{output}'

rule coverage_spmr_all_reads_by_stage:
    input:
        pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_stage', '.bw', 'atac728'),
    output:
        pf('atac_{sample}', 'coverage_spmr_all_reads_by_stage', '.bw', 'atac_geo'),
    shell:
        'ln -s `pwd`/{input} `pwd`/{output}'

rule coverage_spmr_q10_reads_by_stage:
    input:
        pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_stage', '.bw', 'atac728'),
    output:
        pf('atac_{sample}', 'coverage_spmr_q10_reads_by_stage', '.bw', 'atac_geo'),
    shell:
        'ln -s `pwd`/{input} `pwd`/{output}'

"""
rule atac728_qc:
    input:
        expand(pf('{sample}', '{step}', '.txt', 'atac728'), sample=config['stages_atac_by_rep'],
            step=[
                'c_r1', # Total reads
                'tg_se.bwa_se.c',
                'tg_se.bwa_se.rm_unmapped.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.c',
                'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c',
                ]),
    output:
        'atac728/atac_qc_counts.tsv',
        'atac728/atac_qc_passed.tsv',
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

rule sfnorm_atac728_wt_glp1:
    input:
        pf('atac728_wt_glp1', '{step}.sfnorm', '_counts.tsv', 'atac728'), # computed in atac_sizefactor.ipynb so not in snakemake
    output:
        pf('atac728_wt_glp1', '{step}.sfnorm', '_sizefactors.tsv', 'atac728'),
    script:
        'sizefactors_atac728_wt_glp1.R'

rule sfnorm_atac728_wt_glp1_scale:
    input:
        pf('atac728_{sample}', '{step}', '_treat_pileup.bw', 'atac728'),
        pf('atac728_wt_glp1', '{step}.sfnorm', '_sizefactors.tsv', 'atac728'),
    output:
        temp(pf('atac728_{sample}', '{step}.sfnorm', '.chroms', 'atac728')),
        temp(pf('atac728_{sample}', '{step}.sfnorm', '.bedGraph', 'atac728')),
        pf('atac728_{sample}', '{step}.sfnorm', '.bw', 'atac728'),
    shell:
        '''
        sample="{wildcards.sample}"
        scalefactor=$(cat {input[1]} | awk -F'\\t' -v OFS='\\t' -v sample=$sample '($1 == sample) {{print $2}}')
        bigWigInfo -chroms {input[0]} | grep '^[[:space:]]' | awk -v OFS='\\t' '{{print $1,$3}}' > {output[0]}
        sort -k1,1 -o {output[0]} {output[0]}
        wiggletools write_bg {output[1]} scale $scalefactor {input[0]}
        sort -k1,1 -k2,2n -o {output[1]} {output[1]}
        bedGraphToBigWig {output[1]} {output[0]} {output[2]}
        '''

rule atac728:
    input:
        expand(pf('{sample}', 'c_r1', '.txt', 'atac728'), sample=config['atac728'].keys()),
        expand(pf('{sample}', 'c_r2', '.txt', 'atac728'), sample=config['atac728_pe'].keys()),
        expand(pf('{sample}', 'tg_pe.bwa_pe', '.bam', 'atac728'), sample=config['atac728_pe'].keys()),
        expand(pf('{sample}', 'tg_se.bwa_se', '.bam', 'atac728'), sample=config['atac728'].keys()),
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.fsizes', '.txt', 'atac728'), sample=config['atac728_pe'].keys()),
        expand(pf('{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac728'), sample=config['atac728_pe'].keys()),
        expand(pf('{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac728'), sample=config['atac728'].keys()),
        expand(pf('{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_treat_pileup.bw', 'atac728'), sample=config['atac728'].keys()),
        expand(pf('{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac728'), sample=config['atac728'].keys()),
        #pf('atac728_wt', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR.sfnorm', '_sizefactors.tsv', 'atac728'),
        expand(pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR.sfnorm', '.bw', 'atac728'), sample=config['stages_wt_by_rep'] + config['stages_glp1_by_rep']),
        #expand(pf('atac728_{sample}_mean', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR.sfnorm.mean_by_stage', '.bw', 'atac728'), sample=config['stages_wt']),
        #expand(pf('{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize1_keepdup_all', '_treat_pileup.bw', 'atac728'), sample=config['atac728'].keys()),
        #expand(pf('atac728_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_stage', '.bw', 'atac728'), stage=atac728_stage),
        #expand(pf('atac728_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_stage', '.bw', 'atac728'), stage=atac728_stage),
        #expand(pf('atac728_{sample}', 'coverage_spmr_all_reads_by_rep', '.bw', 'atac728_geo'), sample=config['stages_atac_by_rep']),
        #expand(pf('atac728_{sample}', 'coverage_spmr_q10_reads_by_rep', '.bw', 'atac728_geo'), sample=config['stages_atac_by_rep']),
        #expand(pf('atac728_{sample}', 'coverage_spmr_all_reads_by_stage', '.bw', 'atac728_geo'), sample=config['stages_atac']),
        #expand(pf('atac728_{sample}', 'coverage_spmr_q10_reads_by_stage', '.bw', 'atac728_geo'), sample=config['stages_atac']),
        #expand(pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize1_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac728'), sample=config['stages_atac_by_rep']),
        #expand(pf('atac728_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize1_keepdup_all', '_treat_pileup.bw', 'atac728'), sample=config['stages_atac_by_rep']),
        #tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.c',
        # Trim / align / generate coverage (all samples)
        #expand('samples/{bid}.r1.fq.gz'), bid=atac_samples()),
        #expand('samples/{bid}.r2.fq.gz'), bid=atac_samples()),
        #expand(pf('{bid}', 'fastq_r1_count', '.txt', 'atac'), bid=atac_samples()),
        #expand(pf('{bid}', 'fastq_r2_count', '.txt', 'atac'), bid=atac_samples()),
        #expand(pf('{bid}', 'tg_pe.bwa_pe.p10_keep', '.bam', 'atac'), bid=atac_samples()),
        #expand(pf('{bid}', 'tg_pe.bwa_pe.p10_keep.c_r1', '.txt', 'atac'), bid=atac_samples()),
        #expand(pf('{bid}', 'tg_pe.bwa_pe.p10_keep.c_r2', '.txt', 'atac'), bid=atac_samples()),
        #expand(pf('{bid}', 'tg_pe.bwa_pe.p10_keep.{macs2_step}', '_treat_pileup.bw', 'atac'), bid=atac_samples(), macs2_step=l_macs2_step),
        # Generate replicates & pseudoreplicates (dm)
        #expand(pf('atac_{stage}_pep1', 'tg_pe.bwa_pe.p10_keep', '.bam', 'atac'), stage=l_stage),
        #expand(pf('atac_{stage}_pep2', 'tg_pe.bwa_pe.p10_keep', '.bam', 'atac'), stage=l_stage),
        #expand(pf('atac_{stage}_rep1', 'tg_pe.bwa_pe.p10_keep.c', '.txt', 'atac'), stage=l_stage),
        #expand(pf('atac_{stage}_rep2', 'tg_pe.bwa_pe.p10_keep.c', '.txt', 'atac'), stage=l_stage),
        #expand(pf('atac_{stage}_pep1', 'tg_pe.bwa_pe.p10_keep.c', '.txt', 'atac'), stage=l_stage),
        #expand(pf('atac_{stage}_pep2', 'tg_pe.bwa_pe.p10_keep.c', '.txt', 'atac'), stage=l_stage),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr', '.tsv', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr_raw', '_mean10.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr_raw', '_mean25.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr_raw', '_mean50.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr_raw', '_mean100.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr_raw', '_mean150.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.dm_idr_raw', '_D2_s3x150_1E6.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('atac_{stage}_rep1,atac_{stage}_rep2', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.concave_idr', '_all.bed', 'atac'), stage=l_stage, macs2_step=l_macs2_step),
        #expand(pf('dm_sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.concave_idr.dm_idr', '_mean.bw', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('dm_sites', 'tg_pe.bwa_pe.p10_keep.{macs2_step}.concave_idr.dm_idr', '_emb.bed', 'atac'), macs2_step=l_macs2_step),
        #expand(pf('{bid}', 'tg_pe.bwa_pe.p10_keep.macs2_se_extsize150_shiftm75_keepdup_auto', '_treat_pileup.bw', 'atac_am'), bid=atac_samples()),

rule atac728_geo:
    input:
        expand(pf('atac_{sample}', 'alignment_q10_reads_by_rep', '.bam', 'atac_geo'), sample=config['stages_atac_by_rep']),
        expand(pf('atac_{sample}', 'alignment_raw_reads_by_rep', '.bam', 'atac_geo'), sample=config['stages_atac_by_rep']),
        expand(pf('atac_{sample}', 'alignment_q10_reads_by_rep', '.bam.bai', 'atac_geo'), sample=config['stages_atac_by_rep']),
        expand(pf('atac_{sample}', 'alignment_raw_reads_by_rep', '.bam.bai', 'atac_geo'), sample=config['stages_atac_by_rep']),

rule atac749:
    input:
        expand('atac749/reads/atac_{sample}.read1.fastq.gz', sample=config['atac749_pe']),
        expand('atac749/reads/atac_{sample}.read2.fastq.gz', sample=config['atac749_pe']),
        expand('atac749/reads/atac_{sample}.fastq.gz', sample=config['atac749_se']),
        expand('atac749/tracks/atac_{sample}.bw', sample=config['stages_wt'] + config['stages_glp1']),
"""

def df_atac_make():
    df_atac = pd.read_csv('annot/S1_accessible_sites/S1a_accessible_sites.tsv' % locals(), sep='\t')
    return df_atac[yp.NAMES_BED3]

rule atac808_nanmax:
    input:
        pf('atac808_{sample}', '{step}', '_treat_pileup.bw', 'atac808'),
    output:
        pf('atac808_{sample}', '{step}.atac_nanmax', '.tsv', 'atac808'),
        pf('atac808_{sample}', '{step}.atac_nanmax', '.bed', 'atac808'),
    run:
        col_ = 'atac_%s_nanmax' % (wildcards.sample,)
        df_sites = df_atac_make()
        df_sites[col_] = list(map(lambda c: int(np.nanmax(c)), yp.read_regions(input[0], df_sites.chrom.tolist(), df_sites.start.tolist(), df_sites.end.tolist())))
        df_sites[[col_]].to_csv(output[0], sep='\t', index=False)
        df_sites[yp.NAMES_BED3 + [col_]].to_csv(output[1], sep='\t', index=False, header=False)

rule atac814_mean_by_stage_treat_pileup:
    input:
        pf('atac814_{stage}_rep1', '{step}', '_treat_pileup.bw', 'atac814'),
        pf('atac814_{stage}_rep2', '{step}', '_treat_pileup.bw', 'atac814'),
    output:
        pf('atac814_{stage}', '{step}.mean_by_stage', '_treat_pileup.bw', 'atac814'),
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} mean {input[0]} {input[1]}
    '''

rule atac808_mean_by_series_wt:
    input:
        expand(pf('atac808_{sample}', '{{step}}', '_treat_pileup.bw', 'atac808'), sample=config['stages_wt_rep']),
    output:
        pf('atac808_wt', '{step}.mean_by_series', '_treat_pileup.bw', 'atac808'),
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} mean {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]}
    '''

rule atac808_mean_by_series_glp1:
    input:
        expand(pf('atac808_{sample}', '{{step}}', '_treat_pileup.bw', 'atac808'), sample=config['stages_glp1_rep']),
    output:
        pf('atac808_glp1', '{step}.mean_by_series', '_treat_pileup.bw', 'atac808'),
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} mean {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]}
    '''

rule atac814_read1: #https://bitbucket.org/snakemake/snakemake/issues/397/unable-to-set-utime-on-symlink-your-python
    input:
        'samples/atac814_{sample}.r1.fq.gz'
    output:
        'atac814_geo/reads/atac_{sample}.read1.fastq.gz'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule atac814_read2: #https://bitbucket.org/snakemake/snakemake/issues/397/unable-to-set-utime-on-symlink-your-python
    input:
        'samples/atac814_{sample}.r2.fq.gz'
    output:
        'atac814_geo/reads/atac_{sample}.read2.fastq.gz'
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule atac814_tracks:
    input:
        pf('atac814_{sample}_rep1', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'),
        pf('atac814_{sample}_rep2', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'),
    output:
        'atac814_geo/tracks/atac_{sample}.bw',
    shell: '''
        scripts/bigWiggleTools.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}
        '''

rule atac814_alignments:
    input:
        pf('atac814_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10', '.bam', 'atac814'),
    output:
        'atac814_geo/alignments/atac_{sample}.bam',
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule atac814:
    input:
        # techreps separately:
        #expand(pf('atac808_{sample}', 'c_r1', '.txt', 'atac808'), sample=config['atac808'].keys()),
        #expand(pf('atac808_{sample}', 'c_r2', '.txt', 'atac808'), sample=config['atac808_pe'].keys()),
        #expand(pf('atac808_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac808'), sample=config['atac808'].keys()),
        #expand(pf('atac808_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac808'), sample=config['atac808_pe'].keys()),
        # techreps pooled by biological sample:
        #expand(pf('atac808_{sample}', 'c_r1', '.txt', 'atac808'), sample=config['stages_rep']),
        #expand(pf('atac808_{sample}', 'c_r2', '.txt', 'atac808'), sample=config['stages_wt_rep']),
        #expand(pf('atac808_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac808'), sample=config['stages_rep']),
        #expand(pf('atac808_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_treat_pileup.bw', 'atac808'), sample=config['stages_rep']),
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac814'), sample=config['stages_wt_rep']),
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac814'), sample=config['stages_wt_rep']),
        #expand(pf('atac808_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200.mean_by_stage', '_treat_pileup.bw', 'atac808'), sample=config['stages_wt']),
        expand(pf('atac814_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'), sample=techreps_collapse(config['atac814'].keys())),
        expand(pf('atac814_{stage}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_stage', '_treat_pileup.bw', 'atac814'), stage=config['stages']),
        # raw sample-based peak heights for differential accessibility tests
        #expand(pf('atac808_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR.atac_nanmax', '.tsv', 'atac808'), sample=config['stages_rep']),
        # normalised stage-specific peak heights for clustering & other downstream analyses
        #expand(pf('atac808_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_stage.atac_nanmax', '.tsv', 'atac808'), sample=config['stages']),
        # pooled by series
        #pf('atac808_wt', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200.mean_by_series', '_treat_pileup.bw', 'atac808'),
        #pf('atac808_wt', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_series', '_treat_pileup.bw', 'atac808'),
        #pf('atac808_glp1', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_series', '_treat_pileup.bw', 'atac808'),
        # GEO submission -- reads
        expand('atac814_geo/reads/atac_{sample}.read1.fastq.gz', sample=config['atac814'].keys()),
        expand('atac814_geo/reads/atac_{sample}.read2.fastq.gz', sample=config['atac814_pe'].keys()),
        # GEO submission -- coverage tracks
        expand('atac814_geo/tracks/atac_{sample}.bw', sample=config['stages']),
        # alignments (to match GEO submission)
        expand('atac814_geo/alignments/atac_{sample}.bam', sample=config['stages_rep']),
        # GEO -- md5 checksums of fastq files
        expand(pf('atac814_{sample}', 'md5sum_r1', '.txt', 'atac814'), sample=config['atac814'].keys()),
        expand(pf('atac814_{sample}', 'md5sum_r2', '.txt', 'atac814'), sample=config['atac814_pe'].keys()),
        # GEO -- read lengths
        expand(pf('atac814_{sample}', 'readlen_r1', '.txt', 'atac814'), sample=config['atac814'].keys()),
        expand(pf('atac814_{sample}', 'readlen_r2', '.txt', 'atac814'), sample=config['atac814_pe'].keys()),
        # GEO -- fragment sizes
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.fsizes', '.txt', 'atac814'), sample=config['atac814_pe'].keys()),

rule atac:
    input:
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac'), bid=[*config['atac']]),
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac'), bid=[*config['atac']]),
