
rule atac_bin10_ce10:
    input:
        pf('{dataset}', '{step}', '_treat_pileup.bw', 'atac'),
    output:
        pf('{dataset}', '{step}.bin10_ce10', '.bw', 'atac'),
    shell:
        'scripts/bigWiggleTools_ce10.ipy write {output} scale 0.1 bin 10 {input}'

rule atac_bin10:
    input:
        pf('{dataset}', '{step}', '_treat_pileup.bw', 'atac'),
    output:
        pf('{dataset}', '{step}.bin10', '.bw', 'atac'),
    shell:
        'scripts/bigWiggleTools.ipy write {output} scale 0.1 bin 10 {input}'

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

rule atac824_mean_by_condition:
    input:
        pf('atac824_{condition}_rep1', '{step}', '_treat_pileup.bw', 'atac824'),
        pf('atac824_{condition}_rep2', '{step}', '_treat_pileup.bw', 'atac824'),
    output:
        pf('atac824_{condition}', '{step}.mean_by_condition', '.bw', 'atac824'),
    shell: '''
        scripts/bigWiggleTools.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}
    '''

rule atac814_tracks_ce11:
    input:
        pf('atac814_{sample}_rep1', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'),
        pf('atac814_{sample}_rep2', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'),
    output:
        'atac814_geo/tracks_ce11/atac_{sample}_ce11.bw',
    shell: '''
        scripts/bigWiggleTools_ce11.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}
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

rule atac814_alignments_ce11:
    input:
        pf('atac814_{sample}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_q10', '.bam', 'atac814'),
    output:
        'atac814_geo/alignments_ce11/atac_{sample}_ce11.bam',
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule macs2_daugherty2017:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_daugherty2017', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_daugherty2017', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_daugherty2017', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_daugherty2017', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_daugherty2017', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_daugherty2017', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_daugherty2017', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        # single-base reads were shifted 75bp 5’ to mimic read distributions of a 150 bp fragment of ChIP-seq
        # The following settings were used for MACS: -g 9e7, -q 5e-2, --nomodel, --extsize 150, -B, --keep-dup all, and --call-summits
        macs2_args = '--format BAM --shift -75    --gsize ce -q 5e-2  --nomodel --extsize 150 --bdg --keep-dup all     --call-summits  --SPMR '
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_daugherty2017 --name {wildcards.bid}.{wildcards.step}.macs2_daugherty2017
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
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
        # Rev3Q3 -- macs2_daugherty2017
        expand(pf('atac814_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_daugherty2017', '_treat_pileup.bw', 'atac814'), sample=[* config['stages_rep'] ]),
        # Rev3Q1/2 -- noSPMR peak height counts
        expand(pf('atac814_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_treat_pileup.bw', 'atac814'), sample=[* config['stages_rep'] ]),

rule atac814_ce11:
    input:
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe_ce11.rm_unmapped_pe.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac814'), sample=config['stages_wt_rep']),
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe_ce11.rm_unmapped_pe.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac814'), sample=config['stages_wt_rep']),
        expand(pf('atac814_{sample}', 'tg_se.bwa_se_ce11.rm_unmapped.rm_chrM.rm_blacklist_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'), sample=config['stages_rep']),
        # ce11 stage-specific coverage tracks
        expand('atac814_geo/tracks_ce11/atac_{sample}_ce11.bw', sample=config['stages']),
        expand('atac814_geo/alignments_ce11/atac_{sample}_ce11.bam', sample=config['stages_rep']),

rule atac814_mapq0:
    input:
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.macs2_pe_lt200', '_treat_pileup.bw', 'atac814'), sample=config['stages_wt_rep']),
        expand(pf('atac814_{sample}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.macs2_pe_lt300', '_treat_pileup.bw', 'atac814'), sample=config['stages_wt_rep']),
        expand(pf('atac814_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac814'), sample=config['stages_rep']),

rule atac824:
    input:
        expand(pf('atac824_{sample}', 'c_r1', '.txt', 'atac824'), sample=techreps_collapse(config['atac824'].keys())),
        expand(pf('atac824_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac824'), sample=techreps_collapse(config['atac824'].keys())),
        expand(pf('atac824_{sample}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.sample_prp.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', 'atac824'), sample=techreps_collapse(config['atac824'].keys())),
        expand(pf('atac824_{condition}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all.mean_by_condition', '.bw', 'atac824'), condition=config['atac824_tissues']),

rule atac:
    input:
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt200', '_treat_pileup.bw', 'atac'), bid=['HS298_JA26_N2_atac_S1']),#[*config['atac']]),
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'atac'), bid=['HS298_JA26_N2_atac_S1']),#[*config['atac']]),
