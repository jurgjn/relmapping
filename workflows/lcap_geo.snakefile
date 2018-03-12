
rule lcap_alignment_raw_reads_by_rep:
    input:
        pf('lcap728_{bid}', 'trim20.bwa_pe', '.bam', 'lcap728'),
    output:
        pf('lcap_{bid}', 'alignment_raw_reads_by_rep', '.bam', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule lcap_alignment_all_reads_by_rep:
    input:
        pf('lcap728_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist', '.bam', 'lcap728'),
    output:
        pf('lcap_{bid}', 'alignment_all_reads_by_rep', '.bam', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule lcap_alignment_q10_reads_by_rep:
    input:
        pf('lcap728_{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10', '.bam', 'lcap728'),
    output:
        pf('lcap_{bid}', 'alignment_q10_reads_by_rep', '.bam', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule coverage_all_reads_by_rep_fwd:
    input:
        pf('lcap728_{sample}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_fwd', '.bw', 'lcap728'),
    output:
        pf('lcap_{sample}', 'coverage_all_reads_by_rep_fwd', '.bw', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule coverage_all_reads_by_rep_rev:
    input:
        pf('lcap728_{sample}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.filled_rev', '.bw', 'lcap728'),
    output:
        pf('lcap_{sample}', 'coverage_all_reads_by_rep_rev', '.bw', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule coverage_q10_reads_by_rep_fwd:
    input:
        pf('lcap728_{sample}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap728'),
    output:
        pf('lcap_{sample}', 'coverage_q10_reads_by_rep_fwd', '.bw', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule coverage_q10_reads_by_rep_rev:
    input:
        pf('lcap728_{sample}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev', '.bw', 'lcap728'),
    output:
        pf('lcap_{sample}', 'coverage_q10_reads_by_rep_rev', '.bw', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule coverage_q10_reads_by_stage_fwd:
    input:
        pf('lcap_{sample}_rep1', 'coverage_q10_reads_by_rep_fwd', '.bw', 'lcap_geo'),
        pf('lcap_{sample}_rep2', 'coverage_q10_reads_by_rep_fwd', '.bw', 'lcap_geo'),
    output:
        temp(pf('lcap_{sample}', 'coverage_q10_reads_by_stage_fwd', '.bdg', 'lcap_geo')),
        pf('lcap_{sample}', 'coverage_q10_reads_by_stage_fwd', '.bw', 'lcap_geo'),
    shell: '''
        wiggletools write_bg {output[0]} sum {input[0]} {input[1]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
    '''

rule coverage_q10_reads_by_stage_rev:
    input:
        pf('lcap_{sample}_rep1', 'coverage_q10_reads_by_rep_rev', '.bw', 'lcap_geo'),
        pf('lcap_{sample}_rep2', 'coverage_q10_reads_by_rep_rev', '.bw', 'lcap_geo'),
    output:
        temp(pf('lcap_{sample}', 'coverage_q10_reads_by_stage_rev', '.bdg', 'lcap_geo')),
        pf('lcap_{sample}', 'coverage_q10_reads_by_stage_rev', '.bw', 'lcap_geo'),
    shell: '''
        wiggletools write_bg {output[0]} sum {input[0]} {input[1]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
    '''

rule startbp_q10_reads_by_rep_fwd:
    input:
        pf('lcap728_{sample}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_fwd', '.bw', 'lcap728'),
    output:
        pf('lcap_{sample}', 'startbp_q10_reads_by_rep_fwd', '.bw', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule startbp_q10_reads_by_rep_rev:
    input:
        pf('lcap728_{sample}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.startbp_rev', '.bw', 'lcap728'),
    output:
        pf('lcap_{sample}', 'startbp_q10_reads_by_rep_rev', '.bw', 'lcap_geo'),
    shell:
        'cp {input} {output}'

rule lcap_geo:
    input:
        expand(pf('lcap_{sample}', 'alignment_raw_reads_by_rep', '.bam', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'alignment_all_reads_by_rep', '.bam', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'alignment_q10_reads_by_rep', '.bam', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'alignment_raw_reads_by_rep', '.bam.bai', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'alignment_all_reads_by_rep', '.bam.bai', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'alignment_q10_reads_by_rep', '.bam.bai', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'coverage_all_reads_by_rep_fwd', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'coverage_all_reads_by_rep_rev', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'coverage_q10_reads_by_rep_fwd', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'coverage_q10_reads_by_rep_rev', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'coverage_q10_reads_by_stage_fwd', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'coverage_q10_reads_by_stage_rev', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'startbp_q10_reads_by_rep_fwd', '.bw', 'lcap_geo'), sample=config['stages_rep']),
        expand(pf('lcap_{sample}', 'startbp_q10_reads_by_rep_rev', '.bw', 'lcap_geo'), sample=config['stages_rep']),

rule lcap749_read1:
    input:
        'samples/lcap728_{sample}.r1.fq.gz'
    output:
        'lcap749/reads/lcap_{sample}.read1.fastq.gz'
    shell:
        'cp {input} {output}'

rule lcap749_read2:
    input:
        'samples/lcap728_{sample}.r2.fq.gz'
    output:
        'lcap749/reads/lcap_{sample}.read2.fastq.gz'
    shell:
        'cp {input} {output}'

rule lcap749_fwd:
    input:
        pf('lcap728_{sample}_rep1', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap728'),
        pf('lcap728_{sample}_rep2', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap728'),
    output:
        'lcap749/tracks_fwd/lcap_{sample}_fwd.bw'
    shell:
        'scripts/bigWiggleTools.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}'

rule lcap749_rev:
    input:
        pf('lcap728_{sample}_rep1', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap728'),
        pf('lcap728_{sample}_rep2', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap728'),
    output:
        'lcap749/tracks_rev/lcap_{sample}_rev.bw'
    shell:
        'scripts/bigWiggleTools.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}'

rule lcap749_rev_neg:
    input:
        pf('lcap728_{sample}_rep1', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap728'),
        pf('lcap728_{sample}_rep2', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev_sfnorm', '.bw', 'lcap728'),
    output:
        'lcap749/tracks_rev_neg/lcap_{sample}_rev_neg.bw'
    shell:
        'scripts/bigWiggleTools.ipy write {output[0]} "scale -0.1" bin 10 mean {input[0]} {input[1]}'

rule lcap749:
    input:
        expand('lcap749/reads/lcap_{sample}.read1.fastq.gz', sample=config['lcap749'].keys()),
        expand('lcap749/reads/lcap_{sample}.read2.fastq.gz', sample=config['lcap749'].keys()),
        expand('lcap749/tracks_fwd/lcap_{sample}_fwd.bw', sample=config['stages_rep']),
        expand('lcap749/tracks_rev/lcap_{sample}_rev.bw', sample=config['stages_rep']),
        expand('lcap749/tracks_rev_neg/lcap_{sample}_rev_neg.bw', sample=config['stages_rep']),

