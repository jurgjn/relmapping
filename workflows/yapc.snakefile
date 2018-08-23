# Call peaks using the yapc peak caller, either on biological replicates, or 
# pseudoreplicate (sample_prp1, sample_prp2) / self-pseudoreplicate (sample_spr1, sample_spr2) data.

rule sample_prp1: # sample pseudoreplicate #1 from two biological replicates
    input:
        pf('{sample}_rep1', '{step}', '.bam', '{prefix}'),
        pf('{sample}_rep2', '{step}', '.bam', '{prefix}'),
    output:
        pf('{sample}_rep1', '{step}.sample_prp', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools merge -u --threads {threads} - {input[0]} {input[1]} | samtools view -h -S --threads {threads} - \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (rand() < .5)' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule sample_prp2: # sample pseudoreplicate #2 from two biological replicates
    input:
        pf('{sample}_rep1', '{step}', '.bam', '{prefix}'),
        pf('{sample}_rep2', '{step}', '.bam', '{prefix}'),
    output:
        pf('{sample}_rep2', '{step}.sample_prp', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools merge -u --threads {threads} - {input[0]} {input[1]} | samtools view -h -S --threads {threads} - \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (.5 <= rand())' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule sample_spr1: # sample self-pseudoreplicate #1 from a single sample
    input:
        pf('{sample}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{sample}_rep1', '{step}.sample_spr', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools view -h -S --threads {threads} {input[0]} \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (rand() < .5)' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule sample_spr2: # sample self-pseudoreplicate #2 from a single sample
    input:
        pf('{sample}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{sample}_rep2', '{step}.sample_spr', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools view -h -S --threads {threads} {input[0]} \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (.5 <= rand())' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule yapc_spr_pe:
    input:
        pf('{sample}_rep1', '{step}', '_treat_pileup.bw', '{prefix}'),
        pf('{sample}_rep2', '{step}', '_treat_pileup.bw', '{prefix}'),
    output:
        pf('{sample}', '{step}.yapc_spr_pe', '_coverage.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_pe', '_d2smooth.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_pe', '.tsv', '{prefix}'),
    params:
        yapc_args = '--smoothing-window-width 75'
    shell:
        '~/relmapping/yapc/yapc {wildcards.prefix}/{wildcards.step}.yapc_spr_pe/{wildcards.sample}.{wildcards.step}.yapc_spr_pe {wildcards.sample} {input[0]} {input[1]} {params.yapc_args}'

rule yapc_spr_se:
    input:
        pf('{sample}_rep1', '{step}', '_treat_pileup.bw', '{prefix}'),
        pf('{sample}_rep2', '{step}', '_treat_pileup.bw', '{prefix}'),
    output:
        pf('{sample}', '{step}.yapc_spr_se', '_coverage.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_se', '_d2smooth.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_se', '.tsv', '{prefix}'),
    params:
        yapc_args = '--smoothing-window-width 150'
    shell:
        '~/relmapping/yapc/yapc {wildcards.prefix}/{wildcards.step}.yapc_spr_se/{wildcards.sample}.{wildcards.step}.yapc_spr_se {wildcards.sample} {input[0]} {input[1]} {params.yapc_args}'
