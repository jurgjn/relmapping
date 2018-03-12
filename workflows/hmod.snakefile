
rule hmod_reads_ce10:
    input:
        lambda wildcards: config['hmod'][wildcards.sample]['reads_ce10'],
    output:
        'hmod749/reads/{sample}.fastq.gz',
    shell:
        'cp {input} {output}'

rule hmod_log2_by_rep_ce10:
    input:
        lambda wildcards: config['hmod'][wildcards.sample]['log2_ce10'],
    output:
        'hmod749/log2_by_rep/{sample}.bw',
    shell:
        'cp {input} {output}'

rule hmod_log2:
    input:
        'hmod749/log2_by_rep/{sample}_rep1.bw',
        'hmod749/log2_by_rep/{sample}_rep2.bw',
    output:
        'hmod749/log2/{sample}.bw',
    shell: '''
        scripts/bigWiggleTools.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}
        '''

"""
rule hmod_tracks_by_rep_ce10:
    input:
        lambda wildcards: config['hmod'][wildcards.sample]['track_ce10'],
    output:
        'hmod749/tracks_by_rep/{sample}.bw',
    shell:
        'cp {input} {output}'

rule hmod_tracks:
    input:
        'hmod749/tracks_by_rep/{sample}_rep1.bw',
        'hmod749/tracks_by_rep/{sample}_rep2.bw',
    output:
        'hmod749/tracks/{sample}.bw',
    shell: '''
        scripts/bigWiggleTools.ipy write {output[0]} scale 0.1 bin 10 mean {input[0]} {input[1]}
        '''
"""

rule hmod749:
    input:
        #expand('hmod749/reads/{sample}.fastq.gz', sample=config['hmod'].keys()),
        expand('hmod749/log2_by_rep/{sample}.bw', sample=config['hmod'].keys()),
        expand('hmod749/log2/{sample}.bw', sample=config['hmod749']),
