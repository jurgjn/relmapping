
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

rule hmod_reads_ce10_md5sum:
    input: 
        'hmod749/reads/{sample}.fastq.gz',
    output:
        'hmod749/reads_md5sum/{sample}.reads_md5sum.txt'
    shell: 'md5sum {input} > {output}'


rule hmod_reads_ce10_readlen:
    input: 
        'hmod749/reads/{sample}.fastq.gz',
    output:
        'hmod749/reads_readlen/{sample}.reads_readlen.txt'
    run:
        (fp_inp, fp_out) = (str(input), str(output))
        n_reads = int(1e6) # estimate read length from the first 1M reads
        counts = collections.Counter(map(len, itertools.islice(hts.FastqReader(fp_inp), n_reads)))
        print(max(counts.keys()))
        with open(fp_out, 'w') as fh_out:
            #fh_out.write('%d-%d' % (min(counts.keys()), max(counts.keys())))
            fh_out.write('%d' % (max(counts.keys()),))

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
        #expand('hmod749/log2_by_rep/{sample}.bw', sample=config['hmod'].keys()),
        #expand('hmod749/log2/{sample}.bw', sample=config['hmod749']),
        expand('hmod749/reads_md5sum/{sample}.reads_md5sum.txt', sample=config['hmod'].keys()),
        expand('hmod749/reads_readlen/{sample}.reads_readlen.txt', sample=config['hmod'].keys()),
