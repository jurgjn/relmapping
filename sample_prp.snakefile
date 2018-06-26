
configfile: 'workflows/config.yaml'

include: 'workflows/shared.snakefile'
include: 'workflows/shared_dm.snakefile'

sys.path.append(os.path.expanduser('~/relmapping/scripts/yarp'))
import yarp as yp

rule spr1: # sample self-pseudoreplicate #1 from a single sample
    input:
        pf('{sample}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{sample}_spr1', '{step}', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools view -h -S --threads {threads} {input[0]} \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (rand() < .5)' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule spr2: # sample self-pseudoreplicate #2 from a single sample
    input:
        pf('{sample}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{sample}_spr2', '{step}', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools view -h -S --threads {threads} {input[0]} \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (.5 <= rand())' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule atac_cb3_spr: # snakemake --use-conda --cores 12 -s sample_prp.snakefile atac_cb3_spr --dryrun
    input:
        expand(pf('{sample}_spr1', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10', '.bam', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}_spr2', 'tg_pe.bwa_pe_cb3.rm_unmapped_pe.rm_chrM.rm_q10', '.bam', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}_spr1', 'tg_se.bwa_se_cb3.rm_unmapped.rm_chrM.rm_q10', '.bam', 'atac'), sample=config['atac_cb3']),
        expand(pf('{sample}_spr2', 'tg_se.bwa_se_cb3.rm_unmapped.rm_chrM.rm_q10', '.bam', 'atac'), sample=config['atac_cb3']),

"""
rule prp1: # sample pseudoreplicate #1 from two biological replicates
    input:
        pf('{cond}_rep1', '{step}', '.bam', '{prefix}'),
        pf('{cond}_rep2', '{step}', '.bam', '{prefix}'),
    output:
        pf('{cond}_prp1', '{step}', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools merge -u --threads {threads} - {input[0]} {input[1]} | samtools view -h -S --threads {threads} - \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (rand() < .5)' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''

rule prp2: # sample pseudoreplicate #2 from two biological replicates
    input:
        pf('{cond}_rep1', '{step}', '.bam', '{prefix}'),
        pf('{cond}_rep2', '{step}', '.bam', '{prefix}'),
    output:
        pf('{cond}_prp2', '{step}', '.bam', '{prefix}'),
    threads: 4
    shell: '''
        samtools merge -u --threads {threads} - {input[0]} {input[1]} | samtools view -h -S --threads {threads} - \
        | awk 'BEGIN{{srand(42);}} substr($1,1,1) == "@" || (.5 <= rand())' \
        | samtools view -b --threads {threads} - > {output[0]}
    '''
"""
