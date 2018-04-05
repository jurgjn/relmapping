# snakemake  --use-conda --cores 12 -s workflows/scap814_pool.snakefile scap814 --dryrun
configfile: 'workflows/config.yaml'
include: 'shared.snakefile'
include: 'shared_dm.snakefile'

def scap814_bam_input(wildcards):
    l_ = [pf(bid, wildcards.step, '.bam', 'scap') for bid in techreps_retrieve(wildcards.sample, config['scap814'])]
    assert len(l_) > 0
    return l_

rule scap814_bam:
    input:
        scap814_bam_input,
    output:
        pf('scap814_{sample}', '{step}', '.bam', 'scap814'),
    run:
        print(input, output)
        if os.path.isfile(str(input)):
            shell("ln -s `pwd`/{input} `pwd`/{output}")
            shell("touch -h `pwd`/{output}")
        else:
            arg_input = ' '.join(input)
            shell('samtools merge %(output)s %(arg_input)s' % locals())

rule scap814:
    input:
        expand(pf('scap814_{sample}', 'tg_se.bwa_se', '.bam', 'scap814'), sample=techreps_collapse(config['scap814'].keys())),

# snakemake  --use-conda --cores 12 -s workflows/scap814_pool.snakefile atac814 --dryrun
def atac814_r1_input(wildcards):
    l_ = ['samples/%s.r1.fq.gz' % (bid,) for bid in techreps_retrieve(wildcards.sample, config['atac814'])]
    assert len(l_) > 0
    return l_

rule atac814_r1_:
    input:
        atac814_r1_input,
    output:
        'samples/atac814_{sample}.r1.fq.gz'
    run:
        print(input, output)
        if os.path.isfile(str(input)):
            shell("ln -s `pwd`/{input} `pwd`/{output}")
            shell("touch -h `pwd`/{output}")
        else:
            arg_input = ' '.join(input)
            shell('cat %(arg_input)s > %(output)s' % locals())

rule atac814_r1:
    input:
        expand('samples/atac814_{sample}.r1.fq.gz', sample=techreps_collapse(config['atac814'].keys())),
