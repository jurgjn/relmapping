# Call yapc peaks on a single data set via self-pseudoreplicates
rule yapc_spr_pe:
    input:
        pf('{sample}_spr1', '{step}', '_treat_pileup.bw', '{prefix}'),
        pf('{sample}_spr2', '{step}', '_treat_pileup.bw', '{prefix}'),
    output:
        pf('{sample}', '{step}.yapc_spr_pe', '_coverage.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_pe', '_d2smooth.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_pe', '.tsv', '{prefix}'),
    params:
        yapc_args = '--smoothing-window-width 75'
    shell:
        '~/repos/yapc/yapc {wildcards.prefix}/{wildcards.step}.yapc_spr_pe/{wildcards.sample}.{wildcards.step}.yapc_spr_pe {wildcards.sample} {input[0]} {input[1]} {params.yapc_args}'

rule yapc_spr_se:
    input:
        pf('{sample}_spr1', '{step}', '_treat_pileup.bw', '{prefix}'),
        pf('{sample}_spr2', '{step}', '_treat_pileup.bw', '{prefix}'),
    output:
        pf('{sample}', '{step}.yapc_spr_se', '_coverage.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_se', '_d2smooth.bw', '{prefix}'),
        pf('{sample}', '{step}.yapc_spr_se', '.tsv', '{prefix}'),
    params:
        yapc_args = '--smoothing-window-width 150'
    shell:
        '~/repos/yapc/yapc {wildcards.prefix}/{wildcards.step}.yapc_spr_se/{wildcards.sample}.{wildcards.step}.yapc_spr_se {wildcards.sample} {input[0]} {input[1]} {params.yapc_args}'
