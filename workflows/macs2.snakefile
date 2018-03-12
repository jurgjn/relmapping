
rule macs2_pe_lt100:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_pe_lt100', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_pe_lt100', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_pe_lt100', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_pe_lt100', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_pe_lt100', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_pe_lt100', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_pe_lt100', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --bdg --SPMR --gsize ce --nolambda'
    threads: 10
    shell: '''
        samtools view -f 3 --threads {threads} -h {input[0]} | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 100*100)' | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} \
            --treatment {output[5]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_pe_lt100 --name {wildcards.bid}.{wildcards.step}.macs2_pe_lt100
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_pe_lt150:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_pe_lt150', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_pe_lt150', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_pe_lt150', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_pe_lt150', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_pe_lt150', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_pe_lt150', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_pe_lt150', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --bdg --SPMR --gsize ce --nolambda'
    threads: 10
    shell: '''
        samtools view -f 3 --threads {threads} -h {input[0]} | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 150*150)' | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} \
            --treatment {output[5]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_pe_lt150 --name {wildcards.bid}.{wildcards.step}.macs2_pe_lt150
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_pe_lt200:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_pe_lt200', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_pe_lt200', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_pe_lt200', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_pe_lt200', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_pe_lt200', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_pe_lt200', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_pe_lt200', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --bdg --SPMR --gsize ce --nolambda'
    threads: 1
    shell: '''
        samtools view -f 3 --threads {threads} -h {input[0]} | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 200*200)' | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} \
            --treatment {output[5]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_pe_lt200 --name {wildcards.bid}.{wildcards.step}.macs2_pe_lt200
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_pe_lt300:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_pe_lt300', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_pe_lt300', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_pe_lt300', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_pe_lt300', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_pe_lt300', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_pe_lt300', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_pe_lt300', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --bdg --SPMR --gsize ce --nolambda'
    threads: 1
    shell: '''
        samtools view -f 3 --threads {threads} -h {input[0]} | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 300*300)' | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} \
            --treatment {output[5]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_pe_lt300 --name {wildcards.bid}.{wildcards.step}.macs2_pe_lt300
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_pe_lt300_keepdup_all:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_pe_lt300_keepdup_all', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --bdg --SPMR --gsize ce --nolambda --keep-dup all'
    threads: 1
    shell: '''
        samtools view -f 3 --threads {threads} -h {input[0]} | awk 'substr($1,1,1) == "@" || (int($9)*int($9) < 300*300)' | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} \
            --treatment {output[5]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_pe_lt300_keepdup_all --name {wildcards.bid}.{wildcards.step}.macs2_pe_lt300_keepdup_all
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize200:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize200', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize200', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize200', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize200', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize200', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize200', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize200', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 200'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize200 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize200
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize200_keepdup_all:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize200_keepdup_all', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 200 --keep-dup all'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize200_keepdup_all --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize200_keepdup_all
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize150:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize150', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize150', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize150', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize150', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize150', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize150', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize150', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize150 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize150
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize125:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize125', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize125', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize125', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize125', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize125', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize125', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize125', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 125'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize125 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize125
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize100:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize100', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize100', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize100', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize100', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize100', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize100', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize100', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 100'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize100 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize100
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize200_shiftm100:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize200_shiftm100', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 200 --shift -100'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize200_shiftm100 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize200_shiftm100
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize100_shiftm50:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize100_shiftm50', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 100 --shift -50'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize100_shiftm50 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize100_shiftm50
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize150_shiftm75:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150 --shift -75'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize150_shiftm75 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize150_shiftm75
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize1_shiftm75_keepdup_all: # Bigwig tracks of ATAC-seq cut sites, useful for counting tags hitting specific peaks
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize1_shiftm75_keepdup_all', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --gsize ce --nolambda --nomodel --extsize 1 --shift -75 --keep-dup all' #--SPMR no normalisation as DESeq works with raw counts
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize1_shiftm75_keepdup_all --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize1_shiftm75_keepdup_all
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize1_keepdup_all: # Bigwig tracks of ATAC-seq cut sites, useful for counting tags hitting specific peaks
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize1_keepdup_all', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --gsize ce --nolambda --nomodel --extsize 1 --keep-dup all' #--SPMR no normalisation as DESeq works with raw counts
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize1_keepdup_all --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize1_keepdup_all
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize150_shiftm75_keepdup_all:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150 --shift -75  --keep-dup all'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_all --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_all
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize150_shiftm75_keepdup_all_noSPMR:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --gsize ce --nolambda --nomodel --extsize 150 --shift -75  --keep-dup all'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_all_noSPMR
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize150_shiftm75_keepdup_auto:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_auto', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150 --shift -75  --keep-dup auto'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_auto --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_auto
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize150_shiftm75_keepdup_1:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize150_shiftm75_keepdup_1', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150 --shift -75  --keep-dup 1'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_1 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize150_shiftm75_keepdup_1
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize50_shiftm25:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize50_shiftm25', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 50 --shift -25'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize50_shiftm25 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize50_shiftm25
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''

rule macs2_se_extsize20_shiftm10:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '_peaks.narrowPeak', '{prefix}'), # [0]
        pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '_peaks.xls', '{prefix}'), # [1]
        pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '_summits.bed', '{prefix}'), # [2]
        temp(pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '_treat_pileup.bdg', '{prefix}')), # [3]
        pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '_treat_pileup.bw', '{prefix}'), # [4]
        temp(pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '_treat.bam', '{prefix}')), # [5]
        temp(pf('{bid}', '{step}.macs2_se_extsize20_shiftm10', '.chroms', '{prefix}')), # [6]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 20 --shift -10'
    shell: '''
        samtools view --threads {threads} -h {input[0]} | samtools view -b - > {output[5]}
        macs2 callpeak {params.macs2_args} --treatment {output[5]} --outdir {wildcards.prefix}/{wildcards.step}.macs2_se_extsize20_shiftm10 --name {wildcards.bid}.{wildcards.step}.macs2_se_extsize20_shiftm10
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[6]}
        sort -k1,1 -k2,2n -o {output[3]} {output[3]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
        bedGraphToBigWig {output[3]} {output[6]} {output[4]}
    '''
