rule bwa_pe: # rm -f ... gets rid of temporary files left from a previously crashed run...
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', '{step}', '.r2.fq.gz', '{prefix}'),
        'shared/ce10.fa',
        'shared/ce10.fa.amb',
        'shared/ce10.fa.ann',
        'shared/ce10.fa.bwt',
        'shared/ce10.fa.pac',
        'shared/ce10.fa.sa',
    output:
        pf('{bid}', '{step}.bwa_pe', '.bam', '{prefix}'),
    threads: 4
    params: genome_fa='shared/ce10.fa'
    shell: '''
        rm -f {output}_*
        bwa sampe {params.genome_fa} \
            <(bwa aln -t {threads} {params.genome_fa} {input[0]}) \
            <(bwa aln -t {threads} {params.genome_fa} {input[1]}) \
            {input[0]} \
            {input[1]} \
        | samtools view -b -@ {threads} - \
        | samtools sort -@ {threads} -O sam -T {output}_tmpA -n - \
        | samtools fixmate -O sam - - \
        | samtools sort -@ {threads} -O bam -T {output}_tmpB - \
        > {output}
        '''

rule bwa_pe_ce11: # rm -f ... gets rid of temporary files left from a previously crashed run...
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', '{step}', '.r2.fq.gz', '{prefix}'),
        'shared/ce11.fa',
        'shared/ce11.fa.amb',
        'shared/ce11.fa.ann',
        'shared/ce11.fa.bwt',
        'shared/ce11.fa.pac',
        'shared/ce11.fa.sa',
    output:
        pf('{bid}', '{step}.bwa_pe_ce11', '.bam', '{prefix}'),
    threads: 4
    params: genome_fa='shared/ce11.fa'
    shell: '''
        rm -f {output}_*
        bwa sampe {params.genome_fa} \
            <(bwa aln -t {threads} {params.genome_fa} {input[0]}) \
            <(bwa aln -t {threads} {params.genome_fa} {input[1]}) \
            {input[0]} \
            {input[1]} \
        | samtools view -b -@ {threads} - \
        | samtools sort -@ {threads} -O sam -T {output}_tmpA -n - \
        | samtools fixmate -O sam - - \
        | samtools sort -@ {threads} -O bam -T {output}_tmpB - \
        > {output}
        '''

rule bwa_pe_cb3: # rm -f ... gets rid of temporary files left from a previously crashed run...
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', '{step}', '.r2.fq.gz', '{prefix}'),
        'shared/cb3.fa',
        'shared/cb3.fa.amb',
        'shared/cb3.fa.ann',
        'shared/cb3.fa.bwt',
        'shared/cb3.fa.pac',
        'shared/cb3.fa.sa',
    output:
        pf('{bid}', '{step}.bwa_pe_cb3', '.bam', '{prefix}'),
    threads: 4
    params: genome_fa='shared/cb3.fa'
    shell: '''
        rm -f {output}_*
        bwa sampe {params.genome_fa} \
            <(bwa aln -t {threads} {params.genome_fa} {input[0]}) \
            <(bwa aln -t {threads} {params.genome_fa} {input[1]}) \
            {input[0]} \
            {input[1]} \
        | samtools view -b -@ {threads} - \
        | samtools sort -@ {threads} -O sam -T {output}_tmpA -n - \
        | samtools fixmate -O sam - - \
        | samtools sort -@ {threads} -O bam -T {output}_tmpB - \
        > {output}
        '''

rule bwa_pe_CB4: # rm -f ... gets rid of temporary files left from a previously crashed run...
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', '{step}', '.r2.fq.gz', '{prefix}'),
        'shared/CB4.fa',
        'shared/CB4.chroms',
        'shared/CB4.fa.amb',
        'shared/CB4.fa.ann',
        'shared/CB4.fa.bwt',
        'shared/CB4.fa.pac',
        'shared/CB4.fa.sa',
    output:
        pf('{bid}', '{step}.bwa_pe_CB4', '.bam', '{prefix}'),
    threads: 4
    params: genome_fa='shared/CB4.fa'
    shell: '''
        rm -f {output}_*
        bwa sampe {params.genome_fa} \
            <(bwa aln -t {threads} {params.genome_fa} {input[0]}) \
            <(bwa aln -t {threads} {params.genome_fa} {input[1]}) \
            {input[0]} \
            {input[1]} \
        | samtools view -b -@ {threads} - \
        | samtools sort -@ {threads} -O sam -T {output}_tmpA -n - \
        | samtools fixmate -O sam - - \
        | samtools sort -@ {threads} -O bam -T {output}_tmpB - \
        > {output}
        '''

rule bwa_pe_ecoli:
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', '{step}', '.r2.fq.gz', '{prefix}'),
    output:
        pf('{bid}', '{step}.bwa_pe_ecoli', '.bam', '{prefix}'),
    threads: 4
    params: genome_fa='shared/ecoli.fa'
    shell: '''
        bwa sampe {params.genome_fa} \
            <(bwa aln -t {threads} {params.genome_fa} {input[0]}) \
            <(bwa aln -t {threads} {params.genome_fa} {input[1]}) \
            {input[0]} \
            {input[1]} \
        | samtools view -b -@ {threads} - \
        | samtools sort -@ {threads} -O sam -T {output}_tmpA -n - \
        | samtools fixmate -O sam - - \
        | samtools sort -@ {threads} -O bam -T {output}_tmpB - \
        > {output}
        '''

rule bwa_se:
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        'shared/ce10.fa',
        'shared/ce10.fa.amb',
        'shared/ce10.fa.ann',
        'shared/ce10.fa.bwt',
        'shared/ce10.fa.pac',
        'shared/ce10.fa.sa',
    output: pf('{bid}', '{step}.bwa_se', '.bam', '{prefix}')
    threads: 4
    params: genome_fa='shared/ce10.fa'
    shell: 'bwa samse {params.genome_fa} <(bwa aln -t {threads} {params.genome_fa} {input[0]}) {input[0]} | samtools sort -@ {threads} -O bam -T {output}_tmp -o {output} -'

rule bwa_se_ce11:
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        'shared/ce11.fa',
        'shared/ce11.fa.amb',
        'shared/ce11.fa.ann',
        'shared/ce11.fa.bwt',
        'shared/ce11.fa.pac',
        'shared/ce11.fa.sa',
    output: pf('{bid}', '{step}.bwa_se_ce11', '.bam', '{prefix}')
    threads: 4
    params: genome_fa='shared/ce11.fa'
    shell: 'bwa samse {params.genome_fa} <(bwa aln -t {threads} {params.genome_fa} {input[0]}) {input[0]} | samtools sort -@ {threads} -O bam -T {output}_tmp -o {output} -'

rule bwa_se_cb3:
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        'shared/cb3.fa',
        'shared/cb3.fa.amb',
        'shared/cb3.fa.ann',
        'shared/cb3.fa.bwt',
        'shared/cb3.fa.pac',
        'shared/cb3.fa.sa',
    output: pf('{bid}', '{step}.bwa_se_cb3', '.bam', '{prefix}')
    threads: 4
    params: genome_fa='shared/cb3.fa'
    shell: 'bwa samse {params.genome_fa} <(bwa aln -t {threads} {params.genome_fa} {input[0]}) {input[0]} | samtools sort -@ {threads} -O bam -T {output}_tmp -o {output} -'

rule bwa_se_CB4:
    input:
        pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
        'shared/CB4.fa',
        'shared/CB4.fa.amb',
        'shared/CB4.fa.ann',
        'shared/CB4.fa.bwt',
        'shared/CB4.fa.pac',
        'shared/CB4.fa.sa',
    output: pf('{bid}', '{step}.bwa_se_CB4', '.bam', '{prefix}')
    threads: 4
    params: genome_fa='shared/CB4.fa'
    shell: 'bwa samse {params.genome_fa} <(bwa aln -t {threads} {params.genome_fa} {input[0]}) {input[0]} | samtools sort -@ {threads} -O bam -T {output}_tmp -o {output} -'

rule bwa_se_ecoli:
    input: pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}'),
    output: pf('{bid}', '{step}.bwa_se_ecoli', '.bam', '{prefix}')
    threads: 4
    params: genome_fa='shared/ecoli.fa'
    shell: 'bwa samse {params.genome_fa} <(bwa aln -t {threads} {params.genome_fa} {input}) {input} | samtools sort -@ {threads} -O bam -T {output}_tmp -o {output} -'
