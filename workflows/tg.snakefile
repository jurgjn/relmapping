rule trim20_r1:
    input:
        'samples/{bid}.r1.fq.gz',
    output:
        pf('{bid}', 'trim20', '.r1.fq.gz', '{prefix}'),
    threads: 1
    shell:
        '''
        pigz -c -d -p {threads} {input} | fastx_trimmer -l 20 -Q33 | pigz -c -p {threads} - > {output}
        '''

rule trim20_r2:
    input:
        'samples/{bid}.r2.fq.gz',
    output:
        pf('{bid}', 'trim20', '.r2.fq.gz', '{prefix}'),
    threads: 1
    shell:
        '''
        pigz -c -d -p {threads} {input} | fastx_trimmer -l 20 -Q33 | pigz -c -p {threads} - > {output}
        '''

rule trim20_1M_r1:
    input:
        'samples/{bid}_1M.r1.fq.gz',
    output:
        pf('{bid}', 'trim20_1M', '.r1.fq.gz', '{prefix}'),
    threads: 1
    shell:
        '''
        pigz -c -d -p {threads} {input} | fastx_trimmer -l 20 -Q33 | pigz -c -p {threads} - > {output}
        '''

rule trim20_1M_r2:
    input:
        'samples/{bid}_1M.r2.fq.gz',
    output:
        pf('{bid}', 'trim20_1M', '.r2.fq.gz', '{prefix}'),
    threads: 1
    shell:
        '''
        pigz -c -d -p {threads} {input} | fastx_trimmer -l 20 -Q33 | pigz -c -p {threads} - > {output}
        '''

rule trim14_r1:
    input:
        'samples/{bid}.r1.fq.gz',
    output:
        pf('{bid}', 'trim14', '.r1.fq.gz', '{prefix}'),
    threads: 1
    shell:
        '''
        pigz -c -d -p {threads} {input} | fastx_trimmer -l 14 -Q33 | pigz -c -p {threads} - > {output}
        '''

rule trim14_r2:
    input:
        'samples/{bid}.r2.fq.gz',
    output:
        pf('{bid}', 'trim14', '.r2.fq.gz', '{prefix}'),
    threads: 1
    shell:
        '''
        pigz -c -d -p {threads} {input} | fastx_trimmer -l 14 -Q33 | pigz -c -p {threads} - > {output}
        '''

rule tg_pe_raw:
    input:
        'samples/{bid}.r1.fq.gz',
        'samples/{bid}.r2.fq.gz',
    output:
        '{prefix}/tg_pe/_raw/{bid}.r1_val_1.fq.gz',
        '{prefix}/tg_pe/_raw/{bid}.r2_val_2.fq.gz',
        '{prefix}/tg_pe/_raw/{bid}.r1_val_1_fastqc.html',
        '{prefix}/tg_pe/_raw/{bid}.r2_val_2_fastqc.html',
        '{prefix}/tg_pe/_raw/{bid}.r1_val_1_fastqc.zip',
        '{prefix}/tg_pe/_raw/{bid}.r2_val_2_fastqc.zip',
        '{prefix}/tg_pe/_raw/{bid}.r1.fq.gz_trimming_report.txt',
        '{prefix}/tg_pe/_raw/{bid}.r2.fq.gz_trimming_report.txt'
    shell: 'trim_galore {input[0]} {input[1]} --output_dir {wildcards.prefix}/tg_pe/_raw --paired --fastqc'

rule tg_pe:
    input:
        '{prefix}/tg_pe/_raw/{bid}.r1_val_1.fq.gz',
        '{prefix}/tg_pe/_raw/{bid}.r2_val_2.fq.gz',
    output:
        pf('{bid}', 'tg_pe', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', 'tg_pe', '.r2.fq.gz', '{prefix}')
    shell:
        '''
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        '''

rule tg_pe_raw_q15:
    input:
        'samples/{bid}.r1.fq.gz',
        'samples/{bid}.r2.fq.gz',
    output:
        '{prefix}/tg_pe_q15/_raw/{bid}.r1_val_1.fq.gz',
        '{prefix}/tg_pe_q15/_raw/{bid}.r2_val_2.fq.gz',
        '{prefix}/tg_pe_q15/_raw/{bid}.r1_val_1_fastqc.html',
        '{prefix}/tg_pe_q15/_raw/{bid}.r2_val_2_fastqc.html',
        '{prefix}/tg_pe_q15/_raw/{bid}.r1_val_1_fastqc.zip',
        '{prefix}/tg_pe_q15/_raw/{bid}.r2_val_2_fastqc.zip',
        '{prefix}/tg_pe_q15/_raw/{bid}.r1.fq.gz_trimming_report.txt',
        '{prefix}/tg_pe_q15/_raw/{bid}.r2.fq.gz_trimming_report.txt'
    shell: 'trim_galore {input[0]} {input[1]} --output_dir {prefix}/tg_pe_q15/_raw --paired --fastqc -q 15'

rule tg_pe_q15:
    input:
        '{prefix}/tg_pe_q15/_raw/{bid}.r1_val_1.fq.gz',
        '{prefix}/tg_pe_q15/_raw/{bid}.r2_val_2.fq.gz',
    output:
        pf('{bid}', 'tg_pe_q15', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', 'tg_pe_q15', '.r2.fq.gz', '{prefix}')
    shell:
        '''
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        '''

rule tg_pe_raw_q10:
    input:
        'samples/{bid}.r1.fq.gz',
        'samples/{bid}.r2.fq.gz',
    output:
        '{prefix}/tg_pe_q10/_raw/{bid}.r1_val_1.fq.gz',
        '{prefix}/tg_pe_q10/_raw/{bid}.r2_val_2.fq.gz',
        '{prefix}/tg_pe_q10/_raw/{bid}.r1_val_1_fastqc.html',
        '{prefix}/tg_pe_q10/_raw/{bid}.r2_val_2_fastqc.html',
        '{prefix}/tg_pe_q10/_raw/{bid}.r1_val_1_fastqc.zip',
        '{prefix}/tg_pe_q10/_raw/{bid}.r2_val_2_fastqc.zip',
        '{prefix}/tg_pe_q10/_raw/{bid}.r1.fq.gz_trimming_report.txt',
        '{prefix}/tg_pe_q10/_raw/{bid}.r2.fq.gz_trimming_report.txt'
    shell: 'trim_galore {input[0]} {input[1]} --output_dir {prefix}/tg_pe_q10/_raw --paired --fastqc -q 10'

rule tg_pe_q10:
    input:
        '{prefix}/tg_pe_q10/_raw/{bid}.r1_val_1.fq.gz',
        '{prefix}/tg_pe_q10/_raw/{bid}.r2_val_2.fq.gz',
    output:
        pf('{bid}', 'tg_pe_q10', '.r1.fq.gz', '{prefix}'),
        pf('{bid}', 'tg_pe_q10', '.r2.fq.gz', '{prefix}')
    shell:
        '''
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        '''

rule tg_se_raw:
    input:
        'samples/{bid}.r1.fq.gz'
    output:
        '{prefix}/tg_se/_raw/{bid}.r1_trimmed_fastqc.html',
        '{prefix}/tg_se/_raw/{bid}.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se/_raw/{bid}.r1_trimmed.fq.gz',
        '{prefix}/tg_se/_raw/{bid}.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se/_raw --fastqc'

rule tg_se:
    input: '{prefix}/tg_se/_raw/{bid}.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'

rule tg_se_1M_raw:
    input:
        'samples/{bid}_1M.r1.fq.gz'
    output:
        '{prefix}/tg_se_1M/_raw/{bid}_1M.r1_trimmed_fastqc.html',
        '{prefix}/tg_se_1M/_raw/{bid}_1M.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se_1M/_raw/{bid}_1M.r1_trimmed.fq.gz',
        '{prefix}/tg_se_1M/_raw/{bid}_1M.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se_1M/_raw --fastqc'

rule tg_se_1M:
    input: '{prefix}/tg_se_1M/_raw/{bid}_1M.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se_1M', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'

rule tg_se_raw_q0:
    input:
        'samples/{bid}.r1.fq.gz'
    output:
        '{prefix}/tg_se_q0/_raw/{bid}.r1_trimmed_fastqc.html',
        '{prefix}/tg_se_q0/_raw/{bid}.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se_q0/_raw/{bid}.r1_trimmed.fq.gz',
        '{prefix}/tg_se_q0/_raw/{bid}.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se_q0/_raw --fastqc -q 0'

rule tg_se_q0:
    input: '{prefix}/tg_se_q0/_raw/{bid}.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se_q0', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'

rule tg_se_raw_q10:
    input:
        'samples/{bid}.r1.fq.gz'
    output:
        '{prefix}/tg_se_q10/_raw/{bid}.r1_trimmed_fastqc.html',
        '{prefix}/tg_se_q10/_raw/{bid}.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se_q10/_raw/{bid}.r1_trimmed.fq.gz',
        '{prefix}/tg_se_q10/_raw/{bid}.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se_q10/_raw --fastqc -q 10'

rule tg_se_q10:
    input: '{prefix}/tg_se_q10/_raw/{bid}.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se_q10', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'

rule tg_se_raw_q15:
    input:
        'samples/{bid}.r1.fq.gz'
    output:
        '{prefix}/tg_se_q15/_raw/{bid}.r1_trimmed_fastqc.html',
        '{prefix}/tg_se_q15/_raw/{bid}.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se_q15/_raw/{bid}.r1_trimmed.fq.gz',
        '{prefix}/tg_se_q15/_raw/{bid}.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se_q15/_raw --fastqc -q 15'

rule tg_se_q15:
    input: '{prefix}/tg_se_q15/_raw/{bid}.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se_q15', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'

rule tg_se_raw_q20:
    input:
        'samples/{bid}.r1.fq.gz'
    output:
        '{prefix}/tg_se_q20/_raw/{bid}.r1_trimmed_fastqc.html',
        '{prefix}/tg_se_q20/_raw/{bid}.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se_q20/_raw/{bid}.r1_trimmed.fq.gz',
        '{prefix}/tg_se_q20/_raw/{bid}.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se_q20/_raw --fastqc -q 20'

rule tg_se_q20:
    input: '{prefix}/tg_se_q20/_raw/{bid}.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se_q20', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'

rule tg_se_raw_len16:
    input:
        'samples/{bid}.r1.fq.gz'
    output:
        '{prefix}/tg_se_len16/_raw/{bid}.r1_trimmed_fastqc.html',
        '{prefix}/tg_se_len16/_raw/{bid}.r1_trimmed_fastqc.zip',
        '{prefix}/tg_se_len16/_raw/{bid}.r1_trimmed.fq.gz',
        '{prefix}/tg_se_len16/_raw/{bid}.r1.fq.gz_trimming_report.txt',
    shell:
        'trim_galore {input} --output_dir {wildcards.prefix}/tg_se_len16/_raw --fastqc --length 16'

rule tg_se_len16:
    input: '{prefix}/tg_se_len16/_raw/{bid}.r1_trimmed.fq.gz',
    output: pf('{bid}', 'tg_se_len16', '.r1.fq.gz', '{prefix}'),
    shell: 'cp {input[0]} {output[0]}'
