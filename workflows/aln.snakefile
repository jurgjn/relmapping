rule fqz_c_r1_samples:
    input: 'samples/{bid}.r1.fq.gz'
    output: pf('{bid}', 'c_r1', '.txt', '{prefix}')
    shell: "gunzip -c {input} | wc -l | awk '{{print $1/4}}' > {output}"

rule fqz_c_r2_samples:
    input: 'samples/{bid}.r2.fq.gz'
    output: pf('{bid}', 'c_r2', '.txt', '{prefix}')
    shell: "gunzip -c {input} | wc -l | awk '{{print $1/4}}' > {output}"

rule md5sum_r1_samples:
    input: 'samples/{bid}.r1.fq.gz'
    output: pf('{bid}', 'md5sum_r1', '.txt', '{prefix}')
    shell: 'md5sum {input} > {output}'

rule md5sum_r2_samples:
    input: 'samples/{bid}.r2.fq.gz'
    output: pf('{bid}', 'md5sum_r2', '.txt', '{prefix}')
    shell: 'md5sum {input} > {output}'

rule readlen_r1_samples:
    input: 'samples/{bid}.r1.fq.gz'
    output: pf('{bid}', 'readlen_r1', '.txt', '{prefix}')
    run:
        (fp_inp, fp_out) = (str(input), str(output))
        n_reads = int(1e6) # estimate read length from the first 1M reads
        counts = collections.Counter(map(len, itertools.islice(hts.FastqReader(fp_inp), n_reads)))
        print(max(counts.keys()))
        with open(fp_out, 'w') as fh_out:
            #fh_out.write('%d-%d' % (min(counts.keys()), max(counts.keys())))
            fh_out.write('%d' % (max(counts.keys()),))

rule readlen_r2_samples:
    input: 'samples/{bid}.r2.fq.gz'
    output: pf('{bid}', 'readlen_r2', '.txt', '{prefix}')
    run:
        (fp_inp, fp_out) = (str(input), str(output))
        n_reads = int(1e6) # estimate read length from the first 1M reads
        counts = collections.Counter(map(len, itertools.islice(hts.FastqReader(fp_inp), n_reads)))
        with open(fp_out, 'w') as fh_out:
            #fh_out.write('%d-%d' % (min(counts.keys()), max(counts.keys())))
            fh_out.write('%d' % (max(counts.keys()),))

rule fqz_c_r1:
    input: pf('{bid}', '{step}', '.r1.fq.gz', '{prefix}')
    output: pf('{bid}', '{step}.c_r1', '.txt', '{prefix}')
    shell: "gunzip -c {input} | wc -l | awk '{{print $1/4}}' > {output}"

rule fqz_c_r2:
    input: pf('{bid}', '{step}', '.r2.fq.gz', '{prefix}')
    output: pf('{bid}', '{step}.c_r2', '.txt', '{prefix}')
    shell: "gunzip -c {input} | wc -l | awk '{{print $1/4}}' > {output}"

rule c:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c', '.txt', '{prefix}')
    shell: 'samtools view -c {input} > {output}'

rule c_r1:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_r1', '.txt', '{prefix}')
    shell: 'samtools view -f 64 -c {input} > {output}'

rule c_r2:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_r2', '.txt', '{prefix}')
    shell: 'samtools view -f 128 -c {input} > {output}'

rule c_aln: # -F 4 = (not) unmapped
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_aln', '.txt', '{prefix}')
    shell: 'samtools view -c -F 4 {input} > {output}'

rule c_aln_r1: # -f 64 = first in pair; -F 4 = (not) unmapped
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_aln_r1', '.txt', '{prefix}')
    shell: 'samtools view -c -f 64 -F 4 {input} > {output}'

rule c_q10:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_q10', '.txt', '{prefix}')
    shell: 'samtools view -c -q 10 {input} > {output}'

rule c_p10:
    # 67=properly paired, first in pair
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_p10', '.txt', '{prefix}')
    shell: 'samtools view -c -f 67 -F 12 -q 10 {input} > {output}'

rule c_qeq0:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_qeq0', '.txt', '{prefix}')
    shell: "samtools view {input} | awk '$5 == 0' - | wc -l > {output}"

rule c_chrM:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}')
    output: pf('{bid}', '{step}.c_chrM', '.txt', '{prefix}')
    shell: 'samtools view -c -F 4 {input[0]} chrM > {output}'

rule c_blacklist:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.c_blacklist', '.txt', '{prefix}')
    shell:
        'samtools view -c -F 4 -L shared/ce10_blacklist.bed {input[0]} > {output}'

rule c_atac:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.c_atac', '.txt', '{prefix}'),
    shell:
        'samtools view -c -L /mnt/home1/ahringer/jj374/lab/HTSProcessing/dm160615/dm160615.peaks_only.bed {input[0]} > {output}'

rule q10_keep:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.q10_keep', '.bam', '{prefix}'),
        pf('{bid}', '{step}.q10_drop', '.bam', '{prefix}'),
    shell:
        # 4 = read unmapped
        'samtools view -b -F 4 -q 10 -L shared/ce10_keep.bed -U {output[1]} {input} > {output[0]}'

rule p10_keep:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.p10_keep', '.bam', '{prefix}'),
        pf('{bid}', '{step}.p10_drop', '.bam', '{prefix}'),
    shell:
        #  3 = paired, mapped in proper pair
        # 12 = read unmapped, mate unmapped
        'samtools view -b -f 3 -F 12 -q 10 -L shared/ce10_keep.bed -U {output[1]} {input} > {output[0]}'

rule c_fr: # Count reads in FR-orientation
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_fr', '.txt', '{prefix}')
    shell: '''
        samtools view -h {input} \
            | awk '(substr($1,1,1) == "@") || (int($9) > 0)' \
            | samtools view -c - > {output}
    '''

rule c_fr_rmdup: # Count non-duplicate reads (pre-filtering for in FR-orientation)
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.c_fr_rmdup', '.txt', '{prefix}')
    shell: '''
        samtools view -h {input} \
            | awk '(substr($1,1,1) == "@") || (int($9) > 0)' \
            | samtools rmdup - - \
            | samtools view -c - > {output}
    '''

rule fsizes:
    input: pf('{bid}', '{step}', '.bam', '{prefix}')
    output: pf('{bid}', '{step}.fsizes', '.txt', '{prefix}')
    shell: '''
        samtools view {input} \
        | awk '(int($9) > 0) {{print $9}}' \
        | sort | uniq -c | sort -k2,2n | awk -v OFS='\t' '{{print $2,$1}}' \
        > {output}
    '''

rule rm_unmapped:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_unmapped', '.bam', '{prefix}'),
    shell:
        'samtools view -b -F 4 {input[0]} > {output}'

rule rm_unmapped_pe:
    '''
    -f 3 = read mapped & in a proper pair
    -F 12 = read unmapped | mate unmapped
    '''
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_unmapped_pe', '.bam', '{prefix}'),
    shell:
        'samtools view -b -f 3 -F 12 {input[0]} > {output}'

#rule rm_rRNA:
#    input:
#        pf('{bid}', '{step}', '.bam', '{prefix}'),
#        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
#    output:
#        pf('{bid}', '{step}.rm_rRNA', '.bam', '{prefix}'),
#    shell:
#        'samtools view -u -F 4 {input[0]} | samtools view -b -L shared/WS253_ce10.rRNA.bed -U {output} - > /dev/null'

rule rm_rRNA_broad:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_rRNA_broad', '.bam', '{prefix}'),
    shell:
        'samtools view -u {input[0]} | samtools view -b -L shared/WS253_ce10.rRNA_broad.bed -U {output} - > /dev/null'

rule rm_blacklist:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_blacklist', '.bam', '{prefix}'),
    shell:
        'samtools view -b -L shared/ce10_blacklist.bed -U {output} {input[0]} > /dev/null'

rule rm_blacklist_ce11:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_blacklist_ce11', '.bam', '{prefix}'),
    shell:
        'samtools view -b -L shared/ce11_blacklist.bed -U {output} {input[0]} > /dev/null'

rule rm_chrM:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}')
    output:
        pf('{bid}', '{step}.rm_chrM', '.bam', '{prefix}'),
    shell:
        'samtools view -b {input[0]} chrI chrII chrIII chrIV chrV chrX > {output}'

rule rm_q10:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_q10', '.bam', '{prefix}'),
    shell:
        'samtools view -b -q 10 {input} > {output}'

rule rm_q10_reversed:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_q10_reversed', '.bam', '{prefix}'),
    shell:
        'samtools view -b -q 10 {input} -U {output} > /dev/null'
