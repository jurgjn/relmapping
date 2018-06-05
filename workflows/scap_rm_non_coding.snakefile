rule scap_rm_non_coding:
    """
    https://doi.org/10.1101/gr.153668.112
    > Finally, clusters overlapping rRNA, tRNA, miRNA, snRNA, snoRNA, or snlRNA genes (from Ensembl release 61/WS220) were excluded from the set. 
    
    Distribution of WS260_ce10 "gene records" by gene_biotype; asterisk indicates inclusion in non_coding filtering set.
    protein_coding    20210
    piRNA             15364
    ncRNA              7741
    pseudogene         1791
    tRNA                612*
    snoRNA              345*
    miRNA               257*
    lincRNA             172
    snRNA               130*
    antisense           100
    rRNA                 20*
    """
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_non_coding', '.bam', '{prefix}'),
    shell:
        'samtools view -b -L WS260_ce10/WS260_ce10.transcripts.non_coding.bed -U {output} {input[0]} > /dev/null'

rule scap_rm_non_coding_ce11:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}'),
    output:
        pf('{bid}', '{step}.rm_non_coding_ce11', '.bam', '{prefix}'),
    shell:
        'samtools view -b -L annot_ce11/canonical_geneset/WS260_ce11.transcripts.non_coding.bed -U {output} {input[0]} > /dev/null'

rule scap_rm_non_coding_browse:
    input:
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_fwd', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_rev', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
