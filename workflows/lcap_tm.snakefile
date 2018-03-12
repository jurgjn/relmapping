# Replace with *actual* sample names, this example (when run using `smj 4 lcap_tm`) would process:
#     samples/HS569_lcRNA_glp1_P22_D14.r1.fq.gz
#     samples/HS569_lcRNA_glp1_P22_D14.r2.fq.gz
#     samples/HS569_lcRNA_glp1_P39_D14.r1.fq.gz
#     samples/HS569_lcRNA_glp1_P39_D14.r2.fq.gz
#
l_lcap_tm = ['HS569_lcRNA_glp1_P22_D14', 'HS569_lcRNA_glp1_P22_D14']

rule lcap_tm:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'lcap_tm'), bid=l_lcap_tm),
        expand(pf('{bid}', 'c_r2', '.txt', 'lcap_tm'), bid=l_lcap_tm),
        # Coverage tracks, q10
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_fwd', '.bw', 'lcap_tm'), bid=l_lcap_tm),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.filled_rev', '.bw', 'lcap_tm'), bid=l_lcap_tm),
