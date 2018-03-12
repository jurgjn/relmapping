"""
https://doi.org/10.1016/j.molcel.2012.01.017
> peak calling was done using MACS (Zhang et al., 2008). For DRIP-seq, peaks were called using all mapped reads enforcing a greater than 10-fold enrichment above both input and RNase H-pre-treated control data sets.

https://doi.org/10.1101/gr.158436.113
> Peak calling was first performed by MACS 1.4.2 (Zhang et al. 2008) using the matching input library as control. DRIP peaks were further assigned onto restriction fragments using custom Java and Perl scripts. Regions common to DRIP 1 and DRIP 2 were considered consensus DRIP-seq peaks.

https://doi.org/10.1038/ng.3672
> Read density along the genome was calculated by tiling the genome into 200-bp windows (non-overlapping) and counting the number of sequence fragments within each window, using the qCount function of the QuasR package (see URLs). To compensate for differences in the read depths of the various libraries, we divided each sample by the total number of mapped reads and multiplied by the average library size. Log2 expression levels were calculated after adding a pseudocount of 1 (y = log2(x + 1)). ChIP-seq signals are displayed as average enrichment of IP − input (log2).
> Figure 5: met-2 set-25 worms accumulate RNA:DNA hybrids at repeat elements.
> (f) DRIP-seq example showing the R-loop signal over a RE cluster. The IP signal was normalized to the input and the RNase H control values were subtracted.
> Supplementary Figure 5: R-loop accumulation on repeat elements (REs).
> (b) DRIP (log2) - RNaseH
> (e) DRIP-seq examples showing the R-loop signal over two RE clusters. The ChIP signal from antibody S9.6, which is specific for RNA:DNA hybrids, was normalized to input, and the RNase H control values were subtracted.
"""

rule macs2_DRIPseq_vs_input:
    input:
        pf('DRIPseq_wt', '{step}', '.bam', '{prefix}'),
        pf('DRIPseq_input_wt', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_control_lambda.bdg', '{prefix}')), # [0]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_control_lambda.bw', '{prefix}'), # [1]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_peaks.narrowPeak', '{prefix}'), # [2]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_peaks.xls', '{prefix}'), # [3]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_summits.bed', '{prefix}'), # [4]
        temp(pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_treat_pileup.bdg', '{prefix}')), # [5]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_treat_pileup.bw', '{prefix}'), # [6]
        temp(pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_FE.bdg', '{prefix}')), # [7]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_FE.bw', '{prefix}'), # [8]
        temp(pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_logLR.bdg', '{prefix}')), # [9]
        pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '_logLR.bw', '{prefix}'), # [10]
        temp(pf('DRIPseq_vs_input', '{step}.macs2_DRIPseq', '.chroms', '{prefix}')), # [11]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --gsize ce --bdg --SPMR'
    threads: 1
    shell: '''
        macs2 callpeak {params.macs2_args} \
            --treatment {input[0]} \
            --control {input[1]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_DRIPseq --name DRIPseq_vs_input.{wildcards.step}.macs2_DRIPseq
        macs2 bdgcmp -t {output[5]} -c {output[0]} -o {output[7]} -m FE
        macs2 bdgcmp -t {output[5]} -c {output[0]} -o {output[9]} -m logLR
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[11]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        sort -k1,1 -k2,2n -o {output[5]} {output[5]}
        sort -k1,1 -k2,2n -o {output[7]} {output[7]}
        sort -k1,1 -k2,2n -o {output[9]} {output[9]}
        bedGraphToBigWig {output[0]} {output[11]} {output[1]}
        bedGraphToBigWig {output[5]} {output[11]} {output[6]}
        bedGraphToBigWig {output[7]} {output[11]} {output[8]}
        bedGraphToBigWig {output[9]} {output[11]} {output[10]}
    '''

rule macs2_DRIPseq_vs_RNAseH:
    input:
        pf('DRIPseq_wt', '{step}', '.bam', '{prefix}'),
        pf('DRIPseq_RNAseH_wt', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_control_lambda.bdg', '{prefix}')), # [0]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_control_lambda.bw', '{prefix}'), # [1]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_peaks.narrowPeak', '{prefix}'), # [2]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_peaks.xls', '{prefix}'), # [3]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_summits.bed', '{prefix}'), # [4]
        temp(pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_treat_pileup.bdg', '{prefix}')), # [5]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_treat_pileup.bw', '{prefix}'), # [6]
        temp(pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_FE.bdg', '{prefix}')), # [7]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_FE.bw', '{prefix}'), # [8]
        temp(pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_logLR.bdg', '{prefix}')), # [9]
        pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '_logLR.bw', '{prefix}'), # [10]
        temp(pf('DRIPseq_vs_RNAseH', '{step}.macs2_DRIPseq', '.chroms', '{prefix}')), # [11]
    conda:
        'envs/py2.yaml'
    params:
        macs2_args = '--format BAMPE --gsize ce --bdg --SPMR'
    threads: 1
    shell: '''
        macs2 callpeak {params.macs2_args} \
            --treatment {input[0]} \
            --control {input[1]} \
            --outdir {wildcards.prefix}/{wildcards.step}.macs2_DRIPseq --name DRIPseq_vs_RNAseH.{wildcards.step}.macs2_DRIPseq
        macs2 bdgcmp -t {output[5]} -c {output[0]} -o {output[7]} -m FE
        macs2 bdgcmp -t {output[5]} -c {output[0]} -o {output[9]} -m logLR
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[11]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        sort -k1,1 -k2,2n -o {output[5]} {output[5]}
        sort -k1,1 -k2,2n -o {output[7]} {output[7]}
        sort -k1,1 -k2,2n -o {output[9]} {output[9]}
        bedGraphToBigWig {output[0]} {output[11]} {output[1]}
        bedGraphToBigWig {output[5]} {output[11]} {output[6]}
        bedGraphToBigWig {output[7]} {output[11]} {output[8]}
        bedGraphToBigWig {output[9]} {output[11]} {output[10]}
    '''

l_drip = ['DRIPseq_input_wt', 'DRIPseq_wt', 'DRIPseq_RNAseH_wt']

rule centroids_lt300:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
    output:
        temp(pf('{bid}', '{step}.centroids_lt300', '.chroms', '{prefix}')),
        temp(pf('{bid}', '{step}.centroids_lt300', '.bed', '{prefix}')),
        temp(pf('{bid}', '{step}.centroids_lt300', '.bdg', '{prefix}')),
        pf('{bid}', '{step}.centroids_lt300', '.bw', '{prefix}'),
    threads: 1
    shell: '''
        # create .chroms file from alignment
        samtools view -H {input[0]} | grep '@SQ' | awk -F'\\t' -v OFS='\\t' '{{print substr($2, 4), substr($3, 4)}}' > {output[0]}

        # select paired-end, properly paired reads
        # FR-orientation records, size-select <300bp, convert to centroids
        samtools view -f 3 --threads {threads} {input[0]} \
        | awk -F'\\t' -v OFS='\\t' '((int($9) > 0) && (int($9) <= 300)) {{print $3,int(($4-1+$4+$9-1)/2),int(($4-1+$4+$9-1)/2)+1}}' \
        > {output[1]}

        # sort centroid records by coordinate:
        sort -k1,1 -k2,2n -o {output[1]} {output[1]}

        # convert centroids to coverage:
        bedtools genomecov -bg -i {output[1]} -g {output[0]} > {output[2]}

        # convert coverage from bedGraph to bigWig:
        bedGraphToBigWig {output[2]} {output[0]} {output[3]}
    '''

l_dnase = ['HS089+090+092+099-01', 'HS089+090+092+099-02', 'HS089+090+092+099-03', 'HS089+090+092+099-04', 'HS089+090+092+099-05', 'HS089+090+092+099-06', 'HS089+090+092+099-07', 'HS089+090+092+099-08', 'HS089+090+092+099-09', 'HS089+090+092+099-10',]

rule adhoc:
    input:
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'adhoc'), bid=l_drip),
        pf('DRIPseq_vs_input', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.macs2_DRIPseq', '_peaks.narrowPeak', 'adhoc'),
        pf('DRIPseq_vs_RNAseH', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.macs2_DRIPseq', '_peaks.narrowPeak', 'adhoc'),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.macs2_pe_lt300', '_peaks.narrowPeak', 'adhoc'), bid=l_dnase),
        expand(pf('{bid}', 'trim20.bwa_pe.rm_unmapped_pe.rm_chrM.rm_rRNA_broad.rm_blacklist.rm_q10.centroids_lt300', '.bw', 'adhoc'), bid=l_dnase),
