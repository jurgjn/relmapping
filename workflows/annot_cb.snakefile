df_annot_cb_data = pd.DataFrame()

for k, v in config['annot_cb']['atac_samples'].items():
    sample_rep1 = 'annot_cb_atac_%(k)s_rep1' % locals()
    sample_rep2 = 'annot_cb_atac_%(k)s_rep2' % locals()
    l_fp_rep1 = v['rep1']
    l_fp_rep2 = v['rep2']

    df_annot_cb_data = df_annot_cb_data.append(pd.Series({
        'sample': sample_rep1,
        'condition': k,
        'assay': 'atac',
        'r1_raw': l_fp_rep1,
        'r1_pooled': 'samples/%(sample_rep1)s.r1.fq.gz' % locals(),
    }), ignore_index=True)
    df_annot_cb_data = df_annot_cb_data.append(pd.Series({
        'sample': sample_rep2,
        'condition': k,
        'assay': 'atac',
        'r1_raw': l_fp_rep2,
        'r1_pooled': 'samples/%(sample_rep2)s.r1.fq.gz' % locals(),
    }), ignore_index=True)


for k, v in config['annot_cb']['lcap_samples'].items():
    sample_rep1 = 'annot_cb_lcap_%(k)s_rep1' % locals()
    sample_rep2 = 'annot_cb_lcap_%(k)s_rep2' % locals()
    l_fp_rep1_read1 = v['rep1_read1']
    l_fp_rep1_read2 = v['rep1_read2']
    l_fp_rep2_read1 = v['rep2_read1']
    l_fp_rep2_read2 = v['rep2_read2']

    df_annot_cb_data = df_annot_cb_data.append(pd.Series({
        'sample': sample_rep1,
        'condition': k,
        'assay': 'lcap',
        'r1_raw': l_fp_rep1_read1,
        'r2_raw': l_fp_rep1_read2,
        'r1_pooled': 'samples/%(sample_rep1)s.r1.fq.gz' % locals(),
        'r2_pooled': 'samples/%(sample_rep1)s.r2.fq.gz' % locals(),
    }), ignore_index=True)
    df_annot_cb_data = df_annot_cb_data.append(pd.Series({
        'sample': sample_rep2,
        'condition': k,
        'assay': 'lcap',
        'r1_raw': l_fp_rep2_read1,
        'r2_raw': l_fp_rep2_read2,
        'r1_pooled': 'samples/%(sample_rep2)s.r1.fq.gz' % locals(),
        'r2_pooled': 'samples/%(sample_rep2)s.r2.fq.gz' % locals(),
    }), ignore_index=True)

df_annot_cb_data = df_annot_cb_data[['sample', 'assay', 'condition', 'r1_pooled', 'r2_pooled', 'r1_raw', 'r2_raw']]
df_annot_cb_data = df_annot_cb_data.set_index('sample')

def samples_cb_r1_input(wildcards):
    return df_annot_cb_data.query('index == "%s"' % (wildcards.sample,))['r1_raw'][0]

def samples_cb_r2_input(wildcards):
    return df_annot_cb_data.query('index == "%s"' % (wildcards.sample,))['r2_raw'][0]

rule samples_cb_r1:
    input:
        samples_cb_r1_input,
    output:
        'samples/{sample}.r1.fq.gz'
    run:
        arg_input = ' '.join(input)
        shell('cat %(arg_input)s > %(output)s' % locals())

rule samples_cb_r2:
    input:
        samples_cb_r2_input,
    output:
        'samples/{sample}.r2.fq.gz'
    run:
        arg_input = ' '.join(input)
        shell('cat %(arg_input)s > %(output)s' % locals())

rule annot_cb:
    input:
        #expand('samples/{sample}.r1.fq.gz', sample=[*df_annot_cb_data.query('assay == "atac"').index]),
        #expand('samples/{sample}.r1.fq.gz', sample=[*df_annot_cb_data.query('assay == "lcap"').index]),
        #expand('samples/{sample}.r2.fq.gz', sample=[*df_annot_cb_data.query('assay == "lcap"').index]),
        expand(pf('{sample}', 'tg_se.bwa_se_CB4', '.bam', 'atac'), sample=[*df_annot_cb_data.query('assay == "atac"').index]),
        expand(pf('{sample}', 'tg_pe.bwa_pe_CB4', '.bam', 'lcap'), sample=[*df_annot_cb_data.query('assay == "lcap"').index]),
