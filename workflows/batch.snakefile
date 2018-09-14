
def files_from_dataset(r):
    # k = final processing label, v = processing internal details
    for k, v in config['batch'][r['Assay']][r['Organism']].items():
        scale_ = v['Scale']
        genome_ = v['Genome']
        res_ = v['Resolution']
        filetype_ = v['Filetype']
        file_ = '%s.%s.%s' % (r['Dataset'], k, v['Filetype'])
        raw_ = pf(r['Dataset'], v['_step'], v['_suffix'], v['_prefix'])

        r1_isfile = os.path.isfile('samples/%s.r1.fq.gz' % (r['Dataset'],))
        r2_isfile = os.path.isfile('samples/%s.r2.fq.gz' % (r['Dataset'],))
        if (v['_input_type'] == 'paired-end') and r1_isfile and r2_isfile:
            yield(r['Dataset'], k, scale_, genome_, res_, filetype_, file_, raw_)

        elif (v['_input_type'] == 'single-end') and r1_isfile:
            yield(r['Dataset'], k, scale_, genome_, res_, filetype_, file_, raw_)

df_batch_datasets = pd.read_csv('batch/datasets.txt', delim_whitespace=True)
cols_batch_files = ['Dataset', 'Processing', 'Scale', 'Genome', 'Resolution', 'Filetype', 'File', 'Raw']
gen_batch_files = (pd.DataFrame.from_records(files_from_dataset(r), columns=cols_batch_files) for i, r in df_batch_datasets.iterrows())
df_batch_files = pd.concat([* gen_batch_files ], ignore_index=True, axis=0)
d_batch_file_raw = {file: raw for file, raw in zip(df_batch_files['File'], df_batch_files['Raw'])}

rule batch_processed_files:
    input:
        'batch/datasets.txt',
    output:
        'batch/processed_files.txt',
    run:
        assert input[0] == 'batch/datasets.txt'
        df_batch_files[['Dataset', 'Processing', 'Scale', 'Genome', 'Resolution', 'Filetype', 'File']].to_csv(output[0], sep='\t', index=False)

rule batch_processed_file:
    input:
        lambda wildcards: d_batch_file_raw[wildcards.file],
    output:
        'batch/processed_files/{file}',
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule batch:
    input:
        ['batch/processed_files/%s' % (file_,) for file_ in df_batch_files['File']] \
        + ['batch/processed_files.txt']

rule batch_debug:
    input:
        #pf('HS645_ATAC_TG_CB2', 'tg_se.bwa_se_CB4.rm_unmapped', '.bam', 'atac'),
        pf('HS645_ATAC_TG_CB2', 'tg_pe.bwa_pe_CB4.rm_unmapped_pe.rm_q10.sample_spr.macs2_pe_lt200.yapc_spr_pe', '.tsv', 'atac'),
        pf('HS645_ATAC_TG_CB2', 'tg_se.bwa_se_CB4.rm_unmapped.rm_q10.sample_spr.macs2_se_extsize150_shiftm75_keepdup_all.yapc_spr_se', '.tsv', 'atac'),

rule batch_single: # Special rule to process one file at a time for jadb integration
    input:
        'batch/processed_files/%s.%s.%s' % (config.get('Dataset', 'NA'), config.get('Processing', 'NA'), config.get('Filetype', 'NA'))
