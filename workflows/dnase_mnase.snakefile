def make_dnase_mnase_samples():
    fp_inp = 'workflows/Worm Accessibility Data Sets - accessibility.tsv'
    fp_out = 'dnase_mnase/dnase_mnase_samples.tsv'
    df_ = pd.read_csv(fp_inp, sep='\t').query('(Strain == "N2") & (Enzyme == "dnase" | Enzyme == "mnase")').sort_values(['Enzyme', 'Stage'])
    df_ = df_.loc[df_['Bioinformatics ID(s)'] == df_['Bioinformatics ID(s)']] # All sequenced libraries
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA0-N2-dnase-test')] # Early test data, partial digestion course
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA11-N2-dnase-S2.1xtraHS082')] # Pippen/beads tests
    df_ = df_.loc[~(df_['Library series ID'] == 'L2-JAX-N2-dnase-SX')] # Pippen/beads tests
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA12-N2-mnase-S2')] # Djem's time course
    df_.reset_index(drop=True).to_csv(fp_out, sep='\t', index=None)
    return df_

def l_dnase_mnase_samples():
    if not(os.path.isfile('dnase_mnase/dnase_mnase_samples.tsv')):
        make_dnase_mnase_samples()
    l_bid_dnase_mnase = []
    for l_bid_ in pd.read_csv('dnase_mnase/dnase_mnase_samples.tsv', sep='\t')['Bioinformatics ID(s)']:
        for bid in l_bid_.split(';'):
            l_bid_dnase_mnase.append(bid)
    return l_bid_dnase_mnase

def dnase_mnase_label_bid_make():
    dnase_mnase_label_bid = {}
    df_ = make_dnase_mnase_samples()#.head(5)
    for i, r in df_.iterrows():
        libid = r['Library series ID'].rstrip()
        l_bid = r['Bioinformatics ID(s)'].rstrip().split(';')
        l_annot = r['Digestion conditions (homogenized)'].rstrip().split(';')
        for bid, annot in zip(l_bid, l_annot):
            label = '%s/%s_%s_%s' % (libid, libid, annot, bid)
            #print(label, bid, libid, annot)
            #print(pf(label, 'spmr_lt300', '.bw', 'dnase_mnase_labelled'))
            dnase_mnase_label_bid[label] = bid
    return dnase_mnase_label_bid

d_dnase_mnase_label_bid = dnase_mnase_label_bid_make()

rule dnase_mnase_labelled_spmr_lt300:
    input:
        lambda wildcards: pf('%s', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'dnase_mnase') % (d_dnase_mnase_label_bid[wildcards.label],),
    output:
        pf('{label}', 'spmr_lt300', '.bw', 'dnase_mnase_labelled'),
    shell:
        '''
        ln -s `pwd`/{input} `pwd`/{output}
        touch -h `pwd`/{output}
        '''

rule dnase_mnase:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'c_r2', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'tg_pe.bwa_pe', '.bam', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'dnase_mnase'), bid=l_dnase_mnase_samples()),

rule dnase_mnase_labelled:
    input:
        expand(pf('{label}', 'spmr_lt300', '.bw', 'dnase_mnase_labelled'), label=d_dnase_mnase_label_bid.keys()),
