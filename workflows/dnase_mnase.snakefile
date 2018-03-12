def make_dnase_mnase_samples():
    fp_inp = 'workflows/Worm Accessibility Data Sets - accessibility.tsv',
    fp_out = 'dnase_mnase/dnase_mnase_samples.tsv',
    df_ = pd.read_csv(fp_inp, sep='\t').query('(Strain == "N2") & (Enzyme == "dnase" | Enzyme == "mnase")').sort_values(['Enzyme', 'Stage'])
    df_ = df_.loc[df_['Bioinformatics ID(s)'] == df_['Bioinformatics ID(s)']] # All sequenced libraries
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA0-N2-dnase-test')] # Early test data, partial digestion course
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA11-N2-dnase-S2.1xtraHS082')] # Pippen/beads tests
    df_ = df_.loc[~(df_['Library series ID'] == 'L2-JAX-N2-dnase-SX')] # Pippen/beads tests
    df_ = df_.loc[~(df_['Library series ID'] == 'EMB-JA12-N2-mnase-S2')] # Djem's time course
    df_.reset_index(drop=True).to_csv(fp_out, sep='\t', index=None)

def l_dnase_mnase_samples():
    if not(os.path.isfile('dnase_mnase/dnase_mnase_samples.tsv')):
        make_dnase_mnase_samples()
    l_bid_dnase_mnase = []
    for l_bid_ in pd.read_csv('dnase_mnase/dnase_mnase_samples.tsv', sep='\t')['Bioinformatics ID(s)']:
        for bid in l_bid_.split(';'):
            l_bid_dnase_mnase.append(bid)
    return l_bid_dnase_mnase

rule dnase_mnase:
    input:
        expand(pf('{bid}', 'c_r1', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'c_r2', '.txt', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'tg_pe.bwa_pe', '.bam', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
        expand(pf('{bid}', 'tg_pe.bwa_pe.rm_unmapped_pe.rm_chrM.rm_blacklist.rm_q10.macs2_pe_lt300', '_treat_pileup.bw', 'dnase_mnase'), bid=l_dnase_mnase_samples()),
