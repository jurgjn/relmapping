def scap_rm_exonic_step(df_exons_gtf_, fp_inp, fp_out, fc_th, v=False, window_width=50, exon_flank=5):
    """
    Adhoc filtering of short cap read coverage to discard likely false-positive signal seen at (some) highly transcribed exons:
    1. Flatten list of gtf-formatted exon annotations by taking the union of all exons (as provided by <df_exons_gtf_>).
        - exon boundaries are expanded <exon_flank>bp on both sides, as false positive short cap signal can also extend slightly beyond the annotated exon boundaries
    2. Create "strides" of <window_width>bp windows, falling within the union of all exons ("unionised exons").
    3. For every stride, calculate total short cap signal (BigWig file specified by <fp_inp>) within the stride.
    4. For every gene, determine expected stride background coverage by taking the average of all strides completely overlapping an exon of particular gene.
    5. Keep exonic signal only within strides where the fold change -- ratio of signal of the particular stride divided by the gene-specific background -- is above fc_th.
    6. Store output in <fp_out> as a "bedGraph"-style BigWig
    """
    # Union of all exons
    df_exon = df_exons_gtf_.copy()
    df_exon['start'] = df_exon['start'] - 1 - exon_flank # gtf-to-bed coordinate change
    df_exon['end'] = df_exon['end'] + exon_flank
    bt_uxon = BedTool.from_dataframe(df_exon[['chrom', 'start', 'end', 'gene_id']]).merge(c='4', o='distinct')
    print(bt_uxon.count(), 'number of "union exon records"')

    # df_strides <- exon-overlapping strides
    names_ = ['exon_chrom', 'exon_start', 'exon_end', 'gene_id', 'stride_chrom', 'stride_start', 'stride_end']
    df_uxon = pd.read_csv(bt_uxon.fn, sep='\t', names=names_[:4])
    l_exon_chrom = []
    l_exon_start = []
    l_exon_end = []
    l_gene_id = []
    l_stride_chrom = []
    l_stride_start = []
    l_stride_end = []
    for i, r in df_uxon.iterrows():
        for stride_start in range(r['exon_start'], r['exon_end'] - window_width + 1):
            stride_end = stride_start + window_width
            l_exon_chrom.append(r['exon_chrom'])
            l_exon_start.append(r['exon_start'])
            l_exon_end.append(r['exon_end'])
            l_gene_id.append(r['gene_id'])
            l_stride_chrom.append(r['exon_chrom'])
            l_stride_start.append(stride_start)
            l_stride_end.append(stride_end)
    df_strides = pd.DataFrame()
    df_strides['exon_chrom'] = l_exon_chrom
    df_strides['exon_start'] = l_exon_start
    df_strides['exon_end'] = l_exon_end
    df_strides['gene_id'] = l_gene_id
    df_strides['stride_chrom'] = l_stride_chrom
    df_strides['stride_start'] = l_stride_start
    df_strides['stride_end'] = l_stride_end
    print(len(df_strides), 'exon-overlapping strides')

    # ga_inp <- unfiltered input coverage
    # ga_out <- output array initially full copy of input; later NA-d/restored based on individual strides and ga_inp
    fh_inp = pyBigWig.open(fp_inp)
    chroms = collections.OrderedDict()
    for chrom in sorted(fh_inp.chroms().keys()):
        chroms[chrom] = fh_inp.chroms()[chrom]
    ga_inp = collections.OrderedDict()
    ga_out = collections.OrderedDict()
    for chrom, size in chroms.items():
        print(chrom, size)
        ga_inp[chrom] = np.array(fh_inp.values(chrom, 0, size))
        ga_out[chrom] = np.copy(ga_inp[chrom])
    fh_inp.close()
    
    # scap_nansum <- read count within current stride
    df_strides['scap_nansum'] = list(map(lambda chrom, start, end: np.nansum(ga_inp[chrom][start:end]), 
         df_strides['stride_chrom'], df_strides['stride_start'], df_strides['stride_end']))

    # Browse stride-to-gene assignments
    if v: df_strides[['stride_chrom', 'stride_start', 'stride_end', 'gene_id', 'scap_nansum']].to_csv(fp_out + '_df_strides.bed', header=False, index=False, sep='\t')

    # scap_nansum_mean <- mean stride read count by gene_id
    df_strides_bg = df_strides.merge(df_strides.groupby('gene_id').agg({'scap_nansum': np.mean}),
                 left_on='gene_id', right_index=True, suffixes=('', '_mean'))
    if v: df_strides_bg.to_csv(fp_out + '_df_strides_bg.tsv', header=True, index=False, sep='\t')

    # scap_nansum_fc <- filtering fold change
    df_strides_bg['scap_nansum_fc'] = df_strides_bg['scap_nansum'] / df_strides_bg['scap_nansum_mean']

    # Query for sub-threshold ("failing") strides / above-threshold ("passing") strides
    q_neg = 'scap_nansum_fc <= %s' % (fc_th,)
    q_pos = 'scap_nansum_fc > %s' % (fc_th,)
    
    # Browse above/below-threshold strides
    if v: df_strides_bg.query(q_neg)[['stride_chrom', 'stride_start', 'stride_end']].to_csv(fp_out + '_q_neg.bed', header=False, index=False, sep='\t')
    if v: df_strides_bg.query(q_pos)[['stride_chrom', 'stride_start', 'stride_end']].to_csv(fp_out + '_q_pos.bed', header=False, index=False, sep='\t')

    print('%d of %d below-threshold strides: set to NA...' % (len(df_strides_bg.query(q_neg)), len(df_strides_bg)))
    for i, r in df_strides_bg.query(q_neg).iterrows():
        # failing strides set to NA (not zero), as true txn initiation is masked by exonic signal; hence unknown
        ga_out[r['stride_chrom']][r['stride_start']:r['stride_end']] = float('nan')

    print('%d of %d above-threshold strides: restore to input signal...' % (len(df_strides_bg.query(q_pos)), len(df_strides_bg)))
    for i, r in df_strides_bg.query(q_pos).iterrows():
        ga_out[r['stride_chrom']][r['stride_start']:r['stride_end']] = ga_inp[r['stride_chrom']][r['stride_start']:r['stride_end']]

    # Write output file of interval-based bigWig based on array signal stored in ga_out
    fh_out = pyBigWig.open(fp_out, 'w')
    fh_out.addHeader([(chrom, size) for chrom, size in chroms.items()])
    for chrom, size in chroms.items():
        #fh_out.addEntries(chrom, 0, values=list(ga_out[chrom]), span=1, step=1)
        #^ NA values in bedGraph (vs fixedStep) do not cause scale to be set to 'inf' in IGV...
        l_chrom = []
        l_start = []
        l_ends = []
        l_vals = []
        pos_iv = 0
        val_iv = float('nan')
        for pos_i, val_i in enumerate(list(ga_out[chrom]) + [float('nan')]):
            # change in value
            if (val_iv != val_i):
                if (val_iv == val_iv):
                    l_chrom.append(chrom)
                    l_start.append(pos_iv)
                    l_ends.append(pos_i)
                    l_vals.append(val_iv)
                # reset entry
                pos_iv = pos_i
                val_iv = val_i
        fh_out.addEntries(l_chrom, l_start, ends=l_ends, values=l_vals)
    fh_out.close()

rule scap_fwd_rm_exonic_fc0:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic', '_fc0.bw', 'scap'),
    run:
        fc_th = 0
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_rev_rm_exonic_fc0:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic', '_fc0.bw', 'scap'),
    run:
        fc_th = 0
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_fwd_rm_exonic_fc3:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic', '_fc3.bw', 'scap'),
    run:
        fc_th = 3
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_rev_rm_exonic_fc3:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic', '_fc3.bw', 'scap'),
    run:
        fc_th = 3
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_fwd_rm_exonic_fc5:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic', '_fc5.bw', 'scap'),
    run:
        fc_th = 5
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_rev_rm_exonic_fc5:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic', '_fc5.bw', 'scap'),
    run:
        fc_th = 5
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_fwd_rm_exonic_fc10:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic', '_fc10.bw', 'scap'),
    run:
        fc_th = 10
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_rev_rm_exonic_fc10:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic', '_fc10.bw', 'scap'),
    run:
        fc_th = 10
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_fwd_rm_exonic_fc20:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic', '_fc20.bw', 'scap'),
    run:
        fc_th = 20
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_rev_rm_exonic_fc20:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic', '_fc20.bw', 'scap'),
    run:
        fc_th = 20
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th)

rule scap_fwd_rm_exonic_w25_fc0:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic_w25', '_fc0.bw', 'scap'),
    run:
        fc_th = 0
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_rev_rm_exonic_w25_fc0:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic_w25', '_fc0.bw', 'scap'),
    run:
        fc_th = 0
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_fwd_rm_exonic_w25_fc3:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic_w25', '_fc3.bw', 'scap'),
    run:
        fc_th = 3
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_rev_rm_exonic_w25_fc3:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic_w25', '_fc3.bw', 'scap'),
    run:
        fc_th = 3
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_fwd_rm_exonic_w25_fc5:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic_w25', '_fc5.bw', 'scap'),
    run:
        fc_th = 5
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_rev_rm_exonic_w25_fc5:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic_w25', '_fc5.bw', 'scap'),
    run:
        fc_th = 5
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_fwd_rm_exonic_w25_fc10:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic_w25', '_fc10.bw', 'scap'),
    run:
        fc_th = 10
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_rev_rm_exonic_w25_fc10:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic_w25', '_fc10.bw', 'scap'),
    run:
        fc_th = 10
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_fwd_rm_exonic_w25_fc20:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic_w25', '_fc20.bw', 'scap'),
    run:
        fc_th = 20
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_rev_rm_exonic_w25_fc20:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic_w25', '_fc20.bw', 'scap'),
    run:
        fc_th = 20
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), fc_th=fc_th, window_width=25)

rule scap_fwd_rm_exonic:
    input:
        pf('{bid}', '{step}.firstbp_fwd', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_fwd.rm_exonic', '.bw', 'scap'),
    run:
        strand = "+"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), window_width=50, fc_th=5)

rule scap_rev_rm_exonic:
    input:
        pf('{bid}', '{step}.firstbp_rev', '.bw', 'scap'),
    output:
        pf('{bid}', '{step}.firstbp_rev.rm_exonic', '.bw', 'scap'),
    run:
        strand = "-"
        df_raw = yp.read_wbgtf('WS260_ce10/WS260_ce10.transcripts.gtf.gz', parse_attr=False)
        df_exon = yp.df_gfftags_unpack(df_raw, name='attribute')\
            .query('(feature == "exon") & (gene_biotype == "protein_coding")')\
            .query('strand == "%s"' % (strand,))
        scap_rm_exonic_step(df_exon, str(input), str(output), window_width=50, fc_th=5)

rule scap_rm_exonic_browse:
    input:
        # scap_rm_exonic defaults: window_size=50, fc=5, scap541_emb_l3_ya ("nearly-TAP-only" data)
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_fwd.rm_exonic', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.rm_non_coding.firstbp_rev.rm_exonic', '.bw', 'scap'), bid=['scap541_emb_l3_ya']),
        # window size 50bp
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic', '_fc0.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic', '_fc0.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic', '_fc3.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic', '_fc3.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic', '_fc5.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic', '_fc5.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic', '_fc10.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic', '_fc10.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic', '_fc20.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic', '_fc20.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        # window size 25bp
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic_w25', '_fc0.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic_w25', '_fc0.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic_w25', '_fc3.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic_w25', '_fc3.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic_w25', '_fc5.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic_w25', '_fc5.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic_w25', '_fc10.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic_w25', '_fc10.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_fwd.rm_exonic_w25', '_fc20.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
        #expand(pf('{bid}', 'tg_se.bwa_se.rm_unmapped.rm_chrM.rm_blacklist.rm_q10.firstbp_rev.rm_exonic_w25', '_fc20.bw', 'scap'), bid=['scap541_emb_l3_ya', 'scap541_emb_to_ya']),
