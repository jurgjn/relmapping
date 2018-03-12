# fig_chip_stages(): no of modENCODE ChIP data sets per stage
#dp_finalPk = lp('.wget/encodedcc.stanford.edu/ftp/modENCODE_VS_ENCODE/Regulation/Worm/peakCalls/finalPk')
#dp_foldChange = lp('.wget/encodedcc.stanford.edu/ftp/modENCODE_VS_ENCODE/Regulation/Worm/signal/foldChange')

# list all data; extract stage name
def stage_from_fn(fn):
    # Calculate stage from file name
    stage = []
    if "_EE_" in fn: stage.append("EE") # early embryo
    if "_EM_" in fn: stage.append("EM") # mixed embryo?
    if "_LE_" in fn: stage.append("LE") # late embryo
    if "_L1_" in fn: stage.append("L1")
    if "_L2_" in fn: stage.append("L2")
    if "_L3_" in fn: stage.append("L3")
    if "_L4_" in fn: stage.append("L4")
    if "_YA_" in fn: stage.append("YA")
    if "_LY_" in fn: stage.append("LY") # ???
    if "_S1_" in fn: stage.append("S1") # L1 arrest?
    if "_D4_" in fn: stage.append("D4") # dauer larvae?
    if len(stage) != 1: # sanity check
        raise Exception("Unknown stage: stage=%(stage)s fn=%(fn)s" % locals())
    return stage[0]

def finalPk_from_index(index):
    #print l_finalPk
    r = [finalPk for finalPk in l_finalPk if index in finalPk]
    assert len(r) <= 1
    if len(r) == 1:
        return r[0]
    else:
        return '.'

#l_finalPk = list(glob.glob('%(dp_finalPk)s/*.narrowPeak.gz' % locals()))
#l_foldChange = list(glob.glob('%(dp_foldChange)s/*.bw' % locals()))

#print(len(l_finalPk), len(l_foldChange))

# remove RNA Pol-II/Pol-III and TBP-1 (a Pol-II subunit)
def is_factor(s):
    r = not('POL2' in s or 'POLIII' in s or 'TBP1' in s)
    if not r:
        print('excluding: %(s)s' % locals())
    return r

#df_chip = pd.DataFrame(index=[os.path.basename(fp).rstrip('.fc.signal.bw') for fp in l_foldChange])
#df_chip['stage_modencode'] = list(map(stage_from_fn, df_chip.index))
#df_chip['foldChange'] = l_foldChange
#df_chip['finalPk'] = list(map(finalPk_from_index, df_chip.index))
#df_chip = df_chip.query('~(finalPk == ".")')
#print(len(df_chip))
#df_chip = df_chip[list(map(is_factor, df_chip.index))]
#print(len(df_chip))
#df_chip = df_chip.query('stage_modencode in ["EM", "LE", "L1", "L2", "L3", "L4", "YA"]')
#print(len(df_chip))
#df_chip['stage'] = list(map(lambda s: {'EM': 'em', 'LE': 'em', 'L1': 'l1', 'L2': 'l2', 'L3': 'l3', 'L4': 'l4', 'YA': 'ya'}[s], df_chip['stage_modencode']))
#print(df_chip.head())

rule chip_distribution:
    input:
        lambda wildcards: lp(df_chip.ix[wildcards.chip_id]['foldChange']),
    output:
        pf('{chip_id}', 'chip_distr', '.png', 'chip_distribution'),
    run:
        # Load regulatory elements
        fp_sites = 'annot/regulatory_elements.tsv' % locals()
        df_sites = pd.read_csv(fp_sites, sep='\t')

        # Annotate with short cap
        fp_fwd = lp(pf('scap_pooled', 'bwa062_se30.q10_firstbp', '_fwd.bw'))
        fp_rev = lp(pf('scap_pooled', 'bwa062_se30.q10_firstbp', '_rev.bw'))
        flank_len = 50
        df_sites['scap_mode_fwd'] = df_sites['start'] - flank_len + list(yp.read_regions(fp_fwd, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanargmax))
        df_sites['scap_mode_rev'] = df_sites['start'] - flank_len + list(yp.read_regions(fp_rev, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanargmax))
        print(np.median(df_sites['scap_mode_fwd'] - df_sites['scap_mode_rev']), 'median mode')
        df_sites['scap_mode_count_fwd'] = list(yp.read_regions(fp_fwd, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanmax))
        df_sites['scap_mode_count_rev'] = list(yp.read_regions(fp_rev, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanmax))

        gdf_srev = yp.GenomicDataFrame(df_sites.query('scap_mode_count_rev > 5'), pos_column='scap_mode_rev')
        gdf_mode = yp.GenomicDataFrame(df_sites, pos_column='atac_mode')
        gdf_sfwd = yp.GenomicDataFrame(df_sites.query('scap_mode_count_fwd > 5'), pos_column='scap_mode_fwd')
        gdf_srev.add_track('factor', str(input), flank_len=300, bin_size=1)
        gdf_mode.add_track('factor', str(input), flank_len=300, bin_size=1)
        gdf_sfwd.add_track('factor', str(input), flank_len=300, bin_size=1)

        fig = plt.figure(figsize=(16,20))
        fig.suptitle(os.path.basename(str(input)))

        grid = yp.ImageGrid(fig, 111, nrows_ncols=(5,3), aspect=False, axes_pad=0.15, )
        grid[0].set_title('Rev short cap\n(%d with >5 tags at mode)' % (len(gdf_srev.r),))
        grid[1].set_title('Peak accessibility\n(all %d sites)' % (len(gdf_mode.r),))
        grid[2].set_title('Fwd short cap\n(%d with >5 tags at mode)' % (len(gdf_sfwd.r),))
        grid[0].set_ylabel('All sites')
        grid[3].set_ylabel('Bidirectional promoters')
        grid[6].set_ylabel('Fwd promoters')
        grid[9].set_ylabel('Rev promoters')
        grid[12].set_ylabel('Not promoters')

        gdf_srev.t['factor'].plot(color=yp.RED, label='Fwd motif', ax=grid[0])#f=np.sum)
        gdf_mode.t['factor'].plot(color=yp.RED, label='Fwd motif', ax=grid[1])#f=np.sum)
        gdf_sfwd.t['factor'].plot(color=yp.RED, label='Fwd motif', ax=grid[2])#f=np.sum)
        grid[2].legend(loc='upper right')

        q_ = "(annot == 'coding_promoter') & (strand == '.')"
        gdf_srev.query(q_).t['factor'].plot(color=yp.RED, ax=grid[3])#f=np.sum)
        gdf_mode.query(q_).t['factor'].plot(color=yp.RED, ax=grid[4])#f=np.sum)
        gdf_sfwd.query(q_).t['factor'].plot(color=yp.RED, ax=grid[5])#f=np.sum)

        q_ = "(annot == 'coding_promoter') & (strand == '+')"
        gdf_srev.query(q_).t['factor'].plot(color=yp.RED, ax=grid[6])#f=np.sum)
        gdf_mode.query(q_).t['factor'].plot(color=yp.RED, ax=grid[7])#f=np.sum)
        gdf_sfwd.query(q_).t['factor'].plot(color=yp.RED, ax=grid[8])#f=np.sum)

        q_ = "(annot == 'coding_promoter') & (strand == '-')"
        gdf_srev.query(q_).t['factor'].plot(color=yp.RED, ax=grid[9])#f=np.sum)
        gdf_mode.query(q_).t['factor'].plot(color=yp.RED, ax=grid[10])#f=np.sum)
        gdf_sfwd.query(q_).t['factor'].plot(color=yp.RED, ax=grid[11])#f=np.sum)

        q_ = "(annot != 'coding_promoter')"
        gdf_srev.query(q_).t['factor'].plot(color=yp.RED, ax=grid[12])#f=np.sum)
        gdf_mode.query(q_).t['factor'].plot(color=yp.RED, ax=grid[13])#f=np.sum)
        gdf_sfwd.query(q_).t['factor'].plot(color=yp.RED, ax=grid[14])#f=np.sum)

        for i, ax in enumerate(grid):
            ax.axvline(0, linestyle='dashed', color='k', linewidth=1)

        plt.gca().legend(loc='upper right')
        plt.savefig(str(output), bbox_inches='tight')

#rule chip_process:
#    input:
#        expand(pf('{chip_id}', 'chip_distr', '.png', 'chip_distribution'), chip_id=df_chip.index.tolist()),
