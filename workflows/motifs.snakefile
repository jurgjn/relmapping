
def df_cisbp_make():
    fp_cisbp = htp('.wget/cisbp.ccbr.utoronto.ca/tmp/Caenorhabditis_elegans_2017_01_13_6-32_am/TF_Information.txt')
    df_cisbp = pd.read_csv(fp_cisbp, sep='\t')
    print(len(df_cisbp), 'raw records in Cis-BP')
    print(len(df_cisbp.query('Motif_ID != "."')), 'with Motif_ID set')
    print(len(set(df_cisbp.query('Motif_ID != "."')['Motif_ID'])), 'unique Motif_ID values')
    fp_motif = htp('.raw/wget/cisbp.ccbr.utoronto.ca/tmp/Caenorhabditis_elegans_2017_01_13_6-32_am/pwms_all_motifs/%s.txt')
    df_cisbp['Motif_ID_st_size'] = list(map(lambda motif_id: os.stat(fp_motif % (motif_id,)).st_size if motif_id != '.' else 0, df_cisbp['Motif_ID']))
    print(len(df_cisbp.query('Motif_ID_st_size == 12')), 'motifs with length 0')
    return df_cisbp

def motifs_all():
    return list(sorted(set(df_cisbp_make().query('Motif_ID_st_size > 12')['Motif_ID'])) + ['Chen13_INR', 'Chen13_SP1_Q6', 'Chen13_TATA_01'])

rule cisbp_fimo_inp:
    input: htp('.raw/wget/cisbp.ccbr.utoronto.ca/tmp/Caenorhabditis_elegans_2017_01_13_6-32_am/pwms_all_motifs/{motif_id}.txt')
    output: pf('{motif_id}', 'fimo_inp', '.meme', 'motifs')
    run:
        mat = np.loadtxt(str(input), skiprows=1, usecols=(1,2,3,4))
        with open(str(output), 'w') as fh:
            print('MEME version 4', file=fh)
            print('MOTIF cisbp_fimo_inp' % locals(), file=fh)
            print('letter-probability matrix: alength=4 w=%d' % (mat.shape[0]), file=fh)
            for row in mat:
                print(' ' + '\t'.join(map(str, row)), file=fh)

#rule cisbp_fimo_inp_run:
#    input:
#        expand(pf('{motif_id}', 'fimo_inp', '.meme', 'motifs'), motif_id=motifs_all())

rule fimo_raw:
    input: pf('{motif_id}', 'fimo_inp', '.meme', 'motifs')
    output:
        pf('{motif_id}', 'fimo_raw', '/', 'motifs'),
        pf('{motif_id}', 'fimo_raw', '/fimo.gff', 'motifs')
    params:
        ce10_fa = 'shared/ce10.fa'
    conda:
        'envs/meme.yaml'
    shell: 'fimo --oc {output[0]} {input} {params.ce10_fa}' # --oc required as snakemake automatically creates output directory

rule fimo:
    input:  pf('{motif_id}', 'fimo_raw', '/fimo.gff', 'motifs')
    output: pf('{motif_id}', 'fimo', '.gff', 'motifs')
    shell: 'cp {input} {output}'

rule plot_motif_distribution:
    input:
        pf('{motif_id}', 'fimo', '.gff', 'motifs'),
        pf('{motif_id}', 'fimo_inp', '.meme', 'motifs'),
    output:
        pf('{motif_id}', 'fimo.plot_motif_distribution', '.png', 'motifs')
    run:
        df_motif = pd.read_csv(str(input[0]), sep='\t', comment='#', names=yp.NAMES_GTF)
        df_motif['start'] = df_motif['start'] - 1

        ga_motif_fwd = hts.GenomicArray(chroms=yp.chroms_ce10, stranded=False, typecode='i')
        for i, r in df_motif.query("strand == '+'").iterrows():
            iv = hts.GenomicInterval(r['chrom'], r['start'], r['end'], r['strand'])
            ga_motif_fwd[iv] += 1

        ga_motif_rev = hts.GenomicArray(chroms=yp.chroms_ce10, stranded=False, typecode='i')
        for i, r in df_motif.query("strand == '-'").iterrows():
            iv = hts.GenomicInterval(r['chrom'], r['start'], r['end'], r['strand'])
            ga_motif_rev[iv] += 1

        # Load regulatory elements
        fp_sites = 'annot/regulatory_elements.tsv' % locals()
        df_sites = pd.read_csv(fp_sites, sep='\t')

        # Annotate with short cap modes
        fp_fwd = lp(pf('scap_pooled', 'bwa062_se30.q10_firstbp', '_fwd.bw'))
        fp_rev = lp(pf('scap_pooled', 'bwa062_se30.q10_firstbp', '_rev.bw'))
        flank_len = 50
        df_sites['scap_mode_fwd'] = df_sites['start'] - flank_len + list(yp.read_regions(fp_fwd, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanargmax))
        df_sites['scap_mode_rev'] = df_sites['start'] - flank_len + list(yp.read_regions(fp_rev, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanargmax))
        print(np.median(df_sites['scap_mode_fwd'] - df_sites['scap_mode_rev']), 'median mode')
        df_sites['scap_mode_count_fwd'] = list(yp.read_regions(fp_fwd, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanmax))
        df_sites['scap_mode_count_rev'] = list(yp.read_regions(fp_rev, df_sites['chrom'], df_sites['start'] - flank_len, df_sites['end'] + flank_len, f=np.nanmax))

        # Load motif data centered at rev/fwd short cap and peak accessibility
        gdf_srev = yp.GenomicDataFrame(df_sites.query('scap_mode_count_rev > 5'), pos_column='scap_mode_rev')
        gdf_mode = yp.GenomicDataFrame(df_sites, pos_column='atac_mode')
        gdf_sfwd = yp.GenomicDataFrame(df_sites.query('scap_mode_count_fwd > 5'), pos_column='scap_mode_fwd')
        gdf_srev.t['_fwd'] = yp.GenomicDataFrameTrack(flank_len=300, bin_size=1)
        gdf_srev.t['_rev'] = yp.GenomicDataFrameTrack(flank_len=300, bin_size=1)
        gdf_mode.t['_fwd'] = yp.GenomicDataFrameTrack(flank_len=300, bin_size=1)
        gdf_mode.t['_rev'] = yp.GenomicDataFrameTrack(flank_len=300, bin_size=1)
        gdf_sfwd.t['_fwd'] = yp.GenomicDataFrameTrack(flank_len=300, bin_size=1)
        gdf_sfwd.t['_rev'] = yp.GenomicDataFrameTrack(flank_len=300, bin_size=1)
        gdf_srev.t['_fwd'].m_from_ga(ga=ga_motif_fwd, chroms=gdf_srev.r['chrom'], starts=gdf_srev.r[gdf_srev.pos_column], strand='+')
        gdf_srev.t['_rev'].m_from_ga(ga=ga_motif_rev, chroms=gdf_srev.r['chrom'], starts=gdf_srev.r[gdf_srev.pos_column], strand='-')
        gdf_mode.t['_fwd'].m_from_ga(ga=ga_motif_fwd, chroms=gdf_mode.r['chrom'], starts=gdf_mode.r[gdf_mode.pos_column], strand='+')
        gdf_mode.t['_rev'].m_from_ga(ga=ga_motif_rev, chroms=gdf_mode.r['chrom'], starts=gdf_mode.r[gdf_mode.pos_column], strand='-')
        gdf_sfwd.t['_fwd'].m_from_ga(ga=ga_motif_fwd, chroms=gdf_sfwd.r['chrom'], starts=gdf_sfwd.r[gdf_sfwd.pos_column], strand='+')
        gdf_sfwd.t['_rev'].m_from_ga(ga=ga_motif_rev, chroms=gdf_sfwd.r['chrom'], starts=gdf_sfwd.r[gdf_sfwd.pos_column], strand='-')

        # Show all available CisBP metadata when available
        try:
            import textwrap
            title_str_ = '; '.join([
                '%s: %s' % (k, v,) for k, v in list(df_cisbp_make().query('Motif_ID == "%s"' % (wildcards.motif_id)).iterrows())[0][1].items()
                ])
            title_str = '\n'.join(textwrap.wrap(title_str_, width=100))
        except:
            title_str = wildcards.motif_id

        # Plot
        fig = plt.figure(figsize=(16,24))
        fig.suptitle(title_str)
        grid = yp.ImageGrid(fig, 111, nrows_ncols=(7,3), aspect=False, axes_pad=0.15, )
        grid[0].set_title('Rev short cap\n(%d with >5 tags at mode)' % (len(gdf_srev.r),))
        grid[1].set_title('Peak accessibility\n(all %d sites)' % (len(gdf_mode.r),))
        grid[2].set_title('Fwd short cap\n(%d with >5 tags at mode)' % (len(gdf_sfwd.r),))
        grid[0].set_ylabel('All sites')
        grid[3].set_ylabel('Bidirectional promoters')
        grid[6].set_ylabel('Fwd promoters')
        grid[9].set_ylabel('Rev promoters')
        #grid[12].set_ylabel('non_coding_promoter')
        grid[12].set_ylabel('strong_enhancer')
        grid[15].set_ylabel('other_element')
        grid[18].set_ylabel('~coding_promoter')
                
        gdf_srev.t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[0])#f=np.sum)
        gdf_srev.t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[0])#f=np.sum)
        gdf_mode.t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[1])#f=np.sum)
        gdf_mode.t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[1])#f=np.sum)
        gdf_sfwd.t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[2])#f=np.sum)
        gdf_sfwd.t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[2])#f=np.sum)
        grid[2].legend(loc='upper right')

        q_ = "(annot == 'coding_promoter') & (strand == '.')"
        gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[3])#f=np.sum)
        gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[3])#f=np.sum)
        gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[4])#f=np.sum)
        gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[4])#f=np.sum)
        gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[5])#f=np.sum)
        gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[5])#f=np.sum)

        q_ = "(annot == 'coding_promoter') & (strand == '+')"
        gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[6])#f=np.sum)
        gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[6])#f=np.sum)
        gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[7])#f=np.sum)
        gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[7])#f=np.sum)
        gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[8])#f=np.sum)
        gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[8])#f=np.sum)

        q_ = "(annot == 'coding_promoter') & (strand == '-')"
        gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[9])#f=np.sum)
        gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[9])#f=np.sum)
        gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[10])#f=np.sum)
        gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[10])#f=np.sum)
        gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[11])#f=np.sum)
        gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[11])#f=np.sum)

        # (Too few regions)
        #q_ = "(annot == 'non_coding_promoter')"
        #gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[12])#f=np.sum)
        #gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[12])#f=np.sum)
        #gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[13])#f=np.sum)
        #gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[13])#f=np.sum)
        #gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[14])#f=np.sum)
        #gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[14])#f=np.sum)

        q_ = "(annot == 'strong_enhancer')"
        gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[12])#f=np.sum)
        gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[12])#f=np.sum)
        gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[13])#f=np.sum)
        gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[13])#f=np.sum)
        gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[14])#f=np.sum)
        gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[14])#f=np.sum)

        q_ = "(annot == 'other_element')"
        gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[15])#f=np.sum)
        gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[15])#f=np.sum)
        gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[16])#f=np.sum)
        gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[16])#f=np.sum)
        gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[17])#f=np.sum)
        gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[17])#f=np.sum)

        q_ = "(annot != 'coding_promoter')"
        gdf_srev.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[18])#f=np.sum)
        gdf_srev.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[18])#f=np.sum)
        gdf_mode.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[19])#f=np.sum)
        gdf_mode.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[19])#f=np.sum)
        gdf_sfwd.query(q_).t['_fwd'].plot(color=yp.RED, label='Fwd motif', ax=grid[20])#f=np.sum)
        gdf_sfwd.query(q_).t['_rev'].plot(color=yp.BLUE, label='Rev motif', ax=grid[20])#f=np.sum)

        for i, ax in enumerate(grid):
            ax.axvline(0, linestyle='dashed', color='k', linewidth=1)

        #plt.gca().legend(loc='upper right')

        # Add sequence logo
        def motif_logo(counts):
            data = weblogolib.LogoData.from_counts(
                alphabet = weblogolib.std_alphabets['dna'],
                counts = counts,
            )
            options = weblogolib.LogoOptions()
            options.resolution = 300
            #options.color_scheme = weblogolib.base_pairing
            options.color_scheme = weblogolib.classic
            format = weblogolib.LogoFormat(data, options)
            png = weblogolib.png_formatter( data, format)
            img = matplotlib.image.imread(io.BytesIO(png))
            return img

        # Very naive PWM extraction via np.loadtxt...
        # ...though still more functional than the parser in biopython: https://github.com/biopython/biopython/issues/232
        ascii_letters_ = 'abcdfghijklmnopqrstuvwxyz' + 'ABCDFGHIJKLMNOPQRSTUVWXYZ'
        counts_=np.loadtxt(str(input[1]), comments=list(ascii_letters_))
        #comments=list(string.ascii_letters) would incorrectly cut into scientific notation:
        #  0.000604784281206	7.90469533067e-07	3.29460651483e-05	0.999361479184

        img = motif_logo(counts=counts_)
        ax_weblogo = fig.add_axes([0.02, 0.9, 0.2, 0.1])
        ax_weblogo.imshow(img)
        ax_weblogo.axis('off')

        # Save file
        plt.savefig(str(output), bbox_inches='tight')

rule motifs_counts:
    input:
        pf('{motif_id}', 'fimo', '.gff', 'motifs'),
    output:
        pf('{motif_id}', 'fimo.counts', '.tsv', 'motifs')
    run:
        def count_motif(fp_motif, strand):
            df_motif = pd.read_csv(fp_motif, sep='\t', comment='#', names=yp.NAMES_GTF)
            df_motif['start'] = df_motif['start'] - 1

            ga_motif = hts.GenomicArray(chroms=yp.chroms_ce10, stranded=False, typecode='i')
            for i, r in df_motif.query("strand == '%s'" % (strand,)).iterrows():
                pos = (r['start'] + r['end']) // 2 # "One base pair per motif"
                iv = hts.GenomicInterval(r['chrom'], pos, pos + 1, r['strand'])
                ga_motif[iv] += 1

            fp_sites = 'annot/regulatory_elements.tsv' % locals()
            df_sites = pd.read_csv(fp_sites, sep='\t')
            gdf_mode = yp.GenomicDataFrame(df_sites, pos_column='atac_mode')
            gdf_mode.t['_motif'] = yp.GenomicDataFrameTrack(flank_len=100, bin_size=1)
            gdf_mode.t['_motif'].m_from_ga(ga=ga_motif, chroms=gdf_mode.r['chrom'], starts=gdf_mode.r[gdf_mode.pos_column], strand='+')
            return gdf_mode.t['_motif'].m.sum(axis=1)

        df_motif_counts = pd.DataFrame()
        df_motif_counts[wildcards.motif_id + '_fwd'] = list(map(int, count_motif(str(input), '+')))
        df_motif_counts[wildcards.motif_id + '_rev'] = list(map(int, count_motif(str(input), '-')))
        df_motif_counts.to_csv(str(output), sep='\t', index=False)

#rule motifs_counts_agg:
#    input:
#        expand(pf('{motif_id}', 'fimo.counts', '.tsv', 'motifs'), motif_id=motifs_all()),
#    output:
#        pf('_all', 'fimo.counts', '.tsv', 'motifs'),
#    run:
#        pd.concat([pd.read_csv(str(input_), sep='\t') for input_ in input], axis=1).to_csv(str(output), sep='\t', index=False)

rule scanMotifGenomeWide:  # scanMotifGenomeWide outputs 1-based start coordinate in a .bed-file(!??) => fix via awk
    input:
        pf('{motif_id}', 'homer_inp', '.motif', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_raw.bed', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    shell:
        '''
        scanMotifGenomeWide.pl {input} ce10 -bed -p 1 > {output[0]}
        awk -F'\\t' -v OFS='\\t' '{{print $1,$2-1,$3,$4,$5,$6}}' {output[0]} > {output[1]}
        '''

rule scanMotifGenomeWideCoverage_0fwd:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_0fwd.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_0fwd.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "+") {{print $1,$2,$2+1}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_0rev:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_0rev.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_0rev.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "-") {{print $1,$3-1,$3}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_2fwd:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_2fwd.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_2fwd.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "+") {{print $1,$2+2,$2+3}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_2rev:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_2rev.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_2rev.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "-") {{print $1,$3-3,$3-2}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_5fwd:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_5fwd.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_5fwd.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "+") {{print $1,$2+5,$2+6}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_5rev:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_5rev.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_5rev.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "-") {{print $1,$3-6,$3-5}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_6fwd:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_6fwd.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_6fwd.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "+") {{print $1,$2+6,$2+7}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule scanMotifGenomeWideCoverage_6rev:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'),
    output:
        temp(pf('{motif_id}', 'scanMotifGenomeWide', '_6rev.bedGraph', 'motifs')),
        pf('{motif_id}', 'scanMotifGenomeWide', '_6rev.bw', 'motifs'),
    shell:
        '''
        awk -F'\\t' -v OFS='\\t' '($6 == "-") {{print $1,$3-7,$3-6}}' {input} | bedtools genomecov -bga -i stdin -g shared/ce10.chroms > {output[0]}
        sort -k1,1 -k2,2n -o {output[0]} {output[0]}
        bedGraphToBigWig {output[0]} shared/ce10.chroms {output[1]}
        '''

rule motifs_chen2013:
    input:
        #expand(pf('{motif_id}', 'fimo', '.gff', 'motifs'), motif_id=motifs_all()),
        #expand(pf('{motif_id}', 'fimo.plot_motif_distribution', '.png', 'motifs'), motif_id=motifs_all()),
        #pf('_all', 'fimo.counts', '.tsv', 'motifs'), # motis_counts_agg
        #expand(pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'), motif_id=['AATAAA', 'AATAAA1']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_0fwd.bw', 'motifs'), motif_id=['AATAAA', 'AATAAA1']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_0rev.bw', 'motifs'), motif_id=['AATAAA', 'AATAAA1']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_0fwd.bw', 'motifs'), motif_id=['Chen13_TATA_01_5', 'Chen13_TATA_01_6', 'Chen13_TATA_01_7', 'Chen13_TATA_01_8', 'Chen13_TATA_01_9', 'TATAAAA_0', 'TATAAAA_1']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_0rev.bw', 'motifs'), motif_id=['Chen13_TATA_01_5', 'Chen13_TATA_01_6', 'Chen13_TATA_01_7', 'Chen13_TATA_01_8', 'Chen13_TATA_01_9', 'TATAAAA_0', 'TATAAAA_1']),
        # INR, positioned at 'A' of 'TCA' as in (Chen et al., 2013); 8/9 barf as motif search leads to zero hits...
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_5fwd.bw', 'motifs'), motif_id=['Chen13_INR_2', 'Chen13_INR_3', 'Chen13_INR_4', 'Chen13_INR_5', 'Chen13_INR_6', 'Chen13_INR_7']),#, 'Chen13_INR_8', 'Chen13_INR_9']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_5rev.bw', 'motifs'), motif_id=['Chen13_INR_2', 'Chen13_INR_3', 'Chen13_INR_4', 'Chen13_INR_5', 'Chen13_INR_6', 'Chen13_INR_7']),#, 'Chen13_INR_8', 'Chen13_INR_9']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_6fwd.bw', 'motifs'), motif_id=['Chen13_SP1_Q6_6']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_6rev.bw', 'motifs'), motif_id=['Chen13_SP1_Q6_6']),

rule homer_inp_iupac: # HOMER motif input based on consensus sequence using IUPAC coding
    output:
        pf('{motif_id}', 'homer_inp', '.motif', 'motifs'),
    params:
        iupac = lambda wildcards: config['motifs_iupac'][wildcards.motif_id]['iupac']
    shell:
        'seq2profile.pl {params.iupac} 0 {wildcards.motif_id} > {output}'

rule plot_at_tss:
    input:
        pf('{motif_id}', 'scanMotifGenomeWide', '_0fwd.bw', 'motifs'),
    output:
        pf('{motif_id}', 'scanMotifGenomeWide.plot_at_tss', '.png', 'motifs'),
    run:
        (fp_inp, fp_out) = (str(input), str(output))
        fp_ = 'annot/S2_regulatory_annotation/S2_regulatory_annotation.tsv'
        df_regl = pd.read_csv(fp_, sep='\t')
        gdf = yp.GenomicDataFrame(df_regl, pos_column='tss_fwd')
        gdf_p = gdf.query('(annot_fwd == "coding_promoter")')
        gdf_e = gdf.query('(annot == "putative_enhancer")')

        for gdf_ in [gdf_p, gdf_e]:
            assert(os.path.isfile(fp_inp))
            gdf_.add_track(wildcards.motif_id, fp_inp, flank_len=50, bin_size=1, f_bin=np.nanmean)

        fig = plt.figure(figsize=(10,8))
        ax1 = plt.subplot(2, 2, 1)
        plt.suptitle('%s (consensus sequence: %s)' % (wildcards.motif_id, config['motifs_iupac'][wildcards.motif_id]['iupac']), fontsize=10)
        gdf_p.t[wildcards.motif_id].plot()#q(f=np.nanmean, q=3)
        #plt.gca().set_ylim(0, 0.2)
        plt.gca().set_ylabel('coding_promoter')
        plt.subplot(2, 2, 3, sharey=ax1)
        gdf_e.t[wildcards.motif_id].plot()#q(f=np.nanmean, q=3)
        #plt.gca().set_ylim(0, 0.2)
        plt.gca().set_ylabel('putative_enhancer')

        for gdf_ in [gdf_p, gdf_e]:
            assert(os.path.isfile(fp_inp))
            gdf_.add_track(wildcards.motif_id, fp_inp, flank_len=200, bin_size=1, f_bin=np.nanmean)

        plt.subplot(2, 2, 2, sharey=ax1)
        gdf_p.t[wildcards.motif_id].plot()
        plt.subplot(2, 2, 4, sharey=ax1)
        gdf_e.t[wildcards.motif_id].plot()

        plt.savefig(fp_out, bbox_inches='tight', dpi=100)

rule motifs_iupac:
    input:
        expand(pf('{motif_id}', 'homer_inp', '.motif', 'motifs'), motif_id=list(config['motifs_iupac'].keys())),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '.bed', 'motifs'), motif_id=list(config['motifs_iupac'].keys())),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_0fwd.bw', 'motifs'), motif_id=list(config['motifs_iupac'].keys())),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_0rev.bw', 'motifs'), motif_id=list(config['motifs_iupac'].keys())),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_2fwd.bw', 'motifs'), motif_id=['Sloutskin2015_mammalian_Initiator']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide', '_2rev.bw', 'motifs'), motif_id=['Sloutskin2015_mammalian_Initiator']),
        expand(pf('{motif_id}', 'scanMotifGenomeWide.plot_at_tss', '.png', 'motifs'), motif_id=list(config['motifs_iupac'].keys())),
