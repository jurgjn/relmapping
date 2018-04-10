# -*- coding: utf-8; tab-width: 4; -*-
"""
Yet Another Region Plotter
"""

import collections
import copy as copy_module # Yes, this might make someone spin in their graves...
import csv
import datetime
import os
import random
import subprocess
import sys
import hashlib
#import cPickle
#import cPickle as pkl
import pprint
import tempfile
import atexit
#import commands
import datetime
import math
import os
import re
import pickle
import shutil
import subprocess
import sys
import inspect
import itertools

#import StringIO

import __main__
import datetime
import inspect
import os
import pickle
import pprint
import shutil
import subprocess
import sys
import gzip

import numpy as np
import scipy as sp
import pandas as pd
#pd.core.format.set_option('notebook_repr_html', True)

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt # ...instead of from pylab import * -- see http://carreau.github.io/posts/10-No-PyLab-Thanks.ipynb.html
import mpl_toolkits.axes_grid1
from mpl_toolkits.axes_grid1 import ImageGrid

import sklearn
import sklearn.cluster
import sklearn.decomposition
import sklearn.preprocessing

import matplotlib_venn
import pybedtools
import pysam

import pyBigWig

import twobitreader

from IPython.core.display import HTML

import HTSeq as hts

NAMES_BED9 = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']
NAMES_BED3 = NAMES_BED9[:3]
NAMES_BED6 = NAMES_BED9[:6]
NAMES_NARROWPEAK = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak',]
NAMES_GTF = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
NAMES_MSCORE = ['seq_id', 'dist', 'consensus', 'name', 'strand', 'mscore']

TO_GTF_KWARGS = {'index': False, 'header': False, 'quoting': csv.QUOTE_NONE, 'sep': '\t'}
TO_GTF_GZ_KWARGS = {'index': False, 'header': False, 'quoting': csv.QUOTE_NONE, 'sep': '\t', 'compression': 'gzip'}

'''
https://genome.ucsc.edu/FAQ/FAQformat#format9
> The following definition is used for extended gene prediction tables. In alternative-splicing situations, each transcript has a row in this table. The refGene table is an example of the genePredExt format.

http://software.broadinstitute.org/software/igv/genepred
> These files have the same format. Columns in the file match the columns in the table, as described in the “Gene Predictions (Extended)” section of the genePred table format FAQ. Note: The first column of this file holds an integer, which is not documented in the FAQ and is ignored by IGV.
'''
NAMES_REFGENE = ('undocumented_int', 'name', 'chrom', 'strand',
    'txStart', 'txEnd',
    'cdsStart', 'cdsEnd',
    'exonCount', 'exonStarts', 'exonEnds',
    'score', 'name2',
    'cdsStartStat', 'cdsEndStat', 'exonFrames',)

# Color-blindness safe colors from http://www.nature.com/nmeth/journal/v8/n6/full/nmeth.1618.html
# http://jfly.iam.u-tokyo.ac.jp/color/
BLACK   = '#000000' # 0,0,0
ORANGE  = '#e69f00' # 230,159,0
SKYBLUE = '#56b4e9' # 86,180,233
GREEN   = '#009e73' # 0,158,115
YELLOW  = '#f0e442' # 240,228,66
BLUE    = '#0072b2' # 0,114,178
RED     = '#d55e00' # 213,94,0
PURPLE  = '#cc79a7' # 204,121,167

TO_GTF_KWARGS = {'index': False, 'header': False, 'quoting': csv.QUOTE_NONE, 'sep': '\t'}
TO_GTF_GZ_KWARGS = {
    'index': False, 
    'header': False, 
    'quoting': csv.QUOTE_NONE, 
    'sep': '\t', 
    'compression': 'gzip',
}

chroms_ce10 = {
    'chrI': 15072423,
    'chrII': 15279345,
    'chrIII': 13783700,
    'chrIV': 17493793,
    'chrM': 13794,
    'chrV': 20924149,
    'chrX': 17718866
}
d_ensembl_ucsc = {
    'I':     'chrI',
    'II':    'chrII',
    'III':   'chrIII',
    'IV':    'chrIV',
    'V':     'chrV',
    'X':     'chrX',
    'MtDNA': 'chrM',
}
d_ucsc_ensembl = {
    'chrI':   'I',
    'chrII':  'II',
    'chrIII': 'III',
    'chrIV':  'IV',
    'chrV':   'V',
    'chrX':   'X',
    'chrM':   'MtDNA',
}

def nanmean_pc(*args, **kwargs): return np.nanmean(*args, **kwargs) * 100.0

def f_uk(x): return '{:,}'.format(x)
def float_uk(s): return float(''.join(filter(lambda x: x in '0123456789.', s)))


def rolling_mean_kernel(width): return np.ones(width) / float(width)
f_mean50 = lambda y: np.convolve(y, rolling_mean_kernel(50), mode='same')

def nanmax(c):
    nanmax_ = np.nanmax(c)
    return nanmax_ if nanmax_ == nanmax_ else 0

"""
def read_regions(fp, chroms, starts, ends, f=None):
    fh = pyBigWig.open(fp)
    for chrom, start, end in zip(chroms, starts, ends):
        if f is None:
            yield fh.values(chrom, start, end)
        else:
            yield f(np.array(fh.values(chrom, start, end)))
    fh.close()
"""
def read_regions(fp, chroms, starts, ends, v=False):
    assert os.path.isfile(fp)
    n_clip = 0
    fh = pyBigWig.open(fp)
    for chrom, start_, end_ in zip(chroms, starts, ends):
        # fixes pyBigWig -- RuntimeError: You must supply a chromosome, start and end position.
        start = int(start_)
        end = int(end_)

        # Clip region if necessary
        if (0 <= start) and (end < fh.chroms(chrom)):
            r = fh.values(chrom, start, end)#, numpy=True)
        else:
            r = np.zeros(end - start) # Should get the partial signal
            n_clip += 1

        yield np.array(r)
    fh.close()
    if v: print(n_clip, 'clipped regions')

def mread_regions(fp, chroms, starts, ends, rememoize=False, v=True):
    """
    Naive memoization for querying a set of regions; speeds up representative situations (worm genome, 20k regions) by about an order of magnitude
    (Possible alternative: https://pythonhosted.org/joblib/memory.html)
    """
    assert os.path.isfile(fp)
    # Calculate hash of coordinates and properties of the signal file
    (l_chrom, l_start, l_end) = (list(chroms), list(starts), list(ends))
    regions_hash = hashlib.md5()
    for chrom, start, end in zip(chroms, starts, ends):
        regions_hash.update(('%(chrom)s:%(start)d-%(end)d' % locals()).encode('utf-8'))
    st_fp = os.stat(fp) # Raw os.stat also contains *access* time...
    regions_hash.update((','.join(map(str, [st_fp.st_size, st_fp.st_ctime, st_fp.st_mtime]))).encode('utf-8'))
    regions_hashid = regions_hash.hexdigest()

    # Clear memoization if necessary
    fp_memoized = os.path.join(os.path.expanduser('~/relmapping/tmp'), os.path.split('%(fp)s.mread_regions_%(regions_hashid)s.tmp' % locals())[1])
    if rememoize and os.path.isfile(fp_memoized):
        os.remove(fp_memoized)

    # If necessary, create memoization and store to disk
    if not os.path.isfile(fp_memoized):
        memoized=False
        with open(fp_memoized, 'wb') as fh_memoized:
            for data in read_regions(fp, l_chrom, l_start, l_end):
                np.asarray(data).tofile(fh_memoized)
    else:
        memoized=True

    # Read and yield back memoization
    with open(fp_memoized, 'rb') as fh_memoized:
        data_all = np.fromfile(fh_memoized, count=-1)
    offset = 0
    for (chrom, start, end) in zip(l_chrom, l_start, l_end):
        yield data_all[offset: offset + (end - start)]
        offset += end - start
    if v:
        print('mread_regions %s %s' % (memoized, os.path.basename(fp_memoized)))

def pairwise(iterable):
    """
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    source: https://docs.python.org/2/library/itertools.html
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def is_int(x):
    try:
        if int(x) == x:
            return True
    except:
        return False

class OrderedCounter(collections.Counter, collections.OrderedDict):
    """
    Counter that remembers the order elements are first encountered'
    source: https://docs.python.org/2/library/collections.html#collections.OrderedDict
    """
    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, collections.OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (collections.OrderedDict(self),)

def permilles(a):
    """
    equivalent to (but faster): 10 * sp.stats.percentileofscore(l_tag_count, tag_count, kind="weak")
    print permilles([-0.5, 0, 9, 2])
    [0.25, 0.5, 1.0, 0.75]
    """
    a_sort = np.sort(a)
    return [np.searchsorted(a_sort, a_i, side='right') / float(a_sort.shape[0]) for a_i in a]

def parse_igvstr(s="II:7,641,761-7,642,754"):
    (chrom, coords) = s.split(':')
    (start_s, end_s) = coords.split('-')
    start = int(start_s.replace(',', '')) #Remove commas, if necessary
    start = start - 1
    end = int(end_s.replace(',', ''))
    return [chrom, start, end]

def subplot_row_col(row_max, col_max, row, col):
    return plt.subplot(row_max, col_max, (row - 1) * col_max + col)

def plot_a(data, *args, **kwargs):
    """
    # Plot the average profile of data
    TODO visualise uncertainty with:
        http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.fill_between
    e.g. add "uncertainty regions for .5, .8. .95" (all with different colours?)

    TODO in principle, all aggregate plot data *is* discrete (so plot style=step...)
    """
    data_mean_profile = np.mean(data, 0)
    plot(data_mean_profile, *args, **kwargs)

def quantile_mid(i, q):
    """
    print quantile_mid(0, 5)
    print quantile_mid(1, 5)
    print quantile_mid(2, 5)
    print quantile_mid(3, 5)
    print quantile_mid(4, 5)
    0.1
    0.3
    0.5
    0.7
    0.9
    """
    return (1.0 / q) * (i + 0.5)

## Various helper methods for alignments
def l_iv_from_df(df, col_chrom='chrom', col_start='start', col_end='end'):
    return [htw.GenomicInterval(chrom, start, end) for (chrom, start, end) in zip(df[col_chrom], df[col_start], df[col_end])]

def l_gp_from_df(df, col_chrom='chrom', col_start='start'):
    return [htw.GenomicInterval(chrom, start, start + 1) for (chrom, start) in zip(df[col_chrom], df[col_start])]

def l_ga_iv_counts(l_iv, ga):
    return [sum(list(ga[iv])) for iv in l_iv]

def l_ga_iv_mode(l_iv, ga):
    l_mode = []
    for iv in l_iv:
        offset = np.argmax(list(ga[iv]))
        mode = iv.start + offset
        l_mode.append(mode)
    return l_mode

def l_bw_iv_mode(l_iv, fp_bw):
    # add option: flank_len extends interval for mode search to [-flank_len;+flank_len]
    fh = bbi.BigWigFile(fp_bw)
    l_mode = []
    for iv in l_iv:
        offset = np.argmax(fh.read_np(iv.chrom, iv.start, iv.end))
        mode = iv.start + offset
        l_mode.append(mode)
    fh.close()
    return l_mode

bed_names = ('chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb')

def df_from_csv(fp, ext=None, names=None, *args, **kwargs):
    assert os.path.isfile(fp)
    # TODO: how to over-ride sep/names when necessary?
    if ext is None: ext = os.path.splitext(fp)[1]
    d_ext_names = {
        '.bed3': ('chrom', 'start', 'end'),
        '.bed9': ('chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'),
        '.narrowPeak': ('chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'),
        '.gtf':  ('chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
    }
    if names is None:
        names = d_ext_names[ext]
    return pd.read_csv(fp, sep='\t', names=names, *args, **kwargs)

"""
def df_from_merge(df_a, df_b, label_a='label_a', label_b='label_b'):
    df_a_merge = df_a[['chrom', 'start', 'end']].copy()
    df_b_merge = df_b[['chrom', 'start', 'end']].copy()
    df_a_merge['source'] = label_a
    df_b_merge['source'] = label_b
    df_ab = pd.concat([df_a_merge, df_b_merge])
    df_ab.sort(['chrom', 'start', 'end'], inplace=True)
    # bedtools arguments changed; need to update this eventually
    ps = subprocess.Popen('bedtools merge -i stdin -nms', stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    out, err = ps.communicate(df_ab.to_csv(sep='\t', header=False, index=False))
    assert ps.returncode == 0
    df_ab_merge = pd.read_csv(StringIO.StringIO(out), sep='\t', names=['chrom', 'start', 'end', 'source'])
    df_ab_merge['source'] = df_ab_merge['source'].apply(lambda s: ';'.join(sorted(set(s.split(';')))))
    return df_ab_merge
"""

def df_from_org(s):
    df = pd.read_csv(StringIO.StringIO(s), header=0, sep='\s*\|\s*', engine='python', comment='#')
    return df.drop(df.columns[[0,-1]], 1)

def df_sample(df, k=5):
    rows = random.sample(df.index.tolist(), k)
    return df.ix[rows]

def print_iv(df):
    for i, r in df.iterrows():
        print('%s:%s-%s' % (r['chrom'], r['start'], r['end']))

def df_shuffle_rows(df):
    return df.reindex(np.random.permutation(df.index))

def df_browse_iv(df):
    print(len(df), 'regions:')
    print_iv(df_sample(df))

def df_pack_gfftags(df):
    return [";".join([("%s=%s" % (k, v)).replace(" ", "") for k,v in r.items()]) for i, r in df.iterrows()]

def df_gfftags_pack(df, gfftags_col='name', keep_cols=['chrom', 'start', 'end', 'name', 'score', 'strand']):
    # Warning: barfs when name doesn't exist
    pack_cols = [column for column in df.columns if not (column in keep_cols)]
    def pack_row(r): return htw.gff_attr_string(collections.OrderedDict([(col, r[col]) for col in pack_cols]))
    se_pack = df.apply(pack_row, axis=1, reduce=False, raw=True)
    df_out = df[keep_cols].copy()
    if not(gfftags_col in df_out.columns):
        df_out[gfftags_col] = se_pack
    else:
        df_out[gfftags_col] = df_out[gfftags_col] + ';' + se_pack
    return df_out

def df_gfftags_unpack(df, name='name'):
    df_out = df.drop(name, 1)
    df_name = pd.DataFrame(df[name].apply(hts.parse_GFF_attribute_string).tolist())
    df_name = df_name.convert_objects(convert_numeric=True)
    for col in df_name.columns:
        df_out[col] = df_name[col]
    return df_out

def df_reorder_columns(df, l_head):
    l_rest = list(df.columns)
    for i in l_head: l_rest.remove(i)
    return df[l_head + l_rest]

def parse_gtfattr(attr, names):
    d_col_list = collections.OrderedDict((name, []) for name in names)
    for attr_i in attr:
        d = collections.OrderedDict([(attr_ij.split()[0], attr_ij.split()[1].lstrip('"').rstrip('"')) for attr_ij in attr_i.split(';')])
        for (k,v) in d.items():
            d_col_list[k].append(v)
    return pd.DataFrame(d_col_list)

def dict_from_attr(attr):
    d = {}
    for attr_i in attr.split(';'):
        d[attr_i.split()[0]] = attr_i.split()[1].lstrip('"').rstrip('"')
    return d

def dict_from_attr_bed(attr):
    d = {}
    for attr_i in attr.split(';'):
        try:
            (k, v) = attr_i.split('=')
            d[k] = v.lstrip('"').rstrip('"')
        except ValueError:
            if attr_i == '.':
                continue
            else:
                print(attr_i)
                raise
    return d

def df_from_l_dict(l_dict): return pd.DataFrame.from_dict(collections.OrderedDict(enumerate(l_dict)), orient='index').convert_objects(convert_numeric=True)

def enr_(obs_):
    s = obs_.sum(axis=0).sum(axis=0)
    r = obs_.sum(axis=0) / s
    c = obs_.sum(axis=1) / s
    exp_ = np.outer(c, r) * s
    return np.log2(obs_.divide(exp_))

def imshow_enr(crosstab_, ax=None, row_labels=None, col_labels=None, vmin=-1, vmax=+1, rotation=70, *args, **kwargs):
    if ax is None: ax = plt.gca()
    if row_labels is None: row_labels = crosstab_.index
    if col_labels is None: col_labels = crosstab_.columns
    ax.imshow(enr_(crosstab_), cmap='coolwarm', interpolation='nearest', vmin=vmin, vmax=vmax, *args, **kwargs)
    ax.xaxis.tick_top()
    ax.set_xticks(range(len(col_labels)))
    ax.set_yticks(range(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=rotation)
    ax.set_yticklabels(row_labels)
    for (y, x), c in np.ndenumerate(crosstab_):
        ax.text(x, y, '%d' % (c,), color='k', horizontalalignment='center', verticalalignment='center')

def imshow_rpc(crosstab_, ax=None, cmap='viridis', row_labels=None, col_labels=None, vmin=0, vmax=+100, rotation=70, 
               boundaries=np.linspace(0, 30, 7, endpoint=True), *args, **kwargs):
    # plot a heatmap of a cross-tabulation coloured by row percentages
    r_sum = crosstab_.sum(axis=1)
    rpc = 100 * crosstab_.divide(r_sum, axis='rows')
    if ax is None: ax = plt.gca()
    if row_labels is None: row_labels = crosstab_.index
    if col_labels is None: col_labels = crosstab_.columns
    im = ax.imshow(rpc, cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax, *args, **kwargs)

    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    plt.colorbar(im, cax=cax, boundaries=boundaries)

    ax.xaxis.tick_top()
    ax.set_xticks(range(len(col_labels)))
    ax.set_yticks(range(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=rotation)
    ax.set_yticklabels(row_labels)
    for ((y, x), c), ((yy, xx), p) in zip(np.ndenumerate(crosstab_), np.ndenumerate(rpc)):
        ax.text(x, y, '%d\n(%.1f%%)' % (c,p,), color='k', horizontalalignment='center', verticalalignment='center')

def read_wbgtf(fp, parse_attr=[], *args, **kwargs):
    df = pd.read_csv(fp, sep='\t', names=NAMES_GTF, comment='#', *args, **kwargs)
    if parse_attr:
        return df_gfftags_unpack(df, name='attribute') # Does not preserve attribute order...
    else:
        return df
    #df_attr_col = df_from_l_dict(map(dict_from_attr, df['attribute']))
    #df_attr = pd.concat([df.drop('attribute', axis=1), df_attr_col], axis=1)
    #return df_reorder_columns(df_attr, list(NAMES_GTF[:8]) + parse_attr)

def to_wbgtf(df, fp, **kwargs):
    l_attr = [column for column in df.columns if not(column in NAMES_GTF)]
    def gtf_attr_string(r):
        return(" ; ".join([("%s %s" % (key, '"%(val)s"' % locals() if type(val) is str else val))\
                         for key, val in zip(l_attr, r) if val == val]))
    df_out = df[list(NAMES_GTF[:8])].copy()
    df_out['attribute'] = df[l_attr].apply(gtf_attr_string, axis=1, reduce=False, raw=True)
    df_out.sort_values(['chrom', 'start', 'end', 'strand']).to_csv(fp, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, **kwargs)

def to_csv_gz(fp, df, track_line=b'##displayName=locus_id\n', *args, **kwargs):
    fh = gzip.open(fp, 'wb')
    fh.write(track_line)
    fh.write(df.to_csv(*args, **kwargs).encode())
    fh.close()

def to_gffbed(df_inp, fp, bed_cols=None, trackline='#track gffTags=on', v=False):
    """
    trackline='#track gffTags=on useScore=1'
    http://software.broadinstitute.org/software/igv/BED
    > The GFF Name property will become the display name of the feature.
    > You must URL encode spaces and other whitespace (e.g. replace space with %20).  This is not a requirement
    of gff3, rather required because bed files are whitespace delimited.
    """
    if bed_cols is None:
        bed_cols = []
        for col in NAMES_BED9:
            if col in df_inp.columns: bed_cols.append(col)
            else: break

    if v: print('bed_cols =', bed_cols)

    df_out = df_inp[bed_cols].copy()
    df_attr = df_inp.drop(bed_cols, axis=1)
    if len(bed_cols) > 3:
        df_attr['Name'] = df_inp[bed_cols[3]]
        df_attr = df_attr[df_attr.columns[[range(-1, len(df_attr.columns) - 1)]]]

    if len(bed_cols) == 3: bed_cols.append('name')
    def pack_row(r): return (";".join([("%s=%s" % (k, v)).replace(" ", "%20") for k, v in zip(r.index, r)]))
    df_out[bed_cols[3]] = df_attr.apply(pack_row, axis=1, reduce=False, raw=True)

    with open(fp, 'w') as fh:
        print(trackline, file=fh)
        df_out.to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

    if v: return df_out

def peak_mode(m):
    # http://stackoverflow.com/questions/16202348/numpy-divide-row-by-row-sum
    m_rownorm = m / m.sum(axis=1)[:,None]
    m_rowmode = (m_rownorm * np.arange(1, m.shape[1] + 1)).sum(axis=1)
    return m_rowmode

def sort_clusters(l_l, l_v, summary=np.median):
    # sort cluster labels l_l summarized values in l_v
    l = list(set(l_l))
    v = [summary(l_v[l_l == l_i]) for l_i in l]
    r = sp.stats.rankdata(v, method='dense')#np.argsort(v)
    return [r[l_l_i] for l_l_i in l_l]

def iv_centered_at(chrom, start_, end_, flank_len=500, strand='.'):
    mid = int((start_ + end_ - 1) / 2.0)
    start  = mid - flank_len
    end    = mid + flank_len + 1
    return hts.GenomicInterval(chrom, start, end, strand)

def seq_skew(seq, d_skew):
    w = np.array(map(lambda seq_i: d_skew.get(seq_i, 0), seq))
    return np.sum(w) / float(np.sum(np.abs(w)))

def mfc(x, a): # "moderated log2 fold change for skew calculations, similar to Chen at al., 2013 long cap directionality
    x_p = (x > 0).sum(axis=a)
    x_n = (x < 0).sum(axis=a)    
    return np.log2((x_p + 1) / (x_n + 1))

def read_seq(chroms, starts, ends, fp='shared/ce10.2bit'):
    fh = twobitreader.TwoBitFile(fp)
    for chrom, start, end in itertools.islice(zip(chroms, starts, ends), None):
        yield str(fh[chrom][start:end])

def read_seq_fwd(l_chrom, l_pos, uflank, dflank, *args, **kwargs):
    for seq in read_seq(l_chrom, l_pos - uflank, l_pos + dflank + 1, *args, **kwargs):
        yield seq

def read_seq_rev(l_chrom, l_pos, uflank, dflank, *args, **kwargs):
    for seq in read_seq(l_chrom, l_pos - uflank, l_pos + dflank + 1, *args, **kwargs):
        yield str(hts.Sequence(seq.encode('ascii')).get_reverse_complement())

def to_fa(fp_out_fa, df_fwd, df_rev, seq_col):
    with open(fp_out_fa, 'w') as fh_out_fa:
        for i, r in list(itertools.islice(df_fwd.iterrows(), None)):
            print('>%s:%s:%s:+' % (r['chrom'], r['start'], r['end']), file=fh_out_fa)
            print(r[seq_col], file=fh_out_fa)
        for i, r in list(itertools.islice(df_rev.iterrows(), None)):
            print('>%s:%s:%s:-' % (r['chrom'], r['start'], r['end']), file=fh_out_fa)
            print(r[seq_col], file=fh_out_fa)
    #!grep '^>' {fp_out_fa} | wc -l

# https://doi.org/10.1101/gr.153668.112
# For analyses where a single TSS position was required, we considered the distribution of cap 5' ends 
# within the TIC, and selected the position with the most tags (the mode). In the case of a tie (two or more 
# positions with the same number of tags), we selected the mode closest to the median of the TIC.
def nanargmax_median(a):
    a_max = np.nanmax(a)
    a_max_indices = np.flatnonzero(a == a_max)
    if len(a_max_indices) == 0:
        return 0
    return a_max_indices[len(a_max_indices) // 2]

assert nanargmax_median([1,2,3,3,3,2,1]) == 3
assert nanargmax_median([1,2,2,1]) == 2
assert nanargmax_median([1]) == 0
#assert nanargmax_median([float('nan'), float('nan'), float('nan')]) == 0

"""
mfc(np.array([
    [+1,-1,-1,-1,-1],
    [+1,+1,-1,-1,-1],
    [+1,+1,+1,-1,-1],
    [+1,+1,+1,+1,-1]]), 1)
"""

def skew(x, axis=None):
    x_p = np.nansum((x > 0), axis=axis)
    x_n = np.nansum((x < 0), axis=axis)    
    return (x_p - x_n) / (x_p + x_n)

def askew(x, axis=None):
    x_p = np.nansum((x > 0), axis=axis)
    x_n = np.nansum((x < 0), axis=axis)    
    return -(x_p - x_n) / (x_p + x_n)

def imshow_df(df, d_title_kwargs=None, cluster_col=None, figsize=None, kwargs_grid={}, *args, **kwargs):
    l_cluster_i = df[cluster_col]
    ll_cluster_i = [l_cluster_i[l_cluster_i == v].index.tolist() for v in OrderedCounter(l_cluster_i).keys()]

    kwargs_grid_ = {
        'nrows_ncols': (1, len(d_title_kwargs)),
        'label_mode': 'L',
        'aspect': False,
        'axes_pad': 0.1,
        'cbar_location': 'bottom',
        'cbar_mode': 'edge',
        'cbar_size': '2%',
        'cbar_pad': 0.1,
    }
    kwargs_grid_.update(kwargs_grid)

    fig = plt.figure(figsize=figsize)
    grid = ImageGrid(fig, 111, **kwargs_grid_)
    for i, (ax, (title, kwargs_col)) in enumerate(zip(grid, d_title_kwargs.items())):
        data = [np.mean(df.loc[cluster_i, title]) for cluster_i in ll_cluster_i]
        #ax.axis('off')
        im = ax.imshow(np.array([data]).T, aspect = 'auto', interpolation='none', **kwargs_col)
        ax.set_title(title)#, rotation=45)
        ax.set_xticks([])
        ax.set_yticks(range(len(ll_cluster_i)))
        ax.set_yticklabels([len(l_cluster_i) for l_cluster_i in ll_cluster_i])
        grid.cbar_axes[i].colorbar(im, ticks=[kwargs_col['vmin'], kwargs_col['vmax']])
        grid.cbar_axes[i].set_xticklabels([kwargs_col['vmin'], kwargs_col['vmax']], rotation=45)
        #grid.cbar_axes[i].axis('off')

def pie_by_enr(l_obs, l_exp, l_label, 
            cm=matplotlib.cm.get_cmap('coolwarm'),
            norm=matplotlib.colors.Normalize(vmin=-2.,vmax=2.),
            **kwargs):
    """
    Pie chart of l_obs, sectors coloured by (expected) distribution based on l_exp
    """
    l_enr = []
    l_col = []
    for (obs, exp) in zip(l_obs, l_exp):
        enr = np.log2((obs/sum(l_obs)) / (exp/sum(l_exp)))
        col = cm(norm(enr))
        l_enr.append(enr)
        l_col.append(col)
    plt.pie(l_obs, colors=l_col, labels=l_label, counterclock=False, autopct='%.2f%%', **kwargs)


"""
Choosing a good (default) colour map...
http://matplotlib.org/users/colormaps.html
http://web.stanford.edu/~mwaskom/software/seaborn/tutorial/color_palettes.html

http://www.sandia.gov/~kmorel/documents/ColorMaps/ColorMapsExpanded.pdf
> Fig. 1. The rainbow color map. Know thy enemy.
> Fig. 15. A continuous diverging color map well suited to scientific visualization.

http://matplotlib.org/examples/color/colormaps_reference.html
> coolwarm
Looks similar to Fig. 15.

http://earthobservatory.nasa.gov/blogs/elegantfigures/2013/08/05/subtleties-of-color-part-1-of-6/
http://earthobservatory.nasa.gov/blogs/elegantfigures/2013/08/12/subtleties-of-color-part-3-of-6/

http://blog.yhathq.com/posts/ggplot-for-python.html
http://blog.olgabotvinnik.com/prettyplotlib/

https://github.com/matplotlib/matplotlib/issues/875

http://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality
http://www.sandia.gov/~kmorel/documents/ColorMaps/
"""
class GenomicDataFrameTrack(object):
    def __init__(self, flank_len, bin_size):
        self.flank_len = int(flank_len)
        self.bin_size = int(bin_size)
        self.imshow_kwargs = {
            'aspect': 'auto',
            'cmap': matplotlib.cm.get_cmap("coolwarm"),#, N=1024), # A failed attempt at setting quantization levels
            'interpolation': 'nearest', # Really has an effect on whether you can see faint signals
            'vmin': int(-3),
            'vmax': int(+3),
        }
        self.smooth = None #lambda m: m

    def get_subset(self, l_i):
        """
        Create a new GR based on indices of the old. Potential uses: one cluster, TOP20%, etc
        """
        subset = copy_module.deepcopy(self)
        subset.m = self.m[l_i,]
        return subset

    def flexible_chrom(self, chrom, chroms):
        """
        TODO:
        1) move to bbifile; rename to flexi_chrom
        2) add rule: 'if contains M'
        3) add caching (both for overview as well as speed)
        """
        if chrom in chroms:
            return chrom
        elif chrom.startswith('chr') and chrom[3:] in chroms:
            return chrom[3:]
        elif 'chr' + chrom in chroms:
            return 'chr' + chrom
        else:
            raise Exception('chrom=%(chrom)s not found in %(chroms)s' % locals())

    def m_from_bw(self, chroms, poss, fp, bin_skew=False, f_bin=np.mean, memoized=True, rememoize=False):
        flank_adj = self.flank_len + int(self.bin_size) // 2
        l_chrom = chroms
        l_start = []
        l_end = []
        for pos in poss:
            l_start.append(pos - flank_adj)
            l_end.append(pos + flank_adj + 1)
        t0 = datetime.datetime.now()
        (n_pass, n_fail) = (0, 0)
        self.m = np.zeros([len(l_start), self.n_bins])

        if not bin_skew:
            def bin_signal(raw_signal_at_iv): return f_bin(raw_signal_at_iv.reshape(-1, self.bin_size), 1)
        elif bin_skew:
            def bin_signal(raw_signal_at_iv): return skew(raw_signal_at_iv.reshape(-1, self.bin_size), 1)
        
        if memoized:
            for (i, bin_signal_at_iv) in enumerate(itertools.islice(map(bin_signal, mread_regions(fp, l_chrom, l_start, l_end, rememoize=rememoize)), None)):
                self.m[i,] = bin_signal_at_iv
                n_pass += 1
        else:
            for (i, bin_signal_at_iv) in enumerate(itertools.islice(map(bin_signal, read_regions(fp, l_chrom, l_start, l_end)), None)):
                self.m[i,] = bin_signal_at_iv
                n_pass += 1
        
        # brute-force NaN
        self.m = np.nan_to_num(self.m)
        wct = (datetime.datetime.now() - t0)
        #print "m_from_bw: wct=%(wct)s\tn_pass=%(n_pass)d\tn_fail=%(n_fail)d" % locals()

    def m_from_2bit(self, chroms, poss, fp, d_skew, window_len=100):
        flank_adj = self.flank_len + window_len // 2#+ int(self.bin_size) / 2
        l_chrom = chroms
        l_start = []
        l_end = []
        for pos in poss:
            l_start.append(pos - flank_adj)
            l_end.append(pos + flank_adj)# + 1)
        t0 = datetime.datetime.now()
        (n_pass, n_fail) = (0, 0)

        def rolling_skew(a, window):
            def rolling_window(a, window):
                #http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
                shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
                strides = a.strides + (a.strides[-1],)
                return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
            def skew_window(a):
                return np.sum(a) / float(np.sum(np.abs(a)))
            return list(map(skew_window, rolling_window(np.array(a), window=window)))

        self.m = np.array([
            rolling_skew(np.array(map(lambda seq_i: d_skew.get(seq_i, 0), seq)), window_len)
            for (i, seq) in enumerate(itertools.islice(SeqTrack.read_seq(l_chrom, l_start, l_end, fp), None))
        ])

    def m_strandify(self, strands):
        for i, strand in enumerate(strands):
            if strand == '-':
                self.m[i,:] = self.m[i,:][::-1]

    def m_from_ga(self, ga, chroms, starts, strand='.'):
        t0 = datetime.datetime.now()
        self.m = np.zeros([len(chroms), self.n_bins])
        (n_pass, n_fail) = (0, 0)
        for (i, (chrom, start)) in enumerate(zip(chroms, starts)):
            try:
                # Extract and downsample signal
                flank_adj = self.flank_len + int(self.bin_size) // 2
                iv = hts.GenomicInterval(chrom, start - flank_adj, start + flank_adj + 1, strand=strand)
                raw_signal_at_iv = np.array(list(ga[iv]))
                bin_signal_at_iv = np.mean(raw_signal_at_iv.reshape(-1, self.bin_size), 1)
                # Flip signal at reverse strand regions if desired
                #if stranded and iv.strand == '-':
                #    bin_signal_at_iv = bin_signal_at_iv[::-1]
                # Assign to array
                self.m[i,] = bin_signal_at_iv
                n_pass += 1
            # Obsolete at this point:
            except IndexError:
                # np.NAN would be more appropriate, but sklearn tends to barf on these
                self.m[i,] = 0
                n_fail += 1
        wct = (datetime.datetime.now() - t0)
        #print "m_from_ga: wct=%(wct)s\tn_pass=%(n_pass)d\tn_fail=%(n_fail)d" % locals()

    def m_from_ga_sum(self, ga, chroms, starts, strand='.'):
        t0 = datetime.datetime.now()
        self.m = np.zeros([len(chroms), self.n_bins])
        (n_pass, n_fail) = (0, 0)
        for (i, (chrom, start)) in enumerate(zip(chroms, starts)):
            try:
                # Extract and downsample signal
                flank_adj = self.flank_len + int(self.bin_size) // 2
                iv = hts.GenomicInterval(chrom, start - flank_adj, start + flank_adj + 1, strand=strand)
                raw_signal_at_iv = np.array(list(ga[iv]))
                bin_signal_at_iv = np.sum(raw_signal_at_iv.reshape(-1, self.bin_size), 1)
                # Flip signal at reverse strand regions if desired
                #if stranded and iv.strand == '-':
                #    bin_signal_at_iv = bin_signal_at_iv[::-1]
                # Assign to array
                self.m[i,] = bin_signal_at_iv
                n_pass += 1
            # Obsolete at this point:
            except IndexError:
                # np.NAN would be more appropriate, but sklearn tends to barf on these
                self.m[i,] = 0
                n_fail += 1
        wct = (datetime.datetime.now() - t0)
        #print "m_from_ga: wct=%(wct)s\tn_pass=%(n_pass)d\tn_fail=%(n_fail)d" % locals()

    def m_from_bam(self, r, fp, mode='coverage', d_rchrom_tchrom=None, stranded=True, fsize_min=None, fsize_max=None, fetch_fsize_max=300):
        """
        mode is either coverage, fwd_cuts, rev_cuts
        filter_fsize_min
        max_fsize is to properly query boundaries
        """
        t0 = datetime.datetime.now()
        self.m = np.zeros([len(r), self.n_bins])
        self.d_fsize_count = collections.Counter()
        (n_pass, n_fail) = (0, 0)
        with pysam.AlignmentFile(fp, 'rb') as fh:
            for (iv_i, chrom, start, end) in zip(range(len(r)), r['chrom'], r['start'], r['end']):
                query_start = start - self.flank_len
                query_end = start + self.flank_len
                for (aln_i, aln) in itertools.islice(enumerate(fh.fetch(chrom, query_start - fetch_fsize_max, query_end + fetch_fsize_max)), None):
                    assert aln.template_length > 0
                    f_fsize_min = (fsize_min is None) or (fsize_min <= aln.template_length)
                    f_fsize_max = (fsize_max is None) or (aln.template_length < fsize_max)
                    if f_fsize_min and f_fsize_max:
                        start_i = aln.reference_start - query_start
                        end_i = start_i + aln.template_length

                        if mode == 'coverage':
                            start_i = max(start_i, 0)
                            end_i = min(end_i, self.n_bins)
                            if start_i < end_i: # check if at least partially within query range
                                self.m[iv_i, start_i:end_i] += 1
                                self.d_fsize_count[aln.template_length] += 1
                        elif mode == 'fwd_cuts':
                            if 0 <= start_i and start_i < self.n_bins:
                                self.m[iv_i, start_i] += 1
                        elif mode == 'rev_cuts':
                            if 0 <= end_i and end_i < self.n_bins:
                                self.m[iv_i, end_i] += 1
                        else:
                            assert False, 'mode=%(mode)s unknown' % locals()
                n_pass += 1
        wct = (datetime.datetime.now() - t0)
        print("m_from_bam: wct=%(wct)s\tn_pass=%(n_pass)d\tn_fail=%(n_fail)d" % locals())
        #       except IndexError:
        #           n_fail += 1

    def hist_fsize(self, max_fsize=300, bin_width=1, *args, **kwargs):
        # only generated if mode=coverage at the moment
        bins = range(0, max_fsize + 1, bin_width)
        return plt.hist(x=self.d_fsize_count.keys(), bins=bins, weights=self.d_fsize_count.values(), histtype="step", *args, **kwargs)[-1]

    @property
    def n_raw_bins(self):
        return int(2 * self.flank_len + self.bin_size)

    @property
    def n_bins(self):
        return int((2 * self.flank_len + self.bin_size) // self.bin_size)

    @property
    def xlim(self):
        flank_len_ext = self.flank_len# + self.bin_size / 2
        return (-flank_len_ext, +flank_len_ext)

    @property
    def imshow_extent(self):
        flank_len_ext = self.flank_len + self.bin_size // 2
        return (-flank_len_ext, +flank_len_ext, self.m.shape[0] + 0.5, 0.5)

    def imshow_extent_vfixed(self, vextent):
        return (self.imshow_extent[0], self.imshow_extent[1], 0, vextent)

    def imshow(self, ax=None, set_x_axis=True, vextent=None, nsquashed=None):
        if ax is None: ax = plt.gca()

        if self.smooth is None: # smooth via lambda somehow screws up interpolation...
            m_ = self.m
        else:
            m_ = self.smooth(self.m)

        if not(nsquashed is None):
            m_ = np.array([np.nanmedian(chunk, axis=0) for chunk in np.array_split(m_, nsquashed)])                

        if vextent is None:
            extent_ = self.imshow_extent
        else:
            extent_ = self.imshow_extent_vfixed(vextent)

        r = ax.imshow(m_, extent=extent_, **self.imshow_kwargs)

        if set_x_axis:
            ax.set_xlim([-self.flank_len, self.flank_len])
            ax.set_xticks([-self.flank_len, 0, +self.flank_len])
            ax.set_xticklabels([-self.flank_len, 0, "+%s" % (self.flank_len,)], rotation=45)
        return r

    @property
    def x(self):
        return range(-self.flank_len, +self.flank_len + self.bin_size, self.bin_size)

    def plot(self, ax=None, f=np.mean, set_x_axis=True, f_h=None, *args, **kwargs):
        if ax is None: ax = plt.gca()
        y_ = f(self.m, 0)
        if not(f_h is None): y_ = f_h(y_)
        r = ax.plot(self.x, y_, drawstyle='steps-mid', *args, **kwargs)
        ax.set_xlim(-self.flank_len, +self.flank_len)
        return r
    
    def plot_hs(self, ax=None, f=np.mean, set_x_axis=True, f_h=None, smooth_width=50, *args, **kwargs):
        def rolling_mean_kernel(width): return np.ones(width) / float(width)
        f_h50 = lambda y: np.convolve(y, rolling_mean_kernel(smooth_width), mode='same')
        self.plot(ax=ax, f=f, set_x_axis=set_x_axis, f_h=f_h50, *args, **kwargs)      

    def plotq(self, ax=None, q=5, q_i=None, set_x_axis=True, cmap=matplotlib.cm.get_cmap("viridis"), f=np.median, *args, **kwargs):
        if ax is None: ax = plt.gca()
        for (i, m_i) in zip(range(q), np.array_split(self.m, q)):
            if not(q_i is None) and (q_i != i):
                continue
            if q_i is None:
                #http://stackoverflow.com/questions/15140072/how-to-map-number-to-color-using-matplotlibs-colormap
                color_i = cmap(quantile_mid(i, q))
            else:
                color_i = 'k'
            r = ax.plot(self.x, f(m_i, 0),
                        label="%d of %d" % (i + 1, q),
                        color=color_i,
                        drawstyle='steps-mid',
                        *args, **kwargs)
        ax.set_xlim(-self.flank_len, +self.flank_len)
        return r

class GenomicDataFrameDendrogramTrack(GenomicDataFrameTrack):
    def __init__(self, df, bin_size=1000):
        self.bin_size = int(bin_size)
        #self.imshow_kwargs = {
        #    'aspect': 1/bin_size,# / (float(2 * self.flank_len + 1)) / self.bin_size, # GenomicDataFrameTrack.n_bins
        #    'cmap': matplotlib.cm.get_cmap('coolwarm'),#, N=1024), # A failed attempt at setting quantization levels
        #    'interpolation': 'none',
        #    'vmin': -1,
        #    'vmax': +2,
        #}
        self.imshow_kwargs = {
            'aspect': 'auto',
            'cmap': matplotlib.cm.get_cmap("coolwarm"),#, N=1024), # A failed attempt at setting quantization levels
            'interpolation': 'nearest', # Really has an effect on whether you can see faint signals
            'vmin': int(-1),
            'vmax': int(+1),
        }

        self.df_raw = df.copy()        
        self.m_raw = df.values.copy()
        self.m = self.m_raw

    @property
    def imshow_extent(self):
        return (0, self.bin_size, self.m.shape[0] + 0.5, 0.5)

    def imshow_extent_vfixed(self, vextent):
        return (self.imshow_extent[0], self.imshow_extent[1], 0, vextent)

    def imshow(self, ax=None, set_x_axis=True, vextent=None, nsquashed=None):
        if ax is None: ax = plt.gca()

        m_ = self.m

        if not(nsquashed is None):
            m_ = np.array([chunk.mean(axis=0) for chunk in np.array_split(m_, nsquashed)])                

        if vextent is None:
            extent_ = self.imshow_extent
        else:
            extent_ = self.imshow_extent_vfixed(vextent)

        r = ax.imshow(m_, extent=extent_, **self.imshow_kwargs)

        #if set_x_axis:
        #    ax.set_xlim([-self.flank_len, self.flank_len])
        #    ax.set_xticks([-self.flank_len, 0, +self.flank_len])
        #    ax.set_xticklabels([-self.flank_len, 0, "+%s" % (self.flank_len,)], rotation=45)
        return r

    def imshow_dendrogram(self, ax=None, fig=None, set_x_axis=True, padding_ = -0.035, linkage_kwargs={}, dendrogram_kwargs={}):
        if ax is None: ax = plt.gca()
        if fig is None: fig = plt.gcf()
        
        ax_pos = ax.get_position()
        ax_dendrogram = fig.add_axes([ax_pos.x0, ax_pos.y1+padding_, ax_pos.x1-ax_pos.x0, 1-(ax_pos.y1+padding_)])

        self.linkage_ = sp.cluster.hierarchy.linkage(self.m_raw.T, **linkage_kwargs)
        self.dendrogram_ = sp.cluster.hierarchy.dendrogram(self.linkage_, labels=self.df_raw.columns, ax=ax_dendrogram, **dendrogram_kwargs)
        self.m = self.m_raw[:,self.dendrogram_['leaves']].copy()
        self.labels = self.df_raw.columns[self.dendrogram_['leaves']]

        imshow_kwargs_ = {
            'aspect': 1/self.bin_size,# / (float(2 * self.flank_len + 1)) / self.bin_size, # GenomicDataFrameTrack.n_bins
            'cmap': matplotlib.cm.get_cmap('coolwarm'),#, N=1024), # A failed attempt at setting quantization levels
            'interpolation': 'none',
            'vmin': -1,
            'vmax': +2,
        }

        ax_dendrogram.set_xticks([])
        r = ax.imshow(self.m, **imshow_kwargs_)
        ax.set_xticks(list(range(self.m.shape[1])))
        ax.set_xticklabels(self.labels, rotation=90)
        return r

class GenomicDataFrameScalarTrack(GenomicDataFrameTrack):
    def __init__(self, s, flank_len, bin_size):
        self.flank_len = int(flank_len)
        self.bin_size = int(bin_size)
        self.imshow_kwargs = {
            'aspect': 'auto',# / (float(2 * self.flank_len + 1)) / self.bin_size, # GenomicDataFrameTrack.n_bins
            'cmap': matplotlib.cm.get_cmap('viridis'),#, N=1024), # A failed attempt at setting quantization levels
            'interpolation': 'none',
            'vmin': +0,
            'vmax': +1,
        }
        self.m = np.zeros([len(s), 1])
        for i, s_i in enumerate(s):
            self.m[i,0] = s_i

    def imshow(self, ax=None, set_x_axis=True, vextent=None, nsquashed=None):
        if ax is None: ax = plt.gca()

        m_ = self.m
        if not(nsquashed is None):
            m_ = np.array([np.nanmedian(chunk, axis=0) for chunk in np.array_split(m_, nsquashed)])

        if vextent is None:
            extent_ = self.imshow_extent
        else:
            extent_ = self.imshow_extent_vfixed(vextent)

        r = ax.imshow(m_, extent=extent_, **self.imshow_kwargs)

        if True:#set_x_axis:
            ax.set_xlim([-self.flank_len, self.flank_len])
            ax.set_xticks([])
            ax.set_xticklabels([])
        return r

class GenomicDataFramePeaksTrack(GenomicDataFrameTrack):
    """
    Store overlap data between source peaks (=gdf regions)  and "target" peaks (=peaks in the track)

    Overlap data stored as:
        .p -- raw "target" peaks
        .h -- hits for source peaks in target (DataFrame of objects)
        .m -- "binary" matrix (use colours to indicate whether overlaps?)

    def __init__(ga=None, fp_bw=None, fp_wig=None, ...):
        assert(only one is not None):

    rename classmethods -- single constructor, but:
        m_from_ga # read matrix from ga
        m_from_bw # read matrix from bw
    Visualise a bed-file -- with overlaps
    """
    def __init__(self, flank_len, bin_size):
        self.flank_len = int(flank_len)
        self.bin_size = int(bin_size)
        self.imshow_kwargs = {
            'aspect': 'auto',
            'cmap': matplotlib.cm.get_cmap('RdBu'),#, N=1024), # A failed attempt at setting quantization levels
            'interpolation': 'none',
            'vmin': -1,
            'vmax': +1,
        }
        self.summary_str_color='k'
        self.plot_summary_str = True

    def from_df(self, r, df, column=None, *args, **kwargs):
        """
        This assumes
        1) peaks are unstranded; might want to re-write later as strandedness may be meaningful (e.g. TICs)
        2) non-overlapping features (=e.g. typical ChIP peak calls)
        """
        self.p = df.copy()
        # Build ga of 'target peaks'
        ga_p = htw.GenomicArray(chroms='auto', stranded=False, typecode='O')
        for (i, chrom, start, end) in zip(range(len(self.p)), self.p['chrom'], self.p['start'], self.p['end']):
            ga_p[htw.GenomicInterval(chrom, start, end)] = (chrom, start, end)

        # Build hits list
        self.h = [] # Re-factor as series??
        for (i, chrom, start, end) in zip(range(len(r)), r['chrom'], r['start'], r['end']):
            iv = htw.GenomicInterval(chrom, start, end)
            self.h.append(htw.ga_iv_get_union(ga_p, iv))

        # Build "heatmap array"
        ga_m = htw.GenomicArray(chroms='auto', stranded=False, typecode='d')
        if not(column is None):
            for (i, chrom, start, end, score) in zip(range(len(self.p)), self.p['chrom'], self.p['start'], self.p['end'], self.p[column]):
                ga_m[htw.GenomicInterval(chrom, start, end)] = score
        elif 'strand' in df.columns.values:
            for (i, chrom, start, end, strand) in zip(range(len(self.p)), self.p['chrom'], self.p['start'], self.p['end'], self.p['strand']):
                if strand == '-':
                    ga_m[htw.GenomicInterval(chrom, start, end)] = -1
                else:
                    ga_m[htw.GenomicInterval(chrom, start, end)] = +1
        else:
            for (i, chrom, start, end) in zip(range(len(self.p)), self.p['chrom'], self.p['start'], self.p['end']):
                ga_m[htw.GenomicInterval(chrom, start, end)] = 1

        # Add peak signal to plotter
        self.m_from_ga(r, ga_m)
        #print "add_from_bed: %d features, overlap: %s" % (self.n_bed(title), self.bed_overlap_summary_str(title))

    def get_subset(self, l_i):
        """
        Create a new GR based on indices of the old. Potential uses: one cluster, TOP20%, etc
        """
        #subset = super(GenomicDataFramePeaksTrack, self).get_subset(l_i)
        subset = copy_module.deepcopy(self)
        subset.m = self.m[l_i,]
        subset.p = self.p.copy()
        subset.h = [self.h[i] for i in l_i]
        subset._superset = self # Adhoc way to have enrichment score calculation available; obviously kludge-y
        return subset

    @property
    def n_r(self):
        return len(self.h)

    @property
    def n_p(self):
        return len(self.p)

    @property
    def n_rwp(self):
        return sum([1 if len(h_i) > 0 else 0 for h_i in self.h]) # number of rows that overlap with a BED feature

    @property
    def n_pwr(self):
        return len(set().union(*self.h)) # number of BED features that overlap with a row

    @property
    def sens(self):
        return self.n_rwp / float(self.n_r)

    @property
    def prec(self):
        return self.n_pwr / float(self.n_p)

    def __str__(self):
        return "%d(%d%%)/%d(%d%%)" % (self.n_rwp, 100*self.sens, self.n_pwr, 100*self.prec)

    def imshow(self, *args, **kwargs):
        r = super(GenomicDataFramePeaksTrack, self).imshow(*args, **kwargs)
        ax = kwargs.get('ax', plt.gca())
        try:
            regns_total = self._superset.n_r
            regns_clust = self.n_r
            peaks_total = self._superset.n_p
            peaks_clust = self.n_pwr
            (oddsratio, p_value) = sp.stats.fisher_exact([
                [peaks_clust, peaks_total - peaks_clust],
                [regns_clust, regns_total - regns_clust],
            ], 'greaterq')
            summary_str = "%d(%d%%, %.2fx) /\n%d(%d%%)\n%f %f" % (self.n_rwp, 100*self.sens, self.sens / self._superset.sens, self.n_pwr, 100*self.prec, oddsratio, p_value)
            #print regns_total, regns_clust, peaks_total, peaks_clust, oddsratio, p_value
        except:
            if self._superset.sens > 0:
                summary_str = "%d(%d%%, %.2fx) /\n%d(%d%%)" % (self.n_rwp, 100*self.sens, self.sens / self._superset.sens, self.n_pwr, 100*self.prec)
            else:
                summary_str = ''

        if self.plot_summary_str:
            ax.text(self.n_bins // 2, self.n_r // 2, summary_str, ha='center', va='center', color=self.summary_str_color, fontsize=3)

        #print "summary_str: %(summary_str)s" % locals()
        # TODO: add p-values to enrichment?
        # http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.binom.html
        # Or hypergeometric things like
        #print sp.stats.fisher_exact([[500, 2], [345, 10]])
        return r # To avoid known AttributeError: 'NoneType' object has no attribute 'autoscale_None'

class SeqTrack(object):
    def __init__(self, flank_len, bin_size):
        self.flank_len = int(flank_len)
        self.bin_size = int(bin_size)
        assert bin_size == 1

    def get_subset(self, l_i):
        """
        Create a new GR based on indices of the old. Potential uses: one cluster, TOP20%, etc
        """
        subset = copy_module.deepcopy(self)
        subset.m = [self.m[l_i_j] for l_i_j in l_i]#self.m[l_i,]
        return subset

    @staticmethod
    def read_seq(chroms, starts, ends, fp=None, genome=None):
        if fp is None:
            fp = '.wget/hgdownload.cse.ucsc.edu/goldenPath/%(genome)s/bigZips/%(genome)s.2bit' % locals()
        assert os.path.isfile(fp)
        fh = twobitreader.TwoBitFile(fp)
        for chrom, start, end in itertools.islice(zip(chroms, starts, ends), None):
            yield str(fh[chrom][start:end])

    @classmethod
    def from_2bit(cls, chroms, starts, flank_len, genome):
        self = cls(flank_len=flank_len, bin_size = 1)
        self.m = list(cls.read_seq(chroms, [start - flank_len for start in starts], [start + flank_len + 1 for start in starts], genome=genome))
        return self

    @property
    def n_raw_bins(self):
        return int(2 * self.flank_len + self.bin_size)

    @property
    def n_bins(self):
        return int((2 * self.flank_len + self.bin_size) // self.bin_size)

    @property
    def xlim(self):
        flank_len_ext = self.flank_len# + self.bin_size // 2
        return (-flank_len_ext, +flank_len_ext)

    @property
    def x(self):
        return range(-self.flank_len, +self.flank_len + self.bin_size, self.bin_size)

    @staticmethod
    def seq_skew(s, d = {'C': 1, 'T': 1, 'A': -1, 'G': -1}, min_len=1):
        for k, g in itertools.groupby(s, lambda i: d.get(i, 0)):
            g_n = len(list(g))
            for _ in itertools.repeat(g_n, g_n):
                if g_n >= min_len:
                    yield k * g_n
                else:
                    yield 0

    CG_SKEW = {'C': 1, 'G': -1, 'c': 1, 'g': -1}
    GC_SKEW = {'G': 1, 'C': -1, 'g': 1, 'c': -1}
    TA_SKEW = {'T': 1, 'A': -1, 't': 1, 'a': -1}
    YR_SKEW = {'C': 1, 'T': 1, 'A': -1, 'G': -1, 'c': 1, 't': 1, 'a': -1, 'g': -1}
    def plot_skew(self, d=None, min_len=1, v=None, ax=None, y_pp = None, **kwargs):
        m_skew = np.array(list(map(lambda m_i: list(self.seq_skew(s=m_i, d=d, min_len=min_len)), self.m)))
        mean_diff = np.mean(np.sign(m_skew), 0) / np.mean(np.abs(np.sign(m_skew)), 0)
        if not(y_pp is None):
            mean_diff = y_pp(mean_diff)
        if ax is None: ax = plt.gca()
        r = ax.plot(self.x, mean_diff, drawstyle='steps-mid', **kwargs)
        ax.set_xlim(-self.flank_len, +self.flank_len)
        return r

class GenomicDataFrame(object):
    """
    Store and visualise genomic signal at a chosen set of "interesting regions"

    Generally intended to be used as:
        gh = htw.GenomicRegions(l_gf)
        gh.add_regions(ga, title)
        gh.decompose() # manually include PCA component?
        gh.cluster() # manually set clustering algorithm?
        gh.plot()
            # option to set manual cluster order
            # option to sort by pca weights within the cluster
        gh.plot_subset() # provide a function to select for subsets and plot these?

    TODO: de-couple sub-selecting everything in a cluster, and sub-plotting all clusters
    (Simpler for e.g. aggregate plots etc...)
    """

    def __init__(self, df_r, pos_column=None, strand_column=None, genome='ce10', v=True):
        self.r = df_r.copy() # Stores region location and other metadata; row ordering used for plotting order
        self.t = collections.OrderedDict() # Stores signal at regions from track files

        # Force index to be 0...len-1
        self.r.reset_index(drop=True, inplace=True)

        if pos_column in df_r.columns:
            self.pos_column = pos_column
        else:
            if v: print('gdf: pos_column unspecified, reverting to mean(start, end)')
            self.pos_column = '_mean_start_end'
            self.r[self.pos_column] = list(map(int, self.r[['start', 'end']].mean(axis=1)))

        self.strand_column = strand_column
        self.genome = genome

        assert 'chrom' in self.r.columns
        assert self.pos_column in self.r.columns


    def get_subset(self, l_i, l_j=None): # column indices work on tracks only; great for quick colour maps
        self_subset = GenomicDataFrame(self.r.iloc[l_i], self.pos_column, self.strand_column, self.genome)
        self_subset.r.reset_index(inplace=True) # should ideally be avoided...
        self_subset.genome = self.genome
        for (title, track) in self.t.items():
            if (l_j is None) or (title in l_j):
                self_subset.t[title] = self.t[title].get_subset(l_i)
        return self_subset

    def query(self, expr):
        return self.get_subset(self.r.query(expr).index.tolist())

    def __getitem__(self, q): # Syntactic sugar experiment (TODO? add regex for columns?)
        if type(q) is tuple and len(q) == 2:
            (q_row, l_col) = q[0], q[1]
        else:
            q_row = q
            l_col = None

        if type(q_row) is np.ndarray or type(q_row) is list:
            l_row = q_row
        else:
            l_row = self.r.query(q_row).index.values
        return self.get_subset(l_row, l_col)

    def sort(self, columns=None, ascending=True):
        # For now, leave index as is (so anything that uses index works fluently)
        # Ideally, this should somehow be index-agnostic
        # Ideally, re-written by using get_subset
        self_sorted = GenomicDataFrame(self.r, self.pos_column, self.strand_column, self.genome)
        self_sorted.r.sort_values(by=columns, ascending=ascending, inplace=True)
        for (title, track) in self.t.items():
            #print title, track
            self_sorted.t[title] = self.t[title].get_subset(self_sorted.r.index.tolist())
        self_sorted.r.reset_index(inplace=True, drop=True)
        return self_sorted

    def add_seq(self, flank_len=500):
        self.t['_seq'] = SeqTrack.from_2bit(self.r.chrom, self.r[self.pos_column], flank_len=flank_len, genome=self.genome)

    def add_track(self, title=None, track=None, flank_len=500, bin_size=25, strandify=False, v=False, d_rchrom_tchrom=None, force=False, pos_column=None, *args, **kwargs):

        try:
            is_file = os.path.isfile(track)
            ext = os.path.splitext(track)[1]
            if title is None:
                title = os.path.splitext(os.path.basename(track))[0]
        except:
            is_file = False

        try:
            is_df = 'chrom' in track.columns and 'start' in track.columns and 'end' in track.columns
        except:
            is_df = False

        #if title in self.t.keys() and not force:
        #    if v: print('add_track: skip %(title)s' % locals())
        #    return

        if pos_column is None:
            pos_column = self.pos_column

        if is_file:
            self.t[title] = GenomicDataFrameTrack(flank_len=flank_len, bin_size=bin_size)
            if ext == '.bw':
                self.t[title].m_from_bw(self.r['chrom'], self.r[pos_column], track, *args, **kwargs)
            elif ext == '.2bit':
                self.t[title].m_from_2bit(self.r['chrom'], self.r[pos_column], track, *args, **kwargs)

            #, stranded=self.stranded, d_rchrom_tchrom=d_rchrom_tchrom, *args, **kwargs)
            if strandify:
                if v: print('strandify')
                self.t[title].m_strandify(self.r[self.strand_column])

            # Outdated atm
            #elif ext == '.bam':
            #    self.t[title] = GenomicDataFrameTrack(flank_len=flank_len, bin_size=bin_size)
            #    self.t[title].m_from_bam(self.r, track, d_rchrom_tchrom=d_rchrom_tchrom, *args, **kwargs)

        #elif is_df:
        #    self.t[title] = GenomicDataFramePeaksTrack(flank_len=flank_len, bin_size=bin_size)
        #    self.t[title].from_df(self.r, track, *args, **kwargs)
        #
        # TODO: is_ga
        #elif is_ga
        else:
            assert(False)

    def add_scanMotifGenomeWide(self, title_fwd, title_rev, track, flank_len=500, bin_size=1, ntop=None, pos_column=None, type_=None):
        if pos_column is None:
            pos_column = self.pos_column
        fp_motif = track
        df_motif = pd.read_csv(fp_motif, sep='\t', comment='#', names=NAMES_BED6)\
            .sort_values('score', ascending=False).head(n=ntop)\
            .sort_values(['chrom', 'start', 'end', 'strand'])\
            .reset_index(drop=True)
        df_motif['start'] = df_motif['start'] - 1

        ga_motif_fwd = hts.GenomicArray(chroms=chroms_ce10, stranded=False, typecode='i', storage='ndarray')
        for i, r in df_motif.query("strand == '+'").iterrows():
            if type_ == "5'":
                iv = hts.GenomicInterval(r['chrom'], r['start'], r['start'] + 1, r['strand'])
            elif type_ == "3'":
                iv = hts.GenomicInterval(r['chrom'], r['end'] - 1, r['end'], r['strand'])
            elif is_int(type_):
                iv = hts.GenomicInterval(r['chrom'], r['start'] + type_, r['start'] + type_ + 1, r['strand'])
            else:
                iv = hts.GenomicInterval(r['chrom'], r['start'], r['end'], r['strand'])
            ga_motif_fwd[iv] += 1

        ga_motif_rev = hts.GenomicArray(chroms=chroms_ce10, stranded=False, typecode='i', storage='ndarray')
        for i, r in df_motif.query("strand == '-'").iterrows():
            if type_ == "5'":
                iv = hts.GenomicInterval(r['chrom'], r['end'] - 1, r['end'], r['strand'])
            elif type_ == "3'":
                iv = hts.GenomicInterval(r['chrom'], r['start'], r['start'] + 1, r['strand'])
            elif is_int(type_):
                iv = hts.GenomicInterval(r['chrom'], r['end'] - type_ - 1, r['end'] - type_, r['strand'])
            else:
                iv = hts.GenomicInterval(r['chrom'], r['start'], r['end'], r['strand'])
            ga_motif_rev[iv] += 1
    
        self.t[title_fwd] = GenomicDataFrameTrack(flank_len=flank_len, bin_size=bin_size)
        self.t[title_rev] = GenomicDataFrameTrack(flank_len=flank_len, bin_size=bin_size)
        self.t[title_fwd].m_from_ga_sum(ga=ga_motif_fwd, chroms=self.r['chrom'], starts=self.r[pos_column])
        self.t[title_rev].m_from_ga_sum(ga=ga_motif_rev, chroms=self.r['chrom'], starts=self.r[pos_column])

    def add_fimo(self, title_fwd, title_rev, track, flank_len=500, pos_column=None, type_=None):
        if pos_column is None:
            pos_column = self.pos_column

        fp_motif = track
        df_motif = pd.read_csv(fp_motif, sep='\t', comment='#', names=NAMES_GTF)
        df_motif['start'] = df_motif['start'] - 1

        ga_motif_fwd = hts.GenomicArray(chroms=chroms_ce10, stranded=False, typecode='i', storage='ndarray')
        for i, r in df_motif.query("strand == '+'").iterrows():
            iv = hts.GenomicInterval(r['chrom'], r['start'], r['end'])
            ga_motif_fwd[iv] += 1

        ga_motif_rev = hts.GenomicArray(chroms=chroms_ce10, stranded=False, typecode='i', storage='ndarray')
        for i, r in df_motif.query("strand == '-'").iterrows():
            iv = hts.GenomicInterval(r['chrom'], r['start'], r['end'])
            ga_motif_rev[iv] += 1

        self.t[title_fwd] = GenomicDataFrameTrack(flank_len=flank_len, bin_size=1)
        self.t[title_rev] = GenomicDataFrameTrack(flank_len=flank_len, bin_size=1)
        self.t[title_fwd].m_from_ga(ga=ga_motif_fwd, chroms=self.r['chrom'], starts=self.r[pos_column])
        self.t[title_rev].m_from_ga(ga=ga_motif_rev, chroms=self.r['chrom'], starts=self.r[pos_column])

    @property
    def strands(self):
        return self.r[self.strand_column]

    def add_stranded(self, title_fwd, title_rev, track_fwd, track_rev, reverse_polarity=False, *args, **kwargs):
        self.add_track(title_fwd, track_fwd, *args, **kwargs)
        self.add_track(title_rev, track_rev, *args, **kwargs)
        m_fwd = self.t[title_fwd].m
        m_rev = self.t[title_rev].m
        for i, strand in enumerate(self.strands):
            if strand == '-':
                # temp = m_fwd[i,:]
                # m_fwd[i,:] = m_rev[i,:][::-1]
                # m_rev[i,:] = temp[::-1]
                (m_fwd[i,:], m_rev[i,:]) = (m_rev[i,:][::-1], m_fwd[i,:][::-1])
                if reverse_polarity:
                    m_fwd[i,:] = -m_fwd[i,:]
                    m_rev[i,:] = -m_rev[i,:]

    def add_tracks_from_df(self, df_tracks):
        for index, (title, flank_len, bin_size, flags, fp) in itertools.islice(df_tracks.iterrows(), None):
            print(title, flank_len, bin_size, flags)
            args_add_track = {
                'title': title,
                'track': fp,
                'bin_size': bin_size,
                'flank_len': flank_len,
            }
            if 'ce10' in flags:
                args_add_track['d_rchrom_tchrom'] = d_ensembl_ucsc
            self.add_track(**args_add_track)
            if 'zscore' in flags: self.gdf_zscore(title, fp)
        return df_tracks

    def kmeans(self, l_title=None, normalize=True, *args, **kwargs):
        if l_title is None: l_title = self.t.keys()
        if normalize:
            m = sklearn.preprocessing.normalize(self.get_m_stacked(l_title))
        else:
            m = self.get_m_stacked(l_title)
        km = sklearn.cluster.KMeans(init='k-means++', *args, **kwargs)
        km.fit(m)
        self.r['_'] = km.labels_ # By default, store clustering results as in the '_' column
        return self

    def imshow(self, l_title=None, l_title_cbar=None, cluster_col=None, figsize=None, savefig=False, fontsize=3,
            kwargs_grid={'cbar_size': '5%', 'cbar_pad': 0.5, 'axes_pad': 0.3,}, 
            show_cluster_col=True, row_labels=None, 
            axis_off=False, suptitle=None, 
            vextent=None, nsquashed=None, figB_magic=False, figB_suppl_magic=False, thesis_magic=False, *args, **kwargs):
        """
            #if d_title_imshow_kwargs is None: d_title_imshow_kwargs = collections.defaultdict(dict, {})
            #d_title_imshow_kwargs=None, clust_hline=True, auto_figsize=True, **kwargs):
            Potentially re-write such that every cluster would have its own imshow(), padding, etc...
            TOOD d_title_imshow_kwargs: use regex in title
            TODO Might want to add a warning for more than e.g. 50 clusters?
        """
        if cluster_col is None:
            ll_cluster_i = [self.r.index]
        else:
            l_cluster_i = self.r[cluster_col]
            ll_cluster_i = [l_cluster_i[l_cluster_i == v].index.tolist() for v in OrderedCounter(l_cluster_i).keys()]

        if l_title is None: l_title = list(self.t.keys())
        if l_title_cbar is None: l_title_cbar = list(l_title)
        if figsize is None: figsize = (2 * len(l_title), 20)

        kwargs_grid_ = {
            'nrows_ncols': (len(ll_cluster_i), len(l_title)),
            'label_mode': 'L',
            'aspect': False,
            'axes_pad': 0.15,
            'cbar_location': 'bottom',
            'cbar_mode': 'edge',
            'cbar_size': '2%',
            'cbar_pad': 0.2,
        }
        kwargs_grid_.update(kwargs_grid)
        fig = plt.figure(figsize=figsize)
        if suptitle is None:
            fig.suptitle('%d regions total' % (len(self.r)))
        else:
            fig.suptitle(suptitle)
        plt.subplots_adjust(top=0.95) # optimised to save space at the top
        grid = ImageGrid(fig, 111, **kwargs_grid_)

        def set_axis_off_(ax):
            #ax.axis('off')
            # Remove all of axis (ticks, frames, etc), except ylabels
            #https://stackoverflow.com/questions/2553521/setting-axes-linewidth-without-changing-the-rcparams-global-dict
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(0)
            #https://stackoverflow.com/questions/12998430/remove-xticks-in-a-matplotlib-plot
            ax.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')

        for i, ax in enumerate(grid):#title in enumerate(l_title):
            if axis_off: set_axis_off_(ax)
            #print col_i, row_i
            col_i = i % len(l_title)
            row_i = i // len(l_title)
            title_i = l_title[col_i]
            cluster_i = ll_cluster_i[row_i]

            # http://matplotlib.org/examples/axes_grid/demo_axes_grid.html
            im = self.t[title_i].get_subset(cluster_i).imshow(ax=ax, vextent=vextent, nsquashed=nsquashed)
            try:
                if title_i in l_title_cbar:
                    cmap_ = self.t[title_i].imshow_kwargs['cmap']
                    vmin_ = self.t[title_i].imshow_kwargs['vmin']
                    vmax_ = self.t[title_i].imshow_kwargs['vmax']

                    # Special cases that should be abstracted away into analysis notebooks...
                    if thesis_magic and title_i == 'Promoter\nstate':
                        grid.cbar_axes[i].colorbar(im, ticks=[0, 1])
                        grid.cbar_axes[i].set_xticklabels(['other', 'promoter'], rotation=45, fontsize=fontsize)

                    elif thesis_magic and title_i == 'CV\nrank':
                        grid.cbar_axes[i].colorbar(im, ticks=[0, 1])
                        grid.cbar_axes[i].set_xticklabels(['stable', 'regulated'], rotation=45, fontsize=fontsize)
                    
                    elif figB_magic:
                        if title_i.endswith('\nwt_emb') and row_i == 1:
                            if title_i == 'H3K4me3\nwt_emb':
                                ax_ = fig.add_axes([0.141, 0.125, 0.205, 0.017])

                            elif title_i == 'H3K4me1\nwt_emb':
                                ax_ = fig.add_axes([0.387, 0.125, 0.205, 0.017])

                            elif title_i == 'H3K27me3\nwt_emb':
                                ax_ = fig.add_axes([0.635, 0.125, 0.205, 0.017])

                            for axis in ['top','bottom','left','right']:
                                ax_.spines[axis].set_linewidth(0)

                            cb = matplotlib.colorbar.ColorbarBase(ax_, 
                                cmap=cmap_, 
                                norm=matplotlib.colors.Normalize(vmin=vmin_, vmax=vmax_),
                                ticks=[vmin_, 0, vmax_],
                                orientation='horizontal',
                                drawedges=False
                            )

                            ax_.set_xticklabels([vmin_, '0', '+' + str(vmax_)], rotation=45, fontsize=fontsize)

                        elif title_i == 'CV\nrank':
                            grid.cbar_axes[i].colorbar(im, ticks=[vmin_, vmax_])#, drawedges=False)
                            #grid.cbar_axes[i].set_xticklabels([vmin_, '0', '+' + str(vmax_)], rotation=45, fontsize=fontsize)
                            grid.cbar_axes[i].set_xticklabels(['stable', 'regulated'], rotation=45, fontsize=5)

                        else:
                            grid.cbar_axes[i].axis('off') # Needed to completely kill ImageGrid colour bar

                    elif figB_suppl_magic:
                        if title_i.endswith('\nwt_emb') and row_i == 1:
                            width_ = 0.16
                            if title_i == 'H3K4me3\nwt_emb':
                                ax_ = fig.add_axes([0.14, 0.125, width_, 0.017])

                            elif title_i == 'H3K4me1\nwt_emb':
                                ax_ = fig.add_axes([0.325, 0.125, width_, 0.017])

                            elif title_i == 'H3K36me3_gene\nwt_emb':
                                ax_ = fig.add_axes([0.51, 0.125, width_, 0.017])

                            elif title_i == 'H3K27me3\nwt_emb':
                                ax_ = fig.add_axes([0.697, 0.125, width_, 0.017])

                            for axis in ['top','bottom','left','right']:
                                ax_.spines[axis].set_linewidth(0)

                            cb = matplotlib.colorbar.ColorbarBase(ax_, 
                                cmap=cmap_, 
                                norm=matplotlib.colors.Normalize(vmin=vmin_, vmax=vmax_),
                                ticks=[vmin_, 0, vmax_],
                                orientation='horizontal',
                                drawedges=False
                            )

                            ax_.set_xticklabels([vmin_, '0', '+' + str(vmax_)], rotation=45, fontsize=fontsize)

                        elif title_i == 'CV\nrank':
                            grid.cbar_axes[i].colorbar(im, ticks=[vmin_, vmax_])#, drawedges=False)
                            #grid.cbar_axes[i].set_xticklabels([vmin_, '0', '+' + str(vmax_)], rotation=45, fontsize=fontsize)
                            grid.cbar_axes[i].set_xticklabels(['stable', 'regulated'], rotation=45, fontsize=5)

                        else:
                            grid.cbar_axes[i].axis('off') # Needed to completely kill ImageGrid colour bar

                    else:
                        grid.cbar_axes[i].colorbar(im, ticks=[vmin_, 0, vmax_], drawedges=False)
                        grid.cbar_axes[i].set_xticklabels([vmin_, '0', '+' + str(vmax_)], rotation=45, fontsize=fontsize)

                    if axis_off: 
                        for axis in ['top','bottom','left','right']:
                            grid.cbar_axes[i].spines[axis].set_linewidth(0)

                else:
                    grid.cbar_axes[i].axis('off')


            except AttributeError as e: # To ignore known bug
                print('known AttributeError: %(e)s' % locals())

            if figB_magic:
                fig.text(0.20, 1.015, 'H3K4me3')
                fig.text(0.443, 1.015, 'H3K4me1')
                fig.text(0.69, 1.015, 'H3K27me3')

            if figB_suppl_magic:
                fig.text(0.185, 1.015, 'H3K4me3')
                fig.text(0.373, 1.015, 'H3K4me1')
                fig.text(0.554, 1.015, 'H3K36me3')
                fig.text(0.74, 1.015, 'H3K27me3')

            # Reduce font size for long titles (e.g. jadb files names)
            if row_i == 0:
                if figB_magic:
                    d_ = {
                        'H3K4me3\nwt_emb': 'Emb',
                        'H3K4me3\nwt_l1': 'L1',
                        'H3K4me3\nwt_l2': 'L2',
                        'H3K4me3\nwt_l3': 'L3',
                        'H3K4me3\nwt_l4': 'L4',
                        'H3K4me3\nwt_ya': 'YA',
                        'H3K4me1\nwt_emb': 'Emb',
                        'H3K4me1\nwt_l1': 'L1',
                        'H3K4me1\nwt_l2': 'L2',
                        'H3K4me1\nwt_l3': 'L3',
                        'H3K4me1\nwt_l4': 'L4',
                        'H3K4me1\nwt_ya': 'YA',
                        'H3K27me3\nwt_emb': 'Emb',
                        'H3K27me3\nwt_l1': 'L1',
                        'H3K27me3\nwt_l2': 'L2',
                        'H3K27me3\nwt_l3': 'L3',
                        'H3K27me3\nwt_l4': 'L4',
                        'H3K27me3\nwt_ya': 'YA',
                        'CV\nrank': 'CV\nrank',
                    }
                    ax.set_title(d_[title_i], fontsize=fontsize)
                elif figB_suppl_magic:
                    d_ = {
                        'H3K4me3\nwt_emb': 'Emb',
                        'H3K4me3\nwt_l1': 'L1',
                        'H3K4me3\nwt_l2': 'L2',
                        'H3K4me3\nwt_l3': 'L3',
                        'H3K4me3\nwt_l4': 'L4',
                        'H3K4me3\nwt_ya': 'YA',
                        'H3K4me1\nwt_emb': 'Emb',
                        'H3K4me1\nwt_l1': 'L1',
                        'H3K4me1\nwt_l2': 'L2',
                        'H3K4me1\nwt_l3': 'L3',
                        'H3K4me1\nwt_l4': 'L4',
                        'H3K4me1\nwt_ya': 'YA',
                        'H3K36me3_gene\nwt_emb': 'Emb',
                        'H3K36me3_gene\nwt_l1': 'L1',
                        'H3K36me3_gene\nwt_l2': 'L2',
                        'H3K36me3_gene\nwt_l3': 'L3',
                        'H3K36me3_gene\nwt_l4': 'L4',
                        'H3K36me3_gene\nwt_ya': 'YA',
                        'H3K27me3\nwt_emb': 'Emb',
                        'H3K27me3\nwt_l1': 'L1',
                        'H3K27me3\nwt_l2': 'L2',
                        'H3K27me3\nwt_l3': 'L3',
                        'H3K27me3\nwt_l4': 'L4',
                        'H3K27me3\nwt_ya': 'YA',
                        'CV\nrank': 'CV\nrank',
                    }
                    ax.set_title(d_[title_i], fontsize=fontsize)
                else:
                    ax.set_title(title_i, fontsize=fontsize)            
                #if len(title_i) <= 31:
                #    ax.set_title(title_i, fontsize=fontsize)
                #else:
                #    ax.set_title(title_i, fontsize=3)

            if col_i == 0:
                if not (cluster_col is None) and show_cluster_col:
                    cluster_val = list(OrderedCounter(l_cluster_i).keys())[row_i]
                    ax.set_ylabel('%(cluster_col)s=%(cluster_val)s' % locals(), fontsize=fontsize)
                if not(row_labels is None):
                    ax.set_ylabel(row_labels[row_i], fontsize=fontsize)

            if 'skew' in title_i: ax.axvline(0, color='k', linewidth=0.8, alpha=0.5)

        if savefig:
            fp = '_fig/yp%s.pdf' % (datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),)
            print('yp.savefig: %s' % (fp,))
            plt.savefig(fp, bbox_inches='tight')
        
        return grid

    def get_m_stacked(self, l_title):
        """Returns a set of signal matrices 'stacked together' for kmeans/PCA/etc input"""
        return np.column_stack([self.t[title].m for title in l_title])

    def plot(self, l_title=None, cluster_col=None, figsize=None, ylim=None, *args, **kwargs):
        if cluster_col is None:
            ll_cluster_i = [self.r.index]
        else:
            l_cluster_i = self.r[cluster_col]
            ll_cluster_i = [l_cluster_i[l_cluster_i == v].index.tolist() for v in OrderedCounter(l_cluster_i).keys()]

        if l_title is None: l_title = list(self.t.keys())
        kwargs_grid = {
            'nrows_ncols': (len(ll_cluster_i), len(l_title)),
            'aspect': False,
            'axes_pad': 0.15,
        }
        fig = plt.figure(figsize=figsize)
        fig.suptitle('%d regions total' % (len(self.r)))
        grid = ImageGrid(fig, 111, **kwargs_grid)
        for i, ax in enumerate(grid):#title in enumerate(l_title):
            #print col_i, row_i
            col_i = i % len(l_title)
            row_i = i // len(l_title)
            title_i = l_title[col_i]
            cluster_i = ll_cluster_i[row_i]
            self.t[title_i].get_subset(cluster_i).plot(ax=ax, *args, **kwargs)
            if not (ylim is None): ax.set_ylim(*ylim)
            if row_i == 0:
                ax.set_title(l_title[col_i])

            if col_i == 0 and not (cluster_col is None):
                cluster_val = list(OrderedCounter(l_cluster_i).keys())[row_i]
                ax.set_ylabel('%(cluster_col)s=%(cluster_val)s' % locals())
        return grid

    def save_clust(self, cluster_col, decomp_col, figsize=(45,15)):
        pf = tp.PipeFile(tp.ts(), '_yarp_save_clust', None, os.getcwd())
        self.sort([cluster_col, decomp_col]).imshow(cluster_col=cluster_col, figsize=figsize)
        fp_imshow = str(pf.spawn_suffix('.imshow.pdf'))
        plt.savefig(fp_imshow)

        fp_labels = str(pf.spawn_suffix('.labels.bed'))
        self.r[['chrom', 'start', 'end', cluster_col]].to_csv(fp_labels, sep='\t', index=False, header=False)

        fp_sample = str(pf.spawn_suffix('.sample_10.bed'))
        l_clustlabel = OrderedCounter(self.r[cluster_col]).keys()
        with open(fp_sample, 'a') as fh_sample:
            for clustlabel in l_clustlabel:
                df_sub = df_sample(self.r[self.r[cluster_col] == clustlabel], k=10)
                df_sub[['chrom', 'start', 'end', cluster_col]].to_csv(fh_sample, sep='\t', index=False, header=False)

    # Deprecated stuff starts here?
    #def get_m_subset(self, l_title=None):
    #   if l_title is None:
    #        l_title = self.m.keys()
    #    return np.column_stack([self.m[title] for title in l_title])

    @property
    def n_plots(self):
        return len(self.m.keys())

    @property
    def n_regions(self):
        return len(self.l_iv)

    @property
    def l_xbin_mid(self):
        return range(-self.flank_len, self.flank_len + 1, self.bin_size)

    def get_m_title_sorted_if_possible(self, title):
        if self.m_lexsort_keys is None:
            return self.m[title]
        else:
            return self.m[title][np.lexsort(self.m_lexsort_keys),]

    def sort_key_add(self, keys):
        if (self.m_lexsort_keys is None):
            self.m_lexsort_keys = [keys,]
        else:
            self.m_lexsort_keys.insert(0, keys)

    def sort_keys_clear(self):
        self.m_lexsort_keys = None

    #def cluster_sizes(self):
    #    pprint.pprint(collections.Counter(self.clusterer.labels_))

    def m_profiles_sorted(self):
        return self.m_profiles[np.lexsort(self.m_profiles_lexsort),]

    @property
    def l_title(self): return self.m.keys()

    @property
    def cluster_sizes(self):
        #return OrderedCounter(self.clusterer.labels_[np.lexsort(self.m_profiles_lexsort)]).values()
        return OrderedCounter(np.sort(self.m_lexsort_keys[-1])).values() # There's probably a slightly faster way to do this...
        #return np.sort(collections.Counter(gh.m_lexsort_keys[-1]).values()) # And this isn't it...

    @property
    def n_clusters(self):
        return len(set(self.m_lexsort_keys[-1]))

    def l_counts_fractions(self, l_counts):
        sum_counts = sum(l_counts)
        return map(lambda count: "%d\n(%.1f%%)" % (count, 100.0*count/sum_counts), l_counts)

    def cluster_centre(self, i):
        return sum(self.cluster_sizes[:i]) + self.cluster_sizes[i] // 2

    def zero_line(self, title):
        i = self.m.keys().index(title)
        z = (self.m[title].shape[1] - 1) / 2.0
        self.grid[i].axvline(z, color='k')

    def plot_cqa(self, cluster_i, title, *args, **kwargs):
        """
        # plot cluster quantile aggregates
        # TODO set x axis labels properly
        # TODO: use cluster_i instead! Alternatively, is this routine actually necessary - combine from_cluster_i and
        plot_qa instead?
        """
        cluster_start = int(sum(self.cluster_sizes[:cluster_i]))
        cluster_end   = int(sum(self.cluster_sizes[:cluster_i + 1]))
        self.plot_qa(self.get_m_title_sorted_if_possible(title)[cluster_start:cluster_end], *args, **kwargs)

    def cluster_i(self, i):
        """
        get indices of the cluster, as they appear in the current sorted order
        => easy helper routine to generate a new GH only consisting of the current cluster!
        """
        l_i = np.array(range(self.n_regions))#ii = range(np.array(self.m_lexsort_keys).shape[1]
        l_i_sort = l_i[np.lexsort(self.m_lexsort_keys),]
        c_i = np.array(self.m_lexsort_keys[-1])
        c_i_sort = c_i[np.lexsort(self.m_lexsort_keys),]
        return l_i_sort[c_i_sort == i]

    @classmethod
    def from_gr_cluster(cls, gr, c_i, v=False):
        """
        from_gr_cluster(GenomicRegions, gr, 1).imshow()
        """
        return cls.from_gr_subset(gr, gr.cluster_i(c_i), v)

    def plot_qa_sorted(self, q=5, col=1, col_max=2, l_title = None):
        if l_title is None:
            l_title = self.m.keys()
        n_plots = len(l_title)
        for i, title in enumerate(l_title):
            ax = GenomicRegions.subplot_row_col(n_plots, col_max, col, i+1)
            GenomicRegions.plot_qa(self.get_m_title_sorted_if_possible(title), q)
            ax.set_title(title)
            axvline(50, color='k')

    def plot_mean(self, title, *args, **kwargs):
        data = self.get_m_title_sorted_if_possible(title)
        data_mean = np.mean(data, 0)
        plot(self.l_xbin_mid, data_mean, drawstyle='steps-mid', *args, **kwargs)

    def plot_quantiles(self, title, q=5, cmap=matplotlib.cm.get_cmap("Reds"), *args, **kwargs):
        data = self.get_m_title_sorted_if_possible(title)
        for (i, data_i) in zip(range(q), np.array_split(data, q)):
            #http://stackoverflow.com/questions/15140072/how-to-map-number-to-color-using-matplotlibs-colormap
            color_i = cmap(quantile_mid(i, q))
            data_i_mean = np.mean(data_i, 0)
            plot(self.l_xbin_mid, data_i_mean, label="%d of %d" % (i + 1, q), color=color_i, drawstyle='steps-mid', *args, **kwargs)
        #legend()

    def plot_cluster_profiles(self, l_title=None, common_ylim=None, q=5, cmap=matplotlib.cm.get_cmap("Reds")):
        """
        Plot average profiles of clusters
        """
        if l_title is None: l_title = self.m.keys()
        n_tracks = len(l_title)
        n_clusters = self.n_clusters
        fig = plt.figure()
        grid = ImageGrid(fig, 111, nrows_ncols=(n_tracks, n_clusters), aspect=False, share_all=True, axes_pad=0)
        # Iterate over plots; plot data
        for (i, grid_i) in enumerate(grid):
            (row_i, col_i) = (i // n_clusters, i % n_clusters)
            gr_col = GenomicRegions.from_gr_cluster(self, col_i)
            data = gr_col.get_m_title_sorted_if_possible(l_title[row_i])
            for (i, data_i) in zip(range(q), np.array_split(data, q)):
                color_i = cmap(quantile_mid(i, q)) # strip gpl.
                data_i_mean = np.mean(data_i, 0)
                grid_i.plot(self.l_xbin_mid, data_i_mean, label="%d of %d" % (i + 1, q), color=color_i, drawstyle='steps-mid')
            #data_mean = np.mean(data, 0)
            #grid_i.plot(self.l_xbin_mid, data_mean, drawstyle='steps-mid')
            grid_i.set_xlabel('Genomic coordinate (bp)')
            grid_i.set_ylabel(l_title[row_i])
        # Axis limits
        flank_len_ext = self.flank_len# + self.bin_size / 2
        grid_i.set_xlim(-flank_len_ext, +flank_len_ext)
        if not (common_ylim is None): grid_i.set_ylim(common_ylim[0],common_ylim[1])
        # Cluster titles
        for (clust_i, clust_size) in enumerate(self.cluster_sizes):
            grid[clust_i].set_title("%d regions" % (clust_size,))

# source: https://github.com/konstantint/matplotlib-venn/blob/master/matplotlib_venn/_venn2.py
# modified to have labels above the circles
def venn2_fudge_labels_to_top(subsets, set_labels=('A', 'B'), set_colors=('r', 'g'), alpha=0.4, normalize_to=1.0, ax=None, subset_label_formatter=None):
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['10', '01', '11']]
    elif len(subsets) == 2:
        subsets = matplotlib_venn.compute_venn2_subsets(*subsets)

    if subset_label_formatter is None:
        subset_label_formatter = str

    areas = matplotlib_venn._venn2.compute_venn2_areas(subsets, normalize_to)
    centers, radii = matplotlib_venn._venn2.solve_venn2_circles(areas)
    regions = matplotlib_venn._venn2.compute_venn2_regions(centers, radii)
    colors = matplotlib_venn._venn2.compute_venn2_colors(set_colors)

    if ax is None:
        ax = plt.gca()
    matplotlib_venn._venn2.prepare_venn_axes(ax, centers, radii)

    # Create and add patches and subset labels
    patches = [r.make_patch() for r in regions]
    for (p, c) in zip(patches, colors):
        if p is not None:
            p.set_facecolor(c)
            p.set_edgecolor('none')
            p.set_alpha(alpha)
            ax.add_patch(p)
    label_positions = [r.label_position() for r in regions]
    subset_labels = [ax.text(lbl[0], lbl[1], subset_label_formatter(s), va='center', ha='center') if lbl is not None else None for (lbl, s) in zip(label_positions, subsets)]

    # Position set labels
    if set_labels is not None:
        padding = np.mean([r * 0.1 for r in radii])
        label_positions = [centers[0] - np.array([0.0, - radii[0] - padding]),
                           centers[1] - np.array([0.0, - radii[1] - padding])]
        labels = [ax.text(pos[0], pos[1], txt, size='large', ha='right', va='bottom') for (pos, txt) in zip(label_positions, set_labels)]
        labels[1].set_ha('left')
    else:
        labels = None
    return matplotlib_venn._venn2.VennDiagram(patches, subset_labels, labels, centers, radii)

class GenomicVenn2:
    def __init__(self, bt_a, bt_b, label_a=None, label_b=None, v=False):
        self.bt_a = bt_a
        self.bt_b = bt_b
        self.df_a = self.bt_a.to_dataframe()
        self.df_b = self.bt_b.to_dataframe()
        self.n_a = self.bt_a.count()
        self.n_b = self.bt_b.count()
        assert(self.n_a > 0)
        assert(self.n_b > 0)
        self.label_a = label_a if not(label_a is None) else 'label_a'
        self.label_b = label_b if not(label_b is None) else 'label_b'

        try:
            self.df_a_only = self.bt_a.subtract(self.bt_b, A=True).to_dataframe()
        except pd.errors.EmptyDataError: # Biological feature sets can, in fact, have complete overlap....
            self.df_a_only = pd.DataFrame(columns=self.df_a.columns)

        try:
            self.df_b_only = self.bt_b.subtract(self.bt_a, A=True).to_dataframe()
        except pd.errors.EmptyDataError:
            self.df_b_only = pd.DataFrame(columns=self.df_b.columns)

        if v: print('%d len(df_a_only)' % (len(self.df_a_only),))
        if v: print('%d len(df_b_only)' % (len(self.df_b_only),))
        self.df_a_with_b = pd.merge(self.df_a, self.df_a_only, how='outer', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
        self.df_b_with_a = pd.merge(self.df_b, self.df_b_only, how='outer', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
        assert len(self.df_a) == len(self.df_a_only) + len(self.df_a_with_b)
        assert len(self.df_b) == len(self.df_b_only) + len(self.df_b_with_a)
        self.n_a_only = len(self.df_a_only)#self.bt_a.subtract(self.bt_b, A=True).count()
        self.n_b_only = len(self.df_b_only)#self.bt_b.subtract(self.bt_a, A=True).count()

    @property
    def n_a_with_b(self): return self.n_a - self.n_a_only

    @property
    def n_b_with_a(self): return self.n_b - self.n_b_only

    @property
    def f_a_with_b_of_a(self): return self.n_a_with_b / float(self.n_a)

    @property
    def f_b_with_a_of_b(self): return self.n_b_with_a / float(self.n_b)

    def __repr__(self):
        return "%d(%.1f%%) / %d(%.1f%%)" % (self.n_a_with_b, 100*self.f_a_with_b_of_a, self.n_b_with_a, 100*self.f_b_with_a_of_b)

    def plot(self, fudge_labels_to_top=False, style='verbose', *args, **kwargs):
        '''
        Draws a two-way Venn diagram with overlapping vs non-overlapping areas normalised overlapping vs non-overlapping peaks.
        '''
        figsize = kwargs.pop('figsize', None)
        if not(figsize is None):
            plt.figure(figsize=figsize)

        if style == 'verbose':
            set_labels_ = (
                "%s\n(%s regions)" % (self.label_a, f_uk(self.n_a)),
                "%s\n(%s regions)" % (self.label_b, f_uk(self.n_b))
            )
            label_a_only = "%s\n(%.1f%%)" % (f_uk(self.n_a_only), 100*self.n_a_only/float(self.n_a))
            label_b_only = "%s\n(%.1f%%)" % (f_uk(self.n_b_only), 100*self.n_b_only/float(self.n_b))
            label_a_with_b = "%s(%.1f%%)\n/\n%s(%.1f%%)" % (f_uk(self.n_a_with_b), 100*self.f_a_with_b_of_a, f_uk(self.n_b_with_a), 100*self.f_b_with_a_of_b)

        elif style == 'compact':
            set_labels_ = (
                "%s\n(n=%s)" % (self.label_a, f_uk(self.n_a)),
                "%s\n(n=%s)" % (self.label_b, f_uk(self.n_b))
            )
            label_a_only = "%.1f%%\n " % (100*self.n_a_only/float(self.n_a),)
            label_b_only = "\n%.1f%%" % (100*self.n_b_only/float(self.n_b),)
            label_a_with_b = "%.1f%%\n%.1f%%" % (100*self.f_a_with_b_of_a, 100*self.f_b_with_a_of_b)

        area_a = self.n_a_only / float(self.n_a_with_b)
        area_b = self.n_b_only / float(self.n_b_with_a)
        if not fudge_labels_to_top:
            v = matplotlib_venn.venn2(subsets = (area_a, area_b, 1), set_labels = set_labels_, *args, **kwargs)
        else:
            v = venn2_fudge_labels_to_top(subsets = (area_a, area_b, 1), set_labels = set_labels_, *args, **kwargs)

        v.get_label_by_id('10').set_text(label_a_only)
        v.get_label_by_id('01').set_text(label_b_only)
        v.get_label_by_id('11').set_text(label_a_with_b)
        return v

def GenomicVenn3(bt_a, bt_b, bt_c, **kwargs):
    """Three-way Venn diagram to reflect overlap of peak sets. Overlaps are calculated in by counting base pairs in peaks that overlap. (not peaks nor just base pairs)"""
    figsize = kwargs.pop('figsize', None)
    if not(figsize is None): plt.figure(figsize=figsize)
    bt_a_ = pybedtools.BedTool([(r[0], r[1], r[2], 'a') for r in bt_a])
    bt_b_ = pybedtools.BedTool([(r[0], r[1], r[2], 'b') for r in bt_b])
    bt_c_ = pybedtools.BedTool([(r[0], r[1], r[2], 'c') for r in bt_c])
    bt = pybedtools.BedTool(itertools.chain([bt_a_, bt_b_, bt_c_])).sort().merge(c='4', o='collapse')
    s_a = set(['\t'.join([chrom, start, end, name]) for (chrom, start, end, name) in bt if 'a' in name])
    s_b = set(['\t'.join([chrom, start, end, name]) for (chrom, start, end, name) in bt if 'b' in name])
    s_c = set(['\t'.join([chrom, start, end, name]) for (chrom, start, end, name) in bt if 'c' in name])
    v = matplotlib_venn.venn3(subsets = (s_a, s_b, s_c), **kwargs)
    return v
