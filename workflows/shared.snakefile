import collections
import glob
import os
import os.path
import math
import csv
import numpy as np
import scipy as sp
import pandas as pd
import itertools
#import matplotlib
#import matplotlib.pyplot as plt
from pprint import pprint
import yaml
import io
import string
from pprint import pprint

import numpy as np
import scipy as sp
import scipy.signal
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pyBigWig
import weblogolib

import HTSeq as hts

import twobitreader

import pybedtools
from pybedtools import BedTool

chroms_ce10 = {
    'chrI': 15072423,
    'chrII': 15279345,
    'chrIII': 13783700,
    'chrIV': 17493793,
    'chrM': 13794,
    'chrV': 20924149,
    'chrX': 17718866
}

def read_int(fp_inp):
    assert(os.path.isfile(fp_inp))
    #print(fp_inp)
    with open(fp_inp) as fh_inp:
        return int(fh_inp.read())

def pf(id, step, suffix, prefix=''): return os.path.join(prefix, step, '%(id)s.%(step)s%(suffix)s' % locals())

def parse_pf(fp):
    (fd, fn) = os.path.split(fp)
    (prefix, step) = os.path.split(fd) # step: 'deepest' directory; prefix: everything above step directory
    sample = fn[:fn.find(step)].rstrip('.') # sample: find step in file name; subtract '.'
    suffix = fn[fn.find(step) + len(step):] # suffix: everything after step; DO NOT subtract '.'
    return (sample, step, suffix, prefix)

def htp(fp): return os.path.expanduser(os.path.join('~/lab/HTSProcessing', fp))
def lp(fp): return os.path.join(os.path.expanduser('~/lab/HTSProcessing'), fp) # Path relative to old cluster pipeline

class PipeFile:
    def __init__(self, sample, step, suffix, prefix):
        (self.sample, self.step, self.suffix, self.prefix) = sample, step, suffix, prefix

    @classmethod
    def from_fp(cls, fp):
        (fd, fn) = os.path.split(fp)
        (prefix, step) = os.path.split(fd) # step: 'deepest' directory; prefix: everything above step directory
        sample = fn[:fn.find(step)].rstrip('.') # sample: find step in file name; subtract '.'
        suffix = fn[fn.find(step) + len(step):] # suffix: everything after step; DO NOT subtract '.'
        return cls(sample, step, suffix, prefix)

def ls_bid(step, suffix, prefix):
    fd_step = os.path.join(prefix, step)
    l_fp = [os.path.join(step, fn) for fn in os.listdir(fd_step)]
    l_bid = []
    for fp in l_fp:
        l_bid.append(PipeFile.from_fp(fp).sample)
    return l_bid

def read_regions(fp, chroms, starts, ends, f=None):
    fh = pyBigWig.open(fp)
    for chrom, start, end in zip(chroms, starts, ends):
        if f is None:
            yield fh.values(chrom, start, end)
        else:
            yield f(np.array(fh.values(chrom, start, end)))
    fh.close()

"""
rule sample_1M:
    input: htp('.bid/{bid}.r1.fq.gz'), htp('.bid/{bid}.r2.fq.gz')
    output: htp('.bid/{bid}_1M.r1.fq.gz'), htp('.bid/{bid}_1M.r2.fq.gz')
    shell: # trap error code 141 thrown due to broken pipe -- http://stackoverflow.com/questions/22464786/ignoring-bash-pipefail-for-error-code-141
        '''
        gunzip -c {input[0]} | head -n 4000000 | gzip -c - > {output[0]} || if [[ $? -eq 141 ]]; then true; else exit $?; fi
        gunzip -c {input[1]} | head -n 4000000 | gzip -c - > {output[1]} || if [[ $? -eq 141 ]]; then true; else exit $?; fi
        '''

rule sample_1M_test:
    input:
        expand(htp('.bid/{bid}_1M.r1.fq.gz'), bid=['chen13_emb_lcap1', 'chen13_emb_lcap2'])
"""

rule sample_1M_r1:
    input:
        'samples/{bid}.r1.fq.gz',
    output:
        'samples/{bid}_1M.r1.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 1000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_1M_r2:
    input:
        'samples/{bid}.r2.fq.gz',
    output:
        'samples/{bid}_1M.r2.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 1000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_5M_r1:
    input:
        'samples/{bid}.r1.fq.gz',
    output:
        'samples/{bid}_5M.r1.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 5000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_5M_r2:
    input:
        'samples/{bid}.r2.fq.gz',
    output:
        'samples/{bid}_5M.r2.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 5000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_10M_r1:
    input:
        'samples/{bid}.r1.fq.gz',
    output:
        'samples/{bid}_10M.r1.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 10000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_10M_r2:
    input:
        'samples/{bid}.r2.fq.gz',
    output:
        'samples/{bid}_10M.r2.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 10000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_20M_r1:
    input:
        'samples/{bid}.r1.fq.gz',
    output:
        'samples/{bid}_20M.r1.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 20000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule sample_20M_r2:
    input:
        'samples/{bid}.r2.fq.gz',
    output:
        'samples/{bid}_20M.r2.fq.gz',
    params:
        seed = 42, # Tried, but could not resist
        frac_or_number = 20000000,
    shell:
        'seqtk sample -s {params.seed} <(pigz -c -d -p {threads} {input}) {params.frac_or_number} | pigz -c -p {threads} - > {output}'

rule samtools_index:
    input:
        '{inp}.bam'
    output:
        '{inp}.bam.bai'
    shell:
        'samtools index {input}'

rule ce10:
    output:
        'shared/ce10.2bit',
        'shared/ce10.chroms',
        'shared/ce10.fa',
        'shared/ce10.fa.amb',
        'shared/ce10.fa.ann',
        'shared/ce10.fa.bwt',
        'shared/ce10.fa.pac',
        'shared/ce10.fa.sa',
        'shared/ce10.1.ebwt',
    shell: '''
        curl http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit -o shared/ce10.2bit
        curl http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.chrom.sizes -o shared/ce10.chroms
        twoBitToFa shared/ce10.2bit shared/ce10.fa
        bwa index shared/ce10.fa
        bowtie-build shared/ce10.fa shared/ce10
        '''

rule ce11:
    output:
        'shared/ce11.2bit',
        'shared/ce11.chroms',
        'shared/ce11.fa',
        #'shared/ce11.fa.amb',
        #'shared/ce11.fa.ann',
        #'shared/ce11.fa.bwt',
        #'shared/ce11.fa.pac',
        #'shared/ce11.fa.sa',
        #'shared/ce11.1.ebwt',
    shell: '''
        curl http://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit -o shared/ce11.2bit
        curl http://hgdownload.cse.ucsc.edu/goldenPath/ce11/bigZips/ce11.chrom.sizes -o shared/ce11.chroms
        twoBitToFa shared/ce11.2bit shared/ce11.fa
        bwa index shared/ce11.fa
        #bowtie-build shared/ce11.fa shared/ce11
        '''

"""
bwa index shared/ecoli.fa

curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -o shared/hg19.2bit
curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -o shared/hg19.chroms

cat shared/ce10.fa shared/ecoli.fa > shared/ce10_ecoli.fa
bwa index shared/ce10_ecoli.fa
"""

rule cb3:
    output:
        'shared/cb3.2bit',
        'shared/cb3.chroms',
        'shared/cb3.fa',
    shell: '''
        curl http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/cb3.2bit -o shared/cb3.2bit
        curl http://hgdownload.cse.ucsc.edu/goldenPath/cb3/bigZips/cb3.chrom.sizes -o shared/cb3.chroms
        twoBitToFa shared/cb3.2bit shared/cb3.fa
        bwa index shared/cb3.fa
        '''

rule mean10:
    input:
        pf('sites', '{step}', '.bw', '{prefix}'),
    output:
        pf('sites', '{step}', '_mean10.bw', '{prefix}'),
    shell: '''
        scripts/bigWigConvolve {input} mean10 {output}
    '''

rule mean25:
    input:
        pf('sites', '{step}', '.bw', '{prefix}'),
    output:
        pf('sites', '{step}', '_mean25.bw', '{prefix}'),
    shell: '''
        scripts/bigWigConvolve {input} mean25 {output}
    '''

rule mean50:
    input:
        pf('sites', '{step}', '.bw', '{prefix}'),
    output:
        pf('sites', '{step}', '_mean50.bw', '{prefix}'),
    shell: '''
        scripts/bigWigConvolve {input} mean50 {output}
    '''

rule mean100:
    input:
        pf('sites', '{step}', '.bw', '{prefix}'),
    output:
        pf('sites', '{step}', '_mean100.bw', '{prefix}'),
    shell: '''
        scripts/bigWigConvolve {input} mean100 {output}
    '''

rule mean150:
    input:
        pf('sites', '{step}', '.bw', '{prefix}'),
    output:
        pf('sites', '{step}', '_mean150.bw', '{prefix}'),
    shell: '''
        scripts/bigWigConvolve {input} mean150 {output}
    '''

rule D2_s3x150_1E6:
    input:
        pf('sites', '{step}', '.bw', '{prefix}'),
    output:
        pf('sites', '{step}', '_D2_s3x150_1E6.bw', '{prefix}'),
    shell: '''
        scripts/bigWigConvolve {input} D2_s3x150_1E6 {output}
    '''
