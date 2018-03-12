#!/usr/bin/env python
import argparse
import sys

import HTSeq as hts

def inv_dist(i, j):
    """
    Distance between too GenomicInterval-s.
        0 if adjacent or overlapping
        inf if on different chromosomes or strands
    """
    if i.overlaps(j): return(0)
    if (i.chrom == j.chrom) and (i.strand == j.strand):
        if (i.end <= j.start):
            return(j.start - i.end)
        if (j.end <= i.start):
            return(i.start - j.end)
        assert(False)
    else:
        return float('inf')

assert(inv_dist(hts.GenomicInterval("chr1", 1, 5, "+"), hts.GenomicInterval("chr1", 5, 6, "+")) == 0)
assert(inv_dist(hts.GenomicInterval("chr1", 1, 5, "+"), hts.GenomicInterval("chr1", 6, 7, "+")) == 1)
assert(inv_dist(hts.GenomicInterval("chr1", 1, 2, "+"), hts.GenomicInterval("chr1", 8, 9, "+")) == 6)
assert(inv_dist(hts.GenomicInterval("chr1", 5, 6, "+"), hts.GenomicInterval("chr1", 1, 5, "+")) == 0)
assert(inv_dist(hts.GenomicInterval("chr1", 6, 7, "+"), hts.GenomicInterval("chr1", 1, 5, "+")) == 1)
assert(inv_dist(hts.GenomicInterval("chr1", 8, 9, "+"), hts.GenomicInterval("chr1", 1, 2, "+")) == 6)
assert(inv_dist(hts.GenomicInterval("chr1", 1, 2, "+"), hts.GenomicInterval("chr2", 8, 9, "+")) == float('inf'))
assert(inv_dist(hts.GenomicInterval("chr1", 8, 9, "+"), hts.GenomicInterval("chr1", 1, 2, "-")) == float('inf'))
assert(inv_dist(hts.GenomicInterval("chr1", 1, 3, "+"), hts.GenomicInterval("chr1", 2, 4, "+")) == 0)
assert(inv_dist(hts.GenomicInterval("chr1", 1, 3, "+"), hts.GenomicInterval("chr1", 2, 4, "-")) == float('inf'))

def tics_raw_cluster(l_tic, stack_stack_thresh = 50, stack_singleton_thresh=25):
    """
    Generate short cap clusters by merging consecutive (GenomicInterval, val) pairs in l_tic. This *modifies* l_tic!

    This is why we cannot do "seeding" and "extending" in a single pass:

    seq: --------------------------------
    cap:   1 <25nt> 1 <25nt> 8 <26nt> 1
    TIC:   |-----------------|

    A single-pass clustering algorithm could not link the leftmost singleton to the entire cluster.
    """
    i = 1
    while (i < len(l_tic)):
        d = inv_dist(l_tic[i-1][0], l_tic[i][0])
        # Neighbouring stacks check
        f_stack_stack = (d <= stack_stack_thresh) and (min(l_tic[i-1][1], l_tic[i][1]) > 1)
        # Neighbouring stack + singleton cap check
        f_stack_singleton = (d <= stack_singleton_thresh) and (min(l_tic[i-1][1], l_tic[i][1]) == 1) and (max(l_tic[i-1][1], l_tic[i][1]) > 1)
        assert(not(f_stack_stack and f_stack_singleton))
        # Should the pair under consideration be merged?
        if f_stack_stack or f_stack_singleton:
            # Merge l_tic[i-1] and l_tic[i]
            l_tic[i-1][0].start = min(l_tic[i-1][0].start, l_tic[i][0].start)
            l_tic[i-1][0].end = max(l_tic[i-1][0].end, l_tic[i][0].end)
            l_tic[i-1] = (l_tic[i-1][0], l_tic[i-1][1] + l_tic[i][1])
            del l_tic[i]
            i -= 1
        else:
            i += 1

# Straightforward unit tests for cap_cluster:
l_tic_test = [
    (hts.GenomicInterval("chr1", 1, 2, "+"), 1),
    (hts.GenomicInterval("chr1", 27, 27, "+"), 1),
    (hts.GenomicInterval("chr1", 52, 52, "+"), 1),
    (hts.GenomicInterval("chr1", 74, 75, "+"), 8),
    (hts.GenomicInterval("chr1", 156, 157, "+"), 1),
    (hts.GenomicInterval("chr1", 1, 2, "-"), 1),
    (hts.GenomicInterval("chr1", 27, 27, "-"), 1),
    (hts.GenomicInterval("chr1", 52, 52, "-"), 1),
    (hts.GenomicInterval("chr1", 74, 75, "-"), 8),
    (hts.GenomicInterval("chr1", 156, 157, "-"), 1),
]
tics_raw_cluster(l_tic_test)
assert(len(l_tic_test) == 4)
assert((l_tic_test[0][0].start, l_tic_test[0][0].end, l_tic_test[0][1]) == (1, 75, 11))
assert((l_tic_test[1][0].start, l_tic_test[1][0].end, l_tic_test[1][1]) == (156, 157, 1))
assert((l_tic_test[2][0].start, l_tic_test[2][0].end, l_tic_test[2][1]) == (1, 75, 11))
assert((l_tic_test[3][0].start, l_tic_test[3][0].end, l_tic_test[3][1]) == (156, 157, 1))

parser = argparse.ArgumentParser(description='A re-implementation of the transcription initiation clustering strategy from (Chen et al 2013 GR).')
parser.add_argument('--stack_stack_threshold', default=50)
parser.add_argument('--stack_singleton_threshold', default=25)
args = parser.parse_args()

l_tic = []
for line in sys.stdin:
    tokens = line.rstrip().split("\t")
    try:
        chrom = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        val = float(tokens[3])
        l_tic.append((hts.GenomicInterval(chrom, start, end), val * (end - start)))
    except:
        print >>sys.stderr, 'Skipping line: %(line)s' % locals(),

print >>sys.stderr, 'Parsed %d lines from input' % (len(l_tic),)
tics_raw_cluster(l_tic, args.stack_stack_threshold, args.stack_singleton_threshold)
print >>sys.stderr, 'Called %d TICs' % (len(l_tic),)

for (tic, val) in l_tic:
    print '\t'.join(map(str, [tic.chrom, tic.start, tic.end, int(val)]))
