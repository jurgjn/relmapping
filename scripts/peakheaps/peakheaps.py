#!/usr/bin/env python

import itertools
import sys
import numpy as np
import HTSeq as hts

class MergePeaksAcrossDigestions:
	def __init__(self, l_fp):#, min_peak_score=None, min_core_score=None, iv=None):
		# Read peaks
		l_peak_fp = []
		for (i, fp) in enumerate(l_fp):
			fh = hts.BED_Reader(fp)
			l_peak = [gf for gf in fh]
			print >>sys.stderr, "%d peaks found in %s" % (len(l_peak), fp)
			#for peak in l_peak: peak.name = {"source": i}
			l_peak_fp.append(l_peak)

		self.l_peak_fp = l_peak_fp

	def make_adjacency_dict(self, min_adj_score=25):
		"""
		- For every pair "adjacent" tracks:
			- Make GenomicArrayOfSets with thickIv as intervals for A
			- Iterate over features in B
				- If overlap, add to v_e (symmetrically!)
		result gets placed in self.v_e
		"""
		def gf_merge(gf_i, gf_j):
			# Centre of A is contained in core of B, and centre of B is contained in core of A
			gp_i_centre = hts.GenomicPosition(gf_i.iv.chrom, (gf_i.iv.start + gf_i.iv.end) / 2, gf_i.iv.strand)
			gp_j_centre = hts.GenomicPosition(gf_j.iv.chrom, (gf_j.iv.start + gf_j.iv.end) / 2, gf_j.iv.strand)
			return(gp_i_centre.is_contained_in(gf_j.iv) and gp_j_centre.is_contained_in(gf_i.iv))

		v_e = { peak: set([]) for l_peak in self.l_peak_fp for peak in l_peak }
		def adjacency_dict_ij(l_gf_i, l_gf_j):
			# works on v_e in-place!
			ga_j = hts.GenomicArray(chroms="auto", stranded=False, typecode='O')
			for gf_jj in l_gf_j:
				ga_j[gf_jj.iv] = gf_jj
			for gf_ii in l_gf_i:
				for (iv, gf_jj) in ga_j[gf_ii.iv].steps():
					if not (gf_jj is None):
						if gf_merge(gf_ii, gf_jj):
							v_e[gf_ii].add(gf_jj)
							v_e[gf_jj].add(gf_ii)

		def pairwise_with_skip(l, max_diff):
			for (i, j) in itertools.combinations(range(len(l)), 2):
				if (j - i) <= max_diff:
					yield((l[i], l[j]))
				else:
					continue

		# iterate over neighbouring all combinations/pairwise digestions/pairwise-with-skip
		# ...and collect similar peaks into the v_e neighbourhood dictionary
		#for (l_gf_i, l_gf_j) in pairwise(self.l_peak_fp):
		#for (l_gf_i, l_gf_j) in pairwise_with_skip(self.l_peak_fp, max_diff=4):
		for (l_gf_i, l_gf_j) in itertools.combinations(self.l_peak_fp, 2):
			adjacency_dict_ij(l_gf_i, l_gf_j) # this updates v_e in place
		self.v_e = v_e
			
	def find_connected_components(self):
		"""
		Breadth-first walk of the adjacency graph to find connected components
		"""
		v_e = self.v_e
		v_u = set(v_e.keys()) # unvisited nodes
		v_c = [] # a list of connected components
		while v_u:
			v_i = set([v_u.pop()]) # all nodes belonging to the current conncected component
			v_q = set(list(v_i)) # unvisited edges of the current connected component
			while v_q:
				try: # choose arbitrary unvisited neighbour node, lookup neighbours
					v_n = v_e[v_q.pop()]
				except KeyError: # raised when dealing with a 1-length connected component => peak is not in v_e
					v_n = set([]) 
				v_n.difference_update(v_i) # remove neighbours that have already been visited
				v_u.difference_update(v_n) # remove neighbours from list of unvisited nodes
				v_i.update(v_n) # add neighbours to the connected component
				v_q.update(v_n) # add neighbours to the queue
			v_c.append(v_i) # add finalised connected component to results

		#assert(len(l_peak) == sum([len(c) for c in v_c]))
		print >>sys.stderr, "%d connected components " % (len(v_c),)
		#pprint.pprint(v_c[random.choose(v_c.keys(), 3)])
		self.v_c = v_c

	def write_merged_peaks(self):
		"""
		Output connected components as merged peaks
		"""
		def mean_iv(l_iv):
			return hts.GenomicInterval(chrom = l_iv[0].chrom,
									   start = np.mean([iv.start for iv in l_iv]),
									   end = np.mean([iv.end for iv in l_iv]),
									   strand = l_iv[0].strand
									   )
		l_gf_res = []

		for c in self.v_c:
			gf_mean_iv = mean_iv([gf_i.iv for gf_i in c])
			gf_name = ','.join(sorted([gf_i.name for gf_i in c]))
			print '\t'.join(map(str, [gf_mean_iv.chrom, gf_mean_iv.start, gf_mean_iv.end, gf_name, len(c), gf_mean_iv.strand]))

if __name__ == "__main__":
	l_fp = sys.argv[1:]
	mergePeaks = MergePeaksAcrossDigestions(l_fp)
	mergePeaks.make_adjacency_dict()
	mergePeaks.find_connected_components()
	mergePeaks.write_merged_peaks()
