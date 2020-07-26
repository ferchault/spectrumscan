#!/usr/bin/env python
import functools
import numpy as np
import networkx as nx

class IntegerPartitions(object):
    @staticmethod
    @functools.lru_cache(maxsize=64)
    def _do_partition(total, maxelements, around=None, maxdz=None):
        if (around is None) != (maxdz is None):
            raise ValueError("Cannot define center or radius alone.")

        if maxelements == 1:
            if around is not None and maxdz < abs(total - around[-maxelements]):
                return []
            else:
                return [[total]]
        res = []

        # get range to cover
        if around is None:
            first = 0
            last = total
            limit = None
        else:
            first = max(0, around[-maxelements] - maxdz)
            last = min(total, around[-maxelements] + maxdz)
        for x in range(first, last + 1):
            if around is not None:
                limit = maxdz - abs(x - around[-maxelements])
            for p in IntegerPartitions._do_partition(
                total - x, maxelements - 1, around, limit
            ):
                res.append([x] + p)
        return res

    @staticmethod
    def partition(total, maxelements, around=None, maxdz=None):
        if around is not None:
            return IntegerPartitions._do_partition(
                total, maxelements, tuple(around), maxdz
            )
        else:
            return IntegerPartitions._do_partition(total, maxelements)

def treat_molstring(molstring):
	# get CC bond count
	parts = molstring.split()
	parts = [_.split("-") for _ in parts[1:]]
	bonds = [(int(_[0]), int(_[1])) for _ in parts]
	ccbonds = 13
	for a,b in bonds:
		if 7 <= a <= 12 and 13 <= b <= 18:
			mixed = (a, b)
			break
	else:
		ccbonds = 12

	# get basevector
	base = np.zeros(12)
	for a, b in bonds:
		if 7 <= a <= 18:
			base[a-7] += 1
		if 7 <= b <= 18:
			base[b-7] += 1

		# test for terminal O
		a, b = min(a, b), max(a, b)
		if a < 7:
			for c, d in bonds:
				if (c == a and d != b) or (d == a and c != b):
					break
			else:
				base[b-7] += 1				
	print (base)
	
	# find partitions
	attributions = [(0, 1), (1, 2), (2,3), (3,4), (4,5), (5,0),    (6, 7), (7,8), (8,9), (9,10), (10, 11), (11, 6)]
	double_bonds = None
	if ccbonds == 13:
		attributions.append([_-7 for _ in mixed])
	for ndouble in range(6):
		ps = np.array((IntegerPartitions.partition(ndouble, ccbonds)))
		ps = ps[ps.max(axis=1) <= 1]
		additions = np.zeros((len(ps), 12))
		
		for col, attribution in enumerate(attributions):
			additions[:, attribution[0]] += ps[:, col]
			additions[:, attribution[1]] += ps[:, col]
		combined = base + additions
		valid = np.where((combined.min(axis=1) == 4) & (4 == combined.max(axis=1)))[0]
		if len(valid) > 0:
			double_bonds = ps[valid[0]]
			break

	if double_bonds is None:
		raise NotImplementedError("no valid pattern ?!")
	
	# build matching graph. Idea: two sets of nodes, left to right: double bond, right to left: single bond. All simple paths starting left except for the last step are conjugated double bonds
	G = nx.DiGraph()
	for a, b in bonds:
		if a > 18 or b > 18:
			# contains H
			continue
		if a < 7 and b < 7:
			# OO bond
			continue
		a, b = min(a, b), max(a, b)
		if a < 7:
			# OC bond
			for c, d in bonds:
				if (c == a and d != b) or (d == a and c != b):
					break
			else:
				# no other bond for this oxygen, carboxyl group
				G.add_edge("lhs_%d" % a, "rhs_%d" % b)
				G.add_edge("lhs_%d" % b, "rhs_%d" % a)
		else:
			# CC bond done separately
			continue
	for attribution, is_double in zip(attributions, double_bonds):
		a, b = attribution[0]+7, attribution[1]+7
		if is_double == 0:
			# single bond
			G.add_edge("rhs_%d" % a, "lhs_%d" % b)
			G.add_edge("rhs_%d" % b, "lhs_%d" % a)
		else:
			# double bond
			G.add_edge("lhs_%d" % a, "rhs_%d" % b)
			G.add_edge("lhs_%d" % b, "rhs_%d" % a)

	startnodes = [_ for _ in G.nodes() if _.startswith("lhs_")]
	endnodes = [_ for _ in startnodes if G.out_degree(_) == 0]
	found_paths = []
	for startnode in startnodes:
		# include the up to one neighbor that has a single bond
		in_edges = []
		for u,v in G.in_edges(startnode):
			in_edges = [u]
		for path in nx.all_simple_paths(G,startnode,endnodes + in_edges):
			if path[-1].startswith("lhs"):
				path = path[:-1] # strip last element which is single bonded
			for f in found_paths:
				if set(f) == set(path):
					break
			else:
				found_paths.append(path)
	# remove alternating directionality (enforced above)
	found_paths = [[_.split("_")[1] for _ in __] for __ in found_paths]
	# remove contained paths
	accepted = []
	for path in sorted(found_paths, key=len)[::-1]:
		for a in accepted:
			if set(path).issubset(set(a)):
				break
		else:
			accepted.append(path)
			print (path)
		
			
	
	
	
	
#treat_molstring("#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 0-4 1-6 2-4 2-5 3-5 3-6 7-19 8-20 9-21 10-22 11-23 0-12 13-24 14-25 1-15 16-26 17-27 17-28 18-29 18-30")
treat_molstring("#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 10-13 0-15 1-2 2-3 3-4 4-5 5-6 6-9 1-19 8-20 7-21 12-22 11-23 14-24 16-25 16-26 17-27 17-28 18-29 18-30")
