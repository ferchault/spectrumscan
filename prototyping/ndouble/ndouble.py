#!/usr/bin/env python
import functools
import numpy as np

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
	ccbonds = 12
	for a,b in bonds:
		if 7 <= a <= 12 and 13 <= b <= 18:
			break
	else:
		ccbonds = 13

	# get basevector
	base = np.zeros(12)
	for a, b in bonds:
		if 7 <= a <= 18:
			base[a-7] += 1
		if 7 <= b <= 18:
			base[b-7] += 1
	
	# find partitions
	for ndouble in range(6):
		ps = np.array((IntegerPartitions.partition(ndouble, ccbonds)))
		# convert to bonds
	
treat_molstring("#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 0-4 1-6 2-4 2-5 3-5 3-6 7-19 8-20 9-21 10-22 11-23 0-12 13-24 14-25 1-15 16-26 17-27 17-28 18-29 18-30")
