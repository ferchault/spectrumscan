#!/usr/bin/env python
import numpy as np

def read_file():
	q = np.fromfile("features", sep=" ", dtype=np.int)
	q = q.reshape(-1, 12)
	return q

if __name__ == "__main__":
	q = read_file()
	for col in range(12):
		counts = np.bincount(q[:, col])
		for i in range(len(counts)):
			if counts[i] > 0:
				print (col, i, counts[i])
