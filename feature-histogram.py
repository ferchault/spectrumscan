#!/usr/bin/env python

def read_file():
	q = np.fromfile("features", sep=" ", dtype=np.int)
	q = q.reshape(-1, 12)
	return q

if __name__ == "__main__":
	q = read_file()
	for col in range(12):
		print (col, ' '.join(map(str, np.bincount(q[:, col]))))
