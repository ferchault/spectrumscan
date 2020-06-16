#!/usr/bin/env python
def mol2bonds(molstring):
	elements = "O" * 7 + "C" * 12 + "H" * 12
	bonds = molstring.strip().split()[1:]
	counts = {"OO": 0, "OC": 0, "OH": 0, "CC": 0, "CH": 0}
	for bond in bonds:
		parts = bond.split("-")
		a, b = int(parts[0]), int(parts[1])
		counts[elements[a]+elements[b]] += 1
	return counts

if __name__ == "__main__":
	molstring = "#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 1-19 2-20 3-21 4-22 5-23 8-24 9-25 0-10 10-26 1-11 11-27 2-12 12-28 0-14 3-15 4-16 5-16 6-17 18-29 18-30 7-13"
	print(mol2bonds(molstring))
