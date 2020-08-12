#!/usr/bin/env python
import sys

def mol2bonds(molstring):
    bonds = molstring.strip().split()[1:]
    O_used = [0,0,0,0,0,0,0]
    for bond in bonds:
        parts = bond.split("-")
        a, b = int(parts[0]), int(parts[1])
        if a < 7:
            O_used[a] += 1
        if b < 7:
            O_used[b] += 1
    
    carbonyl = len([_ for _ in O_used if _ == 1])
    print (carbonyl)
    
for line in open(sys.argv[1]):
    mol2bonds(line)
