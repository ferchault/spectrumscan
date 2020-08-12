#!/usr/bin/env python
import sys
import networkx as nx

def mol2bonds(molstring):
    carbonyl = []
    hobindungen = []
    cobindungen = []
    oobindungen = []
    elements = "O" * 7 + "C" * 12 + "H" * 12
    bonds = molstring.strip().split()[1:]
    counts = {"OO": 0, "OC": 0, "OH": 0, "CC": 0, "CH": 0}
    g = nx.Graph()
    O_used = [0,0,0,0,0,0,0]
    for bond in bonds:
        parts = bond.split("-")
        a, b = int(parts[0]), int(parts[1])
        if (a < 7) & (b < 7):
            oobindungen.append([a, b])
            g.add_edge(a, b)
        if (a < 7) & (b > 6) & (b < 19):
            cobindungen.append([a, b])
        if (a < 7) & (b > 18):
            hobindungen.append([a, b])
        if a < 7:
            O_used[a] += 1
        if b < 7:
            O_used[b] += 1
    
        counts[elements[a]+elements[b]] += 1

    carbonyl = len([_ for _ in O_used if _ == 1])

    Ochains = [42, 0, 0,0,0,0,0,0]
    # exploits the fact that we exluded non-connected molecules, so no O-O-O-O rings and no O=O bonds. Therefore, connected components = chains
    for component in nx.connected_components(g):
        Ochains[len(component)] += 1
    print(counts['OO'], counts['OC'],counts['OH'],counts['CC'],counts['CH'],Ochains[7], Ochains[6], Ochains[5], Ochains[4], Ochains[3], Ochains[2],carbonyl)
    
for line in open(sys.argv[1]):
    mol2bonds(line)
