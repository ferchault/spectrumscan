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
    
        counts[elements[a]+elements[b]] += 1

    for i in range(len(cobindungen)):
        first=cobindungen[i][0]   # Sauerstoff
        second=cobindungen[i][1]   # Kohlenstoff
        if first in [k[0] for k in hobindungen]:
            dummy=0
        elif first in [k[0] for k in oobindungen]:
            dummy=0
        elif first in [k[0] for k in cobindungen[:i]]+[k[0] for k in cobindungen[i+1:]]:
            dummy=0
        else:
            carbonyl.append([first, second])

    Ochains = [42, 0, 0,0,0,0,0,0]
    # exploits the fact that we exluded non-connected molecules, so no O-O-O-O rings and no O=O bonds. Therefore, connected components = chains
    for component in nx.connected_components(g):
        Ochains[len(component)] += 1
    print(counts['OO'], counts['OC'],counts['OH'],counts['CC'],counts['CH'],Ochains[7], Ochains[6], Ochains[5], Ochains[4], Ochains[3], Ochains[2],len(carbonyl))
    
for line in open(sys.argv[1]):
    mol2bonds(line)
