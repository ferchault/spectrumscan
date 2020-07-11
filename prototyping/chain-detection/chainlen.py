#!/usr/bin/env python
import networkx as nx

def mol2chains(molstring):
    elements = "O" * 7 + "C" * 12 + "H" * 12
    bonds = molstring.strip().split()[1:]
    g = nx.Graph()
    for bond in bonds:
        parts = bond.split("-")
        a, b = int(parts[0]), int(parts[1])
        if (a < 7) & (b < 7):
            g.add_edge(a, b)
    for component in nx.connected_components(g):
        print (len(component))
print ("should be 7")
mol2chains("#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 0-4 1-6 2-4 2-5 3-5 3-6 7-19 8-20 9-21 10-22 11-23 0-12 13-24 14-25 1-15 16-26 17-27 17-28 18-29 18-30")
print ("should be 5-2")
mol2chains("#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 0-3 0-4 1-4 1-5 2-6 5-7 7-19 8-20 6-10 10-21 11-22 11-23 12-24 12-25 13-26 14-27 15-28 15-29 2-16 17-30 3-18 9-16")
