#!/usr/bin/env python
import sys


def mol2bonds(molstring):
    cobindungen = []
    oobindungen = []
    hobindungen = []
    carbonyl = []
    ooo = []
    oooo = []
    ooooo = []
    oooooo = []
    ooooooo = []
    elements = "O" * 7 + "C" * 12 + "H" * 12
    bonds = molstring.strip().split()[1:]
    counts = {"OO": 0, "OC": 0, "OH": 0, "CC": 0, "CH": 0}
    for bond in bonds:
        parts = bond.split("-")
        a, b = int(parts[0]), int(parts[1])
        if (a < 7) & (b < 7):
            oobindungen.append([a, b])
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

    if len(oobindungen) > 1:
        for i in range(len(oobindungen)):
            first=oobindungen[i][0] 
            second=oobindungen[i][1] 
            for j in range(i+1,len(oobindungen)):
                if oobindungen[j][0] == first:
                    ooo.append([oobindungen[j][1],oobindungen[i][0],oobindungen[i][1]])
                if oobindungen[j][0] == second:
                    ooo.append([oobindungen[i][0],oobindungen[i][1],oobindungen[j][1]])
                if oobindungen[j][0] == second:
                    ooo.append([oobindungen[i][0],oobindungen[i][1],oobindungen[j][1]])
                if oobindungen[j][1] == second:
                    ooo.append([oobindungen[i][0],oobindungen[i][1],oobindungen[j][0]])
	    
        if len(ooo) > 1:
            for i in range(len(ooo)):
                first=ooo[i][0] 
                second=ooo[i][1] 
                third=ooo[i][2] 
                for j in range(i+1,len(ooo)):
                    if (ooo[j][0] == second) & (ooo[j][1] == third):
                        oooo.append([first, second,third,ooo[j][2]])
                    if (ooo[j][0] == second) & (ooo[j][1] == first):
                        oooo.append([ooo[j][2],first, second,third])
                    if (second == ooo[j][2]) & (third == ooo[j][1]):
                        oooo.append([first, second,third,ooo[j][0]])
        
            if len(oooo) > 1:
                for i in range(len(oooo)):
                    first=oooo[i][0] 
                    second=oooo[i][1] 
                    third=oooo[i][2] 
                    fourth=oooo[i][3] 
                    for j in range(i+1,len(oooo)):
                        if (first == oooo[j][1]) & (second == oooo[j][2]):
                            ooooo.append([oooo[j][0],first, second,third,fourth])
                        if (second == oooo[j][3]) & (third == oooo[j][2]):
                            ooooo.append([first, second,third,fourth,oooo[j][0]])
                        if (first == oooo[j][2]) & (second == oooo[j][1]):
                            ooooo.append([oooo[j][3],first, second,third,fourth])
            
                if len(ooooo) > 1:
                    for i in range(len(ooooo)):
                        first=ooooo[i][0] 
                        second=ooooo[i][1] 
                        third=ooooo[i][2] 
                        fourth=ooooo[i][3] 
                        fifth=ooooo[i][4] 
                        for j in range(i+1,len(ooooo)):
                            if (first == ooooo[j][3]) & (second == ooooo[j][2]):
                                oooooo.append([ooooo[j][4],first, second,third,fourth,fifth])
                            if (fourth == ooooo[j][2]) & (fifth == ooooo[j][1]):
                                oooooo.append([first, second,third,fourth,fifth,ooooo[j][0]])
                    if len(oooooo) > 1:
                        print('sext:',oooooo)
                        for i in range(len(oooooo)):
                            first=oooooo[i][0] 
                            second=oooooo[i][1] 
                            third=oooooo[i][2] 
                            fourth=oooooo[i][3] 
                            fifth=oooooo[i][4] 
                            sixth=oooooo[i][5] 
                            for j in range(i+1,len(oooooo)):
                                if (first == oooooo[j][1]) & (second == oooooo[j][2]):
                                    ooooooo.append([oooooo[j][0],first, second,third,fourth,fifth,sixth])
                                if (second == oooooo[j][5]) & (third == oooooo[j][4]):
                                    ooooooo.append([first, second,third,fourth,fifth,sixth,oooooo[j][0]])
        
    n_sept= len(ooooooo)
    n_sext= len(oooooo)-2*len(ooooooo)
    n_quint= len(ooooo)-2*n_sext-3*n_sept
    n_quad= len(oooo)-2*n_quint-3*n_sext-4*n_sept
    n_trip= len(ooo)-2*n_quad-3*n_quint-4*n_sext-5*n_sept
    n_dup= len(oobindungen)-2*n_trip-3*n_quad-4*n_quint-5*n_sext-6*n_sept
    print(counts['OO'], counts['OC'],counts['OH'],counts['CC'],counts['CH'],n_sept,n_sext,n_quint,n_quad,n_trip,n_dup,len(carbonyl))
    
for line in open(sys.argv[1]):
    mol2bonds(line)
