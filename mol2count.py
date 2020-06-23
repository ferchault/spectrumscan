#!/usr/bin/env python

from time import time

start=time()

#f = open("/mnt/FOXTROT/PRODUCTION/batches/A/molecules", "r")
#f = open("test.A", "r")
f = open("line.dat", "r")
#f = open("/mnt/FOXTROT/PRODUCTION/batches/B/molecules", "r")
lines = f.readlines()
f.close()

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

    print('OH',hobindungen)
    print('OC',cobindungen)
    print('OO',oobindungen)
    
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

    print('Carbonara',carbonyl)

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
    print(counts,'7 OOOOOOO:',n_sept,'6 OOOOOO:',n_sext,'5 OOOOO:',n_quint,'4 OOOO:',n_quad,'3 OOO:',n_trip,'2 OO:',n_dup,'Carbonyl:',len(carbonyl))
    
    return
    #return counts
    
for line in lines:
    print(line)
    mol2bonds(line)


end=time()
print('TIME:',end-start)

#if __name__ == "__main__":
#	molstring = "#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 1-19 2-20 3-21 4-22 5-23 8-24 9-25 0-10 10-26 1-11 11-27 2-12 12-28 0-14 3-15 4-16 5-16 6-17 18-29 18-30 7-13"




