#!/usr/bin/env python

def aromatic_ring_count(molstring):
    """ Assumes the particular order of both elements and pairs in a molstring."""
    parts = molstring.split()
    parts = [_.split("-") for _ in parts[1:]]
    bonds = [(int(_[0]), int(_[1])) for _ in parts]
    for a,b in bonds:
        if 7 <= a <= 12 and 13 <= b <= 18:
            return 0

    rings = 0
    for site in range(7,13):
        count = 0
        for a, b in bonds:
            if site == a or site == b:
                count += 1
                if a < 7:
                    obondcount = 0
                    for c, d in bonds:
                        if c == a or d == a:
                            obondcount += 1
                    if obondcount == 1:
                        count -= 100
        if count != 3:
            break
    else:
        rings += 1

    for site in range(13, 19):
        count = 0
        for a, b in bonds:
            if site == a or site == b:
                count += 1
                if a < 7:
                    obondcount = 0
                    for c, d in bonds:
                        if c == a or d == a:
                            obondcount += 1
                    if obondcount == 1:
                        count -= 100
        if count != 3:
            break
    else:
        rings += 1

    return rings

if __name__ == "__main__":
    #molstring0 = "#OUT 7-13"
    #molstring0_2 = "#OUT 7-6"
    #molstring1 = "#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 0-4 1-6 2-4 2-5 3-5 3-6 7-19 8-20 9-21 10-22 11-23 0-12 13-24 14-25 1-15 16-26 17-27 17-28 18-29 18-30"
    #molstring2 = "#OUT 7-8 13-14 8-9 14-15 9-10 15-16 10-11 16-17 11-12 17-18 7-12 13-18 0-4 1-6 2-4 2-5 3-5 3-6 7-19 8-20 9-21 10-22 11-23 0-12 13-24 14-25 1-15 16-26 17-27 18-29"
    #print (aromatic_ring_count(molstring0))
    #print (aromatic_ring_count(molstring0_2))
    #print (aromatic_ring_count(molstring1))
    #print (aromatic_ring_count(molstring2))
    import sys
    for line in open(sys.argv[1]):
        print (aromatic_ring_count(line))
