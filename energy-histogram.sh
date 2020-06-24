#!/bin/bash
for column in 6 7 8 9 10 11;
do 
	paste <(cat features | cut -d" " -f $column) <(cat geometries | sed 's/.*"energy": .//;s/".*//') | grep -v "^0" | grep -v ERROR | grep -v -- -- > $column.energies
done
