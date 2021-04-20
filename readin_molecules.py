import random
import os
import numpy as np
import qml

def get_energies(filename, key="cc2"):
    """ Returns a dictionary with excitation energies for each xyz-file.
    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    energies = dict()

    for line in lines:
        tokens = line.split()

        index_name = int(tokens[0])
        index_name = '/research/CNSM-Enrico/enrico/random/random-%05d.xyz' % index_name
        cc2 = float(tokens[2])  #third column S1 excitation energy
        #cc2 = float(tokens[3])  #fourth column S1 osc strength
        #cc2 = float(tokens[4])  #fiths column S2 excitation energy
        #cc2 = float(tokens[5])  #sixth column S2 osc strength
        pbe0 = float(tokens[2])

        if key=="cc2":
            energies[index_name] = cc2

        elif key=="delta":
            energies[index_name] = cc2 - pbe0
        else:
            energies[index_name] = cc2

    return energies

adc2_energy = get_energies("/research/CNSM-Enrico/enrico/random/MACHINE_LEARNING_NEW/energies_TZVP_new.dat", key="cc2")

lines_list = open('/research/CNSM-Enrico/enrico/random/MACHINE_LEARNING_NEW/index_nonegative_new.dat').read().splitlines()    #index_nonegative.dat contains only indices of xyz files that have positive S1 excitation energy
newlist = []
newadc2_energy = dict()
for f in lines_list:
    xxx = '/research/CNSM-Enrico/enrico/random/random-%05d.xyz' % int(f.strip())
    if adc2_energy[xxx] > 0:
         newlist.append(xxx)
         newadc2_energy[xxx] = adc2_energy[xxx]


compounds = [qml.Compound(xyz=f) for f in newlist]

for mol in compounds:
    mol.properties = newadc2_energy[mol.name]

print("Total Number of molecules:", len(compounds))

energy_cc2 = np.array([mol.properties for mol in compounds])



mean_prop = np.mean(energy_cc2)
print('Mittelwert über alle Werte::',mean_prop)

yyy = 0
zzz = 0
#print(energy_cc2)
#print(newadc2_energy)
for f in energy_cc2:
     xxx = abs(mean_prop - f)
     yyy = yyy +xxx
     zzz = zzz +f
     #print(f, xxx, yyy)

yyy = yyy / len(energy_cc2)
print('Nullinie:',yyy)


#print(compounds)

