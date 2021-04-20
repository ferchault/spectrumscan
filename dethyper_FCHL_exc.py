#!/usr/bin/env python
from __future__ import print_function
import numpy as np

import qml
import random
from copy import deepcopy
from qml.kernels import gaussian_kernel
from qml.math import cho_solve

from qml.representations import generate_fchl_acsf

from qml.kernels import get_local_kernel
from qml.kernels import get_local_symmetric_kernel

import tqdm  # progress bar

from tutorial_data import compounds
from tutorial_data import energy_cc2

if __name__ == "__main__":

    # For every compound generate a coulomb matrix
    Qall = []
    print('Generating representations')
    #for mol in tqdm.tqdm(compounds):
    for mol in compounds:

        #mol.generate_coulomb_matrix(size=29, sorting="row-norm")
        mol.representation = generate_fchl_acsf(mol.nuclear_charges, mol.coordinates, gradients=False, pad=33, elements=[1,6,7,8])
        # mol.generate_bob(size=23, asize={"O":3, "C":7, "N":3, "H":16, "S":1})
        Qall.append(mol.nuclear_charges)
            

    # Make a big 2D array with all the 
    X = np.array([mol.representation for mol in compounds])
    # X = np.array([mol.bob for mol in compounds])

    #split into training/validation and final test
    N_train_val = len(compounds)
    N_final_test = len(compounds)
    X_train_val = X[:N_train_val]
    Q_train_val = Qall[:N_train_val]
    Y_train_val = energy_cc2[:N_train_val]
    #Y_train_val = energy_delta[:N_train_val]


    sigma = [0.001*2**i for i in range(3,4)]
    #sigma = [0.001*2**i for i in range(3,8)]
    #sigma = [0.001*2**i for i in range(9,14)]
    #sigma = [0.001*2**i for i in range(15,20)]
    #nbins = 20
    nbins = 5
    xxx = (nbins-1)/nbins
    bins = range(nbins)
    #ntraining =[100, 1000, 3000, 5000, 9000]
    ntraining =[100]
    #for sig in tqdm.tqdm(sigma):
    for sig in sigma:
        print('Calculating Kernel Matrix')
        #calculate K for all N_train_val
        #K = gaussian_kernel(X_train_val, X_train_val, sig)
        K = get_local_kernel(X_train_val, X_train_val, Q_train_val, Q_train_val, sig)
        #lambda1 = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
        #lambda1 = [1e-10, 1e-9, 1e-8, 1e-6]
        lambda1 = [1e-7]
        for train in ntraining:
           for lamb in lambda1:
              maes = []
              for ibin in bins:
                 split = np.array(list(range(N_train_val)))
                 random.shuffle(split)
                 training_index = split[:int(xxx*train)]
                 val_index = split[-int(train/nbins):]
                 # Assign bin for training
                 Y_train = Y_train_val[training_index]
                 # Assign remaining bins to cross validation 
                 Y_val = Y_train_val[val_index]
      
                 C = deepcopy(K[training_index][:,training_index])
                 # Add a small lambda to the diagonal of the kernel matrix
                 C[np.diag_indices_from(C)] += lamb
       
                 # Use the built-in Cholesky-decomposition to solve
                 alpha = cho_solve(C, Y_train) 
                   
                 # Validate 
                 Yss = np.dot((K[training_index][:,val_index]).T, alpha)
                 diff = Yss  - Y_val
                 mae = np.mean(np.abs(diff))
                 maes.append(mae)
       
              # Calculate mean-absolute-error (MAE):
              s = np.std(maes)/np.sqrt(nbins)
              print(sig, lamb, train, sum(maes)/len(maes), s)
              #print(sig, lamb, train, sum(maes)/len(maes), s, file=open("direkt_S1_TZVP_osc_no_neg.dat", "a"))
              #print(sig, lamb, train, sum(maes)/len(maes), s, file=open("direkt_S1_TZVP_no_neg.dat", "a"))

