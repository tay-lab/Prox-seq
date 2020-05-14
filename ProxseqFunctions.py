# -*- coding: utf-8 -*-
"""
Author: Hoang Van temp_change
Address: Pritzker School of Molecular Engineering
         The University of Chicago
         Chicago, IL 60637, USA

This file contains the functions used to analyze PLA data obtained from Prox-seq
"""

# Import packages
import numpy as np
import math
import random
import pandas as pd

import scipy.spatial as spatial
import scipy.stats as stats
import matplotlib.pyplot as plt

import copy
import datetime

import itertools

# =============================================================================
# # Random point generators
# # Return an n-by-3 array, where the columns are x, y and z coordinates
# =============================================================================
def randomPointGen2D(n):
    # Generate n random points on the surface of a unit sphere
    # Ref: http://mathworld.wolfram.com/SpherePointPicking.html
    
    theta = np.random.uniform(0, 2*math.pi, size=(n,))
    z = np.random.uniform(-1,1, size=(n,))
    x = np.sqrt(1-z**2)*np.cos(theta)
    y = np.sqrt(1-z**2)*np.sin(theta)

    return np.array([x,y,z]).T

def randomPointGen3D(n):
    # Generate n random points inside a unit sphere
    # Ref: https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability/87238#87238
    u = np.random.uniform(0,1, size=(n,))
    x = np.random.normal(size=(n,))
    y = np.random.normal(size=(n,))
    z = np.random.normal(size=(n,))
    
    rescale = np.power(u,1/3) / np.sqrt(x*x+y*y+z*z)
    
    return (np.array([x,y,z]) * rescale).T

# =============================================================================
# # Simulate single cell PLA data
# # Assume there are N protein targets, and some of the proteins form pairwise complexes
# # Assume non-complex proteins and protein complexes are randmoly distributed on a sphere
# # Assume saturated antibody binding: all proteins and protein complexes are bound by PLA probes
# # If a pair of PLA probes A and B are within a certain distance (ie, the ligation distance), they are ligated
# # If more than one pair of probes are within the ligation distance, there are 2 options: all of them are ligated, or only one pair is
# =============================================================================
def simulatePLA(num_complex, probeA_ns, probeB_ns, cell_d=10000, PLA_dist=50,
                n_cells=100, ligation_efficiency=1, ligate_all=False,
                cell_variance=True, mode='2D',
                seed_num=None, sep=':'):
    # A function to simulate PLA counts of a cocktail of N targets
    # Input:
    #   num_complex (N-by-N): the number of complexes on each cell (entry [i,j] is the abundance of complex i:j)
    #   probeA_ns (N-by-1): the number of expressed proteins bound by probe A (entry [i] is the abundance of non-complex forming protein i, bound by PLA probe A)
    #   probeB_ns (N-by-1): the number of expressed proteins bound by probe B (entry [j] is the abundance of non-complex forming protein j, bound by PLA probe B)
    #   cell_d: cell diameter in nanometer
    #   PLA_dist: ligation distance in nanometer
    #   n_cells: number of cells to simulate
    #   ligation_efficiency: the chance of a PLA pair being ligated
    #   ligate_all: whether only 1 PLA pair or all pairs are allowed to be ligated
    #   cell_variance: whether to simulate variance due to cell size (log normal distribution)
    #   mode: 2D (PLA probes are on cell surface, default) or 3D (PLA probes are intracellular)
    #   seed_num: seed number for RNG
    #   sep: separator format for PLA product (default is ':')
    # Output:
    #   dge: simulated PLA count data
    #   dge_actual_complexes: simulated complex abundance

    # Seed number
    np.random.seed(seed_num)
    random.seed(seed_num)

    # Number of targets
    num_complex = np.array(num_complex)
    N_targets = num_complex.shape[0]

    # Initialize dge dictionary
    dge = {}
    # Look up dictionary: key is the complex identity, value is its index (or row number in dge matrix)
    complex_ind = {f'{i}{sep}{j}':(i*N_targets+j) for i in range(N_targets) for j in range(N_targets)}

    # Index matrix: to track the identity of each probe A and B
    probeA_ind = np.array([s.split(sep)[0] for s in complex_ind.keys()]) # element ij = i (ie, target of probe A)
    probeB_ind = np.array([s.split(sep)[1] for s in complex_ind.keys()]) # element ij = j (ie, target of probe B)

    # Convert list to numpy array
    probeA_ns = np.array(probeA_ns)
    probeB_ns = np.array(probeB_ns)

    # Data frame to store scaled complex amount of each single cell
    dge_actual_complexes = {}

    # Start simulation
    # Iterate through each single cell
    print(f'{datetime.datetime.now().replace(microsecond=0)}     Start simulation')
    for cell_i in range(n_cells):

        dge[cell_i] = np.zeros((N_targets**2,))

        # Add cell variance with log-normal distribution
        if cell_variance:
            scale_i = np.random.lognormal(mean=0,sigma=0.5,size=1)
            num_complex_i = (num_complex*scale_i).round().astype(int)
            probeA_ns_i = (probeA_ns*scale_i).round().astype(int)
            probeB_ns_i = (probeB_ns*scale_i).round().astype(int)
        else:
            num_complex_i = copy.deepcopy(num_complex).astype(int)
            probeA_ns_i = copy.deepcopy(probeA_ns).astype(int)
            probeB_ns_i = copy.deepcopy(probeB_ns).astype(int)

        # Save the actual complex amount
        dge_actual_complexes[cell_i] = num_complex_i.reshape(-1,)

        # Randomly distribute the protein targets
        if mode == '2D':
            protein_target_i = cell_d/2*randomPointGen2D(num_complex_i.sum())
        elif mode == '3D':
            protein_target_i = cell_d/2*randomPointGen3D(num_complex_i.sum())

        # Probe binding is assumed to be saturated, so all protein targets have probe A and probe B
        probeA_i = copy.deepcopy(protein_target_i)
        probeB_i = copy.deepcopy(protein_target_i)

        # Probe A target
        probeA_target = np.repeat(probeA_ind.flatten(), num_complex_i.flatten())
        probeB_target = np.repeat(probeB_ind.flatten(), num_complex_i.flatten())

        # Add non-specific binding
        if probeA_ns.sum() > 0:
            if mode == '2D':
                probeA_i = np.vstack((probeA_i, cell_d/2*randomPointGen2D(probeA_ns_i.sum())))
            elif mode == '3D':
                probeA_i = np.vstack((probeA_i, cell_d/2*randomPointGen3D(probeA_ns_i.sum())))
            probeA_target = np.concatenate((probeA_target, np.repeat(range(N_targets),probeA_ns_i)))
        if probeB_ns.sum() > 0:
            if mode == '2D':
                probeB_i = np.vstack((probeB_i, cell_d/2*randomPointGen2D(probeB_ns_i.sum())))
            elif mode == '3D':
                probeB_i = np.vstack((probeB_i, cell_d/2*randomPointGen3D(probeB_ns_i.sum())))
            probeB_target = np.concatenate((probeB_target, np.repeat(range(N_targets),probeB_ns_i)))

        # Calculate pairwise euclidean distance
        pairwise_dist = spatial.distance.cdist(probeA_i, probeB_i, metric="euclidean")

        # Ligation


        # Store the index of ligated probe A and B => list of lists
        ligated_probe = [[],[]]
        for i in range(pairwise_dist.shape[0]):

            # Indices of probe B that can be ligated
            proximity_probeB = np.argwhere(pairwise_dist[i,:]<=PLA_dist).flatten()
            if ligate_all:
                # Each PLA probe can be ligated with all proximity PLA probes B
                dge[cell_i][[complex_ind[f"{probeA_target[i]}{sep}{probeB_target[s]}"] for s in proximity_probeB]] += len(proximity_probeB)

            else:
                # Each PLA probe A or B can only be ligated at most once
                # Iterate through each probe A, then check for probe B within the ligation distance
                # If more than 1 probe B is, then chose the partner probe B randomly to be ligated
                # The chosen probe B is excluded from further ligation

                probeB_blacklist = set([]) # index of excluded probes B

                ligation_set = set(proximity_probeB) - probeB_blacklist
                # Random ligation
                if len(ligation_set) > 0:
                    # Will ligation happen
                    if random.random() < ligation_efficiency*np.random.normal(loc=1,scale=0.2,size=1):
                        chosen = random.sample(ligation_set, 1)
                        probeB_blacklist.add(chosen[0])
                        ligated_probe[0].append(i)
                        ligated_probe[1].append(chosen[0])

                        # Save to dge dictionary
                        dge[cell_i][complex_ind[f"{probeA_target[i]}{sep}{probeB_target[chosen[0]]}"]] += 1

        # Keep track of time
        if (cell_i+1) % 10 == 0:
            print(f'{datetime.datetime.now().replace(microsecond=0)}     Processed {cell_i+1:>5} cells')

    # Convert dictionary to pandas data frame
    complex_ind_new = [f'{i+1}{sep}{j+1}' for i in range(N_targets) for j in range(N_targets)] # update protein id from 0 to 1, 1 to 2, etc.
    dge = pd.DataFrame(dge, index=complex_ind_new)
    dge_actual_complexes = pd.DataFrame(dge_actual_complexes, index=complex_ind_new)
    return (dge, dge_actual_complexes)


# =============================================================================
# # Function to calculate the expected count of a PLA product according to the noise model
# # Ei,j = (Xi,. + X.,j)/(X.,.)
# # where Xi,. means sum of Xi,j over all j
# =============================================================================
def calculateExpected(data, PLA_list=None, sep=':'):
    # Input:
    #   data: data frame of dge
    #   PLA_list: list of PLA products for which expected count is calculated
    #             if None, calculate expected count for all PLA products
    # Output:
    #   A data frame of expected count (rows = PLA, columns = single cells)

    # Initialize output
    if PLA_list is None:
        PLA_list = data.index
    output = pd.DataFrame(columns=data.columns, index=PLA_list)

    # Get AB1 and AB2 of each row of data
    AB1 = np.array([s.split(sep)[0] for s in data.index])
    AB2 = np.array([s.split(sep)[1] for s in data.index])
    for i in PLA_list:
        output.loc[i,:] = data.loc[AB1==i.split(sep)[0],:].sum(axis=0)*data.loc[AB2==i.split(sep)[1],:].sum(axis=0)/data.sum(axis=0)

    return output


# =============================================================================
# # Function for estimating complex abundance
# # This function is used to estimate the abundance of true complexes from a digital expression matrix
# # Iteratively solve a system of quadratic equations
# # Goal: find the adjustment value for each complex, such that the adjusted complex abundance lies on the line
# # The adjustment values are the complex true abundance
# # Start with asusming the adjustment values are 0
# # Solve for the adjustment value of each complex
# # If the adjustment values of a complex across single cells do not reject the null hypothesis of a t-test (average higher than mean_cutoff), then the adjustment value of the complex is set to 0
# # The adjustment values from the last iteration are used to update the adjustment values in the next iteration
# # The process stops when the values converge, or when the maximum number of iterations is reached
# =============================================================================
def estimateComplexes(data, non_complexes=[], non_express=[], mean_cutoff=1, sym_weight=0.25, df_guess=None, start_complex=[], nIter=100, tol=5, sep=':'):
    # Input:
    #   data: dge data frame (columns are single cells)
    #   non_complex: list of non-true complex
    #   non_express: list of protein targets that do not form complexes (eg, isotype antibodies)
    #   slope_cutoff: in each iteraction, complexes with (slope - 1) below this cutoff are considered non-complex
    #   sym_weight: the weight factor used to enforce symmetry condition
    #   df_guess: data frame of first guesses of true complex abundance (must be the same shape as df)
    #             if None, use 0 as the first guess
    #   start_complex: a list of complexes to estimate during the first iteration (ie, in the first iteration, only these complexes will be used for calculation)
    #               if empty, perform calculation as normal through all complexes in the first iteration
    #   nIter: max number of iterations to perform
    #   tol: the allowed maximum change in solution between current and last iteration
    #   sep: the separator convention in the PLA complexes (default is ':')
    # Output:
    #   A data frame with the same shape as df, containing estimated complex abundance

    # Convert input data frame into numpy array
    dge = data.copy().values
    # Convert non_complex and non_express lists to sets
    non_complexes = set(non_complexes)
    non_express = set(non_express)
    start_complex = set(start_complex)

    # Get a list of probe A and B targets
    probeA = np.array([s.split(sep)[0] for s in data.index])
    probeB = np.array([s.split(sep)[1] for s in data.index])


    # Initialize a numpy array to store estimated complex amount
    if df_guess is None:
        dge_out = np.zeros(dge.shape)
    else:
        dge_out = df_guess.values

    # Iteration
    loop_num = 0
    max_change = tol + 1
    while (loop_num < nIter) and (max_change > tol):
        # Array to store the change in the complex estimates
        temp_change = np.zeros(dge.shape) + tol + 1

        temp_dge = dge - dge_out
        for i in range(dge.shape[0]):

            # Go through start_complex first
            if len(start_complex) > 0:
                if (loop_num == 0) and (i not in start_complex):
                    temp_change[i,:] = 0
                    continue
            temp_complex = data.index[i]
            temp_probeA, temp_probeB = temp_complex.split(sep) # target of probe A and B

            # Apply the constraints
            if (temp_complex in non_complexes) or (temp_probeA in non_express) or (temp_probeB in non_express):
                temp_change[i,:] = 0
                continue

            temp_expected = temp_dge[probeA==temp_probeA,:].sum(axis=0)*temp_dge[probeB==temp_probeB,:].sum(axis=0)/temp_dge.sum(axis=0)
            temp_diff = dge[i,:] - temp_expected

            # Check to see if the estimated abundance passes the mean_cutoff
            # Ha: sample mean > mean_cutoff
            tval, tp = stats.ttest_1samp(temp_diff, mean_cutoff)
            if (tval < 0) or (tp/2 > 0.01):
                # check for symmetry
                temp_symmetry = dge_out[data.index==f"{temp_probeB}{sep}{temp_probeA}",:]
                if np.mean(temp_symmetry) > mean_cutoff:
                    temp_diff = sym_weight*temp_symmetry
                else:
                    temp_change[i,:] = 0
                    continue

            # Force negative values to be zero <---- should be done after t-test
            temp_diff[temp_diff < 0] = 0

            # Check if observed is 0 but estimated is non 0, then force the estimated to be 0 <---- should be done after t-test
            temp_diff[(temp_diff > 0) & (dge[i,:] == 0)] = 0

            # Store changes in the solutions/estimates
            temp_change[i,:] = temp_diff - dge_out[i,:]

            # Store the new solutions/estimates
            dge_out[i,:] = temp_diff

        # Round the adjustment amount
        dge_out = np.round(dge_out)
        # Save the maximum change in the solution for convergence check
        max_change = abs(temp_change).max()
        
        loop_num += 1

    print(f"estimateComplexes done: Loop number {loop_num}, tolerance {max_change:.2f}")
    return pd.DataFrame(data=dge_out, index=data.index, columns=data.columns)
