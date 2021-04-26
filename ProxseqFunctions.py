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
#
# # Cell variance: Poisson
# =============================================================================
def simulatePLA(num_complex, probeA_ns, probeB_ns, cell_d=10000, PLA_dist=50,
                n_cells=100, ligation_efficiency=1, ligate_all=False,
                cell_variance=False, mode='2D',
                seed_num=None, sep=':'):
    '''
    A function to simulate PLA counts of a cocktail of N targets.
    
    Parameters
    ----------
    num_complex : numpy array 
        An NA-by-NB array containing the number of complexes on each cell
        (NA and NB is the number of targets of probe A and B), and element ij is
        the abundance of complex i:j.
    probeA_ns : numpy array
        An NA-by-1 array containing the number of expressed proteins bound by
        probe A (entry [i] is the abundance of non-complex forming protein i,
                 bound by PLA probe A).
    probeB_ns : numpy array
        An NB-by-1 array containing the number of expressed proteins bound by
        probe B (entry [j] is the abundance of non-complex forming protein j,
                 bound by PLA probe B).
    cell_d : float
        The cell diameter in nanometer.
        Default is 50.
    PLA_dist : float
        The ligation distance in nanometer.
        Default is 10,000.
    n_cells : int
        The number of cells to simulate.
    ligation_efficiency : float
        The chance of a PLA pair being ligated (between 0 and 1).
    ligate_all : boolean
        Whether only 1 PLA pair or all pairs are allowed to be ligated.
        Default is False
    cell_variance : boolean
        Whether to simulate variance due to cell size (log normal distribution).
        Default is True
    mode : string
        '2D' (PLA probes are on cell surface, default) or '3D' (PLA probes are
                                                                intracellular).
    seed_num : float, optional
        The seed number for RNG.
        Default is None.
    sep: string, optional
        The separator format for PLA product.
        Default is ':'.
    
    Returns
    -------
    (dge, dge_complex_true, dge_ns_true)
    dge : pandas data frame
        The simulated PLA count data.
    dge_complex_true : pandas data frame
        The true complex abundance.
    dge_ns_true : pandas data frame
        The true abundance of non-complex forming probes.
        
    '''

    # Seed number
    np.random.seed(seed_num)
    random.seed(seed_num)

    num_complex = np.array(num_complex)
    
    # Number of probe A and B targets
    NA_targets = num_complex.shape[0]
    NB_targets = num_complex.shape[1]

    # Initialize dge dictionary
    dge = {}
    # Look up dictionary: key is the complex identity, value is its index (or row number in dge matrix)
    complex_ind = {f'{i}{sep}{j}':(i*NB_targets+j) for i in range(NA_targets) for j in range(NB_targets)}
    

    # Index matrix: to track the identity of each probe A and B
    probeA_ind = np.array([s.split(sep)[0] for s in complex_ind.keys()]) # element ij = i (ie, target of probe A)
    probeB_ind = np.array([s.split(sep)[1] for s in complex_ind.keys()]) # element ij = j (ie, target of probe B)

    # Data frame to store actual complex abundance of each single cell
    dge_complex_true = {}
    
    # Data frame to store actual non-complexing forming probe abundance of each single cell
    dge_ns_true = {}

    # Cell variance: Poisson
    # variance_probeA_ns = []
    # variance_probeB_ns = []
    # variance_complex = []
    # if cell_variance:
    #     for i in range(len(probeA_ns)):
    #         variance_probeA_ns.append(np.random.poisson(probeA_ns[i], size=n_cells))
    #     for i in range(len(probeB_ns)):
    #         variance_probeB_ns.append(np.random.poisson(probeB_ns[i], size=n_cells))
    #     for i in range(num_complex.shape[0]):
    #         variance_complex.append([])
    #         for j in range(num_complex.shape[1]):
    #             variance_complex[i].append(np.random.poisson(num_complex[i,j], size=n_cells))
    #     variance_probeA_ns = np.array(variance_probeA_ns)
    #     variance_probeB_ns = np.array(variance_probeB_ns)
    #     variance_complex = np.reshape(np.array(variance_complex),
    #                                   (num_complex.shape[0],num_complex.shape[1],n_cells))
    
    # Cell variance: negative binomial: variance = 10*mean
    variance_probeA_ns = []
    variance_probeB_ns = []
    variance_complex = []
    temp_p = 1/10
    if cell_variance:
        for i in range(len(probeA_ns)):
            variance_probeA_ns.append(np.random.negative_binomial(n=probeA_ns[i]*temp_p/(1-temp_p), p=temp_p, size=n_cells))
        for i in range(len(probeB_ns)):
            variance_probeB_ns.append(np.random.negative_binomial(n=probeB_ns[i]*temp_p/(1-temp_p), p=temp_p, size=n_cells))
        for i in range(num_complex.shape[0]):
            variance_complex.append([])
            for j in range(num_complex.shape[1]):
                if num_complex[i,j] == 0:
                    variance_complex[i].append(np.zeros(n_cells).astype(int))
                else:
                    variance_complex[i].append(np.random.negative_binomial(n=num_complex[i,j]*temp_p/(1-temp_p), p=temp_p, size=n_cells))
        variance_probeA_ns = np.array(variance_probeA_ns)
        variance_probeB_ns = np.array(variance_probeB_ns)
        variance_complex = np.reshape(np.array(variance_complex),
                                      (num_complex.shape[0],num_complex.shape[1],n_cells))
    
    # Start simulation
    # Iterate through each single cell
    print(f'{datetime.datetime.now().replace(microsecond=0)}     Start simulation')
    for cell_i in range(n_cells):

        dge[cell_i] = np.zeros((NA_targets*NB_targets,))

        # Add cell variance if required
        if cell_variance:
            probeA_ns_i = variance_probeA_ns[:,cell_i]
            probeB_ns_i = variance_probeB_ns[:,cell_i]
            num_complex_i = variance_complex[:,:,cell_i]
        else:
            num_complex_i = copy.deepcopy(num_complex).astype(int)
            probeA_ns_i = copy.deepcopy(probeA_ns).astype(int)
            probeB_ns_i = copy.deepcopy(probeB_ns).astype(int)

        # Save the true complex abundance
        dge_complex_true[cell_i] = num_complex_i.reshape(-1,)
        
        # Save the true non-complex forming abundance
        dge_ns_true[cell_i] = np.hstack((probeA_ns_i, probeB_ns_i))

        # Randomly distribute the protein complexes
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
            probeA_target = np.concatenate((probeA_target, np.repeat(range(NA_targets),probeA_ns_i)))
        if probeB_ns.sum() > 0:
            if mode == '2D':
                probeB_i = np.vstack((probeB_i, cell_d/2*randomPointGen2D(probeB_ns_i.sum())))
            elif mode == '3D':
                probeB_i = np.vstack((probeB_i, cell_d/2*randomPointGen3D(probeB_ns_i.sum())))
            probeB_target = np.concatenate((probeB_target, np.repeat(range(NB_targets),probeB_ns_i)))

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
    complex_ind_new = [f'{i+1}{sep}{j+1}' for i in range(NA_targets) for j in range(NB_targets)] # update protein id from 0 to 1, 1 to 2, etc.
    return (pd.DataFrame(dge, index=complex_ind_new),
            pd.DataFrame(dge_complex_true, index=complex_ind_new),
            pd.DataFrame(dge_ns_true, index=[f"{i+1}_A" for i in range(NA_targets)]+[f"{i+1}_B" for i in range(NB_targets)]))

# old variance model
def simulatePLAcopy(num_complex, probeA_ns, probeB_ns, cell_d=10000, PLA_dist=50,
                n_cells=100, ligation_efficiency=1, ligate_all=False,
                cell_variance=True, mode='2D',
                seed_num=None, sep=':'):
    '''
    A function to simulate PLA counts of a cocktail of N targets.
    
    Parameters
    ----------
    num_complex : numpy array 
        An NA-by-NB array containing the number of complexes on each cell
        (NA and NB is the number of targets of probe A and B), and element ij is
        the abundance of complex i:j.
    probeA_ns : numpy array
        An NA-by-1 array containing the number of expressed proteins bound by
        probe A (entry [i] is the abundance of non-complex forming protein i,
                 bound by PLA probe A).
    probeB_ns : numpy array
        An NB-by-1 array containing the number of expressed proteins bound by
        probe B (entry [j] is the abundance of non-complex forming protein j,
                 bound by PLA probe B).
    cell_d : float
        The cell diameter in nanometer.
        Default is 50.
    PLA_dist : float
        The ligation distance in nanometer.
        Default is 10,000.
    n_cells : int
        The number of cells to simulate.
    ligation_efficiency : float
        The chance of a PLA pair being ligated (between 0 and 1).
    ligate_all : boolean
        Whether only 1 PLA pair or all pairs are allowed to be ligated.
        Default is False
    cell_variance : boolean
        Whether to simulate variance due to cell size (log normal distribution).
        Default is True
    mode : string
        '2D' (PLA probes are on cell surface, default) or '3D' (PLA probes are
                                                                intracellular).
    seed_num : float, optional
        The seed number for RNG.
        Default is None.
    sep: string, optional
        The separator format for PLA product.
        Default is ':'.
    
    Returns
    -------
    (dge, dge_complex_true, dge_ns_true)
    dge : pandas data frame
        The simulated PLA count data.
    dge_complex_true : pandas data frame
        The true complex abundance.
    dge_ns_true : pandas data frame
        The true abundance of non-complex forming probes.
        
    '''

    # Seed number
    np.random.seed(seed_num)
    random.seed(seed_num)

    num_complex = np.array(num_complex)
    
    # Number of probe A and B targets
    NA_targets = num_complex.shape[0]
    NB_targets = num_complex.shape[1]

    # Initialize dge dictionary
    dge = {}
    # Look up dictionary: key is the complex identity, value is its index (or row number in dge matrix)
    complex_ind = {f'{i}{sep}{j}':(i*NB_targets+j) for i in range(NA_targets) for j in range(NB_targets)}
    

    # Index matrix: to track the identity of each probe A and B
    probeA_ind = np.array([s.split(sep)[0] for s in complex_ind.keys()]) # element ij = i (ie, target of probe A)
    probeB_ind = np.array([s.split(sep)[1] for s in complex_ind.keys()]) # element ij = j (ie, target of probe B)

    # Data frame to store actual complex abundance of each single cell
    dge_complex_true = {}
    
    # Data frame to store actual non-complexing forming probe abundance of each single cell
    dge_ns_true = {}

    # Start simulation
    # Iterate through each single cell
    print(f'{datetime.datetime.now().replace(microsecond=0)}     Start simulation')
    for cell_i in range(n_cells):

        dge[cell_i] = np.zeros((NA_targets*NB_targets,))

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

        # Save the true complex abundance
        dge_complex_true[cell_i] = num_complex_i.reshape(-1,)
        
        # Save the true non-complex forming abundance
        dge_ns_true[cell_i] = np.hstack((probeA_ns_i, probeB_ns_i))

        # Randomly distribute the protein complexes
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
            probeA_target = np.concatenate((probeA_target, np.repeat(range(NA_targets),probeA_ns_i)))
        if probeB_ns.sum() > 0:
            if mode == '2D':
                probeB_i = np.vstack((probeB_i, cell_d/2*randomPointGen2D(probeB_ns_i.sum())))
            elif mode == '3D':
                probeB_i = np.vstack((probeB_i, cell_d/2*randomPointGen3D(probeB_ns_i.sum())))
            probeB_target = np.concatenate((probeB_target, np.repeat(range(NB_targets),probeB_ns_i)))

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
    complex_ind_new = [f'{i+1}{sep}{j+1}' for i in range(NA_targets) for j in range(NB_targets)] # update protein id from 0 to 1, 1 to 2, etc.
    return (pd.DataFrame(dge, index=complex_ind_new),
            pd.DataFrame(dge_complex_true, index=complex_ind_new),
            pd.DataFrame(dge_ns_true, index=[f"{i+1}_A" for i in range(NA_targets)]+[f"{i+1}_B" for i in range(NB_targets)]))


# =============================================================================
# # Function to calculate the total protein abundance
# # Homodimers result in twice the count for the protein
# =============================================================================
def calculateProteinAbundance(data, sep=':'):
    '''
    Calculate protein abundance by summing the counts of each target, regardless of probe A or B.
    

    Parameters
    ----------
    data : pandas data frame
        Columns are cell barcodes, rows are PLA products
    sep : string, optional
        The separator format for PLA product.
        Default is ':'.

    Returns
    -------
    Returns a pandas data frame.
    Each row is the total abundance of proteins 1, 2,...

    '''
    
    # Get AB1 and AB2 of each row of data
    AB1 = np.array([s.split(sep)[0] for s in data.index])
    AB2 = np.array([s.split(sep)[1] for s in data.index])
    
    # Get the unique antibody targets
    AB_unique = np.unique(np.concatenate((AB1,AB2)))
    AB_unique.sort()
    
    # Initialize output dataframes
    output = pd.DataFrame(0, index=AB_unique, columns=data.columns)
    
    for i in output.index:
        output.loc[i,:] = (data.loc[AB1==i,:]).sum(axis=0) + (data.loc[AB2==i,:]).sum(axis=0)
        # output.loc[i,:] = data.loc[(AB1==i) | (AB2==i),:].sum(axis=0)
        # # Add homodimer counts once more
        # if f"{i}{sep}{i}" in data.index:
        #     output.loc[i,:]  = output.loc[i,:] + data.loc[f"{i}{sep}{i}",:]
    
    return output


# =============================================================================
# # Function to calculate the probe A and B abundance
# =============================================================================
def calculateProbeAbundance(data, sep=':'):
    '''
    Calculate probe abundance by summing the counts of probes A and B of each target.

    Parameters
    ----------
    data : pandas data frame
        Columns are cell barcodes, rows are PLA products
    sep : string, optional
        The separator format for PLA product.
        Default is ':'.

    Returns
    -------
    Returns a pandas data frame.
    Each row is the total abundance of probe A1, A2,... and B1, B2,...

    '''
    
    # Get AB1 and AB2 of each row of data
    AB1 = np.array([s.split(sep)[0] for s in data.index])
    AB2 = np.array([s.split(sep)[1] for s in data.index])
    
    # Get the unique AB1 and AB2 probe targets
    AB1_unique = list(set(AB1))
    AB2_unique = list(set(AB2))
    AB1_unique.sort()
    AB2_unique.sort()
    
    # Initialize temporary data frames
    output1 = pd.DataFrame(0, index=AB1_unique, columns=data.columns) # store abundance of all probe A
    output2 = pd.DataFrame(0, index=AB2_unique, columns=data.columns) # store abundance of all probe B
    
    for i in output1.index:
        output1.loc[i,:] = data.loc[AB1==i,:].sum(axis=0)
    for i in output2.index:
        output2.loc[i,:] = data.loc[AB2==i,:].sum(axis=0)
    output1.index = [f"{i}_A" for i in output1.index]
    output2.index = [f"{i}_B" for i in output2.index]
    
    return pd.concat([output1,output2])

# =============================================================================
# # Function to calculate the expected count of a PLA product according to the noise model
# # Ei,j = (Xi,. + X.,j)/(X.,.)
# # where Xi,. means sum of Xi,j over all j
# =============================================================================
def calculateExpected(data, PLA_list=None, sep=':'):
    '''
    Calculate the expected count of a PLA product using marginal probabilities.
    
    Parameters
    ----------
    data : pandas data frame
        Input digital PLA count matrix.
    PLA_list: list, optional
        List of PLA products for which expected count is calculated.
        If None (the default), calculate expected count for all PLA products.
                  
    Returns
    -------
    A data frame of expected count (rows = PLA, columns = single cells).
    
    '''

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
def estimateComplexes(data, non_complex=[], mean_cutoff=1, sym_weight=0.25,
                      df_guess=None, start_complex=[], nIter=100, tol=5, sep=':'):
    '''
    Estimate complex abundance by iteratively solving a system of quadratic equations.
    
    Parameters
    ----------
    data : pandas data frame
        Input digital PLA expression matrix (PLA products x single cells).
    non_complex : list
        List of PLA products or proteins that do no form protein complexes.
        Example: X:Y means X:Y does not form a complex, while X means X does
        not form complexes with any other proteins.
        Default is [].
    non_express : list
        List of protein targets that do not form complexes (eg, isotype antibodies).
        Default is [].
    mean_cutoff : float
        PLA products whose estimated complex abundance at each iteration fails
        the 1-sided t-test sample mean>mean_cutoff is kept as 0.
    sym_weight : float (0 <= sym_weight <= 1).
        The weight factor used to enforce symmetry condition.
    df_guess : pandas data frame
        First guesses of true complex abundance (must be the same shape as data).
        If None (the default), use 0 as the first guess.
    start_complex : list
         List of complexes to estimate during the first iteration (ie, in the
         first iteration, only these complexes will be used for calculation).
    nIter : int
        Max number of iterations to perform.
    tol : float
        If the change in solution between current and last iteration is below
        this value, convergence is reached.
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A data frame with the same shape as df, containing estimated complex abundance
    
    '''

    # Convert input data frame into numpy array
    dge = data.copy().values
    # Convert non_complexes list to sets
    non_complex = set(non_complex)
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
            if (temp_complex in non_complex) or (temp_probeA in non_complex) or (temp_probeB in non_complex):
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

            # Check if observed is 0 but estimated is non 0, then force the estimated to be 0
            # This should only be done after t-test
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

# =============================================================================
# # Function to estimate complex abundance, using non-complex forming probe A and B
# # Suppose probe Ak and Bl do not form any complex
# # The expected count of non-complex PLA product i:j is
# #             Eij = Xil*Xkj/Xkl
# #     where Xij is the observed count of PLA product i:j
# # The complex abundance of i:j is Xij - Eij
# =============================================================================
def estimateComplexes2(data, refA, refB, sep=':'):
    '''
    Estimate complex abundance by using non-complex forming probes

    Parameters
    ----------
    data : pandas data frame
        Input digital PLA expression matrix (PLA products x single cells).
    refA : list
        List of reference non-complex forming probe A
    refB : list
        List of reference non-complex forming probe B
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A data frame with the same shape as df, containing estimated complex abundance

    '''
    out = pd.DataFrame(0, index=data.index, columns=data.columns)
    
    for i in out.index:
        AB1, AB2 = i.split(sep)
        if (AB1 == refA[0]) or (AB2 == refB[0]):
            continue
        
        out.loc[i,:] = data.loc[i,:] - data.loc[f"{AB1}{sep}{refB[0]}",:]*data.loc[f"{refA[0]}{sep}{AB2}",:]/data.loc[f"{refA[0]}{sep}{refB[0]}",:]
    
    return out
    

# =============================================================================
# # Function to estimate probe abundance and complex abundance using the multinomial model
# =============================================================================
def MultinomialObjFun(counts, Ai, Bj, cij, log_barrier=None, t_prime=2):
    '''
    Objective function: negative log-likehood function.

    Parameters
    ----------
    counts : 2D numpy array
        Array of PLA product count of a single cell, where element ij is the
        count of PLA product i:j.
    Ai : 1D numpy array
        Array of current estimated abundance of all probes A, where i-th element
        is the estimated abundance of probe Ai.
    Bj : 1D numpy array
        Array of current estimated abundance of all probes B.
    cij : 2D numpy array
        Array of estimated complex abundance of a single cell, where element ij
        is the count of complex i:j.
    log_barrier: None or float, optional
        The value of gamma used in the log barrier function. If None, the log
        barrier is not applied to the objective function.
        Default is None.
    t_prime: float, optional
        The value used in the log barrier function. Only used when log_barrier
        is not None.
        Default is 2.

    Returns
    -------
    The value of the negative log-likehood function

    '''
    # p: the scaling factor for product Ai*Bj
    # p = 1/(counts.sum())
    p = 1/1000
    
    B_grid, A_grid = np.meshgrid(Bj, Ai)
    prod_AB = A_grid * B_grid
    
    # Deal with zeros in log
    term1 = np.log(p*prod_AB + cij)
    term1[(prod_AB==0) & (cij==0)] = 0
    
    out = -(counts*term1).sum() + counts.sum()*np.log((p*prod_AB+cij).sum())
    
    # Add log-barrier
    if log_barrier is not None:
        if np.any((log_barrier*A_grid-cij < 0) | (log_barrier*B_grid-cij < 0)):
            out = np.inf
        else:
            tempA = np.log(log_barrier*A_grid - cij)
            tempB = np.log(log_barrier*B_grid - cij)
            tempA[(A_grid==0) & (cij==0)] = 0
            tempB[(B_grid==0) & (cij==0)] = 0
            out -= (tempA + tempB).sum()/t_prime

    return out

def GDSolver(data, non_complex=[],
             Ai0=None, Bj0=None, cij0=None,
             log_barrier=None, t_prime=2,
             nIter=100, step0=10,
             tol1=1e-3, tol2=1e-3, tol3=1e-3, sep=':'):
    '''
    Projected gradient descent algorithm to minimize MultinomialObjFun().

    Parameters
    ----------
    data : pandas series
        Array storing the PLA products count of a single cell, ordered by probe
        A target, then probe B target (ie, [1:1, 1:2, 1:3,..., 3:1, 3:2, 3:3,...]).
    non_complex : list, optional
        List of PLA products or proteins that do no form protein complexes.
        Example: X:Y means X:Y does not form a complex, while X means X does
        not form complexes with any other proteins.
        Default is [].
    Ai0 : numpy array, optional
        Starting values for probe A abundance Ai. If None, use the abundance of
        the corresponding PLA probes.
        Default is None.
    Bj0 : numpy array, optional
        Starting values for probe B abundance Bj. If None, use the abundance of
        the corresponding PLA probes.
        Default is None.
    cij0 : numpy array, optional
        Starting values for complex abundance cij. If None, use 0.
        Default is None.
    log_barrier: None or float, optional
        The value of gamma used in the log barrier function. If None, the log
        barrier is not applied to the objective function.
        Default is None.
    t_prime: float, optional
        The value used in the log barrier function. Only used when log_barrier
        is not None.
        Default is 2.
    nIter : int, optional
        The number of maximum iterations to perform.
        Default is 100.
    step0 : float, optional
        Initial step size for gradient descent.
        Default is 10.
    tol1 : float, optional
        Convergence condition 1: if the gradient's norm is below this value,
        convergence is reached.
        Default is 1e-3.
    tol2 : float, optional
        Convergence condition 2: if the change in negative log-likelihood
        function is below this value, convergence is reached.
        Default is 1e-3.
    tol3 : float, optional
        Convergence condition 3: if the percentage change in negative
        log-likelihood function is below this value, convergence is reached.
        Default is 1e-3.
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A dictionary containing the following data:
    Ai : 1D numpy array
        Array of current estimated abundance of all probes A, where i-th element
        is the estimated abundance of probe Ai.
    Bj : 1D numpy array
        Array of current estimated abundance of all probes B.
    cij : 2D numpy array
        Array of estimated complex abundance of a single cell, where element ij
        is the count of complex i:j.
    convergence_info : pandas data frame
        Return the convergence values for the algorithm at each iteration,
        including the objective function's value, the norm of the gradient,
        and the step size.
    d_Ak : 1D array
        Derivative of objective function w.r.t. Ak.
    d_Bk : 1D array
        Derivative of objective function w.r.t. Bk.
    d_ckl : 1D array
        Derivative of objective function w.r.t. ckl.

    '''
    
    # Identity of antibody A and B of each PLA product
    i_targets_all = [s.split(sep)[0] for s in data.index]
    j_targets_all = [s.split(sep)[1] for s in data.index]
    
    # Get the array of unique probe targets
    i_targets = np.sort(np.array(list(set(i_targets_all))))
    j_targets = np.sort(np.array(list(set(j_targets_all))))
    
    # Sort data frame by index alphabetically
    data = data[[f"{i}{sep}{j}" for i in i_targets for j in j_targets]]
    
    # Counts data reshaped into non_complex's format
    # rows are target of antibody A, columns are target of antibody B
    counts = data.to_numpy().reshape((len(i_targets), len(j_targets)))
    
    # # Remove rows and columns that only contain 0
    # i_targets = i_targets[counts.sum(axis=1)>0]
    # j_targets = j_targets[counts.sum(axis=0)>0]
    # counts = counts[counts.sum(axis=1)>0,:]
    # counts = counts[:,counts.sum(axis=0)>0]
    
    # Initializes Ak, Bk and ckl
    if Ai0 is None:
        Ai = counts.sum(axis=1)
    else:
        Ai = Ai0
    if Bj0 is None:
        Bj = counts.sum(axis=0)
    else:
        Bj = Bj0
    # Ak[:] = Ak.min()
    # Bk[:] = Bk.min()
    if cij0 is None:
        cij = np.zeros(counts.shape)
    else:
        cij = cij0
    
    # Scaling factor for product Ai*Bj
    # p = 1/(counts.sum())
    p = 1/1000

    # Convert list of non-complex-forming products into a matrix of indicator: 1 = non-complex forming
    # Same format as cij
    non_complex_indicator = np.zeros(cij.shape)
    if non_complex:  
        for s in non_complex:
            if sep in s:
                i1, i2 = s.split(sep)
                non_complex_indicator[i_targets==i1, j_targets==i2] = 1
            else:
                non_complex_indicator[i_targets==s,:] = 1
                non_complex_indicator[:,j_targets==s] = 1

    # A dictionary to store the value of negative log likelihood function after
    # each iteration, the norm of the gradient, and the corresponding step size
    convergence_info = {'ObjFun':[], 'gradnorm':[],'step':[]} 

    # Gradient descent
    # Projected gradient descent: handles positive bound constraint ()
    # log-barrier function: handles inequality constraint
    
    for _ in range(nIter):
        B_grid, A_grid = np.meshgrid(Bj, Ai)
        prod_AB = A_grid * B_grid
        
        # Current negative log likelihood
        current_ObjFun = MultinomialObjFun(counts, Ai, Bj, cij, log_barrier=log_barrier, t_prime=t_prime)
        
        # Derivative w.r.t. Ak
        d_Ak = - np.nansum(counts*p*B_grid/(p*prod_AB + cij),axis=1) + counts.sum()*p*B_grid.sum(axis=1)/(p*prod_AB+cij).sum()
        
        # Derivative w.r.t. Bk
        d_Bk = - np.nansum(counts*p*A_grid/(p*prod_AB + cij),axis=0) + counts.sum()*p*A_grid.sum(axis=0)/(p*prod_AB+cij).sum()

        # Derivative w.r.t. ckl
        c_const = counts.sum()/(p*prod_AB + cij).sum()
        d_ckl = -counts/(p*prod_AB + cij) + c_const
        
        
        # Due to projection and drop-out, sometimes the derivatives can be zero
        # Deal with these zeros by replacing np.nan and np.inf with 0
        d_Ak[~np.isfinite(d_Ak)] = 0
        d_Bk[~np.isfinite(d_Bk)] = 0
        d_ckl[~np.isfinite(d_ckl)] = 0
        
        # Add the derivative of the log-barrier function
        if log_barrier is not None:
            # Derivatives of log-barrier function
            dphi_Ak = 1/t_prime/(log_barrier*A_grid - cij)
            dphi_Ak[~np.isfinite(dphi_Ak)] = 0
            dphi_Bk = 1/t_prime/(log_barrier*B_grid - cij)
            dphi_Bk[~np.isfinite(dphi_Bk)] = 0
            dphi_ckl = -(dphi_Ak + dphi_Bk)

            # Add the contribution from log-barrier derivatives
            d_Ak -= log_barrier*dphi_Ak.sum(axis=1)
            d_Bk -= log_barrier*dphi_Bk.sum(axis=0)
            d_ckl -= dphi_ckl
        
        
        # Stop condition: sufficiently small gradient
        grad_normsq = (d_Ak**2).sum() + (d_Bk**2).sum() + (d_ckl**2).sum()
        convergence_info['gradnorm'].append(np.sqrt(grad_normsq))
        if np.sqrt(grad_normsq) < tol1:
            break
        
        # Step size: backtracking line search (Convex optimization, B & V, algorithm 9.2)
        # alpha = 0.01, beta = 0.8
        step = step0
        while True: # decrease step size
            temp_Ai = Ai - step*d_Ak
            temp_Bj = Bj - step*d_Bk
            temp_cij = cij - step*d_ckl
            
            # Projection
            temp_Ai[temp_Ai < 0] = 0
            temp_Bj[temp_Bj < 0] = 0
            temp_cij[temp_cij < 0] = 0
            
            # Backtracking line search
            if (MultinomialObjFun(counts, temp_Ai, temp_Bj, temp_cij,
                                  log_barrier=log_barrier, t_prime=t_prime) >
                current_ObjFun - 0.01*step*grad_normsq):
                step *= 0.8
            else:
                # Increment
                Ai = temp_Ai
                Bj = temp_Bj
                cij = temp_cij
                break

        # Force non-complex condition
        cij[non_complex_indicator==1] = 0
        
        # Update convergence information
        convergence_info['step'].append(step)
        new_ObjFun = MultinomialObjFun(counts, Ai, Bj, cij, log_barrier=log_barrier, t_prime=t_prime)
        convergence_info['ObjFun'].append(new_ObjFun)
        
        # Stop condition: negative loglikehood function converges
        if (abs(current_ObjFun - new_ObjFun) < tol2) or (abs(current_ObjFun - new_ObjFun)/current_ObjFun*100 < tol3):
            break
        
    # Create dataframe of probes A and B abundance, and complex abundance
    Ai = pd.DataFrame(Ai, index=i_targets, columns=['Ai'])
    Bj = pd.DataFrame(Bj, index=j_targets, columns=['Bj'])
    cij = pd.DataFrame(cij, index=i_targets, columns=j_targets)
    
    # Create dataframes of derivatives
    d_Ak = pd.DataFrame(d_Ak, index=i_targets, columns=['d_Ak'])
    d_Bk = pd.DataFrame(d_Bk, index=j_targets, columns=['d_Bk'])
    d_ckl = pd.DataFrame(d_ckl, index=i_targets, columns=j_targets)
    
    # For rare cases where all PLA products are forced to be non-complex
    if (not convergence_info['ObjFun']) or (not convergence_info['step']):
        convergence_info['ObjFun'] = [np.nan]
        convergence_info['step'] = [np.nan]
    
    return {'Ai':Ai, 'Bj':Bj, 'cij':cij, 'convergence_info':pd.DataFrame(convergence_info),
            'd_Ak':d_Ak, 'd_Bk':d_Bk, 'd_ckl':d_ckl}

def BarrierSolver(data, log_barrier, t_prime0=2, mu=10,
                  non_complex=[],
                  nIter=100, step0=10, tol=1e-6, sep=':'):
    '''
    Solve the optimization problem with inequality constraint using the barrier
    method.
    
    Parameters
    ----------
    data : pandas series
        Array storing the PLA products count of a single cell, ordered by probe
        A target, then probe B target (ie, [1:1, 1:2, 1:3,..., 3:1, 3:2, 3:3,...]).
    log_barrier : float
        The value of gamma used in the log barrier function. If None, the log
        barrier is not applied to the objective function.
    t_prime0 : int or float, optional
        The starting value of t prime.
        Default is 2.
    mu: int or float, optional
        The increment factor for t prime.
        Default is 10.
    non_complex : list, optional
        List of PLA products or proteins that do no form protein complexes.
        Example: X:Y means X:Y does not form a complex, while X means X does
        not form complexes with any other proteins.
        Default is [].
    nIter : int, optional
        The number of maximum iterations to perform during the centering step.
        Default is 100.
    step0 : float, optional
        Initial step size for GDSolver during the centering step.
        Default is 10.
    tol: float, optional
        The tolerance (duality gap) for stopping the outer iterations.
        Default is 1e-6.
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A dictionary containing the following data:
    Ai : 1D numpy array
        Array of current estimated abundance of all probes A, where i-th element
        is the estimated abundance of probe Ai.
    Bj : 1D numpy array
        Array of current estimated abundance of all probes B.
    cij : 2D numpy array
        Array of estimated complex abundance of a single cell, where element ij
        is the count of complex i:j.
    convergence_info : pandas data frame
        Return the convergence values for the algorithm at each iteration,
        including the objective function's value, the norm of the gradient,
        and the step size.
        
    '''
    
    # Calculate m, the number of inequalities
    i_targets_all = [s.split(sep)[0] for s in data.index]
    j_targets_all = [s.split(sep)[1] for s in data.index]
    i_targets = np.sort(np.array(list(set(i_targets_all))))
    j_targets = np.sort(np.array(list(set(j_targets_all))))
    m = 2*len(i_targets)*len(j_targets)
    
    # Sort data frame by index alphabetically
    data = data[[f"{i}{sep}{j}" for i in i_targets for j in j_targets]]
    
    # Initialize values of Ai, Bj and cij
    counts = data.to_numpy().reshape((len(i_targets), len(j_targets)))
    Ai = counts.sum(axis=1)
    Bj = counts.sum(axis=0)
    cij = np.zeros(counts.shape)
        
    # Dictionary to store convergence info
    convergence_info = {'t_prime':[], 'inner_nIter':[], 'ObjFun':[], 'gradnorm':[]}
    
    # Outer iteration
    t_prime = t_prime0
    while m/t_prime >= tol:
        
        # Centering step
        temp_out = GDSolver(data, Ai0=Ai, Bj0=Bj, cij0=cij,
                            log_barrier=log_barrier, t_prime=t_prime,
                            nIter=nIter, step0=step0)
        Ai = temp_out['Ai'].to_numpy().reshape((-1,))
        Bj = temp_out['Bj'].to_numpy().reshape((-1,))
        cij = temp_out['cij'].to_numpy()
        
        # Update convergence info
        convergence_info['t_prime'].append(t_prime)
        convergence_info['inner_nIter'].append(temp_out['convergence_info'].shape[0])
        convergence_info['ObjFun'].append(temp_out['convergence_info']['ObjFun'].iloc[-1])
        convergence_info['gradnorm'].append(temp_out['convergence_info']['gradnorm'].iloc[-1])
        
        # Update t prime
        t_prime *= mu
    
    # Create dataframe of probes A and B abundance, and complex abundance
    Ai = pd.DataFrame(Ai, index=i_targets, columns=['Ai'])
    Bj = pd.DataFrame(Bj, index=j_targets, columns=['Bj'])
    cij = pd.DataFrame(cij, index=i_targets, columns=j_targets)
    
    return {'Ai':Ai, 'Bj':Bj, 'cij':cij,
            'convergence_info':pd.DataFrame(convergence_info)}

# =============================================================================
# # Function to estimate probe abundance and complex abundance using the two-part model
# =============================================================================
def TwopartObjFun(counts, Ai, Bj, cij, p=(50/10000)**2):
    '''
    Calculate the objective function of the two-part model

    Parameters
    ----------
    counts : 2D numpy array
        Array of PLA product count of a single cell, where element ij is the
        count of PLA product i:j.
    Ai : 1D numpy array
        Array of current estimated abundance of all probes A, where i-th element
        is the estimated abundance of probe Ai.
    Bj : 1D numpy array
        Array of current estimated abundance of all probes B.
    cij : 2D numpy array
        Array of estimated complex abundance of a single cell, where element ij
        is the count of complex i:j.
    p: float, optional
        The probabability of forming a PLA product due to background ligation.
        Default is 2D case for ligation distance=50 nm, cell diameter=10um.

    Returns
    -------
    The value of the least-square objective function.

    '''
        
    B_grid, A_grid = np.meshgrid(Bj, Ai)
    
    out = 1/2*((counts - p*A_grid*B_grid - cij)**2)

    return np.nansum(out)

def GDSolver2(data, non_complex=[],
             Ai0=None, Bj0=None, cij0=None, p=(50/10000)**2,
             form='full', ini0=100, reduced_whitelist=[],
             nIter=100, step0=10, reduced_step0=200,
             tol1=1e-3, tol2=1e-3, tol3=1e-3, sep=':'):
    '''
    Use projected gradient descent to minimize the two-part objective function.

    data : pandas series
        Array storing the PLA products count of a single cell, ordered by probe
        A target, then probe B target (ie, [1:1, 1:2, 1:3,..., 3:1, 3:2, 3:3,...]).
    non_complex : list, optional
        List of PLA products or proteins that do no form protein complexes.
        Example: X:Y means X:Y does not form a complex, while X means X does
        not form complexes with any other proteins.
        Default is [].
    Ai0 : numpy array, optional
        Starting values for probe A abundance Ai. The values must be arranged
        alphabetically by target names.
        If None, set all values to 0.
        Default is None.
    Bj0 : numpy array, optional
        Starting values for probe B abundance Bj. The values must be arranged
        alphabetically by target names.
        If None, set all values to 0.
        Default is None.
    cij0 : numpy array, optional
        Starting values for complex abundance cij. If None, set all values to 0.
        Default is None.
    p: float, optional
        The probabability of forming a PLA product due to background ligation.
        Default is 2D case for ligation distance=50 nm, cell diameter=10um.
    form: string, optional
        The form of the optimization problem to be soved. If 'full', solve the
        full optimization problem. If 'reduced', solve the reduced problem, in
        which c_ij is set to 0.
        Default is 'full'.
    ini0: int or float, optional
        The value used to initialize Ai and Bj for solving the reduced form.
        Default is 100.
    reduced_whitelist: list, optional
        The list of PLA products that the reduced problem will solve for to initialize
        Ai and Bj for the full problem.
        Default is [].
    nIter : int, optional
        The number of maximum iterations to perform.
        Default is 100.
    step0 : float, optional
        Initial step size for gradient descent.
        Default is 10.
    reduced_step0 : float, optional
        Initial step size for gradient descent when solving the reduced problem
        to initialize for the full problem.
        Default is 200.
    tol1 : float, optional
        Convergence condition 1: if the gradient's norm is below this value,
        convergence is reached.
        Default is 1e-3.
    tol2 : float, optional
        Convergence condition 2: if the change in negative log-likelihood
        function is below this value, convergence is reached.
        Default is 1e-3.
    tol3 : float, optional
        Convergence condition 3: if the percentage change in negative
        log-likelihood function is below this value, convergence is reached.
        Default is 1e-3.
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A dictionary containing the following data:
    Ai : 1D numpy array
        Array of current estimated abundance of all probes A, where i-th element
        is the estimated abundance of probe Ai.
    Bj : 1D numpy array
        Array of current estimated abundance of all probes B.
    cij : 2D numpy array
        Array of estimated complex abundance of a single cell, where element ij
        is the count of complex i:j.
    convergence_info : pandas data frame
        Return the convergence values for the algorithm at each iteration,
        including the objective function's value, the norm of the gradient,
        and the step size.
    d_Ak : 1D array
        Derivative of objective function w.r.t. Ak.
    d_Bk : 1D array
        Derivative of objective function w.r.t. Bk.
    d_ckl : 1D array
        Derivative of objective function w.r.t. ckl.
    '''
    
    # Identity of antibody A and B of each PLA product
    i_targets_all = [s.split(sep)[0] for s in data.index]
    j_targets_all = [s.split(sep)[1] for s in data.index]
    
    # Get the array of unique probe targets
    i_targets = np.sort(np.array(list(set(i_targets_all))))
    j_targets = np.sort(np.array(list(set(j_targets_all))))
    
    # If the input data doesn't have exactly i*j PLA products, fill in the missing products with nan
    temp_all_possible_products = [f"{i}{sep}{j}" for i in i_targets for j in j_targets]
    for temp_row in [s for s in temp_all_possible_products if s not in data.index]:
        data.at[temp_row] = np.nan
      
    # Sort data frame by index alphabetically
    data = data.loc[temp_all_possible_products]
    
    # Counts data reshaped into non_complex's format
    # rows are target of antibody A, columns are target of antibody B
    counts = data.to_numpy().reshape((len(i_targets), len(j_targets)))   
    
    ###### Solve the reduced problem to initialize Ai and Bj for the full problem
    if (form == 'full'):
        if ((Ai0 is None) or (Bj0 is None)):
            temp_ini = GDSolver2(data=data, non_complex=non_complex,
                                 Ai0=None, Bj0=None, cij0=None, p=p,
                                 form='reduced', ini0=ini0, nIter=nIter*2, step0=reduced_step0,
                                 tol1=tol1, tol2=tol2, tol3=tol3, sep=sep)
            if Ai0 is None:
                Ai = temp_ini['Ai'].to_numpy().reshape((-1,))
            if Bj0 is None:
                Bj = temp_ini['Bj'].to_numpy().reshape((-1,))
        
        else:
            Ai = Ai0
            Bj = Bj0
        
        if cij0 is None:
            cij = np.zeros(counts.shape)
        else:
            cij = cij0
    
    ###### Initialize Ai and Bj for the reduced problem
    elif form == 'reduced':
        # Initializes Ak, Bk and ckl
        if Ai0 is None:
            Ai = np.zeros(i_targets.shape) + ini0
            Ai[counts.sum(axis=1) == 0] = 0
        else:
            Ai = Ai0
        if Bj0 is None:
            Bj = np.zeros(j_targets.shape) + ini0
            Bj[counts.sum(axis=0) == 0] = 0
        else:
            Bj = Bj0
        if cij0 is None:
            cij = np.zeros(counts.shape)
        else:
            cij = cij0
        
    else:
        raise ValueError("Invalid form input.")
        
    
    # Initialize the derivative w.r.t. c_kl and temp_cij (for form=='reduced')
    d_ckl = np.zeros(counts.shape)
    temp_cij = np.zeros(counts.shape)
    
    # Convert list of non-complex-forming products into a matrix of indicator: 1 = non-complex forming
    # Same format as cij
    non_complex_indicator = np.zeros(counts.shape)
    if non_complex:  
        for s in non_complex:
            if sep in s:
                i1, i2 = s.split(sep)
                non_complex_indicator[i_targets==i1, j_targets==i2] = 1
            else:
                non_complex_indicator[i_targets==s,:] = 1
                non_complex_indicator[:,j_targets==s] = 1

    # A dictionary to store the value of negative log likelihood function after
    # each iteration, the norm of the gradient, and the corresponding step size
    convergence_info = {'ObjFun':[], 'gradnorm':[],'step':[]} 

    # Gradient descent
    # Projected gradient descent: handles positive bound constraint ()
    
    
    for _ in range(nIter):
        B_grid, A_grid = np.meshgrid(Bj, Ai)
        temp_sum = (counts - p*A_grid * B_grid - cij) # calculate the term inside the sum

        # Current negative log likelihood
        current_ObjFun = TwopartObjFun(counts, Ai, Bj, cij, p=p)
        
        # Derivatives
        d_Ak = -p*np.nansum(B_grid*temp_sum, axis=1)
        d_Bk = -p*np.nansum(A_grid*temp_sum, axis=0)
        if form == 'full':
            d_ckl = -temp_sum
        
        # Stop condition: sufficiently small gradient
        grad_normsq = np.nansum(d_Ak**2) + np.nansum(d_Bk**2) + np.nansum(d_ckl**2)
        if np.sqrt(grad_normsq) < tol1:
            break
        
        # Step size: backtracking line search (Convex optimization, B & V, algorithm 9.2)
        # alpha = 0.01, beta = 0.8
        step = step0
        while True: # decrease step size
            temp_Ai = Ai - step*d_Ak
            temp_Bj = Bj - step*d_Bk
            
            # Projection
            temp_Ai[temp_Ai < 0] = 0
            temp_Bj[temp_Bj < 0] = 0
            
            # Update temp_cij only in form='full'
            if form == 'full':
                temp_cij = cij - step*d_ckl
                temp_cij[temp_cij < 0] = 0
            
            # Backtracking line search
            if (TwopartObjFun(counts, temp_Ai, temp_Bj, temp_cij, p=p)
                > current_ObjFun - 0.01*step*grad_normsq):
                step *= 0.8
                if step < 0.01:
                    break
            else:
                # Increment
                Ai = temp_Ai
                Bj = temp_Bj
                cij = temp_cij
                break
            
        # Force non-complex condition
        cij[non_complex_indicator==1] = 0
        
        # Update convergence information
        convergence_info['gradnorm'].append(np.sqrt(grad_normsq))
        convergence_info['step'].append(step)
        new_ObjFun = TwopartObjFun(counts, Ai, Bj, cij, p=p)
        convergence_info['ObjFun'].append(new_ObjFun)
        
        # Stop condition: negative loglikehood function converges
        if (abs(current_ObjFun - new_ObjFun) < tol2) or (abs(current_ObjFun - new_ObjFun)/current_ObjFun*100 < tol3):
            break
        
    # Create dataframe of probes A and B abundance, and complex abundance
    Ai = pd.DataFrame(Ai, index=i_targets, columns=['Ai'])
    Bj = pd.DataFrame(Bj, index=j_targets, columns=['Bj'])
    cij = pd.DataFrame(cij, index=i_targets, columns=j_targets)
    
    # Create dataframes of derivatives
    d_Ak = pd.DataFrame(d_Ak, index=i_targets, columns=['d_Ak'])
    d_Bk = pd.DataFrame(d_Bk, index=j_targets, columns=['d_Bk'])
    d_ckl = pd.DataFrame(d_ckl, index=i_targets, columns=j_targets)
    
    # For rare cases where all PLA products are forced to be non-complex
    if (not convergence_info['ObjFun']) or (not convergence_info['step']):
        convergence_info['ObjFun'] = [np.nan]
        convergence_info['step'] = [np.nan]
    
    # Display some information about convergence
    print(f"{form} problem: converge after {len(convergence_info['step'])} steps.")
    
    # return temp_ini
    return {'Ai':Ai, 'Bj':Bj, 'cij':cij, 'convergence_info':pd.DataFrame(convergence_info),
            'd_Ak':d_Ak, 'd_Bk':d_Bk, 'd_ckl':d_ckl}