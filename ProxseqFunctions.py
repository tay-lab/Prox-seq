# -*- coding: utf-8 -*-
"""
Author: Hoang Van temp_change
Address: Pritzker School of Molecular Engineering
         The University of Chicago
         Chicago, IL 60637, USA

This file contains the functions used to analyze PLA count data obtained from Prox-seq
"""

# Import packages
import numpy as np
import math
import random
import pandas as pd

import scipy.spatial as spatial
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

import copy
import datetime

import itertools

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
        output.loc[i,:] = data.loc[AB1==i.split(sep)[0],:].sum(axis=0).to_numpy()*data.loc[AB2==i.split(sep)[1],:].sum(axis=0).to_numpy()/data.sum(axis=0).to_numpy()

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
def estimateComplexes(data, non_complex=[], mean_cutoff=0, p_cutoff=0.05, p_adjust=True, sym_weight=0.25,
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
        Default is 0.
    p_cutoff : float
        The alpha level to decide if the 1-sided t-test is sinificant.
        Default is 0.05.
    p_adjust : boolean
        Whether to perform FDR correction for the one-sided t-test.
        Default is True.
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


        # List to store the one-sided t-test p-values
        tp_list = []

        temp_dge = dge - dge_out
        # First pass: get all the p-values
        for i in range(dge.shape[0]):

            # Go through start_complex first
            if len(start_complex) > 0:
                if (loop_num == 0) and (i not in start_complex):
                    tp_list.append(1)
                    continue
            temp_complex = data.index[i]
            temp_probeA, temp_probeB = temp_complex.split(sep) # target of probe A and B

            # Apply the constraints
            if (temp_complex in non_complex) or (temp_probeA in non_complex) or (temp_probeB in non_complex):
                tp_list.append(1)
                continue

            temp_expected = temp_dge[probeA==temp_probeA,:].sum(axis=0)*temp_dge[probeB==temp_probeB,:].sum(axis=0)/temp_dge.sum(axis=0)
            temp_diff = dge[i,:] - temp_expected

            # Check to see if the estimated abundance passes the mean_cutoff
            # Ha: sample mean > mean_cutoff
            tval, tp = stats.ttest_1samp(temp_diff, mean_cutoff)
            if (tval > 0):
                tp_list.append(tp/2)
            else:
                tp_list.append(1-tp/2)

        # Multiple comparison correction
        if p_adjust:
            _, tp_adj, _,_ = multipletests(tp_list, alpha=p_cutoff, method='fdr_bh')
        else:
            tp_adj = tp_list

        # Array to store the change in the complex estimates
        temp_change = np.zeros(dge.shape) + tol + 1
        # Second pass: calculate protein complex
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
            if (tp_adj[i] > p_cutoff):
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
