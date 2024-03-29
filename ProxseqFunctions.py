# -*- coding: utf-8 -*-
"""
Author: Hoang Van Phan
Address: Tay Lab
         Pritzker School of Molecular Engineering
         The University of Chicago
         Chicago, IL 60637, USA

This file contains the functions used to analyze PLA product count data obtained
from Prox-seq
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
    Calculate protein abundance by summing the counts of each target, regardless
    of the target of probe A and B.

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

    # Get probeA and probeB of each row of data
    probeA = np.array([s.split(sep)[0] for s in data.index])
    probeB = np.array([s.split(sep)[1] for s in data.index])

    # Get the unique antibody targets
    AB_unique = np.unique(np.concatenate((probeA,probeB)))
    AB_unique.sort()

    # Initialize output dataframes
    output = pd.DataFrame(0, index=AB_unique, columns=data.columns)

    for i in output.index:
        output.loc[i,:] = (data.loc[probeA==i,:]).sum(axis=0) + (data.loc[probeB==i,:]).sum(axis=0)

    return output


# =============================================================================
# # Function to calculate the probe A and B abundance
# =============================================================================
def calculateProbeAbundance(data, sep=':'):
    '''
    Calculate probe abundance by summing the counts of probes A and B of each
    protein target.

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

    # Get probeA and probeB of each row of data
    probeA = np.array([s.split(sep)[0] for s in data.index])
    probeB = np.array([s.split(sep)[1] for s in data.index])

    # Get the unique probeA and probeB probe targets
    probeA_unique = list(set(probeA))
    probeB_unique = list(set(probeB))
    probeA_unique.sort()
    probeB_unique.sort()

    # Initialize temporary data frames
    output1 = pd.DataFrame(0, index=probeA_unique, columns=data.columns) # store abundance of all probe A
    output2 = pd.DataFrame(0, index=probeB_unique, columns=data.columns) # store abundance of all probe B

    for i in output1.index:
        output1.loc[i,:] = data.loc[probeA==i,:].sum(axis=0)
    for i in output2.index:
        output2.loc[i,:] = data.loc[probeB==i,:].sum(axis=0)
    output1.index = [f"{i}_A" for i in output1.index]
    output2.index = [f"{i}_B" for i in output2.index]

    return pd.concat([output1,output2])

# =============================================================================
# # Function to calculate the expected random count of a PLA product if there are
# # no protein interactions in the data
# # Ei,j = (Xi,. + X.,j)/(X.,.)
# # where Xi,. means sum of Xi,j over all j
# =============================================================================
def calculateExpected(data, PLA_list=None, sep=':'):
    '''
    Calculate the expected random count of a PLA product, if no protein interactions
    exist in the data.

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

    # Get probe A and B identity of each row of data
    probeA = np.array([s.split(sep)[0] for s in data.index])
    probeB = np.array([s.split(sep)[1] for s in data.index])
    for i in PLA_list:
        output.loc[i,:] = data.loc[probeA==i.split(sep)[0],:].sum(axis=0).to_numpy()*data.loc[probeB==i.split(sep)[1],:].sum(axis=0).to_numpy()/data.sum(axis=0).to_numpy()

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
def estimateComplexes(data, non_complex=[], mean_cutoff=1, p_cutoff=0.05, p_adjust=True,
                      sym_weight=0.25, df_guess=None, nIter=200, tol=5, sep=':'):
    '''
    Estimate complex abundance by iteratively solving a system of quadratic
    equations. The system of equations is set up based on the expected random
    count of each PLA product.

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
        the 1-sided t-test sample mean > mean_cutoff is kept as 0.
        Default is 1.
    p_cutoff : float
        The alpha level to decide if the 1-sided t-test is sinificant.
        Default is 0.05.
    p_adjust : boolean
        Whether to perform FDR correction for the one-sided t-test.
        Default is True.
    sym_weight : float (0 <= sym_weight <= 1).
        The weight factor used to enforce symmetry condition. 0 means no enforcement.
        Default is 0.25.
    df_guess : pandas data frame
        First guesses of true complex abundance (must be the same shape as data).
        If None (the default), use 0 as the first guess.
    nIter : int
        Max number of iterations to perform.
    tol : float
        If the change in solution between current and last iteration is below
        this value, convergence is reached.
        Default is 5.
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A data frame with the same shape as df, containing estimated complex abundance

    '''

    # The code runs faster if working with numpy array than pandas data frame

    # Convert input data frame into numpy array
    pla = data.to_numpy(copy=True)
    # Convert non_complexes list to sets
    non_complex = set(non_complex)

    # Set of PLA products
    pla_product_set = set(data.index)

    # Get a list of probe A and B targets
    probeA = np.array([s.split(sep)[0] for s in data.index])
    probeB = np.array([s.split(sep)[1] for s in data.index])


    # Initialize a numpy array to store estimated complex amount
    if df_guess is None:
        complex_out = np.zeros(pla.shape)
    else:
        complex_out = df_guess.to_numpy()

    # Iteration
    loop_num = 0
    max_change = tol + 1
    while (loop_num < nIter) and (max_change > tol):

        # Dict to store the one-sided t-test p-values
        tp_all = {}

        # PLA product count minus previous iteration's complex count
        temp_pla = pla - complex_out

        # Calculate the sum of probe A and B
        temp_pla_probeA = {}
        for i in set(probeA):
            temp_pla_probeA[i] = temp_pla[probeA==i,:].sum(axis=0)
        temp_pla_probeB = {}
        for i in set(probeB):
            temp_pla_probeB[i] = temp_pla[probeB==i,:].sum(axis=0)
        temp_pla_sum = temp_pla.sum(axis=0)

        # First pass: get all the p-values
        for i in range(data.shape[0]):

            # if this PLA product is not detected in any cells, skip
            if np.sum(pla[i,:]) == 0:
                continue

            # target of probe A and B
            temp_complex = data.index[i]
            temp_probeA, temp_probeB = temp_complex.split(sep)

            # Apply the constraints
            if (temp_complex in non_complex) or (temp_probeA in non_complex) or (temp_probeB in non_complex):
                continue

            temp_expected = temp_pla_probeA[temp_probeA]*temp_pla_probeB[temp_probeB]/temp_pla_sum
            temp_diff = pla[i,:] - temp_expected

            # Check to see if the estimated abundance passes the mean_cutoff
            # Ha: sample mean > mean_cutoff
            tval, tp = stats.ttest_1samp(temp_diff, mean_cutoff)
            if (tval > 0):
                tp_all[data.index[i]] = tp/2
            else:
                tp_all[data.index[i]] = 1-tp/2

        # Convert p-values dictionary to series
        tp_all = pd.Series(tp_all)
        # Multiple comparison correction
        if p_adjust:
            _, tp_adj, _,_ = multipletests(tp_all.to_numpy(), alpha=p_cutoff, method='fdr_bh')
            tp_adj = pd.Series(tp_adj, index=tp_all.index)
        else:
            tp_adj = tp_all

        # Array to store the change in the complex estimates
        temp_change = np.zeros(pla.shape) + tol + 1
        # Second pass: calculate protein complex
        for i in range(data.shape[0]):
            # if this PLA product is not detected in any cell, skip
            if np.sum(pla[i,:]) == 0:
                temp_change[i,:] = 0
                continue

            # target of probe A and B
            temp_complex = data.index[i]
            temp_probeA, temp_probeB = temp_complex.split(sep)

            # Apply the constraints
            if (temp_complex in non_complex) or (temp_probeA in non_complex) or (temp_probeB in non_complex):
                temp_change[i,:] = 0
                continue

            # Check to see if the estimated abundance passes the mean_cutoff
            # Ha: sample mean > mean_cutoff
            if (tp_adj[data.index[i]] <= p_cutoff):
                temp_expected = temp_pla_probeA[temp_probeA]*temp_pla_probeB[temp_probeB]/temp_pla_sum
                temp_diff = pla[i,:] - temp_expected

            elif (f"{temp_probeB}{sep}{temp_probeA}" in pla_product_set):
                # check for symmetry
                temp_symmetry = complex_out[data.index==f"{temp_probeB}{sep}{temp_probeA}",:]
                if np.mean(temp_symmetry) > mean_cutoff:
                    temp_diff = sym_weight*temp_symmetry
                else:
                    temp_change[i,:] = 0
                    continue
            else:
                temp_change[i,:] = 0
                continue

            # Force negative values to be zero <---- should be done after t-test
            temp_diff[temp_diff < 0] = 0

            # Check if observed is 0 but estimated is non 0, then force the estimated to be 0
            # This should only be done after t-test
            temp_diff[(temp_diff > 0) & (pla[i,:] == 0)] = 0

            # Store changes in the solutions/estimates
            temp_change[i,:] = temp_diff - complex_out[i,:]

            # Store the new solutions/estimates
            complex_out[i,:] = temp_diff

        # Round the adjustment amount
        complex_out = np.round(complex_out)
        # Save the maximum change in the solution for convergence check
        max_change = abs(temp_change).max()

        loop_num += 1

    print(f"estimateComplexes done: Loop number {loop_num}, tolerance {max_change:.2f}")
    return pd.DataFrame(data=complex_out, index=data.index, columns=data.columns)
