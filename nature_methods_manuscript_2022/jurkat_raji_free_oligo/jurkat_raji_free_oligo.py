# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# =============================================================================
# # Simulation of Prox-seq data
# =============================================================================
# Import libraries
import numpy as np
import math
import random
import pandas as pd

import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

import statsmodels.api as sm # linear regression

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import ProxseqFunctions as PF

#*****
mpl.rcdefaults()
# Set font to be arial
mpl.rc('font', **{'sans-serif':'Arial', 'size':12})
mpl.rcParams['mathtext.rm'] = 'sans' # to have non-italic greek letter, use r'$\mathrm{\alpha}$', does NOT work with f-string
mpl.rcParams['axes.titlesize'] = 12
# Set default tick size
mpl.rcParams['xtick.major.size'] = 5.5
mpl.rcParams['ytick.major.size'] = 5.5
mpl.rcParams['xtick.minor.size'] = 2.5
mpl.rcParams['ytick.minor.size'] = 2.5
# Default legend settings
mpl.rcParams['legend.fancybox'] = False
mpl.rcParams['legend.edgecolor'] = 'k'
#*****
# Seed number
np.random.seed(2019)
random.seed(2019)


#%% Import real data
# New data
dge = pd.read_csv("count_matrix.txt.gz",
                  index_col=0, delimiter="\t")

# Plot total number of UMIs
fig, ax = plt.subplots()
ax.hist(dge.sum(axis=0), bins=100)

# Keep cells with > 100 UMIs and below 25k UMIs
dge = dge.loc[:,(dge.sum(axis=0)>100) & (dge.sum(axis=0)<25000)]

# Keep PLA products detected in at least 5 cells
dge = dge.loc[(dge>0).sum(axis=1)>=5,:]

# Cell-type data
dge_T = dge.loc[:,dge.columns.str.contains("Jurkat")]
dge_B = dge.loc[:,dge.columns.str.contains("Raji")]

# Calculate protein: exclude free oligo counts
dge_T_protein = PF.calculateProteinAbundance(dge_T.loc[~dge_T.index.str.contains("free_oligo"),:])
dge_B_protein = PF.calculateProteinAbundance(dge_B.loc[~dge_B.index.str.contains("free_oligo"),:])

# Calculate protein: include free oligo counts
dge_T_protein2 = PF.calculateProteinAbundance(dge_T)
dge_B_protein2 = PF.calculateProteinAbundance(dge_B)

#%% Compare protein abundance vs free oligo
fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(9.7,8.3))
for col, value in enumerate(["CD3","CD28","PD1"]):
    ax[0,col].scatter(dge_T_protein.loc[value,:],
                      dge_T.loc[f"{value}:free_oligo"]+dge_T.loc[f"free_oligo:{value}"],
                      c='k', s=8)
    r = stats.pearsonr(dge_T_protein.loc[value,:],
                       dge_T.loc[f"{value}:free_oligo"]+dge_T.loc[f"free_oligo:{value}"])
    ax[0,col].set_title(f"Jurkat: {value} protein\n{r[0]:.2f}")
    ax[0,col].set_xlabel("PLA product-based")
ax[0,3].axis('off')
ax[0,0].set_ylabel("Free oligo-based")
for col, value in enumerate(["B7","ICAM1","HLADR","PDL1"]):
    ax[1,col].scatter(dge_B_protein.loc[value,:],
                      dge_B.loc[f"{value}:free_oligo"]+dge_B.loc[f"free_oligo:{value}"],
                      c='k', s=8)
    r = stats.pearsonr(dge_B_protein.loc[value,:],
                       dge_B.loc[f"{value}:free_oligo"]+dge_B.loc[f"free_oligo:{value}"])
    ax[1,col].set_title(f"Raji: {value} protein\n{r[0]:.2f}")
    ax[1,col].set_xlabel("PLA product-based")
ax[1,0].set_xticks(range(0,301,100))
ax[1,1].set_xticks(range(0,20001, 10000))
ax[1,2].set_xticks(range(0,15001,5000))
ax[1,0].set_ylabel("Free oligo-based")
ax[2,0].scatter(dge_T_protein.loc["CD147",:],
                dge_T.loc["CD147:free_oligo"]+dge_T.loc["free_oligo:CD147"],
                c='k', s=8)
r = stats.pearsonr(dge_T_protein.loc["CD147",:],dge_T.loc["CD147:free_oligo"]+dge_T.loc["free_oligo:CD147"])
ax[2,0].set_title(f"Jurkat: CD147 protein\n{r[0]:.2f}")
ax[2,0].set_xlabel("PLA product-based")
ax[2,0].set_yticks(range(0,151,50))
ax[2,0].set_ylabel("Free oligo-based")
ax[2,1].scatter(dge_B_protein.loc["CD147",:],
                dge_B.loc["CD147:free_oligo"]+dge_B.loc["free_oligo:CD147"],
                c='k', s=8)
r = stats.pearsonr(dge_B_protein.loc["CD147",:],dge_B.loc["CD147:free_oligo"]+dge_B.loc["free_oligo:CD147"])
ax[2,1].set_title(f"Raji: CD147 protein\n{r[0]:.2f}")
ax[2,1].set_xlabel("PLA product-based")
ax[2,1].set_xticks(range(0,4001,2000))
ax[2,2].axis('off')
ax[2,3].axis('off')
sns.despine(fig=fig)
fig.tight_layout(h_pad=1.7)
fig.savefig(myDir+"figures_revision/free_oligo_vs_protein.svg",
            bbox_inches='tight', pad_inches=0) # Extended figure 4c-e

#%% Ratio of PLA product-based protein to PLA product with free oligo
T_proteins = ["CD3","CD28","PD1","CD147"]
B_proteins = ["B7","ICAM1","HLADR","PDL1","CD147"]

dge_T_protein_ratio = 100*(dge_T_protein/dge_T_protein2).loc[T_proteins,:]
dge_B_protein_ratio = 100*(dge_B_protein/dge_B_protein2).loc[B_proteins,:]

# Convert to long format for boxplot
dge_T_protein_ratio = pd.melt(dge_T_protein_ratio.T.reset_index(),
                              id_vars=["index"], value_vars=dge_T_protein_ratio.index)
dge_B_protein_ratio = pd.melt(dge_B_protein_ratio.T.reset_index(),
                              id_vars=["index"], value_vars=dge_B_protein_ratio.index)

# Plot
fig, ax = plt.subplots(ncols=2, figsize=(4.9,3.2))
# np.random.seed(0)
sns.boxplot(data=dge_T_protein_ratio, x="variable", y="value",
            width=0.5, ax=ax[0], order=T_proteins,
            saturation=1, color="grey", fliersize=2)
# sns.stripplot(data=dge_T_protein_ratio, x="variable", y="value",
#               jitter=0.2, color="k", ax=ax[0],
#               size=2)
ax[0].set_xlabel("Protein")
ax[0].set_ylabel("Jurkat")
ax[0].set_xticklabels(T_proteins, rotation=45)
ax[0].set_ylim(89.5,100.5)
ax[0].set_yticks(np.arange(90,101,5))
sns.boxplot(data=dge_B_protein_ratio, x="variable", y="value",
            width=0.5, ax=ax[1], order=B_proteins,
            saturation=1, color="grey", fliersize=2)
ax[1].set_xlabel("Protein")
ax[1].set_ylabel("Raji")
ax[1].set_xticklabels(B_proteins, rotation=45)
ax[1].set_ylim(89.5,100.5)
ax[1].set_yticks(np.arange(90,101,5))
fig.suptitle("Fraction of PLA product-based\nprotein counts (%)",
             fontsize=12, y=0.93)
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"figures_revision/protein_fraction.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 4a
