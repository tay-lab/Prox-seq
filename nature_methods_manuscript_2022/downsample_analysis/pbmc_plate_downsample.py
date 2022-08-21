# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# Import libraries
import numpy as np
import pandas as pd

# import scipy.spatial as spatial
import scipy.stats as stats

from statsmodels.stats.multitest import multipletests

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker

import ProxseqFunctions as PF

#*****
mpl.rcdefaults()
# Set font to be arial
mpl.rc('font', **{'sans-serif':'Arial', 'size':12})
mpl.rcParams['mathtext.rm'] = 'sans' # to have non-italic greek letter, use r'$\mathrm{\alpha}$', does NOT work with f-string
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['hatch.linewidth'] = 1.5  # previous svg hatch linewidth
# Set default tick size
mpl.rcParams['xtick.major.size'] = 5.5
mpl.rcParams['ytick.major.size'] = 5.5
mpl.rcParams['xtick.minor.size'] = 2.5
mpl.rcParams['ytick.minor.size'] = 2.5
loc = mplticker.MultipleLocator(base=0.5)
# Default legend settings
mpl.rcParams['legend.fancybox'] = False
mpl.rcParams['legend.edgecolor'] = 'k'
# Set default seaborn palette
sns.set_palette('Dark2')
#*****
my_colors = sns.color_palette()

#%% RNA subsampling
my_subsampling = ['0.005','0.01','0.05','0.1','0.25','0.5','0.75','1.0']
pla = pd.DataFrame(0.0, columns=['nReads','nUMIs','nPLAs','nUMIs_protein','nProteins'],
                   index=my_subsampling)
for i in my_subsampling:
    if i == '1.0':
        df = pd.read_csv("SmartPLAy_dge.txt.gz",
                         index_col=0, delimiter="\t")
        pla.at[i,'nReads'] = 37_802_725/df.shape[1]
    else:
        j = i.replace('.','_')
        df = pd.read_csv(f"downsample_{j}.csv",
                         index_col=0) # import number of reads per cell for each downsample fraction
        pla.at[i,'nReads'] = df.sum()/df.shape[0]

        df = pd.read_csv(f"downsample_{j}/PLA_count_matrix.txt.gz",
                         index_col=0, delimiter="\t")

    pla.at[i,'nUMIs'] = df.sum(axis=0).median()
    pla.at[i,'nPLAs'] = (df>0).sum(axis=0).median()

    # proteins
    temp = PF.calculateProteinAbundance(df)
    pla.at[i,'nUMIs_protein'] = temp.sum(axis=0).median()
    pla.at[i,'nProteins'] = (temp>0).sum(axis=0).median()

    del df
    
# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.3,2.8))
ax[0].plot(pla['nReads']/1000, pla['nPLAs'], 'o-', c='k', ms=4)
ax[0].set_ylabel("Median PLA products\nper cell")
ax[0].set_yticks(np.arange(0,61,20))
ax[0].set_title("Detected PLA products", pad=12)
ax[1].plot(pla['nReads']/1000, pla['nUMIs'], 'o-', c='k', ms=4)
ax[1].set_ylabel("Median UMIs per cell")
ax[1].set_yticks(np.arange(0,1001,250))
ax[1].set_title("Detected PLA product UMIs", pad=12)
for i in ax:
    i.axvline(10, c='r', ls='--', lw=1)
    i.set_xticks(np.arange(0,101,25))
    i.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.6)
fig.savefig("pla_downsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended figure 9f

# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.3,2.8))
ax[0].plot(pla['nReads']/1000, pla['nProteins'], 'o-', c='k', ms=4, clip_on=False)
ax[0].set_ylabel("Median proteins\nper cell")
ax[0].set_yticks(np.arange(0,21,5))
ax[0].set_title("Detected proteins", pad=12)
ax[1].plot(pla['nReads']/1000, pla['nUMIs_protein'], 'o-', c='k', ms=4)
ax[1].set_ylabel("Median UMIs per cell")
ax[1].set_yticks(np.arange(0,2001,500))
ax[1].set_title("Detected protein UMIs", pad=12)
for i in ax:
    i.axvline(10, c='r', ls='--', lw=1)
    i.set_xticks(np.arange(0,101,25))
    i.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.6)
fig.savefig("protein_downsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended figure 9g
