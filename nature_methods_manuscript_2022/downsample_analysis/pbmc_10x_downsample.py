# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# Import libraries
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import multipletests

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mplticker
# from mpl_toolkits.mplot3d import Axes3D

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

rna_cluster_id = pd.read_csv("RNA_cluster_id.csv",
                             index_col=0, header=0, delimiter=',')

#%% RNA subsampling

# Median mumber of reads, genes and UMIs
rna = pd.DataFrame(0.0, columns=['nReads','nUMIs','nGenes'],
                   index=['0.1','0.2','0.4','0.6','0.8','1.0'])

# Enter medinum number of UMIs and genes from web_summary.html
rna.at['0.1','nReads'] = 4615
rna.at['0.1','nUMIs'] = 1140
rna.at['0.1','nGenes'] = 565
rna.at['0.2','nReads'] = 8368
rna.at['0.2','nUMIs'] = 1925
rna.at['0.2','nGenes'] = 847
rna.at['0.4','nReads'] = 15659
rna.at['0.4','nUMIs'] = 3002
rna.at['0.4','nGenes'] = 1184
rna.at['0.6','nReads'] = 22995
rna.at['0.6','nUMIs'] = 3722
rna.at['0.6','nGenes'] = 1384
rna.at['0.8','nReads'] = 30334
rna.at['0.8','nUMIs'] = 4217
rna.at['0.8','nGenes'] = 1516
rna.at['1.0','nReads'] = 38314
rna.at['1.0','nUMIs'] = 4623
rna.at['1.0','nGenes'] = 1620

# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.3,2.85))
ax[0].plot(rna['nReads']/1000, rna['nGenes'], 'o-', c='k', ms=4)
ax[0].set_ylabel("Median genes per cell")
ax[0].set_yticks(np.arange(0,1501,500))
ax[0].set_title("Detected genes", pad=12)
ax[1].plot(rna['nReads']/1000, rna['nUMIs'], 'o-', c='k', ms=4)
ax[1].set_ylabel("Median UMIs per cell")
ax[1].set_yticks([0,2000,4000,5000])
ax[1].set_title("Detected mRNA UMIs", pad=12)
for i in ax:
    i.set_xticks(np.arange(0,41,10))
    i.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.6)
fig.savefig("rna_downsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 9a

#%% Check cell type calling accuracy
# Reference singleR labels (import from R)
singleR_ref = pd.read_csv("singleR_labels.csv", index_col=0)

# Accuracy
rna['singleR_accuracy'] = 0.0
rna.at['1.0','singleR_accuracy'] = 1 # reference

for i in rna.index:
    if i == "1.0":
        continue

    j = i.replace(".","_")
    singleR_test = pd.read_csv(f"singleR_labels_{j}.csv", index_col=0)
    temp = singleR_ref.merge(singleR_test, how='inner', left_index=True, right_index=True)

    # Matching labels
    rna.at[i,'singleR_accuracy'] = (temp.iloc[:,0]==temp.iloc[:,1]).sum()/temp.shape[0]

# Plot accuracy
fig, ax = plt.subplots(figsize=(2.7,2.7))
ax.plot(rna['nReads']/1000, rna['singleR_accuracy'],
        'o-', c='k', ms=4, clip_on=False)
ax.set_xticks(np.arange(0,41,10))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.set_ylim(0,1)
ax.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
ax.set_ylabel("Cell type annotation\naccuracy")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("rna_downsample_singleR.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 9d

#%% PLA subsampling
pla_readcounts = {'0_1':9646567, '0_2':19294481, '0_4':38587633,
                  '0_6':57881154, '0_8':77192233, '1_0':96493480}

# Median mumber of reads, genes and UMIs
pla = pd.DataFrame(0.0, columns=['nReads','nUMIs','nPLAs','nUMIs_protein','nProteins'],
                   index=['0_1','0_2','0_4','0_6','0_8','1_0'])


for i in pla.index:
    if i == '1_0':
        dge = pd.read_csv("10x_PLA_count_matrix.txt.gz",
                          index_col=0, delimiter="\t")
    else:
        dge = pd.read_csv(f"count_downsample_{i}/10x_PLA_count_matrix.txt.gz",
                          index_col=0, delimiter="\t")

    pla.at[i,'nReads'] = pla_readcounts[i]/dge.shape[1]
    pla.at[i,'nUMIs'] = dge.sum(axis=0).median()
    pla.at[i,'nPLAs'] = (dge>0).sum(axis=0).median()

    temp = PF.calculateProteinAbundance(dge)
    pla.at[i,'nUMIs_protein'] = temp.sum(axis=0).median()
    pla.at[i,'nProteins'] = (temp>0).sum(axis=0).median()

# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.3,2.8))
ax[0].plot(pla['nReads']/1000, pla['nPLAs'], 'o-', c='k', ms=4)
ax[0].set_ylabel("Median PLA products\nper cell")
ax[0].set_yticks(np.arange(0,101,25))
ax[0].set_title("Detected PLA products", pad=12)
ax[1].plot(pla['nReads']/1000, pla['nUMIs'], 'o-', c='k', ms=4)
ax[1].set_ylabel("Median UMIs per cell")
ax[1].set_yticks([0,200,400,500])
ax[1].set_title("Detected PLA product UMIs", pad=12)
for i in ax:
    i.set_xticks(np.arange(0,11,2.5))
    i.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.6)
fig.savefig("pla_downsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 9b

# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.3,2.8))
ax[0].plot(pla['nReads']/1000, pla['nProteins'], 'o-', c='k', ms=4)
ax[0].set_ylabel("Median proteins\nper cell")
ax[0].set_yticks(np.arange(0,31,10))
ax[0].set_title("Detected proteins", pad=12)
ax[1].plot(pla['nReads']/1000, pla['nUMIs_protein'], 'o-', c='k', ms=4)
ax[1].set_ylabel("Median UMIs per cell")
ax[1].set_yticks(np.arange(0,1001,250))
ax[1].set_title("Detected protein UMIs", pad=12)
for i in ax:
    i.set_xticks(np.arange(0,11,2.5))
    i.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.6)
fig.savefig("protein_downsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 9c

# Number of protein complexes
# Subsample number of reads: rna_clusters 0 and 3
nComplexes_reads = pd.DataFrame(0, index=['0_1','0_2','0_4','0_6','0_8','1_0'],
                                columns=[0,3])
for i in nComplexes_reads.index:
    if i == "1_0":
        temp = pd.read_csv("10x_PLA_count_matrix.txt.gz",
                           index_col=0, delimiter="\t")
    else:
        temp = pd.read_csv(f"count_downsample_{i}/10x_PLA_count_matrix.txt.gz",
                           index_col=0, delimiter="\t")

    # Loop through the rna clusters
    for j in nComplexes_reads.columns:
        # Find the barcodes that are in each cluster
        temp_bc = set(rna_cluster_id.index[rna_cluster_id['x']==j]).intersection(set(temp.columns))
        # Protein complex prediction
        temp_complex = PF.estimateComplexes(temp.loc[:,temp_bc],
                                            nIter=200, tol=5, mean_cutoff=1, p_cutoff=0.05, p_adjust=True)

        nComplexes_reads.at[i,j] = np.sum(temp_complex.sum(axis=1)>0)

    del temp, temp_bc, temp_complex

# Plot
fig, ax = plt.subplots(figsize=(2.7,2.7))
ax.plot(pla["nReads"]/1e3, nComplexes_reads[0], "o-",
        c="#1b9e77", ms=4)
ax.plot(pla["nReads"]/1e3, nComplexes_reads[3], "o-",
        c="#e7298a", ms=4)
ax.set_xlabel("Mean reads per cell\n("+r"$\times$ 1,000)")
ax.set_ylabel("Number of protein\ncomplexes per cluster")
ax.set_xticks(np.arange(0,11,2.5))
ax.set_yticks(np.arange(0,26,5))
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("protein_complex_downsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 9e
