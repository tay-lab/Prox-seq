# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# Import libraries
import numpy as np
import pandas as pd

import scipy.spatial as spatial
import scipy.stats as stats
# from scipy.io import mmread
from scipy.cluster.hierarchy import linkage

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

#%% Import PLA data
dge = pd.read_csv("10x_PLA_count_matrix.txt.gz",
                  index_col=0, delimiter="\t")

# Import cell barcodes from R
good_cb = pd.read_csv("RNA_cell_barcodes.csv",
                      index_col=0, header=0, delimiter=',')

# Import RNA cluster ids from R
rna_cluster_id = pd.read_csv("RNA_cluster_id.csv",
                             index_col=0, header=0, delimiter=',')

# Keep good cell barcodes
dge = dge.loc[:,good_cb['x']]

# Discard rows with only zero
dge = dge.loc[dge.sum(axis=1)>0,:]

# Protein
dge_protein = PF.calculateProteinAbundance(dge)

# Expected
dge_expected = PF.calculateExpected(dge)

# Probe abundance
dge_probeAB = PF.calculateProbeAbundance(dge)

#%% Find CD8:CD9 interactions
# Find CD8 T cells
dge_id = np.array(["Other" for _ in range(dge.shape[1])])
dge_id[(dge_protein.loc["CD4",:]<100) & (dge_protein.loc["CD8",:]>100) & (dge_protein.loc["CD3",:]>100)] = "CD8"
dge_id[(dge_protein.loc["CD8",:]<100) & (dge_protein.loc["CD4",:]>100) & (dge_protein.loc["CD3",:]>100)] = "CD4"

# Find the two subpopulations based on CD9 homodimer
CD9hetero = dge.index.str.contains("CD9:|:CD9$") & (dge.index != "CD9:CD9")

# New id: CD8_CD9 = CD8 T cells with high count of CD9CD9 PLA product
# Color by CD9:CD9 high/low
dge_id2 = np.copy(dge_id).astype("<U25")
temp = (dge.loc["CD9:CD9",:] > dge.loc[CD9hetero,:].sum(axis=0)) & (dge_id != 'Other')
dge_id2[temp] = np.char.add(dge_id2[temp], "_CD9")
dge_id2_export = pd.Series(dge_id2, index=dge.columns)
dge_id2_export.to_csv("CD9_homodimer_cell_types.csv")

# Plot CD8 vs CD4
fig, ax = plt.subplots(figsize=(3,3))
ax.scatter(dge_protein.loc["CD4",dge_id=="CD8"], dge_protein.loc["CD8",dge_id=="CD8"])
ax.scatter(dge_protein.loc["CD4",dge_id=="CD4"], dge_protein.loc["CD8",dge_id=="CD4"])
ax.set_xlabel("CD4")
ax.set_ylabel("CD8")
ax.legend(["CD8","CD4"])

# Plot only CD8 T cells
fig, ax = plt.subplots(ncols=1, figsize=(3,3))
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD8_CD9"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD8_CD9"],
           alpha=0.7, edgecolors="none", s=20)
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD8"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD8"],
           alpha=0.7, edgecolors="none", s=20)
ax.plot([0,3000],[0,3000], color='red')
ax.set_xticks(range(0,3001,1000))
ax.set_ylabel("PLA product count:\nCD9:CD9")
ax.set_xlabel("PLA product count:\nCD9 and non-CD9")
ax.set_title("CD8 T cells")
fig.tight_layout()
sns.despine(fig=fig)

# Plot only CD8 T cells: zoom in
fig, ax = plt.subplots(ncols=1, figsize=(3,3))
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD8_CD9"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD8_CD9"],
           alpha=0.7, edgecolors="none", s=20)
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD8"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD8"],
           alpha=0.7, edgecolors="none", s=20)
ax.plot([0,3000],[0,3000], color='red')
ax.set_xticks(range(0,1001,250))
ax.set_xlim(-50,1000)
ax.set_ylim(-50,1000)
ax.set_ylabel("PLA product count:\nCD9:CD9")
ax.set_xlabel("PLA product count:\nCD9 and non-CD9")
ax.set_title("CD8 T cells")
fig.tight_layout()
sns.despine(fig=fig)
fig.savefig(myDir+"figures/CD9_homo_hetero_zoomed.png", dpi=600,
            bbox_inches='tight', pad_inches=0) # Figure 4e

# Complex count
dge_complex = {}
for i in ["CD8_CD9","CD8","CD4_CD9","CD4"]:
    dge_complex[i] = PF.estimateComplexes(dge.loc[:,dge_id2==i], sym_weight=0.25,
                                          nIter=200, tol=5, mean_cutoff=1, p_cutoff=0.05, p_adjust=True)

# Export CD9:CD8 complex count
# Add CD9:CD8 and CD8:CD9 together
temp1 = pd.DataFrame(dge_complex["CD8_CD9"].loc["CD9:CD8",:]+dge_complex["CD8_CD9"].loc["CD8:CD9",:],
                     columns=["CD9:CD8"])
temp1['type'] = "CD8_CD9"
temp2 = pd.DataFrame(dge_complex["CD8"].loc["CD9:CD8",:]+dge_complex["CD8"].loc["CD8:CD9",:],
                     columns=["CD9:CD8"])
temp2['type'] = "CD8"
pd.concat((temp1,temp2)).to_csv(myDir+"CD9.CD8_complex_count.csv")

# Prepare data
temp1 = {'CD3:CD8':dge_complex['CD8_CD9'].loc["CD3:CD8",:],
         'CD9:CD8':dge_complex['CD8_CD9'].loc["CD9:CD8",:],
         'CD8:CD8':dge_complex['CD8_CD9'].loc["CD8:CD8",:],
         'type':"High"}
temp2 = {'CD3:CD8':dge_complex['CD8'].loc["CD3:CD8",:],
         'CD9:CD8':dge_complex['CD8'].loc["CD9:CD8",:],
         'CD8:CD8':dge_complex['CD8'].loc["CD8:CD8",:],
         'type':"Low"}
for i in temp1:
    if i == 'type':
        continue
    print(i, temp1[i].mean())
for i in temp2:
    if i == 'type':
        continue
    print(i, temp2[i].mean())
print(stats.mannwhitneyu(temp1['CD9:CD8'],temp2['CD9:CD8'], alternative='less'))
temp_violin = pd.concat([pd.DataFrame(temp1), pd.DataFrame(temp2)])



# Violin plot: CD9:CD8 complex only
print(stats.mannwhitneyu(dge_complex["CD8"].loc["CD9:CD8",:],
                         dge_complex["CD8_CD9"].loc["CD9:CD8",:],
                         alternative='greater'))
fig, ax = plt.subplots(ncols=1, figsize=(2.9,3))
np.random.seed(3)
sns.violinplot(data=temp_violin, x='type', y='CD9:CD8', order=["High","Low"], ax=ax, cut=0, inner=None, saturation=1)
sns.stripplot(data=temp_violin, x='type', y='CD9:CD8', order=["High","Low"], ax=ax, color="black",
              jitter=0.3, size=2.5)
ax.set_title("CD9:CD8 complex")
ax.set_yticks(range(0,201,100))
ax.set_xlabel("CD9:CD9 PLA\nproduct level")
ax.set_ylabel("Complex count (UMI)")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"figures/CD9.CD8_complex.png", dpi=600,
            bbox_inches='tight', pad_inches=0) # Figure 4f


#%% Protein complexes of RNA clusters
dge_complex_rna_cluster = {}
for i in range(8):
    dge_complex_rna_cluster[i] = PF.estimateComplexes(dge.loc[:,rna_cluster_id.index[rna_cluster_id['x']==i]],
                                                      nIter=200, tol=5, mean_cutoff=1, p_cutoff=0.05, p_adjust=True)

# Export list of predicted protein complexes per cluster
dge_complex_export = {}
all_complex = set([])
for i in dge_complex_rna_cluster:
    dge_complex_export[i] = pd.Series(dge_complex_rna_cluster[i].index[dge_complex_rna_cluster[i].sum(axis=1)>0], name=i)
    all_complex.update(dge_complex_export[i])
dge_complex_export = pd.DataFrame(dge_complex_export)
# dge_complex_export.to_csv(myDir+"RNA_cluster_predicted_complexes.csv")

# Number of unique complexes across all clusters
all_unique_complex = set([])
for i in all_complex:
    tempA, tempB = i.split(':')
    if (i not in all_unique_complex) and (f"{tempB}:{tempA}" not in all_unique_complex):
        all_unique_complex.add(i)


# Get list of proteins expressed in each cluster
dge_protein_rna_cluster = pd.DataFrame(0, index=dge_protein.index, columns=range(8))
for i in range(8):
    dge_protein_rna_cluster.loc[:,i] = dge_protein.loc[:,rna_cluster_id.index[rna_cluster_id['x']==i]].median(axis=1)

# Expressed proteins: median at least 5
dge_protein_rna_cluster = {}
for i in range(8):
    dge_protein_rna_cluster[i] = pd.Series(dge_protein.index[dge_protein.loc[:,rna_cluster_id.index[rna_cluster_id['x']==i]].median(axis=1) >= 5], name=i)
dge_protein_rna_cluster = pd.DataFrame(dge_protein_rna_cluster)
# dge_protein_rna_cluster.to_csv(myDir+"RNA_cluster_expressed_proteins.csv")

# Plot number of protein complexes per cluster
fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(9.7,5))
# temp = {} # for exporting to Source data
np.random.seed(1)
for i in dge_complex_rna_cluster:
    # print(dge_complex_rna_cluster[i].shape[1])
    row = i // 4
    col = i % 4
    # temp[i] = (dge_complex_rna_cluster[i]>0).sum(axis=0) # for exporting to Source data
    sns.violinplot(y=(dge_complex_rna_cluster[i]>0).sum(axis=0), ax=ax[row,col],
                   saturation=1, cut=0, inner=None, color='silver')
    sns.stripplot(y=(dge_complex_rna_cluster[i]>0).sum(axis=0), ax=ax[row,col],
                  jitter=0.3, color='k', size=2.5)
    ax[row,col].set_xticks([0])
    ax[row,col].set_xticklabels([f"Cluster {i}"])
    ax[row,col].axhline(y=sum(dge_complex_rna_cluster[i].sum(axis=1)>0), c='red')
    ax[row,col].set_ylim(-1,32)
ax[0,0].set_ylabel("Number of predicted\ncomplexes per cell")
ax[1,0].set_ylabel("Number of predicted\ncomplexes per cell")
fig.tight_layout(w_pad=1.7, h_pad=1.7)
sns.despine(fig=fig)
fig.savefig(myDir+"figures/RNA_cluster_nComplexes.png", bbox_inches='tight',
            pad_inches=0, dpi=600) # Extended figure 7b

# Subsample number of cells
n_cells = [10,50,100,500,1000,2000,3000]
temp_dict = {}
for i in range(8):
    # value of the dict is a dataframe: number of complexes detected
    temp_dict[i] = pd.DataFrame(np.nan, columns=['nCells','nComplexes'], index=range(len(n_cells)))
    np.random.seed(0)
    for j in temp_dict[i].index:
        # Number of cells
        temp_dict[i].at[j,'nCells'] = min(n_cells[j],sum(rna_cluster_id['x']==i))

        # Skip if the number of cells doesn't change
        if j > 0:
            if temp_dict[i].at[j,'nCells'] == temp_dict[i].at[j-1,'nCells']:
                continue

        temp_index = np.random.choice(rna_cluster_id.index[rna_cluster_id['x']==i],
                                      size=int(temp_dict[i].at[j,'nCells']), replace=False)
        temp = PF.estimateComplexes(dge.loc[:,temp_index],
                   nIter=200, tol=5, mean_cutoff=1, p_cutoff=0.05, p_adjust=True)
        temp_dict[i].at[j,'nComplexes'] = sum(temp.sum(axis=1)>0)

fig, ax = plt.subplots(figsize=(4.9,3.3))
for i in temp_dict:
    ax.plot(temp_dict[i]['nCells'], temp_dict[i]['nComplexes'],
            label=i, clip_on=False)
ax.set_xlabel("Number of cells\n per cluster")
ax.set_ylabel("Number of predicted\nprotein complexes")
ax.set_xticks(np.arange(0,3001,1000))
# ax.set_ylim(0,35)
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5),
          title="Cluster", frameon=False, ncol=2)
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"figures/RNA_cluster_nComplexes_cells_subsample.svg",
            bbox_inches='tight', pad_inches=0) # Extended figure 7c
