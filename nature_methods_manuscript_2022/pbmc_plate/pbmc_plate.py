# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# Import libraries
import numpy as np
import pandas as pd

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
# Set default seaborn palette
sns.set_palette('Dark2')
#*****
dark2_colors = sns.color_palette()

#%% Import live PLA data
dge = pd.read_csv("/PLA_dge.txt.gz", index_col=0, delimiter="\t")

# Keep cells with at least 10 UMIs
dge = dge.loc[:,dge.sum(axis=0)>=10]

# Remove empty rows
dge = dge.loc[dge.sum(axis=1)>0,:]

#%% Process cells
# Calculate protein abundance
dge_protein = PF.calculateProteinAbundance(dge)

# Get CD4 and CD8 cells
fig, ax = plt.subplots(figsize=(4,4), tight_layout=True)
ax.scatter(dge_protein.loc["CD4",:].to_numpy(), dge_protein.loc["CD8",:].to_numpy())
ax.set_xlabel("CD4")
ax.set_ylabel("CD8")

# CD4 cells: CD3 > 100, CD8 count < 100 and CD4 count > 100
dge_id = np.array(["Other" for _ in range(dge.shape[1])])
dge_id[(dge_protein.loc["CD8",:]<100) & (dge_protein.loc["CD4",:]>100) & (dge_protein.loc["CD3",:]>100)] = "CD4"
dge_id[(dge_protein.loc["CD4",:]<100) & (dge_protein.loc["CD8",:]>100) & (dge_protein.loc["CD3",:]>100)] = "CD8"

# New id: CD9:CD9 high cells = more CD9:CD9 count than CD9 heterodimer PLA count
# CD4_CD9 = CD9:CD9 PLA product high CD4 T cells
# Color by CD9:CD9 PLA product high/low
dge_id2 = np.copy(dge_id).astype("<U25")
temp = (dge.loc["CD9:CD9",:] > dge.loc[CD9hetero,:].sum(axis=0)) & (dge_id != 'Other')
dge_id2[temp] = np.char.add(dge_id2[temp], "_CD9")

# Plot CD8 T cells only
fig, ax = plt.subplots(figsize=(3.1,3.2), tight_layout=True)
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD8_CD9"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD8_CD9"],
           edgecolors="none", alpha=0.7, s=20)
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD8"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD8"],
           edgecolors="none", alpha=0.7, s=20)
ax.plot([0,1200],[0,1200], color='red')
ax.set_xlabel("PLA product count:\nCD9 and non-CD9")
ax.set_ylabel("PLA product count:\nCD9:CD9")
ax.set_title("CD8 T cells")
ax.set_yticks([0,500,1000,1500,2000])
sns.despine(fig=fig)
fig.savefig(myDir+"data_analysis/figures/CD9_homo_hetero_CD8Tcells.svg",
            bbox_inches="tight", pad_inches=0) # Figure 3e

fig, ax = plt.subplots(figsize=(4.6,3), tight_layout=True)
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD4_CD9"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD4_CD9"],
           edgecolors="none", alpha=0.7, s=20, label="CD9:CD9 PLA\nproduct high")
ax.scatter(dge.loc[CD9hetero,dge_id2=="CD4"].sum(axis=0),
           dge.loc["CD9:CD9",dge_id2=="CD4"],
           edgecolors="none", alpha=0.7, s=20, label="CD9:CD9 PLA\n product low")
ax.plot([0,1000],[0,1000], color='red')
ax.set_xlabel("PLA product count:\nCD9 and non-CD9")
ax.set_ylabel("PLA product count:\nCD9:CD9")
ax.set_title("CD4 T cells")
ax.legend(frameon=False,
          loc=6, bbox_to_anchor=(1.02, 0.5))
sns.despine(fig=fig)
fig.savefig(myDir+"data_analysis/figures/CD9_homo_hetero_CD4Tcells.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Extended figure 6a

#%% Separate T cells out to make the figure for protein complex prediction schematic
dgeT = dge.loc[:,dge_protein.loc["CD3",:]>100].copy()

# Calculate expected
dgeT_expected = PF.calculateExpected(dgeT)

# Plot difference in observed and expected
fig, ax = plt.subplots(figsize=(3,3.1))
ax.hist(dgeT.loc["CD28:CD28",:]-dgeT_expected.loc["CD28:CD28",:],
        color='silver', edgecolor='k')
ax.axvline(x=1, c='red')
ax.set_xticks(range(0,81,20))
ax.set_xlabel("Difference between observed\nand expected random PLA count\n of CD28:CD28 (UMI count)")
ax.set_ylabel("Number of cells")
ax.set_title("Iteration 0 (before adjustment)")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"data_analysis/figures/CD28.CD28_complex_example.svg",
            bbox_inches="tight", pad_inches=0) # Supplementary figure 11c

#%% Complex estimate
dge_complex = {}
for i in ["CD4","CD4_CD9","CD8","CD8_CD9"]:
    dge_complex[i] = PF.estimateComplexes(dge.loc[:,dge_id2==i], tol=1, nIter=200, mean_cutoff=1, p_adjust=True)
dge_complex_mean = {}
for i in dge_complex.keys():
    dge_complex_mean[i] = dge_complex[i].mean(axis=1)
dge_complex_mean = pd.DataFrame(dge_complex_mean)

# Plot heatmap showing mean log10(UMI+1) of PLA products, and complex detected
probes_to_plot = ["CD127","CD28","CD3","CD38","CD4","CD8","CD9","HLADR"]
dge_temp1 = dge.loc[:,(dge_id2=="CD4")|(dge_id2=="CD4_CD9")|(dge_id2=="CD8")|(dge_id2=="CD8_CD9")]
dge_temp1 = pd.DataFrame(np.log10(dge_temp1 + 1).mean(axis=1))
dge_temp1.loc[:,"probeA"] = [s.split(':')[0] for s in dge_temp1.index]
dge_temp1.loc[:,"probeB"] = [s.split(':')[1] for s in dge_temp1.index]
dge_temp1 = dge_temp1.pivot(index="probeA", columns="probeB", values=0)
dge_temp1.fillna(value=0, inplace=True)
# dge_temp2 = dge.loc[:,(dge_id2=="CD4")|(dge_id2=="CD4_CD9")|(dge_id2=="CD8")|(dge_id2=="CD8_CD9")]
dge_temp2 = pd.concat(dge_complex.values(), axis=1)
dge_temp2 = pd.DataFrame(np.log10(dge_temp2 + 1).mean(axis=1))
dge_temp2.loc[:,"probeA"] = [s.split(':')[0] for s in dge_temp2.index]
dge_temp2.loc[:,"probeB"] = [s.split(':')[1] for s in dge_temp2.index]
dge_temp2 = dge_temp2.pivot(index="probeA", columns="probeB", values=0)
dge_temp2.fillna(value=0, inplace=True)
fig, ax = plt.subplots(ncols=2, figsize=(8.7,4))
sns.heatmap(dge_temp1, ax=ax[0], linewidths=0.6,
            yticklabels=False, xticklabels=False,
            cbar_kws={'label':'log10(UMI + 1)','ticks':np.arange(0,1.51,0.3)},
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True), vmax=1.5)
ax[0].set_title("Average count of all PLA products")
ax[0].set_xlabel("Probe B target")
ax[0].set_ylabel("Probe A target")
sns.heatmap(dge_temp2.loc[probes_to_plot,probes_to_plot], ax=ax[1], linewidths=1,
            # yticklabels=False, xticklabels=False,
            cbar_kws={'label':'log10(UMI + 1)','ticks':np.arange(0,1.51,0.3)},
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True), vmax=1.5)
ax[1].set_xlabel("Probe B target")
ax[1].set_ylabel("Probe A target")
ax[1].set_title("Average complex count")
ax[1].tick_params(axis='both', length=0)
fig.tight_layout(w_pad=1.75)
fig.savefig(myDir+"data_analysis/figures/complex_heatmap.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Figure 3b

# Heatmap of 2 example single cells
my_cell = "P2-C06" # CD4 T cell
dge_temp1 = pd.DataFrame(np.log10(dge.loc[:,my_cell]+1))
dge_temp1.loc[:,"probeA"] = [s.split(':')[0] for s in dge_temp1.index]
dge_temp1.loc[:,"probeB"] = [s.split(':')[1] for s in dge_temp1.index]
dge_temp1 = dge_temp1.pivot(index="probeA", columns="probeB", values=my_cell)
dge_temp1.fillna(0, inplace=True)
dge_temp2 = pd.DataFrame(np.log10(dge_complex["CD4"].loc[:,my_cell] + 1))
dge_temp2.loc[:,"probeA"] = [s.split(':')[0] for s in dge_temp2.index]
dge_temp2.loc[:,"probeB"] = [s.split(':')[1] for s in dge_temp2.index]
dge_temp2 = dge_temp2.pivot(index="probeA", columns="probeB", values=my_cell)
fig, ax = plt.subplots(ncols=2, figsize=(8.7,4.3)) # (7,3.8)
sns.heatmap(dge_temp1, ax=ax[0], linewidths=0.4, vmax=2.6,
            yticklabels=False, xticklabels=False,
            cbar_kws={'label':'log10(UMI + 1)','ticks':[0,1,2,2.5]},
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True))
ax[0].set_title("Count of all PLA products")
ax[0].set_xlabel("Probe B target")
ax[0].set_ylabel("Probe A target")
sns.heatmap(dge_temp2.loc[probes_to_plot,probes_to_plot], ax=ax[1], linewidths=1, vmax=2.6,
            # square=True,
            # yticklabels=False, xticklabels=False,
            cbar_kws={'label':'log10(UMI + 1)','ticks':[0,1,2,2.5]},
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True))
ax[1].set_xlabel("Probe B target")
ax[1].set_ylabel("Probe A target")
ax[1].set_title("Complex count")
ax[1].tick_params(axis='both', length=0)
fig.suptitle("Single CD4 T cell", y=0.91, fontsize=12)
fig.tight_layout(w_pad=1.5)
fig.savefig(myDir+"data_analysis/figures/CD4Tcell_heatmap.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Figure 3c

# my_cell = "P1-D09" # CD8 T cell
my_cell = "P4-G01" # CD8 T cell
dge_temp1 = pd.DataFrame(np.log10(dge.loc[:,my_cell]+1))
dge_temp1.loc[:,"probeA"] = [s.split(':')[0] for s in dge_temp1.index]
dge_temp1.loc[:,"probeB"] = [s.split(':')[1] for s in dge_temp1.index]
dge_temp1 = dge_temp1.pivot(index="probeA", columns="probeB", values=my_cell)
dge_temp1.fillna(0, inplace=True)
dge_temp2 = pd.DataFrame(np.log10(dge_complex["CD8"].loc[:,my_cell] + 1))
dge_temp2.loc[:,"probeA"] = [s.split(':')[0] for s in dge_temp2.index]
dge_temp2.loc[:,"probeB"] = [s.split(':')[1] for s in dge_temp2.index]
dge_temp2 = dge_temp2.pivot(index="probeA", columns="probeB", values=my_cell)
fig, ax = plt.subplots(ncols=2, figsize=(8.7,4.3)) # figsize=(7,3.8)
sns.heatmap(dge_temp1, ax=ax[0], linewidths=0.4, vmax=2.6,
            yticklabels=False, xticklabels=False,
            cbar_kws={'label':'log10(UMI + 1)','ticks':[0,1,2,2.5]},
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True))
ax[0].set_title("Count of all PLA products")
ax[0].set_xlabel("Probe B target")
ax[0].set_ylabel("Probe A target")
sns.heatmap(dge_temp2.loc[probes_to_plot,probes_to_plot], ax=ax[1], linewidths=1, vmax=2.6,
            # square=True,
            # yticklabels=False, xticklabels=False,
            cbar_kws={'label':'log10(UMI + 1)','ticks':[0,1,2,2.5]},
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True))
ax[1].set_xlabel("Probe B target")
ax[1].set_ylabel("Probe A target")
ax[1].set_title("Complex count")
ax[1].tick_params(axis='both', length=0)
fig.suptitle("Single CD8 T cell", y=0.91, fontsize=12)
fig.tight_layout(w_pad=1.5)
fig.savefig(myDir+"data_analysis/figures/CD8Tcell_heatmap.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Figure 3d


#%% Violin plot protein complex
# Prepare a dataframe for Violin plot
dge_violin = pd.DataFrame(columns=["complex","cell_type","count"])
for i in ["CD3:CD3","CD3:CD4","CD3:CD8","CD4:CD4","CD8:CD8",
          "CD9:CD4","CD9:CD8","CD8:CD9","CD28:CD28","CD9:CD9","CD4:CD38",
          "CD9:CD28","CD28:CD9","CD8:CD28"]:
    for j in ["CD4","CD4_CD9","CD8","CD8_CD9"]:
        temp = pd.DataFrame(columns=["complex","cell_type","count"])
        temp.loc[:,"count"] = dge_complex[j].loc[i,:].to_numpy()
        temp.loc[:,"complex"] = i
        temp.loc[:,"cell_type"] = j
        dge_violin = dge_violin.append(temp, ignore_index=True)

# CD9:CD8 complex
print(stats.mannwhitneyu(dge_complex["CD8"].loc["CD9:CD8",:],
                         dge_complex["CD8_CD9"].loc["CD9:CD8",:],
                         alternative='greater'))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3,3.2))
np.random.seed(1)
sns.violinplot(x="cell_type", y="count",
               data=dge_violin.loc[((dge_violin["cell_type"]=="CD8")|(dge_violin["cell_type"]=="CD8_CD9"))&(dge_violin["complex"]=="CD9:CD8"),:],
               order=["CD8_CD9","CD8"], ax=ax, saturation=1, cut=0, inner=None)
sns.stripplot(x="cell_type", y="count",
              data=dge_violin.loc[((dge_violin["cell_type"]=="CD8")|(dge_violin["cell_type"]=="CD8_CD9"))&(dge_violin["complex"]=="CD9:CD8"),:],
              order=["CD8_CD9","CD8"], ax=ax, color="black", jitter=0.25, size=3)
ax.set_ylabel("Complex count\nper cell (UMI)")
ax.set_xlabel("CD9:CD9 PLA\nproduct level")
ax.set_xticks([0,1])
ax.set_xticklabels(["High","Low"])
ax.set_yticks(list(range(0,151,50)))
ax.set_title("CD9:CD8 complex")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"data_analysis/figures/CD9.CD8_complex_CD8Tcells.svg",
            bbox_inches="tight", pad_inches=0) # Figure 3f

# Plot CD8 and CD9 protein abundance
temp1 = {} # CD8_CD9 subpopulation
temp1["CD3"] = dge_protein.loc["CD3",dge_id2=="CD8_CD9"]
temp1["CD8"] = dge_protein.loc["CD8",dge_id2=="CD8_CD9"]
temp1["CD9"] = dge_protein.loc["CD9",dge_id2=="CD8_CD9"]
temp1['type'] = "High"
temp2 = {} # CD8 subpopulation
temp2["CD3"] = dge_protein.loc["CD3",dge_id2=="CD8"]
temp2["CD8"] = dge_protein.loc["CD8",dge_id2=="CD8"]
temp2["CD9"] = dge_protein.loc["CD9",dge_id2=="CD8"]
temp2['type'] = "Low"
dge_violin = pd.concat((pd.DataFrame(temp1), pd.DataFrame(temp2)))

fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(3.15,8.5))
np.random.seed(1)
for col, i in enumerate(["CD3","CD8","CD9"]):
    sns.violinplot(x="type", y=i, data=dge_violin,
                   ax=ax[col], saturation=1, cut=0, inner=None)
    sns.stripplot(x="type", y=i, data=dge_violin,
                  ax=ax[col], color="black", jitter=0.25, size=3)
    ax[col].set_xlabel("")
    ax[col].set_title(f"{i} protein")
    ax[col].set_ylabel("Protein count\nper cell (UMI)")
ax[-1].set_xlabel("CD9:CD9 PLA\nproduct level")
ax[0].set_yticks(range(0,751,250))
ax[2].set_yticks(range(0,4001,1000))
fig.tight_layout(w_pad=1.7)
sns.despine(fig=fig)
fig.savefig(myDir+"data_analysis/figures/CD8Tcell_subpop_3proteins.svg",
            bbox_inches="tight", pad_inches=0) # Figure 3g

###### CD4 T cells
# CD9:CD8 complex
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3,3))
np.random.seed(1)
sns.violinplot(x="cell_type", y="count",
               data=dge_violin.loc[((dge_violin["cell_type"]=="CD4")|(dge_violin["cell_type"]=="CD4_CD9"))&(dge_violin["complex"]=="CD9:CD4"),:],
               order=["CD4_CD9","CD4"], ax=ax, saturation=1, cut=0, inner=None)
sns.stripplot(x="cell_type", y="count",
              data=dge_violin.loc[((dge_violin["cell_type"]=="CD4")|(dge_violin["cell_type"]=="CD4_CD9"))&(dge_violin["complex"]=="CD9:CD4"),:],
              order=["CD4_CD9","CD4"], ax=ax, color="black", jitter=0.25, size=3)
ax.set_ylabel("Complex count\nper cell (UMI)")
ax.set_xlabel("CD9:CD9 PLA\nproduct level")
ax.set_xticks([0,1])
ax.set_xticklabels(["High","Low"])
# ax.set_yticks(list(range(0,151,50)))
ax.set_ylim(-5,100)
ax.set_title("CD9:CD4 complex")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"data_analysis/figures/CD9.CD4_complex_CD4Tcells.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 6b

# Plot CD4 and CD9 protein abundance
temp1 = {} # CD4_CD9 subpopulation
temp1["CD3"] = dge_protein.loc["CD3",dge_id2=="CD4_CD9"]
temp1["CD4"] = dge_protein.loc["CD4",dge_id2=="CD4_CD9"]
temp1["CD9"] = dge_protein.loc["CD9",dge_id2=="CD4_CD9"]
temp1['type'] = "High"
temp2 = {} # CD4 subpopulation
temp2["CD3"] = dge_protein.loc["CD3",dge_id2=="CD4"]
temp2["CD4"] = dge_protein.loc["CD4",dge_id2=="CD4"]
temp2["CD9"] = dge_protein.loc["CD9",dge_id2=="CD4"]
temp2['type'] = "Low"
dge_violin = pd.concat((pd.DataFrame(temp1), pd.DataFrame(temp2)))

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8,3))
np.random.seed(1)
for col, i in enumerate(["CD3","CD4","CD9"]):
    sns.violinplot(x="type", y=i, data=dge_violin,
                   ax=ax[col], saturation=1, cut=0, inner=None)
    sns.stripplot(x="type", y=i, data=dge_violin,
                  ax=ax[col], color="black", jitter=0.25, size=3)
    ax[col].set_xlabel("CD9:CD9 PLA\nproduct level")
    ax[col].set_ylabel("")
    ax[col].set_title(f"{i} protein")
ax[0].set_ylabel("Protein count\nper cell (UMI)")
ax[0].set_yticks(range(0,1501,500))
ax[1].set_yticks(range(0,1501,500))
ax[2].set_yticks(range(0,3001,1000))
fig.tight_layout(w_pad=1.7)
sns.despine(fig=fig)
fig.savefig(myDir+"data_analysis/figures/CD4Tcell_subpop_3proteins.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 6c
