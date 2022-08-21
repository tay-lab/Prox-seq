# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# Import libraries
import numpy as np
import pandas as pd

# import scipy.spatial as spatial
import scipy.stats as stats
from scipy.io import mmread

# import statsmodels.api as sm
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

#%% Import files
pla = pd.read_csv("10x_PLA_count_matrix.txt.gz",
                  index_col=0, delimiter="\t")

# Check number of UMIs
fig, ax = plt.subplots(ncols=2, figsize=(6,2.8))
ax[0].hist(pla.sum(axis=0), bins=100, color='k')
ax[0].set_xlabel("Number of UMIs")
ax[0].axvline(100, c='r')
ax[1].hist(np.log10(pla.sum(axis=0)), bins=100, color='k')
ax[1].axvline(np.log10(100), c='r')
ax[1].set_xlabel("log10(Number of UMIs)")
sns.despine(fig=fig)
fig.tight_layout()

# Filter cells with > 50 UMI and < 10,000 UMI
pla = pla.loc[:,(pla.sum(axis=0)>50) & (pla.sum(axis=0)<10e3)]
# Keep PLA products detected in at least 1 cells
pla = pla.loc[(pla>0).sum(axis=1)>=1,:]

#%% Import cluster identity from R
pla_cluster = pd.read_csv("PLA_cluster.csv",
                          index_col=0)
pla_cluster = pd.Series(pla_cluster["x"])

# Only keep the cell barcodes in pla_cluster
pla = pla.loc[:,pla_cluster.index]

#%% Calculate counts of non-expressed PLA products
# Probes to exclude in calculation, because they are only detected in the
# full probe panel
full_T_only = ["CD3","CD4","CD45RA","LFA1","IgG","ICAM1","B7"]
full_B_only = ["ICAM1","B7","IgG","CD3","CD4","CD45RA","LFA1"]
full_T_only = "|".join(full_T_only)
full_B_only = "|".join(full_B_only)

# Calculate median of non-specific PLA products
pla_median = {
    "Jurkat_H":pla.loc[~pla.index.str.contains(full_T_only) & pla.index.str.contains("PDL1|HLADR"), pla_cluster=="Jurkat_H"].median(axis=1),
    "Jurkat_F":pla.loc[~pla.index.str.contains(full_T_only) & pla.index.str.contains("PDL1|HLADR"),pla_cluster=="Jurkat_F"].median(axis=1),
    "Raji_H":pla.loc[~pla.index.str.contains(full_B_only) & pla.index.str.contains("PD1|CD28"),pla_cluster=="Raji_H"].median(axis=1),
    "Raji_F":pla.loc[~pla.index.str.contains(full_B_only) & pla.index.str.contains("PD1|CD28"),pla_cluster=="Raji_F"].median(axis=1)
}
pla_median_T = pd.DataFrame({s:pla_median[s] for s in ["Jurkat_H","Jurkat_F"]})
pla_median_B = pd.DataFrame({s:pla_median[s] for s in ["Raji_H","Raji_F"]})

# Plot
fig, ax = plt.subplots(figsize=(3,2.6))
ax.plot([0,1], pla_median_T.to_numpy()[:,:2].T, "o-",
        c="k", ms=2, zorder=5)
pla_median_T["x"] = 0
sns.boxplot(data=pla_median_T, y="Jurkat_H", x="x", ax=ax,
            width=0.5, order=range(5), fliersize=0,
            saturation=1, color="silver")
pla_median_T["x"] = 1
sns.boxplot(data=pla_median_T, y="Jurkat_F", x="x", ax=ax,
            width=0.5, order=range(5), fliersize=0,
            saturation=1, color="silver")
ax.plot([3,4], pla_median_B.to_numpy()[:,:2].T, "o-",
        c="k", ms=2, zorder=5)
pla_median_B["x"] = 3
sns.boxplot(data=pla_median_B, y="Raji_H", x="x", ax=ax,
            width=0.5, order=range(5), fliersize=0,
            saturation=1, color="silver")
pla_median_B["x"] = 4
sns.boxplot(data=pla_median_B, y="Raji_F", x="x", ax=ax,
            width=0.5, order=range(5), fliersize=0,
            saturation=1, color="silver")
ax.set_xlabel("")
ax.set_ylabel("Median count (UMIs)")
ax.set_title("Level of non-specific\nPLA products")
ax.set_xticks([0,1,3,4])
ax.set_xticklabels(["H","F","H","F"])
ax.set_yticks([0,10,20,25])
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("non_specific_pla.svg", bbox_inches='tight', pad_inches=0) # Extended figure 5d

#%% Calculate protein complexes of probes detected in both half and full panels
complex_dict = {
    "Jurkat_H":PF.estimateComplexes(pla.loc[:,pla_cluster=="Jurkat_H"]),
    "Jurkat_F":PF.estimateComplexes(pla.loc[:,pla_cluster=="Jurkat_F"]),
}

# Log-transform, then calculate average
T_markers = ["CD3","CD28","PD1","CD147"]
fig, ax = plt.subplots(ncols=2, figsize=(5.7,3))
for counter, value in enumerate(["Jurkat_H","Jurkat_F"]):
    temp = pd.DataFrame(
        {"mean":np.log10(complex_dict[value]+1).mean(axis=1)})
    temp["probeA"] = [i.split(":")[0] for i in temp.index]
    temp["probeB"] = [i.split(":")[1] for i in temp.index]
    temp = temp.pivot(index="probeA", columns="probeB", values="mean")
    sns.heatmap(temp.loc[T_markers, T_markers], ax=ax[counter],
                vmin=0, vmax=2, square=True, linewidths=0.5,
                cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True),
                cbar_kws={'label':'log10(UMI + 1)','ticks':np.arange(0,2.1,1),'fraction':0.045})
    ax[counter].set_title(value)
    del temp
for i in ax:
    i.set_xlabel("Probe B target")
    i.set_ylabel("Probe A target")
    i.tick_params(axis="both", length=0)
    i.set_xticklabels(i.get_xticklabels(), fontsize=11, rotation=90)
    i.set_yticklabels(i.get_yticklabels(), fontsize=11, rotation=0)
fig.tight_layout()
fig.savefig("Jurkat_complex_probe_panels.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 5g

# Violinplot
fig, ax = plt.subplots(ncols=2, figsize=(4.8,2.7))
np.random.seed(0)
for counter, value in enumerate(["CD28:CD28","CD28:PD1"]):
    temp1 = pd.DataFrame(
        {value:np.log10(complex_dict["Jurkat_H"].loc[value,:]/complex_dict["Jurkat_H"].sum(axis=0)*10e3+1),
         "sample":"Jurkat_H"})
    temp2 = pd.DataFrame(
        {value:np.log10(complex_dict["Jurkat_F"].loc[value,:]/complex_dict["Jurkat_F"].sum(axis=0)*10e3+1),
         "sample":"Jurkat_F"})
    temp = pd.concat([temp1,temp2])

    sns.violinplot(x="sample", y=value, data=temp, ax=ax[counter],
                   cut=0, saturation=1, scale="width", color="darkgrey", inner=None)
    sns.stripplot(x="sample", y=value, data=temp, ax=ax[counter],
                  color="k", jitter=0.25, size=1.5)

    # Label
    ax[counter].set_xlabel("")
    ax[counter].set_ylabel("Normalized level")
    ax[counter].set_title(f"Protein complex {value}")

    # Axes limit
    ax[counter].set_ylim(-0.2,4.5)

    # LogFC
    print(np.log2(temp2[value].mean()/temp1[value].mean()))

    del temp, temp1, temp2
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.8)
fig.savefig("Jurkat_complex_probe_panels_violin.svg",
            bbox_inches='tight', pad_inches=0) # Extended data figure 5h
