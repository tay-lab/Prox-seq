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
my_colors = sns.color_palette('Dark2')
#*****

# Import data
dge = pd.read_csv("Jurkat_Raji_plate_count_matrix.txt.gz",
                  index_col=0, sep="\t")

# Keep cells with at least 100 UMIs
dge = dge.loc[:,dge.sum(axis=0)>=100]

#%% Plot UMI counts and detected PLA products
dge_temp = dge.loc[:,dge.columns.str.contains("separate")].copy()
dge_temp = dge_temp.loc[dge_temp.sum(axis=1)>0,:] # remove empty rows
# Make a dataframe of UMI count and detected PLA products
dge_qc = pd.DataFrame(index=dge_temp.columns)
dge_qc['nPLA'] = (dge_temp>0).sum(axis=0)
dge_qc['nUMI'] = dge_temp.sum(axis=0)

# Plot QC
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6,3))
np.random.seed(2020) # for jitter position
sns.violinplot(y='nPLA', data=dge_qc, cut=0, inner=None, ax=ax[0], saturation=1, color='silver')
sns.stripplot(y='nPLA', data=dge_qc, color='black', jitter=0.3, size=4, ax=ax[0], clip_on=False)
ax[0].set_ylim(-3,90)
ax[0].set_yticks(range(0,101,25))
ax[0].set_ylabel("Number of PLA products")
ax[0].set_title("Detected PLA products")
ax[0].set_xlabel("")
ax[0].set_xticks([0])
ax[0].set_xticklabels(["Plate"])
sns.violinplot(y='nUMI', data=dge_qc, cut=0, inner=None, ax=ax[1], saturation=1, color='silver')
sns.stripplot(y='nUMI', data=dge_qc, color='black', jitter=0.3, size=4, ax=ax[1], clip_on=False)
ax[1].set_ylim(-300,1e4)
ax[1].set_yticks(range(0,10001,2000))
ax[1].set_yticklabels([f"{int(s):,}" for s in ax[1].get_yticks()])
ax[1].set_ylabel("Number of UMIs")
ax[1].set_title("Detected PLA UMIs")
ax[1].set_xlabel("")
ax[1].set_xticks([0])
ax[1].set_xticklabels(["Plate"])
sns.despine(fig=fig)
fig.tight_layout(w_pad=2)

#%% Extract T and B cells
# T cells
dgeT = dge.loc[:,dge.columns.str.contains("Tcells-separate")]
# Remove empty rows
dgeT = dgeT.loc[dgeT.sum(axis=1)>0,:]

# B cells
dgeB = dge.loc[:,dge.columns.str.contains("Bcells-separate")]
# Remove empty rows
dgeB = dgeB.loc[dgeB.sum(axis=1)>0,:]


# Calculate expected
dgeT_expected = PF.calculateExpected(dgeT)

# Plot observed vs expected for 2 example PLA products
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5.7,2.8))
for counter,i in enumerate(["CD3:CD3","CD28:CD28"]):
    ax[counter].scatter(dgeT_expected.loc[i,:], dgeT.loc[i,:],
                        c='black', s=12, clip_on=False)
    _,xmax = ax[counter].get_xlim()
    _,ymax = ax[counter].get_ylim()
    ax_max = max(xmax,ymax)
    ax[counter].set_xlim(-ax_max/20, ax_max)
    ax[counter].set_ylim(-ax_max/20, ax_max)
    ax[counter].plot([0,ax_max], [0,ax_max], c='red', lw=2, clip_on=True)
    ax[counter].set_title(i)
    ax[counter].set_xlabel("Expected random PLA count")
    ax[counter].set_ylabel("Observed PLA count")
    ax[counter].set_aspect('equal')
ax[0].set_xticks([0,200,400,500])
ax[0].set_yticks([0,200,400,500])
ax[1].set_xticks(range(0,151,50))
sns.despine(fig=fig)
fig.tight_layout(w_pad=2)
fig.savefig(myDir+"Tcell_expected_2example.svg",
            bbox_inches="tight", pad_inches=0) # Figure 2b
# temp = dgeT_expected.loc[["CD3:CD3","CD28:CD28"],:]
# temp = dgeT.loc[["CD3:CD3","CD28:CD28"],:]

# Estimate dimer abundance
dgeT_dimers = PF.estimateComplexes(dgeT, nIter=200, non_complex=[], mean_cutoff=1,
                                   tol=1, sym_weight=0.25, p_adjust=True)
dgeT_dimers_avg = dgeT_dimers.mean(axis=1)

# Calculate fraction of PLA count that is estimated as dimer
dgeT_dimers_fraction = dgeT_dimers/dgeT
dgeT_dimers_fraction[dgeT_dimers_fraction.isna()] = 0 # replace nan with 0

#%% Plot changes in expected/observed count with different iterations
def my_set_xyticks(ax, arr):
    # Set same x ticks and y ticks
    ax.set_xticks(arr)
    ax.set_yticks(arr)

fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(11.7,5.8))
for col, i in enumerate([0,1,2,3]):
    # CD3:CD28
    temp = PF.estimateComplexes(dgeT, nIter=i, non_complex=[], mean_cutoff=1,
                                tol=1, sym_weight=0.25, p_adjust=True)
    temp2 = dgeT - temp
    temp2_expected = PF.calculateExpected(temp2)
    ax[0,col].scatter(temp2_expected.loc["CD3:CD28",:], temp2.loc["CD3:CD28",:],
                        c="k", s=12, clip_on=False)
    ax[0,col].plot([0,150], [0,150], c="red", lw=2, clip_on=True)
    ax[0,col].set_title(f"Iteration {i}")
    my_set_xyticks(ax[0,col], np.arange(0,151,50))

    # HLADR:PDL1
    temp = PF.estimateComplexes(dgeB, nIter=i, non_complex=[], mean_cutoff=1,
                                tol=1, sym_weight=0.25, p_adjust=True)
    temp2 = dgeB - temp
    temp2_expected = PF.calculateExpected(temp2)
    ax[1,col].scatter(temp2_expected.loc["HLADR:PDL1",:], temp2.loc["HLADR:PDL1",:],
                        c="k", s=12, clip_on=False)
    ax[1,col].plot([0,100], [0,100], c="red", lw=2, clip_on=True)
    ax[1,col].set_title(f"Iteration {i}")

    my_set_xyticks(ax[1,col], np.arange(0,101,25))
for i in ax.flatten():
    i.set_xlabel("Expected random PLA count")
    i.set_ylabel("Observed PLA count")
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.9, h_pad=2)
fig.savefig("TBcell_expected_example_iterations.svg",
            bbox_inches="tight", pad_inches=0) # Figure 2d, e

#%% Plot heatmap
# Calculate average dimer abundance
dgeT_dimers_avg = pd.DataFrame(dgeT_dimers_avg, columns=["UMI"])
dgeT_dimers_avg["AB1"] = [s.split(":")[0] for s in dgeT_dimers_avg.index]
dgeT_dimers_avg["AB2"] = [s.split(":")[1] for s in dgeT_dimers_avg.index]
dgeT_dimers_avg_square = dgeT_dimers_avg.pivot(index="AB1", columns="AB2", values="UMI")

# Calculate average fraction of dimer abundance
dgeT_dimers_fraction_avg = pd.DataFrame(dgeT_dimers_fraction.mean(axis=1, skipna=True), columns=["fraction"])
dgeT_dimers_fraction_avg["AB1"] = [s.split(":")[0] for s in dgeT_dimers_fraction_avg.index]
dgeT_dimers_fraction_avg["AB2"] = [s.split(":")[1] for s in dgeT_dimers_fraction_avg.index]
dgeT_dimers_fraction_avg_square = dgeT_dimers_fraction_avg.pivot(index="AB1", columns="AB2", values="fraction")

# Reorder rows and columns
my_order = ["CD3","CD28","PD1","CD147","HLADR","ICAM1","PDL1","B7"]
dgeT_dimers_avg_square2 = dgeT_dimers_avg_square.loc[my_order,my_order].copy()
dgeT_dimers_fraction_avg_square2 = dgeT_dimers_fraction_avg_square.loc[my_order,my_order].copy()
# Fill nan with 0
dgeT_dimers_fraction_avg_square2[dgeT_dimers_fraction_avg_square2.isna()] = 0

# Plot heatmap of complex fraction
fig, ax = plt.subplots(figsize=(4.2,3.7))
sns.heatmap(dgeT_dimers_fraction_avg_square2,
            square=True, linewidths=1, ax=ax,
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True),
            vmin=0, vmax=1, cbar_kws={"shrink":0.7})
ax.tick_params(axis="both", length=0)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("Probe B target")
ax.set_ylabel("Probe A target")
fig.tight_layout()
fig.savefig("Tcell_dimer_fraction_heatmap_full.svg",
            bbox_inches="tight", pad_inches=0) # Figure 2f

#%% Plot observed vs expected for T cell markers
dgeT_expected = PF.calculateExpected(dgeT)
AB = ["CD147", "CD28", "CD3", "PD1"]
fig, ax = plt.subplots(nrows=len(AB), ncols=len(AB), figsize=(7.5,7.5))
# temp = {} # for saving source data
for counter_i, i in enumerate(AB):
    for counter_j, j in enumerate(AB):
        ax[counter_i, counter_j].scatter(dgeT_expected.loc[f"{i}:{j}",:], dgeT.loc[f"{i}:{j}",:],
                                         color='k', s=7.5)
        # temp[f"{i}:{j}"] = pd.DataFrame({"expected":dgeT_expected.loc[f"{i}:{j}",:],
        #                                  "observed":dgeT.loc[f"{i}:{j}",:]})
        _, xmax = ax[counter_i, counter_j].get_xlim()
        ax[counter_i, counter_j].plot([0, xmax],[0, xmax], c='red', lw=2)
        if counter_i == len(AB)-1:
            ax[counter_i, counter_j].set_xlabel("Expected\nrandom count")
        if counter_j == 0:
            ax[counter_i, counter_j].set_ylabel("Observed count")
        ax[counter_i, counter_j].set_title(f"{i}:{j}")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("Tcell_expected_Tmarkers.png", bbox_inches="tight", pad_inches=0, dpi=600) # Supplementary figure 4


#%%
# =============================================================================
# # B cells
# =============================================================================
# Calculate random expected count
dgeB_expected = PF.calculateExpected(dgeB)

# Plot 2 example PLA products
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5.7,2.8))
for counter,i in enumerate(["HLADR:HLADR","PDL1:PDL1"]):
    ax[counter].scatter(dgeB_expected.loc[i,:], dgeB.loc[i,:],
                        s=12, c='black', clip_on=False)
    if (counter==0):
        ax_max = 400
    elif (counter==2):
        ax_max = 150
    else:
        _,xmax = ax[counter].get_xlim()
        _,ymax = ax[counter].get_ylim()
        ax_max = max(xmax,ymax)
    ax[counter].set_xlim(-ax_max/20, ax_max)
    ax[counter].set_ylim(-ax_max/20, ax_max)
    ax[counter].plot([0,ax_max], [0,ax_max], c='red', lw=2, clip_on=True)
    # ax[counter].fill_between(x=[0,ax_max], y1=[0,1.5*ax_max], y2=[0,ax_max/1.5],
    #                          facecolor='red', alpha=0.25)
    ax[counter].set_title(i)
    ax[counter].set_xlabel("Expected random PLA count")
    ax[counter].set_ylabel("Observed PLA count")
    ax[counter].set_aspect('equal')
ax[0].set_xticks(range(0,401,100))
ax[0].set_yticks(range(0,401,100))
ax[1].set_xticks(range(0,201,50))
ax[1].set_yticks(range(0,201,50))
sns.despine(fig=fig)
fig.tight_layout(w_pad=2)
fig.savefig(myDir+"Bcell_expected_2example.svg",
            bbox_inches="tight", pad_inches=0) # Figure 2c
# temp = dgeB_expected.loc[["HLADR:HLADR","PDL1:PDL1"],:]
# temp = dgeB.loc[["HLADR:HLADR","PDL1:PDL1"],:]

dgeB_dimers = PF.estimateComplexes(dgeB, nIter=200, non_complex=[], mean_cutoff=1, tol=1, sym_weight=0.25, p_adjust=True)
dgeB_dimers_avg = dgeB_dimers.mean(axis=1)

# Calculate fraction of PLA count that is estimated as dimer
dgeB_dimers_fraction = dgeB_dimers/dgeB
dgeB_dimers_fraction[dgeB_dimers_fraction.isna()] = 0 # replace nan with 0

#%% Plot heatmap of average dimer abundance
dgeB_dimers_avg = pd.DataFrame(dgeB_dimers_avg, columns=["UMI"])
dgeB_dimers_avg["AB1"] = [s.split(":")[0] for s in dgeB_dimers_avg.index]
dgeB_dimers_avg["AB2"] = [s.split(":")[1] for s in dgeB_dimers_avg.index]
dgeB_dimers_avg_square = dgeB_dimers_avg.pivot(index="AB1", columns="AB2", values="UMI")

# Calculate average fraction of dimer abundance
dgeB_dimers_fraction_avg = pd.DataFrame(dgeB_dimers_fraction.mean(axis=1, skipna=True), columns=["fraction"])
dgeB_dimers_fraction_avg["AB1"] = [s.split(":")[0] for s in dgeB_dimers_fraction_avg.index]
dgeB_dimers_fraction_avg["AB2"] = [s.split(":")[1] for s in dgeB_dimers_fraction_avg.index]
# dgeB_dimers_fraction_avg_square = dgeB_dimers_fraction_avg.loc[dgeB_dimers_fraction_avg["AB1"].str.contains("B7|PDL1|ICAM1|CD147|HLADR") & dgeB_dimers_fraction_avg["AB2"].str.contains("B7|PDL1|ICAM1|CD147|HLADR"),:].pivot(index="AB1", columns="AB2", values="fraction")
dgeB_dimers_fraction_avg_square = dgeB_dimers_fraction_avg.pivot(index="AB1", columns="AB2", values="fraction")

# Reorder rows and columns
dgeB_dimers_avg_square2 = dgeB_dimers_avg_square.loc[my_order,my_order].copy()
dgeB_dimers_fraction_avg_square2 = dgeB_dimers_fraction_avg_square.loc[my_order,my_order].copy()
# Fill nan with 0
dgeB_dimers_fraction_avg_square2[dgeB_dimers_fraction_avg_square2.isna()] = 0

# Plot heatmap of complex fraction
fig, ax = plt.subplots(figsize=(4.2,3.7))
sns.heatmap(dgeB_dimers_fraction_avg_square2, ax=ax,
            square=True, linewidths=1,
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True),
            vmin=0, vmax=1, cbar_kws={"shrink":0.7})
ax.tick_params(axis="both", length=0)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("Probe B target")
ax.set_ylabel("Probe A target")
fig.tight_layout()
fig.savefig("Bcell_dimer_fraction_heatmap_full.svg",
            bbox_inches="tight", pad_inches=0) # Figure 2g

# Plot heatmap of complex abundance
fig, ax = plt.subplots(figsize=(5.5,4.5))
sns.heatmap(dgeB_dimers_avg_square2, ax=ax,
            square=True, linewidths=0.75,
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True))
ax.tick_params(axis="both", length=0)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("Probe B target")
ax.set_ylabel("Probe A target")
cbar = ax.collections[0].colorbar
cbar.set_label("Average complex abundance (UMI)")
fig.tight_layout()
fig.savefig("Bcell_dimer_abundance_heatmap_full.svg", bbox_inches="tight")

#%% One-sided Fisher's exact test
def my_fisher_test(pla):
    # Predict complex expression in single-cells using one-sided Fisher's
    # exact test, and return a BH-corrected P-value

    fisher_p = pd.DataFrame(np.nan, index=pla.index, columns=pla.columns)
    probeA = np.array([s.split(':')[0] for s in fisher_p.index])
    probeB = np.array([s.split(':')[1] for s in fisher_p.index])
    for i in fisher_p.index:
        tempA, tempB = i.split(':')
        x00 = (probeA == tempA) & (probeB == tempB)
        x01 = (probeA == tempA) & (probeB != tempB)
        x10 = (probeA != tempA) & (probeB == tempB)
        x11 = (probeA != tempA) & (probeB != tempB)
        temp01 = pla.loc[x01,:].sum(axis=0)
        temp10 = pla.loc[x10,:].sum(axis=0)
        temp11 = pla.loc[x11,:].sum(axis=0)
        for j in fisher_p.columns:
            fisher_p.at[i,j] = stats.fisher_exact(
                [[pla.loc[x00,j], temp01[j]],
                 [temp10[j], temp11[j]]],
                alternative="greater")[1]

    # FDR correction: single cell by single cell
    complex_fisher = pd.DataFrame(np.nan, index=pla.index, columns=pla.columns)
    for i in complex_fisher.columns:
        complex_fisher.loc[:,i] = multipletests(fisher_p.loc[:,i], method="fdr_bh")[1]

    return complex_fisher

T_fisher = my_fisher_test(dgeT)
B_fisher = my_fisher_test(dgeB)

# Expressed complexes mean FDR < 0.05
# Fraction of cells that express a protein complex
T_fisher_frac = pd.DataFrame(
    {"fraction":(T_fisher<0.05).sum(axis=1)/T_fisher.shape[1]})
T_fisher_frac["AB1"] = [s.split(":")[0] for s in T_fisher.index]
T_fisher_frac["AB2"] = [s.split(":")[1] for s in T_fisher.index]
T_fisher_frac_square = T_fisher_frac.pivot(index="AB1", columns="AB2", values="fraction")
T_fisher_frac_square = T_fisher_frac_square.loc[my_order,my_order]
T_fisher_frac_square.fillna(0, inplace=True)

B_fisher_frac = pd.DataFrame(
    {"fraction":(B_fisher<0.05).sum(axis=1)/B_fisher.shape[1]})
B_fisher_frac["AB1"] = [s.split(":")[0] for s in B_fisher.index]
B_fisher_frac["AB2"] = [s.split(":")[1] for s in B_fisher.index]
B_fisher_frac_square = B_fisher_frac.pivot(index="AB1", columns="AB2", values="fraction")
B_fisher_frac_square = B_fisher_frac_square.loc[my_order,my_order]
B_fisher_frac_square.fillna(0, inplace=True)

# Plot heatmap of complex fraction
fig, ax = plt.subplots(ncols=2, figsize=(9,4))
sns.heatmap(T_fisher_frac_square, ax=ax[0],
            square=True, linewidths=1,
            cmap="viridis", vmin=0, vmax=0.8,
            cbar_kws={"label":"Fraction of cells", "shrink":0.7,
                      "ticks":np.arange(0,0.81,0.2)})
sns.heatmap(B_fisher_frac_square, ax=ax[1],
            square=True, linewidths=1,
            cmap="viridis", vmin=0, vmax=0.8,
            cbar_kws={"label":"Fraction of cells", "shrink":0.7,
                      "ticks":np.arange(0,0.81,0.2)})
for i in ax:
    i.tick_params(axis="both", length=0)
    i.xaxis.tick_top()
    i.xaxis.set_label_position('top')
    i.set_xticklabels(i.get_xticklabels(), rotation=90)
    i.set_xlabel("Probe B target")
    i.set_ylabel("Probe A target")
fig.tight_layout(w_pad=1.5)
fig.savefig("fisher_fraction_heatmap.svg",
            bbox_inches="tight", pad_inches=0) # Supplementary figure 6

#%% Plot observed vs expected for B cell markers
dgeB_expected = PF.calculateExpected(dgeB)
AB = ["B7", "CD147", "HLADR", "ICAM1", "PDL1"]
temp = {}
fig, ax = plt.subplots(nrows=len(AB), ncols=len(AB), figsize=(8.5,8.5))
for counter_i, i in enumerate(AB):
    for counter_j, j in enumerate(AB):
        ax[counter_i, counter_j].scatter(dgeB_expected.loc[f"{i}:{j}",:], dgeB.loc[f"{i}:{j}",:],
                                         color='k', s=7.5)
        temp[f"{i}:{j}"] = pd.DataFrame({"expected":dgeB_expected.loc[f"{i}:{j}",:],
                                         "observed":dgeB.loc[f"{i}:{j}",:]})
        _, xmax = ax[counter_i, counter_j].get_xlim()
        ax[counter_i, counter_j].plot([0, xmax],[0, xmax], c='red', lw=2)
        if counter_i == len(AB)-1:
            ax[counter_i, counter_j].set_xlabel("Expected\nrandom count")
        if counter_j == 0:
            ax[counter_i, counter_j].set_ylabel("Observed")
        ax[counter_i, counter_j].set_title(f"{i}:{j}")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("Bcell_expected_Bmarkers.png", bbox_inches="tight", pad_inches=0, dpi=600)

# Use IgG_29:CD147 as the negative control, calculate P values of the homodimers CD3 and CD28
temp2 = dgeB_dimers.loc[["HLADR:HLADR","PDL1:PDL1","IgG_28:ICAM1"],:].copy().T.reset_index()
T_neg_ctrl_complex = pd.melt(temp2, id_vars=['index'],var_name='complex',value_name='UMI')
fig,ax = plt.subplots(1,1, figsize=(3,3.5))
sns.boxplot(x='complex', y='UMI', data=T_neg_ctrl_complex, fliersize=0, saturation=1, ax=ax)
ax.set_xlabel("")
ax.set_ylabel("Estimated complex\nabundance (UMI)")
ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=45, ha='right')
ymin,_ = ax.get_ylim()
ax.set_ylim(ymin, 250)
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("Bcell_neg_ctrl_complex.svg", bbox_inches="tight")

# t-test
print(stats.ttest_ind(dgeB_dimers.loc["HLADR:HLADR",:],dgeB_dimers.loc["IgG_28:ICAM1",:], equal_var=False))
print(stats.ttest_ind(dgeB_dimers.loc["PDL1:PDL1",:],dgeB_dimers.loc["IgG_28:ICAM1",:], equal_var=False))

#%% Compare Prox-seq vs flow for protein
dgeT_protein = PF.calculateProteinAbundance(dgeT)
dgeT_protein_quantile = pd.DataFrame({'median':dgeT_protein.quantile(q=0.5, axis=1),
                                      'q25':dgeT_protein.quantile(q=0.25, axis=1),
                                      'q75':dgeT_protein.quantile(q=0.75, axis=1)})
dgeB_protein = PF.calculateProteinAbundance(dgeB)
dgeB_protein_quantile = pd.DataFrame({'median':dgeB_protein.quantile(q=0.5, axis=1),
                                      'q25':dgeB_protein.quantile(q=0.25, axis=1),
                                      'q75':dgeB_protein.quantile(q=0.75, axis=1)})

# Import summary stats from flow data
dgeTB_protein_flow = pd.read_csv('flow_cytometry.csv', index_col=None)
T_proteins = dgeTB_protein_flow.index[dgeTB_protein_flow["cell_type"]=="Jurkat"].tolist()
B_proteins = dgeTB_protein_flow.index[dgeTB_protein_flow["cell_type"]=="Raji"].tolist()

# Plot median values
fig, ax = plt.subplots(figsize=(3.2,3.2), tight_layout=True)
tempT = dgeTB_protein_flow["protein"][T_proteins].tolist()
tempB = dgeTB_protein_flow["protein"][B_proteins].tolist()
ax.scatter(np.log10(dgeTB_protein_flow.loc[T_proteins,"median"]),
           np.log10(dgeT_protein_quantile.loc[tempT,"median"]),
           s=25)
ax.scatter(np.log10(dgeTB_protein_flow.loc[B_proteins,"median"]),
           np.log10(dgeB_protein_quantile.loc[tempB,"median"]),
           marker='v', s=30)
# Spearman correlation
rho_x = pd.concat((np.log10(dgeTB_protein_flow.loc[T_proteins,"median"]),
                   np.log10(dgeTB_protein_flow.loc[B_proteins,"median"])))
rho_y = pd.concat((np.log10(dgeT_protein_quantile.loc[tempT,"median"]),
                   np.log10(dgeB_protein_quantile.loc[tempB,"median"])))
rho = stats.spearmanr(rho_x,rho_y)
ax.text(0.05,0.90, r"$\mathrm{\rho}$ = "+f"{rho[0]:.2f}", transform=ax.transAxes)
ax.set_xlabel("Flow cytometry\nlog10 (fluorescence intensity)")
ax.set_ylabel("Prox-seq\nlog10 (UMI count)")
ax.set_xlim(2.5,4.5)
ax.set_xticks(np.arange(2.5,4.6,0.5))
ax.set_ylim(1.5,3.7)
ax.legend(["Jurkat","Raji"], loc='lower right', fancybox=False, edgecolor='black')
sns.despine(fig=fig)
fig.savefig("proxseq_vs_flow.svg",
            bbox_inches='tight', pad_inches=0) # Extended Figure 3

stats.spearmanr(dgeTB_protein_flow.loc[T_proteins,"median"],dgeT_protein_quantile.loc[tempT,"median"])
stats.spearmanr(dgeTB_protein_flow.loc[B_proteins,"median"],dgeB_protein_quantile.loc[tempB,"median"])
stats.spearmanr(pd.concat([dgeTB_protein_flow.loc[B_proteins,"median"],dgeTB_protein_flow.loc[B_proteins,"median"]]),
                pd.concat([dgeB_protein_quantile.loc[tempB,"median"],dgeB_protein_quantile.loc[tempB,"median"]]))

#%% Estimate non-specific binding
Tmarkers = ["CD3","CD28","PD1"]
Bmarkers = ["HLADR","ICAM1","PDL1","B7"]
TBmarkers = ["CD147"]

# Combine data frame
dgeT_protein_temp = dgeT_protein.T.copy()
dgeT_protein_temp["type"] = "Jurkat"
dgeB_protein_temp = dgeB_protein.T.copy()
dgeB_protein_temp["type"] = "Raji"
dgeTB_protein = pd.concat((dgeT_protein_temp, dgeB_protein_temp), axis=0)

# Log-transform
dgeTB_protein.iloc[:,:-1] = np.log10(dgeTB_protein.iloc[:,:-1] + 1)

# Plot
my_list = [Tmarkers,Bmarkers,TBmarkers]
p_values = []
fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(11,8.5))
np.random.seed(0)
for row in range(len(my_list)):
    for col in range(4):
        if col >= len(my_list[row]):
            ax[row,col].axis('off')
            continue
        sns.violinplot(data=dgeTB_protein, x="type", y=my_list[row][col], saturation=1, cut=0,
                       inner=None, ax=ax[row,col])
        sns.stripplot(data=dgeTB_protein, x="type", y=my_list[row][col], color='black',
                      jitter=0.25, size=3, ax=ax[row,col])
        p_values.append(stats.ttest_ind(dgeTB_protein.loc[dgeTB_protein["type"]=="Jurkat",my_list[row][col]],
                                        dgeTB_protein.loc[dgeTB_protein["type"]=="Raji",my_list[row][col]],
                                        equal_var=False)[1])
        ax[row,col].set_title(f"{my_list[row][col]} protein")
        ax[row,col].set_xlabel("Cell type")
        ax[row,col].set_ylabel("")
    ax[row,0].set_ylabel("log10(UMI + 1)")
sns.despine(fig=fig)
fig.tight_layout(w_pad=2)
fig.savefig("nonspecific_background.png", bbox_inches='tight',
            pad_inches=0, dpi=600) # Supplementary Figure 3

# Multiple comparison correction
p_values_fdr = multipletests(p_values, method='fdr_bh')

#%% Combine T and B cells, and estimate complex abundance
dgeTB_real = pd.concat([dgeT, dgeB], axis=1)
dgeTB_dimers = PF.estimateComplexes(dgeTB_real, nIter=200, non_complex=[], mean_cutoff=1, tol=1, sym_weight=0.25, p_adjust=True)
dgeTB_dimers_avg = dgeTB_dimers.mean(axis=1)
dgeTB_Tcell = dgeTB_real.columns.str.contains("Tcells") # boolean if cell type is T

# Convert average abundance to square format
dgeTB_dimers_avg = pd.DataFrame(dgeTB_dimers_avg, columns=["UMI"])
dgeTB_dimers_avg["AB1"] = [s.split(":")[0] for s in dgeTB_dimers_avg.index]
dgeTB_dimers_avg["AB2"] = [s.split(":")[1] for s in dgeTB_dimers_avg.index]
dgeTB_dimers_avg_square = dgeTB_dimers_avg.pivot(index="AB1", columns="AB2", values="UMI")
dgeTB_dimers_avg_square2 = dgeTB_dimers_avg_square.loc[my_order,my_order].copy()

# Plot heatmap
fig, ax = plt.subplots(figsize=(5.5,4.5))
sns.heatmap(dgeTB_dimers_avg_square2, ax=ax,
            square=True, linewidths=0.75,
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True), vmin=0, vmax=120)
ax.tick_params(axis="both", length=0)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("Probe B target")
ax.set_ylabel("Probe A target")
cbar = ax.collections[0].colorbar
cbar.set_label("Average complex abundance (UMI)")
fig.tight_layout()
fig.savefig("TBcell_dimer_abundance_heatmap_full.svg", bbox_inches="tight")
