# -*- coding: utf-8 -*-
"""
@author: Hoang Van Phan
"""

# Import libraries
import numpy as np
import pandas as pd

# import scipy.spatial as spatial
import scipy.stats as stats

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

#%% Import PLA data
dge = pd.read_csv("LPS_PAM_count_matrix.txt.gz", index_col=0, delimiter="\t")

# Remove cells below 10 or above 10,000 UMIs
dge = dge.loc[:,(dge.sum(axis=0)>10) & (dge.sum(axis=0)<10000)]

# Remove zero rows
dge = dge.loc[dge.sum(axis=1)>0,:]

#%%
# Ligand
cnd = np.repeat("Control", dge.shape[1])
for i in ["LPS","PAM","Both"]:
    cnd[dge.columns.str.contains(i)] = i
# Time point
tpt = np.repeat("Control", dge.shape[1])
for i in ["5m","2h","12h"]:
    tpt[dge.columns.str.contains(f"-{i}")] = i

# log transform
dge_log = np.log10(dge+1)

#%% Plot proteins
dge_protein = PF.calculateProteinAbundance(dge)

# Calculate median protein abundance of each ligand and time point
dge_temp = pd.DataFrame(0, index=dge_protein.index, columns=range(10))

# Control
dge_temp.rename(columns={0:"Control"}, inplace=True)
dge_temp.loc[:,"Control"] = np.log10(dge_protein.loc[:,cnd=="Control"]+1).mean(axis=1)

# Treatment
for counter1, i in enumerate(["LPS","PAM","Both"]):
    for counter2, j in enumerate(["5m","2h","12h"]):
        dge_temp.rename(columns={(counter1*3+counter2+1):f"{i}-{j}"}, inplace=True)
        dge_temp.loc[:,f"{i}-{j}"] = np.log10(dge_protein.loc[:,(cnd==i)&(tpt==j)]+1).mean(axis=1)

# Calculate z-score
dge_temp_z = dge_temp.apply(stats.zscore, axis=1, result_type='broadcast')

# Heatmap
fig, ax = plt.subplots(figsize=(3.5,4))
sns.heatmap(dge_temp_z, cmap="RdYlGn_r", center=0, vmin=-2, vmax=2,
            linewidths=1, ax=ax, cbar_kws={'label':'z-score\n(Protein level)',
                                           'ticks':[-2,-1,0,1,2]})
ax.tick_params(axis='both', length=0)
ax.set_ylabel("Protein")
fig.savefig("protein_zscore.svg",
            bbox_inches="tight", pad_inches=0) # Figure 5c

# Distribution of proteins
dge_protein_corr = np.log10(dge_protein+1).T.corr(method='pearson')
sns.clustermap(dge_protein_corr)

#%% Plot PLA products
dge_temp = pd.DataFrame(0, index=dge.index, columns=range(10))

# Control
dge_temp.rename(columns={0:"Control"}, inplace=True)
dge_temp.loc[:,"Control"] = np.log10(dge.loc[:,cnd=="Control"]+1).mean(axis=1)

# Treatment
for counter1, i in enumerate(["LPS","PAM","Both"]):
    for counter2, j in enumerate(["5m","2h","12h"]):
        dge_temp.rename(columns={(counter1*3+counter2+1):f"{i}-{j}"}, inplace=True)
        dge_temp.loc[:,f"{i}-{j}"] = np.log10(dge.loc[:,(cnd==i)&(tpt==j)]+1).mean(axis=1)

# Calculate z-score
dge_temp_z = dge_temp.apply(stats.zscore, axis=1, result_type='broadcast')

# Heatmap
f = sns.clustermap(dge_temp_z, col_cluster=False, center=0, cmap='RdBu_r',
                   method='complete', metric='euclidean', figsize=(4.8,10),
                   cbar_pos=(1,0.4,0.03,0.2), cbar_kws={'label':"z-score\n(PLA product level)"})
f.ax_heatmap.tick_params(axis='x', length=0)
f.ax_heatmap.set_ylabel('')
f.ax_row_dendrogram.set_visible(False)
f.savefig("pla_zscore.svg",
          bbox_inches="tight", pad_inches=0) # Figure 5b

#%% Changes in PLA partners
probeA = "TLR2"
dge_temp = pd.DataFrame(0, index=dge.index[[s.split(':')[0]==probeA for s in dge.index]], columns=range(10))

# Control
dge_temp.rename(columns={0:"Control"}, inplace=True)
dge_temp.loc[:,"Control"] = np.log10(dge.loc[dge_temp.index,cnd=="Control"]+1).mean(axis=1)

# Treatment
for counter1, i in enumerate(["LPS","PAM","Both"]):
    for counter2, j in enumerate(["5m","2h","12h"]):
        dge_temp.rename(columns={(counter1*3+counter2+1):f"{i}-{j}"}, inplace=True)
        dge_temp.loc[:,f"{i}-{j}"] = np.log10(dge.loc[[s.split(':')[0]==probeA for s in dge.index],(cnd==i)&(tpt==j)]+1).mean(axis=1)

# Calculate z-score
dge_temp_z = dge_temp.apply(stats.zscore, axis=1, result_type='broadcast')

# Heatmap
fig, ax = plt.subplots(figsize=(3.6,4))
sns.heatmap(dge_temp_z, ax=ax, center=0,
            cmap='RdBu_r', linewidths=1, vmin=-2, vmax=2,
            cbar_kws={'ticks':[-2,-1,0,1,2],'label':"z-score\n(PLA product level)"})
ax.tick_params(axis='both', length=0)
ax.set_ylabel("TLR2-related PLA products")
fig.savefig("TLR2_partners_heatmap.svg",
            bbox_inches="tight", pad_inches=0) # Figure 5d


#%% Mean dynamics plot: PLA products
dge_temp = {}
for i in ["LPS","PAM","Both"]:
    dge_temp[i] = pd.DataFrame(np.nan, index=dge.index, columns=["Control","5m","2h","12h"])

    for j in dge_temp[i].columns:
        if j == "Control":
            dge_temp[i].loc[:,j] = dge.loc[:,cnd==j].mean(axis=1)
        else:
            dge_temp[i].loc[:,j] = dge.loc[:,(cnd==i)&(tpt==j)].mean(axis=1)


    # Only use PLA products detected in at least 10% of cells in control and treatment groups
    PLA_mask1 = (dge.loc[:,(cnd=="Control")|(cnd==i)]>0).sum(axis=1) >= 0.1*sum((cnd=="Control")|(cnd==i))
    dge_temp[i] = dge_temp[i].loc[PLA_mask1,:]

    # Divide by control
    dge_temp[i] = (dge_temp[i]+1).div((dge_temp[i].iloc[:,0]+1), axis=0)


# Plot
fig, ax = plt.subplots(ncols=3, figsize=(6.2,2.4))
for counter1, i in enumerate(["LPS","PAM","Both"]):
    ax[counter1].plot(range(4), dge_temp[i].T, c="gainsboro", zorder=1)
    ax[counter1].plot(range(4), dge_temp[i].mean(axis=0).T, c="red", zorder=9, lw=2.5)
    ax[counter1].fill_between(range(4),
                              dge_temp[i].mean(axis=0)-dge_temp[i].std(axis=0),
                              dge_temp[i].mean(axis=0)+dge_temp[i].std(axis=0),
                              color='red', ec='None', alpha=0.4, zorder=10)
    ax[counter1].plot([0,3],[1,1], c="black", ls='--', zorder=11)
    ax[counter1].set_title(i)
    ax[counter1].set_ylim(-0.2,3)
    ax[counter1].set_yticks(range(0,4,1))
    ax[counter1].set_xticks(range(4))
    ax[counter1].set_xticklabels(["Control     ","5m","2h","12h"])
ax[0].set_ylabel("PLA product\nfold change")
# Highlight PLA product TLR2:TLR2
# For LPS: TLR2:TLR2 increases greatly at 2h
ax[0].plot(range(4), dge_temp['LPS'].loc["TLR2:TLR2",:], zorder=12, c='tab:blue')
ax[1].plot(range(4), dge_temp['PAM'].loc["TLR2:TLR2",:], zorder=12, c='tab:blue')
ax[2].plot(range(4), dge_temp['Both'].loc["TLR2:TLR2",:], zorder=12, c='tab:blue')
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.2)
fig.savefig("mean_dynamics_pla.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Figure 5e


# Mean dynamics plot: proteins
dge_temp = {}
for counter1, i in enumerate(["LPS","PAM","Both"]):
    dge_temp[i] = pd.DataFrame(np.nan, index=dge_protein.index, columns=["Control","5m","2h","12h"])

    for j in dge_temp[i].columns:
        if j == "Control":
            dge_temp[i].loc[:,j] = dge_protein.loc[:,cnd==j].mean(axis=1)
        else:
            dge_temp[i].loc[:,j] = dge_protein.loc[:,(cnd==i)&(tpt==j)].mean(axis=1)

    # Only use proteins detected in at least 10% of cells in control and treatment group
    prot_mask1 = (dge_protein.loc[:,(cnd=="Control")|(cnd==i)]>0).sum(axis=1) >= 0.1*sum((cnd=="Control")|(cnd==i))
    dge_temp[i] = dge_temp[i].loc[prot_mask1,:]

    # Divide by control
    dge_temp[i] = dge_temp[i].div(dge_temp[i].iloc[:,0], axis=0)

# Plot
fig, ax = plt.subplots(ncols=3, figsize=(6.1,2.4))
for counter1, i in enumerate(["LPS","PAM","Both"]):
    ax[counter1].plot(range(4), dge_temp[i].T, c="gainsboro", zorder=1)
    ax[counter1].plot(range(4), dge_temp[i].mean(axis=0).T, c="red", zorder=2, lw=2.5)
    ax[counter1].fill_between(range(4),
                              dge_temp[i].mean(axis=0)-dge_temp[i].std(axis=0),
                              dge_temp[i].mean(axis=0)+dge_temp[i].std(axis=0),
                              color='red', ec='None', alpha=0.4, zorder=10)
    ax[counter1].plot([0,3],[1,1], c="black", ls='--', zorder=11)
    ax[counter1].set_title(i)
    ax[counter1].set_ylim(-0.2,3)
    ax[counter1].set_yticks([0,1,2,3])
    ax[counter1].set_xticks(range(4))
    ax[counter1].set_xticklabels(["Control     ","5m","2h","12h"])
ax[0].set_ylabel("Protein fold change")
# Highlight certain proteins
# For LPS: TLR2 increases at 2h
# For PAM: CD14 increases greatly at 2h
# For Both: TLR6 starts decreasing at 5m
ax[0].plot(range(4), dge_temp['LPS'].loc["TLR2",:], zorder=12, c='tab:blue')
ax[1].plot(range(4), dge_temp['PAM'].loc["CD14",:], zorder=12, c='tab:blue')
ax[2].plot(range(4), dge_temp['Both'].loc["TLR6",:], zorder=12, c='tab:blue')
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.2)
fig.savefig("mean_dynamics_protein.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Figure 5f


# Violin plot the highlighted PLA products
dge_temp = {}
dge_temp['LPS'] = dge.loc[["TLR2:TLR2"],(cnd=="LPS")|(cnd=="Control")].T
dge_temp['LPS']['tpt'] = tpt[(cnd=="LPS")|(cnd=="Control")]
dge_temp['PAM'] = dge.loc[["TLR2:TLR2"],(cnd=="PAM")|(cnd=="Control")].T
dge_temp['PAM']['tpt'] = tpt[(cnd=="PAM")|(cnd=="Control")]
dge_temp['Both'] = dge.loc[["TLR2:TLR2"],(cnd=="Both")|(cnd=="Control")].T
dge_temp['Both']['tpt'] = tpt[(cnd=="Both")|(cnd=="Control")]
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9,8.5))
np.random.seed(1)
for row, i in enumerate(["LPS","PAM","Both"]):
    for col, j in enumerate(dge_temp[i].columns):
        if j == "tpt":
            continue
        sns.stripplot(data=dge_temp[i], x="tpt", y=j, ax=ax[row,col],
                      order=["Control","5m","2h","12h"], jitter=0.3, color='silver', zorder=1)
        sns.pointplot(data=dge_temp[i], x="tpt", y=j, ax=ax[row,col],
                     order=["Control","5m","2h","12h"], ci='sd', color='tab:blue', zorder=10)
fig.tight_layout()

#%% Predict ligand using machine learning and PLA count
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc

# Remove PLA products detected in fewer than 10 cells
dge_ml = dge_log.loc[(dge_log>0).sum(axis=1)>=10,:].copy()
# remove cells with 0 count of PLA products
dge_ml = dge_ml.loc[:,dge_ml.sum(axis=0)>0]

#============ Prediction at 5m only ============#
X_lps = dge_ml.loc[:,(cnd=="LPS") & (tpt=="5m")].T
Y_lps = ["LPS" for _ in range(X_lps.shape[0])]
X_pam = dge_ml.loc[:,(cnd=="PAM") & (tpt=="5m")].T
Y_pam = ["PAM" for _ in range(X_pam.shape[0])]
X_both = dge_ml.loc[:,(cnd=="Both") & (tpt=="5m")].T
Y_both = ["PAM" for _ in range(X_both.shape[0])]
# Combine
X = pd.concat([X_lps, X_pam])
Y = np.array(Y_lps + Y_pam)

# Number of LPS and PAM cells
np.unique(Y, return_counts=True)

# Split data into training and validation sets with 5-fold CV
fpr_5m = np.linspace(0,1,100)
tpr_5m = []
auc_5m = []
fig, ax = plt.subplots(figsize=(2.75,2.55))
kfold_counter = 1
temp = {"fpr":[], "tpr":[]}
for train_idx, test_idx in KFold(n_splits=5, shuffle=True, random_state=1).split(X):
    # print(len(train_idx), len(test_idx))

    # Scale data
    scaler_ml = StandardScaler()
    scaler_ml.fit(X.iloc[train_idx,:])
    X_train = scaler_ml.transform(X.iloc[train_idx,:])
    X_test = scaler_ml.transform(X.iloc[test_idx,:])

    # Train logistic regression
    my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_train, Y[train_idx])
    y_test_lr = my_lr.predict(X_test)
    print(my_lr.score(X_test, Y[test_idx]))
    y_score = my_lr.decision_function(X_test)
    fpr1, tpr1, _ = roc_curve(Y[test_idx], my_lr.predict_proba(X_test)[:,1], pos_label='PAM')

    # Store data for Source data file
    temp["fpr"].append(fpr1)
    temp["tpr"].append(tpr1)

    # 1D interpolation
    tpr1p = np.interp(fpr_5m, fpr1, tpr1)
    tpr1p[0] = 0.0

    tpr_5m.append(tpr1p)
    auc_5m.append(auc(fpr1, tpr1))

    ax.plot(fpr1, tpr1, label=f"Fold {kfold_counter}", lw=2, alpha=0.7)
    kfold_counter += 1
# ax.plot(fpr_5m, np.mean(tpr_5m,axis=0), c='k', lw=3, label='Mean')
ax.plot([0,1], [0,1], c='k', ls='--', lw=1.5)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("5-minute time point")
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))
fig.savefig("logreg_5m_ROC_5fold.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8a

#============ Prediction at 12h only ============#
X_lps = dge_ml.loc[:,(cnd=="LPS") & (tpt=="12h")].T
Y_lps = ["LPS" for _ in range(X_lps.shape[0])]
X_pam = dge_ml.loc[:,(cnd=="PAM") & (tpt=="12h")].T
Y_pam = ["PAM" for _ in range(X_pam.shape[0])]
X_both = dge_ml.loc[:,(cnd=="Both") & (tpt=="12h")].T
Y_both = ["PAM" for _ in range(X_both.shape[0])]
# Combine
X = pd.concat([X_lps, X_pam])
Y = np.array(Y_lps + Y_pam)

# Split data into training and validation sets with 5-fold CV
fpr_12h = np.linspace(0,1,100)
tpr_12h = []
auc_12h = []
fig, ax = plt.subplots(figsize=(2.75,2.55))
kfold_counter = 1
temp = {"fpr":[], "tpr":[]}
for train_idx, test_idx in KFold(n_splits=5, shuffle=True, random_state=1).split(X):
    # print(len(train_idx), len(test_idx))

    # Scale data
    scaler_ml = StandardScaler()
    scaler_ml.fit(X.iloc[train_idx,:])
    X_train = scaler_ml.transform(X.iloc[train_idx,:])
    X_test = scaler_ml.transform(X.iloc[test_idx,:])

    # Train logistic regression
    my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_train, Y[train_idx])
    y_test_lr = my_lr.predict(X_test)
    print(my_lr.score(X_test, Y[test_idx]))
    fpr3, tpr3, _ = roc_curve(Y[test_idx], my_lr.predict_proba(X_test)[:,1], pos_label='PAM')

    # Store data for Source data file
    temp["fpr"].append(fpr3)
    temp["tpr"].append(tpr3)

    # 1D interpolation
    tpr1p = np.interp(fpr_12h, fpr3, tpr3)
    tpr1p[0] = 0.0

    tpr_12h.append(tpr1p)
    auc_12h.append(auc(fpr3, tpr3))

    # Plot to check
    ax.plot(fpr3, tpr3, label=f"Fold {kfold_counter}", lw=2, alpha=0.7)
    kfold_counter += 1
# ax.plot(fpr_12h, np.mean(tpr_12h,axis=0), c='k', lw=3, label='Mean')
ax.plot([0,1], [0,1], c='k', ls='--', lw=1.5)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("12-hour time point")
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))
fig.savefig("logreg_12h_ROC_5fold.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8c

#============ Prediction at 2h only ============#
# Best performance
X_lps = dge_ml.loc[:,(cnd=="LPS") & (tpt=="2h")].T
Y_lps = ["LPS" for _ in range(X_lps.shape[0])]
X_pam = dge_ml.loc[:,(cnd=="PAM") & (tpt=="2h")].T
Y_pam = ["PAM" for _ in range(X_pam.shape[0])]
X_both = dge_ml.loc[:,(cnd=="Both") & (tpt=="2h")].T
# Y_both = ["PAM" for _ in range(X_both.shape[0])]
# Combine
X = pd.concat([X_lps, X_pam])
Y = np.array(Y_lps + Y_pam)

# Split data into training and validation sets with 5-fold CV
fpr_2h = np.linspace(0,1,100)
tpr_2h = []
auc_2h = []
fig, ax = plt.subplots(figsize=(2.75,2.55))
kfold_counter = 1
lr_coef_2h = {}
temp = {"fpr":[], "tpr":[]}
for train_idx, test_idx in KFold(n_splits=5, shuffle=True, random_state=1).split(X):
    # print(len(train_idx), len(test_idx))

    # Scale data
    scaler_ml = StandardScaler()
    scaler_ml.fit(X.iloc[train_idx,:])
    X_train = scaler_ml.transform(X.iloc[train_idx,:])
    X_test = scaler_ml.transform(X.iloc[test_idx,:])

    # Train logistic regression
    my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_train, Y[train_idx])
    y_test_lr = my_lr.predict(X_test)
    print(my_lr.score(X_test, Y[test_idx]))
    fpr2, tpr2, _ = roc_curve(Y[test_idx], my_lr.predict_proba(X_test)[:,1], pos_label='PAM')

    # Store data for Source data file
    temp["fpr"].append(fpr2)
    temp["tpr"].append(tpr2)

    # 1D interpolation
    tpr1p = np.interp(fpr_2h, fpr2, tpr2)
    tpr1p[0] = 0.0

    tpr_2h.append(tpr1p)
    auc_2h.append(auc(fpr2,tpr2))

    # Store coefficients that are common across all 5 folds
    lr_coef_2h[kfold_counter] = pd.Series(my_lr.coef_.reshape(-1,))
    lr_coef_2h[kfold_counter].index = X.columns.copy()

    # Plot to check
    ax.plot(fpr2, tpr2, label=f"Fold {kfold_counter}", lw=2, alpha=0.7)
    kfold_counter += 1
# ax.plot(fpr_2h, np.mean(tpr_2h,axis=0), c='k', lw=3, label='Mean')
ax.plot([0,1], [0,1], c='k', ls='--', lw=1.5)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("2-hour time point")
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))
fig.savefig("logreg_2h_ROC_5fold.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8b

# Concatenate the coefficients across all 5 folds
lr_coef_2h_df = pd.DataFrame(lr_coef_2h)

# Average and plot
lr_coef_2h_mean = lr_coef_2h_df.mean(axis=1).sort_values(ascending=True)
lr_coef_2h_sem = lr_coef_2h_df.std(axis=1)/np.sqrt(5)

# Plot ROC
fig, ax = plt.subplots(figsize=(2.9,4))
ax.plot(fpr_5m, np.mean(tpr_5m, axis=0), lw=2, clip_on=False, zorder=10)
ax.plot(fpr_2h, np.mean(tpr_2h, axis=0), lw=2, clip_on=False, zorder=9)
ax.plot(fpr_12h, np.mean(tpr_12h, axis=0), lw=2, clip_on=False, zorder=9)
ax.plot([0,1], [0,1], ls='--', c='k', clip_on=False)
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("LPS vs PAM prediction", pad=8)
ax.legend([f"5m (AUC = {np.mean(auc_5m):.2f} "+r"$\pm$"+ f" {np.std(auc_5m):.2f})",
           f"2h (AUC = {np.mean(auc_2h):.2f} "+r"$\pm$"+ f" {np.std(auc_2h):.2f})",
           f"12h (AUC = {np.mean(auc_12h):.2f} "+r"$\pm$"+ f" {np.std(auc_12h):.2f})"],
          loc='upper center', bbox_to_anchor=(0.5,-0.25), frameon=False)
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("logreg_AUC_separate_timepoint.svg",
            bbox_inches="tight", pad_inches=0) # Figure 6a

# Only plot coefficients with absolute values above 0.2
coef_to_plot = lr_coef_2h_mean.index[lr_coef_2h_mean.abs()>0.2]

# Prepare a dataframe to plot raw points on bar graph
lr_coef_2h_temp = lr_coef_2h_df.loc[coef_to_plot,:].copy()
lr_coef_2h_temp = pd.melt(lr_coef_2h_temp.reset_index(),
                          id_vars='PLA_pair',
                          var_name='kfold', value_name='level')

# Plot coefficients
fig, ax = plt.subplots(figsize=(4,3.4))
y_temp = np.arange(len(coef_to_plot))
ax.barh(y=coef_to_plot,
        width=lr_coef_2h_mean[coef_to_plot],
        height=0.7, color=my_colors[1], zorder=3)
ax.errorbar(x=lr_coef_2h_mean[coef_to_plot], y=y_temp,
            xerr=lr_coef_2h_sem[coef_to_plot],
            fmt='k', ls='',
            capsize=4, capthick=1.5, zorder=10)
ax.axvline(x=0, color='k', lw=1, zorder=4)
# ax.axhline(y=0.2, color='k', lw=1, ls='--')
# ax.axhline(y=-0.2, color='k', lw=1, ls='--')
np.random.seed(1)
sns.stripplot(data=lr_coef_2h_temp, y='PLA_pair', x='level',
              order=coef_to_plot, size=4, jitter=0.2,
              color="white", edgecolor="black", linewidth=1,
              ax=ax, zorder=7)
ax.set_ylabel("PLA product")
ax.set_xlabel("Logistic regression\ncoefficient's value")
ax.set_xticks(np.arange(-2,2.1,1))
ax.set_xlim(-2,2)
ax.grid(True, axis='both', zorder=0)
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("logreg_coefficients.svg",
            bbox_inches="tight", pad_inches=0) # Figure 6b

# Plot the top 3 PLA product markers
pla_to_plot = ["TLR2:TLR2","MD2:TGFBR1","MD2:CD14"]
dge_to_plot = dge_log.loc[pla_to_plot, (cnd=="Control") | (tpt=="2h")].T
dge_to_plot['cnd'] = "Control"
for i in ["LPS","PAM","Both"]:
    dge_to_plot.loc[dge_to_plot.index.str.contains(i),'cnd'] = i

fig, ax = plt.subplots(nrows=3, figsize=(3,8.5))
np.random.seed(1)
for counter, i in enumerate(pla_to_plot):
    sns.violinplot(x='cnd', y=i, data=dge_to_plot, ax=ax[counter],
                   order=["Control","LPS","PAM","Both"], saturation=1,
                   cut=0, inner=None)
    sns.stripplot(x='cnd', y=i, data=dge_to_plot, ax=ax[counter],
                  order=["Control","LPS","PAM","Both"],
                  color='k', jitter=0.28, size=3)
    ax[counter].set_xlabel("")
    ax[counter].set_ylabel("")
    ax[counter].set_title(f"{i} PLA product", pad=25)
    # 2-sided t-test
    print(i, stats.ttest_ind(dge_to_plot.loc[dge_to_plot["cnd"]=="LPS",i],
                             dge_to_plot.loc[dge_to_plot["cnd"]=="PAM",i],
                             equal_var=False))
    ax[counter].set_ylabel("log10(UMI + 1)")
ax[1].set_yticks(np.arange(0,1.1,0.2))
ax[2].set_yticks(np.arange(0,1.26,0.25))
sns.despine(fig=fig)
fig.tight_layout(h_pad=1.75)
fig.savefig("logreg_top3markers.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8d

#============ Prediction of Both treatment at 2h ============#
# Train on the whole 2h data set
scaler_ml = StandardScaler()
scaler_ml.fit(X)
X_scaled = scaler_ml.transform(X)
X_both_scaled = scaler_ml.transform(X_both)

# Train logistic regression
my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_scaled, Y)
Y_both_pred_proba = my_lr.predict_proba(X_both_scaled)
# Set 70% threshold to be called LPS or PAM, otherwise set to mixed
Y_both_pred = np.array(["Mixed" for _ in range(X_both.shape[0])])
Y_both_pred[Y_both_pred_proba[:,1]<=0.3] = "LPS"
Y_both_pred[Y_both_pred_proba[:,1]>=0.7] = "PAM"
Y_both_pred = pd.Series(Y_both_pred, index=X_both.index)

# Plot stacked bar to show proportion of each response
values, counts = np.unique(Y_both_pred, return_counts=True)
fig, ax = plt.subplots(figsize=(4.1,2.6))
colors_temp = sns.color_palette('Set1')
for i in range(3):
    if values[i] != "Mixed":
        values[i] = values[i] + "-like"
    ax.bar(x=0, height=counts[i]/sum(counts),
           bottom=sum(counts[:i])/sum(counts),
           color=colors_temp[i],
           width=0.7, label=f"{values[i]}: n = {counts[i]}")
ax.legend(title="Predicted response", frameon=False,
          loc='center left', bbox_to_anchor=(1.05,0.5))
ax.set_xlim(-1,1)
ax.set_xticks([0])
ax.set_xticklabels(["Both treatment\n2-hour time point"])
ax.set_ylabel("Proportion of treated cells")
fig.tight_layout()
fig.savefig("logreg_both_2h_proportion.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8f

# Plot the top 3 PLA product markers for the both treatment group
pla_to_plot = ["TLR2:TLR2","MD2:TGFBR1","MD2:CD14"]
dge_to_plot = dge_log.loc[pla_to_plot, Y_both_pred.index].T
dge_to_plot['pred'] = Y_both_pred



# Plot level of all top coefficients
pla_to_plot = coef_to_plot
dge_to_plot = dge_log.loc[pla_to_plot, Y_both_pred.index].T
dge_to_plot['pred'] = Y_both_pred

fig, ax = plt.subplots(ncols=5, nrows=2, figsize=(11.6,4.8))
np.random.seed(1)
for counter, i in enumerate(pla_to_plot):
    row = counter // 5
    col = counter % 5
    sns.violinplot(x='pred', y=i, data=dge_to_plot, ax=ax[row,col],
                   order=["LPS","Mixed","PAM"], saturation=1, palette="Set1",
                   cut=0, inner=None)
    sns.stripplot(x='pred', y=i, data=dge_to_plot, ax=ax[row,col],
                  order=["LPS","Mixed","PAM"],
                  color='k', jitter=0.28, size=3)
    if row == 1:
        ax[row,col].set_xlabel("Predicted response")
    else:
        ax[row,col].set_xlabel("")
    ax[row,col].set_ylabel("log10(UMI + 1)")
    ax[row,col].set_title(i)
ax[0,0].set_yticks(np.arange(0,1.6,0.5))
ax[0,3].set_yticks(np.arange(0,1.1,0.5))
ax[1,0].set_yticks([0,0.2,0.3])
ax[1,2].set_yticks(np.arange(0,1.1,0.5))
ax[1,3].set_yticks(np.arange(0,1.1,0.5))
sns.despine(fig=fig)
fig.suptitle("PLA product level", y=0.95, fontsize=12)
fig.tight_layout(w_pad=1.5)
fig.savefig("logreg_coefficients_both.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8g

# Plot cartoon of logistic regression
x_cartoon = np.linspace(-5,5,100)
y_cartoon = 1/(1 + np.exp(-x_cartoon))
fig, ax = plt.subplots(figsize=(3,2), tight_layout=True)
ax.plot(x_cartoon, y_cartoon, c='k', lw=2, clip_on=True, zorder=10)
x_cartoon_03 = np.log(0.3/0.7)
x_cartoon_07 = np.log(0.7/0.3)
ax.plot([-5,x_cartoon_03], [0.3,0.3], c='#e41a1c', lw=1.5)
ax.plot([x_cartoon_03,x_cartoon_03], [0,0.3], c='#e41a1c', lw=1.5)
ax.plot([-5,x_cartoon_07], [0.7,0.7], c='#4daf4a', lw=1.5)
ax.plot([x_cartoon_07,x_cartoon_07], [0,0.7], c='#4daf4a', lw=1.5)
ax.set_xlim(-5,5)
ax.set_ylim(0,1)
ax.set_yticks([0,0.3,0.7,1])
ax.set_ylabel("Logistic regression\nprobability")
ax.tick_params(axis='y', length=0)
ax.xaxis.set_visible(False)
sns.despine(fig=fig)
ax.get_figure().savefig("logreg_cartoon.svg",
                        bbox_inches="tight", pad_inches=0) # Extended figure 8e

# Plot all PLA products of predicted cells, averaged by predicted response
# Remove PLA products that were not detected in any cell
dge_temp = X_both.loc[:,X_both.sum(axis=0)>0]
# Average by group
temp_mean = {}
for i in ['LPS','Mixed','PAM']:
    temp_mean[i] = dge_temp.loc[Y_both_pred==i,:].mean(axis=0)

# Calculate z-score
dge_temp_z = pd.DataFrame(temp_mean)
dge_temp_z = dge_temp_z.apply(stats.zscore, axis=1, result_type='broadcast').T

# Make list of colors to pass to clustermap
my_color_dict = {'LPS':sns.color_palette('Set1')[0],
                 'Mixed':sns.color_palette('Set1')[1],
                 'PAM':sns.color_palette('Set1')[2],}
my_row_colors = [my_color_dict[s] for s in dge_temp_z.index]

# Clustermap of PLA products
f = sns.clustermap(dge_temp_z, center=0, cmap='RdBu_r',
                   method='complete', metric='euclidean', figsize=(10.3,3.5),
                   row_colors=my_row_colors,
                   cbar_pos=(1.02,0.54,0.03,0.3),
                   cbar_kws={'label':"z-score\n(PLA product level)"})
f.ax_heatmap.set_xlabel("PLA product")
f.ax_heatmap.set_ylabel('')
yticks = f.ax_heatmap.get_yticks()
f.ax_heatmap.set_yticks(yticks)
f.ax_heatmap.set_yticklabels(f.ax_heatmap.get_yticklabels(), rotation=0)
f.ax_heatmap.tick_params(axis='y', length=0)
f.savefig("logreg_predicted_pla_clustermap.png",
          bbox_inches="tight", pad_inches=0, dpi=600) # Extended figure 8h

# Plot all proteins of predicted cells, averaged by predicted response
dge_temp = np.log10(dge_protein.loc[:,Y_both_pred.index] + 1)
# Average by group
temp_mean = {}
for i in ['LPS','Mixed','PAM']:
    temp_mean[i] = dge_temp.loc[:,Y_both_pred==i].mean(axis=1)

# Calculate z-score
dge_temp_z = pd.DataFrame(temp_mean).T
dge_temp_z = dge_temp_z.apply(stats.zscore, axis=0, result_type='broadcast')

# Clustermap of proteins
f = sns.clustermap(dge_temp_z, center=0, cmap='RdYlGn_r',
                   method='complete', metric='euclidean', figsize=(5.5,2.8),
                   row_colors=my_row_colors, linewidths=1.2,
                   cbar_pos=(1.02,0.43,0.03,0.3), colors_ratio=0.065,
                   cbar_kws={'label':"z-score\n(Protein level)"})
f.ax_heatmap.set_xlabel("Protein")
f.ax_heatmap.set_ylabel('')
yticks = f.ax_heatmap.get_yticks()
f.ax_heatmap.set_yticks(yticks)
f.ax_heatmap.set_yticklabels(f.ax_heatmap.get_yticklabels(), rotation=0)
f.ax_heatmap.tick_params(axis='y', length=0)
f.savefig("logreg_predicted_protein_clustermap.png",
          bbox_inches="tight", pad_inches=0, dpi=600) # Extended figure 8i

# Are the mixed cells bad quality?
dge_temp = pd.DataFrame({'nUMI':dge.loc[:,X_both.index].sum(axis=0)})
dge_temp['pred'] = Y_both_pred[dge_temp.index].copy()
dge_temp["nPLA"] = (dge>0).loc[:,dge_temp.index].sum(axis=0)
# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.5,2.7))
np.random.seed(1)
sns.boxplot(data=dge_temp, x="pred", y="nUMI", ax=ax[0],
            order=["LPS","Mixed","PAM"], palette="Set1", saturation=1,
            width=0.6, fliersize=0)
sns.stripplot(data=dge_temp, x="pred", y="nUMI", ax=ax[0], jitter=0.2,
              order=["LPS","Mixed","PAM"], color='white',
              size=4, edgecolor='k', linewidth=0.8)
ax[0].set_yticks([0,1000,2000,2500])
ax[0].set_xlabel("Predicted response")
ax[0].set_ylabel("Number of UMIs")
ax[0].set_title("Detected PLA product UMI")
sns.boxplot(data=dge_temp, x="pred", y="nPLA", ax=ax[1],
            order=["LPS","Mixed","PAM"], palette="Set1", saturation=1,
            width=0.6, fliersize=0)
sns.stripplot(data=dge_temp, x="pred", y="nPLA", ax=ax[1], jitter=0.2,
              order=["LPS","Mixed","PAM"], color='white',
              size=4, edgecolor='k', linewidth=0.8)
ax[1].set_yticks(range(0,151,50))
ax[1].set_xlabel("Predicted response")
ax[1].set_ylabel("Number of PLA products")
ax[1].set_title("Detected PLA products")
sns.despine(fig=fig)
fig.tight_layout(w_pad=2.5)
fig.savefig("QC_2h_both_logreg.svg",
            bbox_inches="tight", pad_inches=0) # Supplementary figure 8a, b

# Compare to LPS and PAM samples
dge_temp1 = dge.loc[:,((cnd=="LPS")|(cnd=="PAM")) & (tpt=="2h")]
dge_temp = pd.DataFrame({'nUMI':dge_temp1.sum(axis=0)})
dge_temp['nPLA'] = (dge_temp1>0).sum(axis=0)
dge_temp['Treatment'] = "LPS"
dge_temp.loc[dge_temp.index.str.contains("PAM"),'Treatment'] = "PAM"
# Plot
fig, ax = plt.subplots(ncols=2, figsize=(5.5,2.7))
np.random.seed(1)
sns.boxplot(data=dge_temp, x="Treatment", y="nUMI", ax=ax[0],
            saturation=1,
            width=0.6, fliersize=0)
sns.stripplot(data=dge_temp, x="Treatment", y="nUMI", ax=ax[0], jitter=0.2,
              color='white',
              size=4, edgecolor='k', linewidth=0.8)
ax[0].set_yticks([0,1000,2000,4000])
ax[0].set_ylabel("Number of UMIs")
ax[0].set_title("Detected PLA product UMI")
sns.boxplot(data=dge_temp, x="Treatment", y="nPLA", ax=ax[1],
            saturation=1,
            width=0.6, fliersize=0)
sns.stripplot(data=dge_temp, x="Treatment", y="nPLA", ax=ax[1], jitter=0.2,
              color='white',
              size=4, edgecolor='k', linewidth=0.8)
ax[1].set_yticks(range(0,151,50))
ax[1].set_ylabel("Number of PLA products")
ax[1].set_title("Detected PLA products")
sns.despine(fig=fig)
fig.tight_layout(w_pad=2.5)
fig.savefig("QC_2h_LPS_PAM_logreg.svg",
            bbox_inches="tight", pad_inches=0) # Supplementary figure 8c, d

#%% Logistic regression on proteins
# Remove PLA products detected in fewer than 10 cells
dge_ml = np.log10(dge_protein+1).copy()

#============ Prediction at 5m only ============#
X_lps = dge_ml.loc[:,(cnd=="LPS") & (tpt=="5m")].T
Y_lps = ["LPS" for _ in range(X_lps.shape[0])]
X_pam = dge_ml.loc[:,(cnd=="PAM") & (tpt=="5m")].T
Y_pam = ["PAM" for _ in range(X_pam.shape[0])]
# Combine
X = pd.concat([X_lps, X_pam])
Y = np.array(Y_lps + Y_pam)

# Number of LPS and PAM cells
np.unique(Y, return_counts=True)

# Split data into training and validation sets with 5-fold CV
fpr_5m = np.linspace(0,1,100)
tpr_5m = []
auc_5m = []
fig, ax = plt.subplots(figsize=(2.75,2.55))
kfold_counter = 1
for train_idx, test_idx in KFold(n_splits=5, shuffle=True, random_state=1).split(X):
    # print(len(train_idx), len(test_idx))

    # Scale data
    scaler_ml = StandardScaler()
    scaler_ml.fit(X.iloc[train_idx,:])
    X_train = scaler_ml.transform(X.iloc[train_idx,:])
    X_test = scaler_ml.transform(X.iloc[test_idx,:])

    # Train logistic regression
    my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_train, Y[train_idx])
    y_test_lr = my_lr.predict(X_test)
    print(my_lr.score(X_test, Y[test_idx]))
    y_score = my_lr.decision_function(X_test)
    fpr1, tpr1, _ = roc_curve(Y[test_idx], my_lr.predict_proba(X_test)[:,1], pos_label='PAM')

    # 1D interpolation
    tpr1p = np.interp(fpr_5m, fpr1, tpr1)
    tpr1p[0] = 0.0

    tpr_5m.append(tpr1p)
    auc_5m.append(auc(fpr1, tpr1))

    ax.plot(fpr1, tpr1, label=f"Fold {kfold_counter}", lw=1, alpha=0.7)
    kfold_counter += 1
ax.plot(fpr_5m, np.mean(tpr_5m,axis=0), c='k', lw=3, label='Mean')
ax.plot([0,1], [0,1], c='r', ls='--', lw=2)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("5-minute time point")
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))

#============ Prediction at 12h only ============#
X_lps = dge_ml.loc[:,(cnd=="LPS") & (tpt=="12h")].T
Y_lps = ["LPS" for _ in range(X_lps.shape[0])]
X_pam = dge_ml.loc[:,(cnd=="PAM") & (tpt=="12h")].T
Y_pam = ["PAM" for _ in range(X_pam.shape[0])]
# Combine
X = pd.concat([X_lps, X_pam])
Y = np.array(Y_lps + Y_pam)

# Split data into training and validation sets with 5-fold CV
fpr_12h = np.linspace(0,1,100)
tpr_12h = []
auc_12h = []
fig, ax = plt.subplots(figsize=(2.75,2.55))
kfold_counter = 1
for train_idx, test_idx in KFold(n_splits=5, shuffle=True, random_state=1).split(X):
    # print(len(train_idx), len(test_idx))

    # Scale data
    scaler_ml = StandardScaler()
    scaler_ml.fit(X.iloc[train_idx,:])
    X_train = scaler_ml.transform(X.iloc[train_idx,:])
    X_test = scaler_ml.transform(X.iloc[test_idx,:])

    # Train logistic regression
    my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_train, Y[train_idx])
    y_test_lr = my_lr.predict(X_test)
    print(my_lr.score(X_test, Y[test_idx]))
    fpr3, tpr3, _ = roc_curve(Y[test_idx], my_lr.predict_proba(X_test)[:,1], pos_label='PAM')

    # 1D interpolation
    tpr1p = np.interp(fpr_12h, fpr3, tpr3)
    tpr1p[0] = 0.0

    tpr_12h.append(tpr1p)
    auc_12h.append(auc(fpr3, tpr3))

    # Plot to check
    ax.plot(fpr3, tpr3, label=f"Fold {kfold_counter}", lw=1, alpha=0.7)
    kfold_counter += 1
ax.plot(fpr_12h, np.mean(tpr_12h,axis=0), c='k', lw=3, label='Mean')
ax.plot([0,1], [0,1], c='r', ls='--', lw=2)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("12-hour time point")
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))

#============ Prediction at 2h only ============#
# Best performance
X_lps = dge_ml.loc[:,(cnd=="LPS") & (tpt=="2h")].T
Y_lps = ["LPS" for _ in range(X_lps.shape[0])]
X_pam = dge_ml.loc[:,(cnd=="PAM") & (tpt=="2h")].T
Y_pam = ["PAM" for _ in range(X_pam.shape[0])]
X_both = dge_ml.loc[:,(cnd=="Both") & (tpt=="2h")].T
# Y_both = ["PAM" for _ in range(X_both.shape[0])]
# Combine
X = pd.concat([X_lps, X_pam])
Y = np.array(Y_lps + Y_pam)

# Split data into training and validation sets with 5-fold CV
fpr_2h = np.linspace(0,1,100)
tpr_2h = []
auc_2h = []
fig, ax = plt.subplots(figsize=(2.75,2.55))
kfold_counter = 1
lr_coef_2h = {}
for train_idx, test_idx in KFold(n_splits=5, shuffle=True, random_state=1).split(X):
    # print(len(train_idx), len(test_idx))

    # Scale data
    scaler_ml = StandardScaler()
    scaler_ml.fit(X.iloc[train_idx,:])
    X_train = scaler_ml.transform(X.iloc[train_idx,:])
    X_test = scaler_ml.transform(X.iloc[test_idx,:])

    # Train logistic regression
    my_lr = LogisticRegression(penalty='l1', solver='liblinear', random_state=1).fit(X_train, Y[train_idx])
    y_test_lr = my_lr.predict(X_test)
    print(my_lr.score(X_test, Y[test_idx]))
    fpr2, tpr2, _ = roc_curve(Y[test_idx], my_lr.predict_proba(X_test)[:,1], pos_label='PAM')

    # 1D interpolation
    tpr1p = np.interp(fpr_2h, fpr2, tpr2)
    tpr1p[0] = 0.0

    tpr_2h.append(tpr1p)
    auc_2h.append(auc(fpr2,tpr2))

    # Store coefficients that are common across all 5 folds
    lr_coef_2h[kfold_counter] = pd.Series(my_lr.coef_.reshape(-1,))
    lr_coef_2h[kfold_counter].index = X.columns.copy()

    # Plot to check
    ax.plot(fpr2, tpr2, label=f"Fold {kfold_counter}", lw=1, alpha=0.7)
    kfold_counter += 1
ax.plot(fpr_2h, np.mean(tpr_2h,axis=0), c='k', lw=3, label='Mean')
ax.plot([0,1], [0,1], c='r', ls='--', lw=2)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("2-hour time point")
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))

# Concatenate the coefficients across all 5 folds
lr_coef_2h_df = pd.DataFrame(lr_coef_2h)

# Average and plot
lr_coef_2h_mean = lr_coef_2h_df.mean(axis=1).sort_values(ascending=True)
lr_coef_2h_sem = lr_coef_2h_df.std(axis=1)/np.sqrt(5)

# Plot ROC
fig, ax = plt.subplots(figsize=(5.3,2.7))
ax.plot(fpr_5m, np.mean(tpr_5m, axis=0), lw=2, clip_on=False, zorder=10)
ax.plot(fpr_2h, np.mean(tpr_2h, axis=0), lw=2, clip_on=False, zorder=9)
ax.plot(fpr_12h, np.mean(tpr_12h, axis=0), lw=2, clip_on=False, zorder=9)
ax.plot([0,1], [0,1], ls='--', c='k', clip_on=False)
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_xticks(np.arange(0,1.1,0.2))
ax.set_yticks(np.arange(0,1.1,0.2))
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("LPS vs PAM prediction", pad=8)
ax.legend([f"5m (AUC = {np.mean(auc_5m):.2f} "+r"$\pm$"+ f" {np.std(auc_5m):.2f})",
           f"2h (AUC = {np.mean(auc_2h):.2f} "+r"$\pm$"+ f" {np.std(auc_2h):.2f})",
           f"12h (AUC = {np.mean(auc_12h):.2f} "+r"$\pm$"+ f" {np.std(auc_12h):.2f})"],
          loc='center left', bbox_to_anchor=(1.05,0.5), frameon=False, title="Time point")
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig("logreg_protein_AUC_separate_timepoint.svg",
            bbox_inches="tight", pad_inches=0) # Extended figure 8j

#%% Quantify protein complex
# Estimate complex abundance for each group
dge_complex_all = {}
dge_temp = dge.loc[:,cnd=="Control"]
dge_complex_all["Control"] = PF.estimateComplexes(dge_temp, nIter=200, mean_cutoff=1, p_cutoff=0.05, p_adjust=True) # control
for i in ["LPS","PAM","Both"]:
    for j in ["5m","2h","12h"]:
        dge_temp = dge.loc[:,(cnd==i)&(tpt==j)]
        dge_complex_all[f"{i}-{j}"] = PF.estimateComplexes(dge_temp, nIter=200, mean_cutoff=1, p_cutoff=0.05, p_adjust=True)

# Mean of complex by sample
dge_complex_all_mean = {}
for i in dge_complex_all:
    dge_complex_all_mean[i] = dge_complex_all[i].mean(axis=1)
dge_complex_all_mean = pd.DataFrame(dge_complex_all_mean)

# Prepare data for violin plot
my_pc = ["TLR2:TLR2","MMP16:CD36","MD2:CD36"] # interesting protein complexes
temp = []
for i in dge_complex_all:
    temp.append(dge_complex_all[i].loc[my_pc,:].T)
    if i == "Control":
        a = "Control"
        b = "Control"
    else:
        a, b = i.split('-')
    temp[-1]['cnd'] = a
    temp[-1]['tpt'] = b
temp_violin = pd.concat(temp, ignore_index=True)

# Plot mean and sem
temp_violin_mean = temp_violin.groupby(['cnd','tpt']).mean()
temp_violin_sem = temp_violin.groupby(['cnd','tpt']).sem(ddof=1)
x = ["Control","5m","2h","12h"]
fig, ax = plt.subplots(ncols=3, figsize=(9,2.8))
for col, i in enumerate(my_pc):
    for counter, j in enumerate(["LPS","PAM","Both"]):
        my_mean = [temp_violin_mean.loc[("Control","Control"),i]]
        my_sem = [temp_violin_sem.loc[("Control","Control"),i]]
        my_mean += temp_violin_mean.loc[(j,x[1:]),i].tolist()
        my_sem += temp_violin_sem.loc[(j,x[1:]),i].tolist()

        ax[col].errorbar(x=range(4), y=my_mean, yerr=my_sem,
                         capsize=3, capthick=2,
                         ls='-', lw=2, marker='.', ms=10)
        ax[col].set_xticks(range(4))
        ax[col].set_xticklabels(x)
ax[0].set_title("TLR2:TLR2 complex")
ax[1].set_title("MMP16:CD36 complex")
ax[2].set_title("MD2:CD36 complex")
ax[0].set_ylim(-0.35,8)
ax[0].set_ylabel("Complex count (UMI)")
ax[-1].legend(["LPS","PAM","Both"], loc='center left', title="Treatment",
              bbox_to_anchor=(1.05,0.5), frameon=False)
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.75)
fig.savefig("complex_dynamics.svg",
            bbox_inches="tight", pad_inches=0) # Figure 6c

#%% Analyze the stochasticity of signaling
# Mean-variance: log-transformed
fig, ax = plt.subplots(figsize=(2.65,2.8))
ax.scatter(dge_log.loc[:,cnd=="Control"].mean(axis=1),
           dge_log.loc[:,cnd=="Control"].var(axis=1), s=20, alpha=0.7, ec="none", zorder=10)
ax.scatter(dge_log.loc[:,(cnd=="LPS") & (tpt=="5m")].mean(axis=1),
           dge_log.loc[:,(cnd=="LPS") & (tpt=="5m")].var(axis=1), s=20, alpha=0.7, ec="none", zorder=10)
ax.axvline(x=1, ls='--', c='k', zorder=1, lw=1.25)
ax.set_xlabel("Mean")
ax.set_ylabel("Variance")
ax.set_xticks(np.arange(0,2.1,0.5))
ax.set_yticks(np.arange(0,0.61,0.2))
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(myDir+"data_analysis/figures/stochasticity_example_mean_var.svg",
            bbox_inches="tight", pad_inches=0) # Figure 6d

# Plot all condition
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(8.5,8.5))
# temp = {} # for saving source data
for row, i in enumerate(["LPS","PAM","Both"]):
    for col, j in enumerate(["5m","2h","12h"]):
        if row == col == 0:
            continue
        ax[row,col].scatter(dge_log.loc[:,(cnd==i) & (tpt==j)].mean(axis=1),
                            dge_log.loc[:,(cnd==i) & (tpt==j)].var(axis=1),
                            c='grey', s=20, alpha=0.7, ec="none", zorder=10)
        # temp[f"{i}-{j}"] = pd.DataFrame({'mean':dge_log.loc[:,(cnd==i) & (tpt==j)].mean(axis=1),
        #                                  'var':dge_log.loc[:,(cnd==i) & (tpt==j)].var(axis=1)})
        ax[row,col].set_xticks(np.arange(0,2.1,0.5))
        ax[row,col].set_yticks(np.arange(0,6.1,0.2))
        ax[row,col].set_ylim(-0.01,0.6)
        ax[row,col].set_title(f"{i} {j}")
        ax[row,col].set_xlabel("Mean")
        ax[row,col].set_ylabel("Variance")
ax[0,0].axis('off')
fig.tight_layout(h_pad=1.3, w_pad=1.3)
sns.despine(fig=fig)
fig.savefig("stochasticity_all_mean_var.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Supplementary figure 9

# PLA products with high mean
temp1 = dge_log.index[dge_log.loc[:,(cnd=="Control")].mean(axis=1)>=1]

# Plot histogram of each of the 9 CD36-related PLA products
print(f"Max count: {dge_log.loc[temp1,:].max().max()}")
my_bins = np.arange(0,3.30,0.25)
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(7.5,7.5))
for i in range(len(temp1)):
    row = i // 3
    col = i % 3
    ax[row,col].hist(dge_log.loc[temp1[i],cnd=="Control"], alpha=0.7,
                     bins=my_bins, density=True)
    ax[row,col].hist(dge_log.loc[temp1[i],(cnd=="LPS")&(tpt=="5m")], alpha=0.5,
                     bins=my_bins, density=True)
    ax[row,col].axvline(x=1, ls='--', c='k', zorder=1, lw=1.25)
    ax[row,col].set_title(f"{temp1[i]}\nPLA product")
    ax[row,col].set_xticks(range(4))
for i in range(3):
    ax[-1,i].set_xlabel("log10(UMI + 1)")
    ax[i,0].set_ylabel("Density")
fig.tight_layout(w_pad=1.5)
sns.despine(fig=fig)
fig.savefig("stochasticity_CD36_PLA_hist.png",
            bbox_inches="tight", pad_inches=0, dpi=600) # Supplementary figure 10

# Average the histogram of all 9 PLA products by treatment vs control
my_bins = np.arange(0,3.30,0.25)
my_hist = {'Control':[], 'Treatment':[]}
for i in range(len(temp1)):
    temp_hist,_ = np.histogram(dge_log.loc[temp1[i],cnd=="Control"], bins=my_bins)
    my_hist['Control'].append(temp_hist/sum(temp_hist))
    temp_hist,_ = np.histogram(dge_log.loc[temp1[i],cnd!="Control"], bins=my_bins)
    my_hist['Treatment'].append(temp_hist/sum(temp_hist))
my_hist['Control'] = np.array(my_hist['Control'])
my_hist['Treatment'] = np.array(my_hist['Treatment'])

fig, ax = plt.subplots(ncols=2, figsize=(5.1,3))
ax[0].plot((my_bins[1:]+my_bins[:-1])/2, np.mean(my_hist['Control'], axis=0),
           color=my_colors[0], lw=2, zorder=10)
ax[0].fill_between((my_bins[1:]+my_bins[:-1])/2,
                   np.mean(my_hist['Control'], axis=0)-np.std(my_hist['Control'], axis=0, ddof=1),
                   np.mean(my_hist['Control'], axis=0)+np.std(my_hist['Control'], axis=0, ddof=1),
                   color=my_colors[0], alpha=0.4, ec='None', zorder=9)
ax[0].set_title("Control")
ax[1].plot((my_bins[1:]+my_bins[:-1])/2, np.mean(my_hist['Treatment'], axis=0),
           color=my_colors[1], lw=2, zorder=10)
ax[1].fill_between((my_bins[1:]+my_bins[:-1])/2,
                   np.mean(my_hist['Treatment'], axis=0)-np.std(my_hist['Treatment'], axis=0, ddof=1),
                   np.mean(my_hist['Treatment'], axis=0)+np.std(my_hist['Treatment'], axis=0, ddof=1),
                   color=my_colors[1], alpha=0.4, ec='None', zorder=9)
ax[1].set_title("Treatment")
for i in range(2):
    ax[i].set_xlabel("PLA product count\nlog10(UMI + 1)")
    ax[i].set_xlim(0,3.25)
    ax[i].set_ylim(-0.03,0.3)
    ax[i].set_xticks(np.arange(0,3.1,1))
    ax[i].set_yticks(np.arange(0,0.31,0.1))
    ax[i].set_ylabel("Density")
    ax[i].grid(True, axis='both', zorder=0)
sns.despine(fig=fig)
fig.tight_layout(w_pad=1.5)
fig.savefig("stochasticity_hist_treatment_grouped.svg",
            bbox_inches="tight", pad_inches=0) # Figure 6f

# Quantify % population with mean of histogram below 1 log10(UMI+1)
below_ctrl = np.sum(np.mean(my_hist['Control'], axis=0)[:4])
below_trt = np.sum(np.mean(my_hist['Treatment'], axis=0)[:4])
print(below_ctrl, below_trt)

# Plot example to show correlated reduction in variance
fig, ax = plt.subplots(ncols=1, figsize=(2.75,2.9), tight_layout=True)
ax.scatter(dge_log.loc["IL-10R:CD36",cnd=="Control"],
           dge_log.loc["MD2:CD36",cnd=="Control"],
           edgecolors='none', zorder=5, s=25)
ax.scatter(dge_log.loc["IL-10R:CD36",(cnd=="LPS")&(tpt=="5m")],
           dge_log.loc["MD2:CD36",(cnd=="LPS")&(tpt=="5m")],
           edgecolors='none', zorder=5, s=25)
ax.axvline(x=1, ls='--', lw=1.5, c='k', zorder=1)
ax.axhline(y=1, ls='--', lw=1.5, c='k', zorder=1)
ax.set_xlabel("IL-10R:CD36 count\nlog10(UMI + 1)")
ax.set_ylabel("MD2:CD36 count\nlog10(UMI + 1)")
ax.set_xticks(np.arange(0,3.1,1))
ax.set_xlim(-0.3,3)
ax.set_ylim(-0.3,3)
sns.despine(fig=fig)
fig.savefig("stochasticity_example_scatter.svg",
            bbox_inches="tight", pad_inches=0) # Figure 6e
