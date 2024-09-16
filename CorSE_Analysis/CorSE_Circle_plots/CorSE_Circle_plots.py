"""
Housekeeping
"""
#laoding packages
import pandas as pd
import numpy as np
import mne
import matplotlib.pyplot as plt
import seaborn as sns
from mne.viz import circular_layout
from mne_connectivity.viz import plot_connectivity_circle
import os

#change output directory
bd = r"C:\Users\Fraser\Desktop\SHANK2 Analysis\CorSE_Analysis\CorSE_Circle_plots"
os.chdir(bd)

#change default plotting dpi
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


"""
setup circle plot layout
"""
#current electrode layout
label_names = [11, 12, 13, 14, 15, 16, 17, 18,
                 21, 22, 23, 24, 25, 26, 27, 28,
                 31, 32, 33, 34, 35, 36, 37, 38,
                 41, 42, 43, 44, 45, 46, 47, 48,
                 51, 52, 53, 54, 55, 56, 57, 58,
                 61, 62, 63, 64, 65, 66, 67, 68,
                 71, 72, 73, 74, 75, 76, 77, 78,
                 81, 82, 83, 84, 85, 86, 87, 88]

#desired electrode grouping
channel_groups = [11, 12, 21, 22,
                  13, 14, 23, 24,
                  15, 16, 25, 26,
                  17, 18, 27, 28,
                  31, 32, 41, 42,
                  33, 34, 43, 44,
                  35, 36, 45, 46,
                  37, 38, 47, 48,
                  51, 52, 61, 62,
                  53, 54, 63, 64,
                  55, 56, 65, 66,
                  57, 58, 67, 68,
                  71, 72, 81, 82,
                  73, 74, 83, 84,
                  75, 76, 85, 86,
                  77, 78, 87, 88]

#convert to string
label_names = [str(x) for x in label_names]
channel_groups = [str(x) for x in channel_groups]

#create layout
bounds = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60]
node_angles = circular_layout(channel_groups, channel_groups, start_pos = 90, group_boundaries = bounds)

#indicies for re-mapping order of correlation matricies
new_idx = [label_names.index(x) for x in channel_groups]


#please excuse the incredibly lazy coding below. I was writing my thesis and didn't have the energy to make this a loop.
"""
load CorSE data
"""
#load R841X data
mutdat = pd.read_csv("G:\Synchrony Analysis\CorSE Analysis 2023\Analyzed Data\Plate 37\Plate-37_D54_A3_R841X_CorrData.csv", header = None)
mutdat = mutdat.to_numpy()
mutdat = mutdat[:, new_idx]
mutdat = mutdat[new_idx,:]
mutdat = np.where(mutdat<0,0,mutdat)

#load R841X-C data
cordat = pd.read_csv("G:\Synchrony Analysis\CorSE Analysis 2023\Analyzed Data\Plate 39\Plate-39_D62_B1_R841X-C_CorrData.csv", header = None)
cordat = cordat.to_numpy()
cordat = cordat[:, new_idx]
cordat = cordat[new_idx,:]
cordat = np.where(cordat<0,0,cordat)

#load KO data
kodat = pd.read_csv("G:\Synchrony Analysis\CorSE Analysis 2023\Analyzed Data\Plate 41\Plate-41_D50_B4_KO_CorrData.csv", header = None)
kodat = kodat.to_numpy()
kodat = kodat[:, new_idx]
kodat = kodat[new_idx,:]
kodat = np.where(kodat<0,0,kodat)

#load CTRL data
ctrldat = pd.read_csv("G:\Synchrony Analysis\CorSE Analysis 2023\Analyzed Data\Plate 41\Plate-41_D50_A1_CTRL_CorrData.csv", header = None)
ctrldat = ctrldat.to_numpy()
ctrldat = ctrldat[:, new_idx]
ctrldat = ctrldat[new_idx,:]
ctrldat = np.where(ctrldat<0,0,ctrldat)

#load DHPG data
dhpgdat = pd.read_csv("G:\Synchrony Analysis\CorSE Analysis 2023\Analyzed Data\Plate 43\Plate-43_D50_A3_R841X DHPG_CorrData.csv", header = None)
dhpgdat = dhpgdat.to_numpy()
dhpgdat = dhpgdat[:, new_idx]
dhpgdat = dhpgdat[new_idx,:]
dhpgdat = np.where(dhpgdat<0,0,dhpgdat)


"""
Plot R841X
"""
#plot top 150 connections
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    mutdat,
    channel_groups,
    n_lines=150,
    colormap= "inferno",
    vmin = 0.5,
    vmax = 1,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("R841X_map_150.pdf")

#plot all connections with strength > 0.5
ncon = (mutdat >= 0.5).sum()
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    mutdat,
    channel_groups,
    n_lines=ncon,
    colormap= "inferno",
    vmin = 0.5,
    vmax = 1,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("R841X_map_0.5.pdf")

#plot heatmap
mask = np.triu(np.ones_like(mutdat, dtype=bool))
f, ax = plt.subplots(figsize=(10, 10))
cmap = sns.color_palette("inferno", as_cmap=True)
sns.heatmap(mutdat, #mask=mask,
            cmap=cmap, 
            vmax=1, 
            square=True, linewidths=0, cbar_kws={"shrink": .5})
f.savefig("R841X_heatmap.pdf")



"""
Plot R841X-C
"""
#plot top 150 connections
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    cordat,
    channel_groups,
    n_lines=150,
    colormap= "inferno",
    #vmin = 0.4,
    #vmax = 0.8,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("R841X-C_map_150.pdf")

#plot all connections with strength > 0.5
ncon = (cordat >= 0.5).sum()
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    cordat,
    channel_groups,
    n_lines=ncon,
    colormap= "inferno",
    vmin = 0.5,
    vmax = 1,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("R841X-C_map_0.5.pdf")

#plot heatmap
mask = np.triu(np.ones_like(cordat, dtype=bool))
f, ax = plt.subplots(figsize=(10, 10))
cmap = sns.color_palette("inferno", as_cmap=True)
sns.heatmap(cordat, 
            #mask=mask,
            cmap=cmap, 
            vmax=1, 
            square=True, linewidths=0, cbar_kws={"shrink": .5})
f.savefig("R841X-C_heatmap.pdf")




"""
Plot KO data
"""
#plot top 150 connections
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    kodat,
    channel_groups,
    n_lines=150,
    colormap= "inferno",
    #vmin = 0.4,
    #vmax = 0.8,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("KO_map_150.pdf")

#plot all connections with strength > 0.5
ncon = (kodat >= 0.5).sum()
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    kodat,
    channel_groups,
    n_lines=ncon,
    colormap= "inferno",
    vmin = 0.5,
    vmax = 1,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("KO_map_0.5.pdf")

#plot heatmap
mask = np.triu(np.ones_like(kodat, dtype=bool))
f, ax = plt.subplots(figsize=(10, 10))
cmap = sns.color_palette("inferno", as_cmap=True)
sns.heatmap(kodat, 
            #mask=mask,
            cmap=cmap, 
            vmax=1, 
            square=True, linewidths=0, cbar_kws={"shrink": .5})
f.savefig("KO_heatmap.pdf")



"""
Plot CTRL data
"""
#plot top 150 connections
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    ctrldat,
    channel_groups,
    n_lines=150,
    colormap= "inferno",
    #vmin = 0.4,
    #vmax = 0.8,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("CTRL_map_150.pdf")

#plot all connections with strength > 0.5
ncon = (ctrldat >= 0.5).sum()
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    ctrldat,
    channel_groups,
    n_lines=ncon,
    colormap= "inferno",
    vmin = 0.5,
    vmax = 1,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("CTRL_map_0.5.pdf")

#plot heatmap
mask = np.triu(np.ones_like(ctrldat, dtype=bool))
f, ax = plt.subplots(figsize=(10, 10))
cmap = sns.color_palette("inferno", as_cmap=True)
sns.heatmap(ctrldat, 
            #mask=mask,
            cmap=cmap, 
            vmax=1, 
            square=True, linewidths=0, cbar_kws={"shrink": .5})
f.savefig("CTRL_heatmap.pdf")



"""
Plot CTRL data
"""
#plot top 150 connections
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    dhpgdat,
    channel_groups,
    n_lines=150,
    colormap= "inferno",
    #vmin = 0.4,
    #vmax = 0.8,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("DHPG_map_150.pdf")

#plot all connections with strength > 0.5
ncon = (dhpgdat >= 0.5).sum()
fig, ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection="polar"))
plot_connectivity_circle(
    dhpgdat,
    channel_groups,
    n_lines=ncon,
    colormap= "inferno",
    vmin = 0.5,
    vmax = 1,
    facecolor = "white",
    textcolor = "black",
    node_angles=node_angles,
    ax=ax)
fig.savefig("DHPG_map_0.5.pdf")

#plot heatmap
mask = np.triu(np.ones_like(dhpgdat, dtype=bool))
f, ax = plt.subplots(figsize=(10, 10))
cmap = sns.color_palette("inferno", as_cmap=True)
sns.heatmap(dhpgdat, 
            #mask=mask,
            cmap=cmap, 
            vmax=1, 
            square=True, linewidths=0, cbar_kws={"shrink": .5})
f.savefig("DHPG_heatmap.pdf")
