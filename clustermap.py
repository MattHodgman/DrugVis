from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import fastcluster
import numpy
sns.set_theme()

distances_long = pd.read_csv('all_drug_distances_melted.tsv', delimiter='\t') # load data
distances = distances_long.pivot('Drug1', 'Drug2', 'Distance') # pivot
g = sns.clustermap(distances, norm=LogNorm(), xticklabels=False, yticklabels=False, linewidths=0.0, cmap='viridis') # cluster and plot
g.savefig('clustermap4.png')