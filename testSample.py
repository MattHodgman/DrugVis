'''
Process:

1. get list of ChEMBL ID's for 3 ATC classifications level2 --> diabete (48), antidiarrheals (44), antithrombotics (72)
2. get those drug's pref_name
3. get those drugs' rows in distance table
4. transform with different sigma values
5. clustermap

'''

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from numpy.core.arrayprint import printoptions
import pandas as pd
import argparse
from math import exp
from numpy import nan
import seaborn as sns
sns.set_theme()


'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Using three classes of drugs, test different sigma values for RBF kernel transformation on semantic distances for best clustering.')
    parser.add_argument('-i', '--input', help='input drug lists', nargs='*', action='store', dest='input', required=False)
    parser.add_argument('-m', '--matrix', help='matrix of distances', type=str, required=True)
    args = parser.parse_args()
    return args


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1] # get filename from end of input path
    dataName = fileName[:fileName.rfind('.')] # get data name by removing extension from file name
    return dataName


'''
Read file rows into list. Each line should just contain a single drugs ChEMBL preferred name.
'''
def readFile(file):
    drugs = [] # list to store drug names

    # parse file
    with open(file) as f:
        drugs = f.readlines()
    f.close()

    drugs = [d.strip() for d in drugs] # remove trailing whitespace

    return drugs


'''
RBF kernel transformation. 
'd' is the semantic distance/similarity between drugs. 
Set sigma as appropriate.
'''
def rbfKernel(d, sigma):
    try:
        k = exp(-d / (2 * sigma^2))
        return k
    except ZeroDivisionError:
        return nan


'''
Make a seaborn clustermap
'''
def clusterMap():
    print(f'Clustering and plotting (sigma={s})...')
    classes = all_distances.pop('Class')
    consensus_clusters = pd.read_csv('../data/3_classes_drug_distances_RBF_kernel_sigma_20-consensus.tsv', delimiter='\t', index_col='Drug')

    distances_s = distances.pivot('Drug1', 'Drug2', s) # pivot

    # do some stuff to color classes
    lut = dict(zip(classes.unique(), ['darkgrey','grey','tan']))
    row_colors = classes.map(lut)
    row_colors.index = distances_s.index
    handles = [Patch(facecolor=lut[name]) for name in lut]

    lut2 = {0:'lightgrey',1:'red',2:'orange',3:'green',4:'blue',5:'purple'}
    row_colors2 = consensus_clusters['Cluster'].map(lut2)
    row_colors2.index = distances_s.index
    handles2 = [Patch(facecolor=lut2[name]) for name in lut2]
    
    sns.set(font_scale=0.6)
    g = sns.clustermap(distances_s, linewidths=0.0, figsize=(20,20), cmap='viridis', row_colors=[row_colors,row_colors2], col_colors=[row_colors,row_colors2]) # cluster and plot

    legend1 = plt.legend(handles, lut, title='ATC Classes level 2',
       bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc=1)
    legend2 = plt.legend(handles2, lut2, title='Consensus Clusters',
       bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc=2)

    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)

    g.savefig(f'clustermap_rbf_kernel_{s}-clustered-consensus.png')


'''
Plot the similarity scores on a histogram.
'''
def histogram(s):
    print('plotting similarities...')
    similarities = distances[s].tolist()
    plt.hist(similarities, bins=100)
    plt.title(s)
    plt.savefig(f'similarities_histo_{s}.png')


'''
Main.
'''
if __name__ == '__main__':

    args = parseArgs() # parse arguments

    # classes = {} # dict where key = drug class and value = list of drugs

    # read input files
    # print('Reading drug lists...')
    # for file in args.input:
    #     classes[getDataName(file)] = readFile(file)

    # extract rows from distance table
    # print('Loading in distance matrix...')
    # all_distances = pd.read_csv(args.matrix, delimiter='\t') # load matrix

    # distances = (all_distances.rename(columns={"Drug": "Drug1"})
    #     .melt("Drug1", var_name="Drug2", value_name="Distance")
    #     .sort_values(by="Drug1")
    # ) # melt data

    # load in melted distances
    print('loading data...')
    distances = pd.read_csv(args.matrix, delimiter='\t') # load matrix

    # calculate RBF kernel using different sigma values
    print('Running sigma tests...')
    sigmas = [2,5,10,15,20,25,35]
    for s in sigmas:
        print(f'Calculating RBF kernel (sigma={s})...')
        distances[s] = distances.apply(lambda row : rbfKernel(row['Distance'], s), axis=1)

        histogram(s)

        # distances.to_csv(f'3_classes_drug_distances_RBF_kernel_sigma_{s}.tsv', index=False, sep='\t') # save to file
        # clusterMap() # make a cluster map