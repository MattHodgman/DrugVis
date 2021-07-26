'''
Process:

1. get list of ChEMBL ID's for 3 ATC classifications level2 --> diabete (48), antidiarrheals (44), antithrombotics (72)
2. get those drug's pref_name
3. get those drugs' rows in distance table
4. transform with different sigma values
5. clustermap

'''


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
    parser.add_argument('-i', '--input', help='input drug lists', nargs='*', action="store", dest="input", required=True)
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
Main.
'''
if __name__ == '__main__':

    args = parseArgs() # parse arguments

    classes = {} # dict where key = drug class and value = list of drugs

    # read input files
    print('Reading drug lists...')
    for file in args.input:
        classes[getDataName(file)] = readFile(file)


    # extract rows from distance table
    print('Loading in distance data...')
    all_distances = pd.read_csv('~/Harvard/sms/data/all_drug_distances_melted.tsv', delimiter='\t')
    all_distances = all_distances.drop('Log10_Distance', axis=1) # drop unnecessary column

    distances = pd.DataFrame(columns = list(all_distances.columns))

    print('Extracting relevant rows...')
    for c,drugs in classes.items():
        rows1 = all_distances['Drug1'].isin(drugs)
        rows2 = all_distances['Drug2'].isin(drugs)

        distances = distances.append(rows1)
        distances = distances.append(rows2)


    # calculate RBF kernel using different sigma values
    print('Running sigma tests...')
    sigmas = [5,10] # [2,4,6,8,10,15,20]
    for s in sigmas:
        print(f'Calculating RBF kernel (sigma={s})...')
        distances[str(s)] = distances.apply(lambda row : rbfKernel(row['Distance']), axis=1)
        
        # clustermap
        print(f'Clustering and plotting (sigma={s})...')
        distances_s = distances.pivot('Drug1', 'Drug2', str(s)) # pivot
        g = sns.clustermap(distances_s, xticklabels=False, yticklabels=False, linewidths=0.0, cmap='viridis') # cluster and plot
        g.savefig(f'clustermap_rbf_kernel_{s}.png')