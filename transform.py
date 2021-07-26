import pandas as pd
from math import exp


'''
RBF kernel transformation. 
'd' is the semantic distance/similarity between drugs. 
Set sigma as appropriate.
'''
def rbfKernel(d):
    sigma = 6 # change this as needed

    try:
        k = exp(-d / (2 * sigma^2))
    except ZeroDivisionError:
        print(d)

    return k

    


'''
Main.
'''
if __name__ == '__main__':

    # get data
    print('getting data')
    distances = pd.read_csv('../data/all_drug_distances_melted.tsv', delimiter='\t') # load in data
    distances = distances.drop('Log10_Distance', axis=1) # drop log10 distance

    print('applying transformation')
    # apply transformation
    distances['RBF_kernel'] = distances.apply(lambda row : rbfKernel(row['Distance']), axis=1)

    print('writing results')
    # write to file
    distances.to_csv('all_drug_distances_melted.tsv', sep='\t')