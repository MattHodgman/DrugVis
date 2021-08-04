import pandas as pd
from math import exp
import sys

# parameters: input melted file name, sigma, output melted transformed file name


'''
RBF kernel transformation. 
'd' is the semantic distance/similarity between drugs. 
Set sigma as appropriate.
'''
def rbfKernel(d):
    sigma = sys.args[2] # change this as needed

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
    distances = pd.read_csv(sys.args[1], delimiter='\t') # load in data
    distances = distances.drop('Log10_Distance', axis=1) # drop log10 distance

    print('applying transformation')
    # apply transformation
    distances['RBF_kernel'] = distances.apply(lambda row : rbfKernel(row['Distance']), axis=1)

    print('writing results')
    # write to file
    distances.to_csv(sys.args[3], sep='\t')