import numpy as np
import pandas as pd
import os
import re 
import matplotlib
from matplotlib import pyplot as plt
from ete3 import NCBITaxa
from  scipy.stats import entropy


#package s ettings
matplotlib.use('Agg')
ncbi = NCBITaxa()



def ComputeMetric(resultsfolder, output, weightsfile):

    """
    Compute a "goodness" metric for the parameters used in the grid search. 
    Generates a barplot of the metric and returns the best identified parameters.
    Return a list of the best parameter set in order [alpha,beta,gamma].

    :param resultsfolder: str, path to PepGM resultsfolder
    :param host: str, host to be excluded from the parameter checked taxa
    :param output: str, name of the output .png file 
    :param weightsfile: str, path to the .csv file with Unipept output
    
    """
    
    #predefine necessary lists
    Params = []
    Metrics = []
    SumProportions = []
    Entropies = []
    WeightCoeffs = []

    #file with weights of taxids
    Weights = pd.read_csv(weightsfile,usecols=['HigherTaxa','scaled_weight'])
    Weights = Weights.groupby(['taxa']).sum().reset_index()
    Maxweight = Weights.max()['weight']
         
    for folders in os.listdir(resultsfolder):
        if os.path.isdir(resultsfolder+ '/' + folders):
            for file in os.listdir(resultsfolder+ '/' + folders):
    
                if file.endswith('.csv'):
               
                    Results = pd.read_csv(resultsfolder+ '/' + folders +'/' + file,names = ['ID','score','type'])
                    TaxIDS = Results.loc[Results['type']=='taxon']
                    TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
                    TaxIDS = TaxIDS.sort_values('score', ascending = False)


                    #compute the metric
                    
                    #what's the weight of the highest scoring taxids?
                    Weight =Weights.loc[Weights['taxa']==int(TaxIDS.ID.head(1).item())]['weight'].head(1).item()
                    WeightCoeff = Weight/Maxweight
                    WeightCoeffs.append(WeightCoeff)
    
                    #compute entropy of the posterior probability distribution
                    Entropy = entropy(TaxIDS['score'])
                    Entropies.append(Entropy)
                  
    
                    Matching = (1/(Entropy))*WeightCoeff
    
                    Metrics.append(Matching)
    
                    #get the corresponding parameters
                    params = re.findall('\d+\.\d+',file)
                    Params.append([params[0],params[1],params[2]])


    zipLists = zip(Metrics,Params)
    SortedPairs = sorted(zipLists,reverse= True)
    
    tuples = zip(*SortedPairs)
    Metrics,Params = [list(tuple) for tuple in tuples]
    
    figure1 = plt.figure(figsize=(30,15), tight_layout=True)
    
    plt.barh(list(range(20)),Metrics[0:20],color = 'mediumvioletred')
    plt.xticks(fontsize =35)
    plt.yticks(list(range(20)),Params[0:20],fontsize=35)
    
    print(output)
    plt.savefig(output)
    plt.close()

    return Params[0]




    




