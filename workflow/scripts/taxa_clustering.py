import re

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sbn
from Bio import SeqIO
from ete3 import NCBITaxa
import argparse


parser = argparse.ArgumentParser(description = 'cluster Taxa based on peptidome similarity and weight attributed')

parser.add_argument('--GraphIN', type = str, help = 'path(s) to the Peptonizer graphml file for which you wish to cluster taxa')
parser.add_argument('--TaxaWeightFile', type = str,help = 'path to file with weighted taxa computed in the WeightTaxa step')
parser.add_argument('--SimilarityThreshold',type =float, help = 'Threshold for the petidome sinilarity at which two taxa should belong to the same cluster')
parser.add_argument('--outSimilarities', type = str, help = 'path to the snilarities iutput csv')
parser.add_argument('--out', type = str, help= 'path to clustered taxa output csv')

args = parser.parse_args()

ncbi = NCBITaxa()



def GetPeptidesperTaxon(Graphin):
    graph = nx.read_graphml(Graphin)
    PeptidomeDict = {}
    for node in graph.nodes(data=True):
        if node[1]['category'] == 'taxon' and node[0]:
            neighbors = graph.neighbors(node[0])
            PeptidomeDict.update({node[0]: [n[:-4] for n in neighbors]})

    return PeptidomeDict

def ComputeDetectedPeptidomeSimilarity(PeptidomeDict):
    SimMatrixMax = []
    SimMatrixMin = []
    Taxa1 = []
    Taxa2 = []
    for taxon1 in PeptidomeDict.keys():
        Taxa1.append(taxon1)
        SimMatrixMaxRow = []
        SimMatrixMinRow = []
        for taxon2 in PeptidomeDict.keys():
            Taxa2.append(taxon2)
            peptides1 = set(PeptidomeDict[taxon1])
            peptides2 = set(PeptidomeDict[taxon2])
            shared = len(peptides1.intersection(peptides2))
            try:
                Sim = shared / ( len(peptides2))
            except:
                Sim = 0

            SimMatrixMaxRow.append(Sim)

        SimMatrixMax.append(SimMatrixMaxRow)

    
    SimilarityFrame = pd.DataFrame(SimMatrixMax, columns = Taxa1, index = Taxa1)

    return SimilarityFrame


def ClusterTaxaBasedOnSimilarity(TaxaWeightFile,SimilarityFrame,SimilarityTreshold):
    WeightSortedTaxa = pd.read_csv(TaxaWeightFile)
    WeightSortedTaxa = WeightSortedTaxa.loc[WeightSortedTaxa.HigherTaxa.isin(SimilarityFrame.index.tolist())]
    ListOfWeightSortedTaxa = WeightSortedTaxa.HigherTaxa.tolist()
    TaxaClusterList = []
    Clusterheads = [WeightSortedTaxa.HigherTaxa[0]]

    while ListOfWeightSortedTaxa:
            taxon1 = ListOfWeightSortedTaxa[0]
            ClusterList = []
            Clusterheads.append(ListOfWeightSortedTaxa[0])

            for taxon2 in WeightSortedTaxa.HigherTaxa:
                
                if SimilarityFrame[str(taxon2)][str(taxon1)]>SimilarityTreshold:
                    check = SimilarityFrame[str(taxon2)][str(taxon1)]
                    ClusterList.append(taxon2)
                    ListOfWeightSortedTaxa.remove(taxon2)

        
            TaxaClusterList.append(ClusterList)
    

    ClusteredWeighsortedTaxa = WeightSortedTaxa.loc[WeightSortedTaxa.HigherTaxa.isin(Clusterheads)]
    ClusteredWeightsortedTaxa2 = WeightSortedTaxa.loc[WeightSortedTaxa.HigherTaxa.isin([clustertaxa[1:] for clustertaxa in TaxaClusterList])]
    ClusteredWeighsortedTaxa.concat(ClusteredWeightsortedTaxa2)
    ClusteredWeighsortedTaxa['Clustermembers'] = TaxaClusterList

    return ClusteredWeighsortedTaxa



Peptidomedict = GetPeptidesperTaxon(args.GraphIN)
similarities = ComputeDetectedPeptidomeSimilarity(Peptidomedict)
similarities.to_csv(args.outSimilarities)

ClusteredTaxa = ClusterTaxaBasedOnSimilarity(args.TaxaWeightFile,similarities,args.SimilarityThreshold)
ClusteredTaxa.to_csv(args.out)



       