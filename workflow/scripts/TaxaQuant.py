import pandas as pd
import argparse


parser = argparse.ArgumentParser(description = 'calculate relative biomass contributions')
parser.add_argument('--PeptonizerResults', type = str, help = 'path to peptonizer results file')    
parser.add_argument('--GraphdataFrame', type = str, help = 'path to graphdataframe file')
parser.add_argument('--PosteriorThreshold', type = float, help = 'threshold for posterior score')
parser.add_argument('--out', type = str, help = 'path to output file')

args = parser.parse_args()

def  DistributePeptides(PeptonizerResults,GraphdataFrame,PosteriorThreshold):
    PeptonizerPosteriors = pd.read_csv(PeptonizerResults)
    Graphdataframe = pd.read_csv(GraphdataFrame)
    PeptonizerPosteriors = PeptonizerPosteriors.loc[PeptonizerPosteriors.score>PosteriorThreshold] #filter out low scoring taxa

    Graphdataframe = Graphdataframe.loc[Graphdataframe.HigherTaxa.isin(PeptonizerPosteriors.ID.tolist())] #filter out taxa not in peptonizer results
    Graphdataframe['peptide_degeneracy']= Graphdataframe.groupby('sequence')['sequence'].transform('count') #count number of taxa per sequence
    Graphdataframe['peptideWeight'] = 1/(Graphdataframe['peptide_degeneracy']) #calculate peptide weight 
    GroupedResults = Graphdataframe.groupby('HigherTaxa')['peptideWeight'].sum().reset_index() #sum number of taxa per higher taxa
    TotalWeight = GroupedResults.peptideWeight.sum() #sum total number of Weights
    GroupedResults['RelativeWeight'] = GroupedResults['peptideWeight']/TotalWeight #calculate relative biomass contribtions based on PSMs
    GroupedResults = GroupedResults.sort_values(by=['RelativeWeight'],ascending=False) #sort by relative weight

    return GroupedResults


GroupedMasscontribution = DistributePeptides(args.PeptonizerResults,args.GraphdataFrame,args.PosteriorThreshold)
GroupedMasscontribution.to_csv(args.out,index=False)
    
    

