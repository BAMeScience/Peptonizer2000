import pandas as pd
import argparse
import matplotlib.pyplot as plt
from ete3 import NCBITaxa
import colorcet as cc
import seaborn as sns

parser = argparse.ArgumentParser(description = 'plot relative weight by higher taxa')
parser.add_argument('--GroupedMassResults', type = str, help = 'path to grouped mass results file')
parser.add_argument('--out', type = str, help = 'path to output file')

args = parser.parse_args()

ncbi = NCBITaxa()

Biomasses = pd.read_csv(args.GroupedMassResults)
Biomasses.drop('peptideWeight', axis=1, inplace=True) # drop the 'PeptideWeights' column
Biomasses.set_index('HigherTaxa', inplace=True) # set HigherTaxa as index
Biomasses.index = Biomasses.index.astype(str) # convert index to string

TaxaNameDict = ncbi.get_taxid_translator(Biomasses.index) # translate taxids to scientific names
TaxaNameDict = {str(key):value for key, value in TaxaNameDict.items()} # convert keys to strings

# Convert the index of the Biomasses dataframe according to the TaxaNameDict dictionary
Biomasses.index = Biomasses.index.map(TaxaNameDict)


Biomasses = Biomasses.T # transpose the dataframe


# Create a horizontal stacked bar chart

palette = sns.color_palette(cc.glasbey, n_colors=25)
Biomasses.plot(kind='barh', stacked=True, color = palette, figsize=(15, 10))

#place legend at bottom of figure, increase text size 
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=5, fontsize=12)

ax = plt.gca() # get current axis
ax.yaxis.set_tick_params(labelleft=False) # remove y-axis labels
ax.set_yticks([])

# Add labels and title
plt.xlabel('Relative Number of PSMs',fontsize = 14)



figure = plt.gcf() # get current figure
figure.tight_layout() # make sure the figure is not cut off

# Show the plot
plt.savefig(args.out)

