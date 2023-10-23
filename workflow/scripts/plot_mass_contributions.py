import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = 'plot relative weight by higher taxa')
parser.add_argument('--GroupedMassResults', type = str, help = 'path to grouped mass results file')
parser.add_argument('--out', type = str, help = 'path to output file')

args = parser.parse_args()

Biomasses = pd.read_csv(args.GroupedMassResults)

# Create a horizontal stacked bar chart
Biomasses.plot(kind='barh', x='relativeweight', y='highertaxa', stacked=True)

# Add labels and title
plt.xlabel('Relative Weight')
plt.ylabel('Higher Taxa')
plt.title('Relative Weight by Higher Taxa')

# Show the plot
plt.savefig(args.out)
