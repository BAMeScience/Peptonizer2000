import pandas as pd

# Load the CSV file
df = pd.read_csv('/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s07_root_new_rescore/CAMPI_SIHUMIx/PeptonizerResults.csv')

# Keep only the 'ID' and 'score' columns and rename them
df = df[['ID', 'score']].rename(columns={'ID': 'taxon_id', 'score': 'probability'})
# Save the modified dataframe to a new CSV file
df.to_csv('/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s07_root_new_rescore/CAMPI_SIHUMIx/PeptonizerResults_reshaped.csv', index=False)
