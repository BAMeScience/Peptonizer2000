import pandas as pd
import os


# find all benchmark subfolders in the results directory
def find_benchmark_subfolders(results_dir):

    subfolders = [f.path for f in os.scandir(results_dir) if f.is_dir()]
    benchmark_subfolders = []
    for subfolder in subfolders:
        subfolders2 = [f.path for f in os.scandir(subfolder) if f.is_dir()]
        for subfolder3 in subfolders2:
            subfolders4 = [f.path for f in os.scandir(subfolder3) if f.is_dir()]
            for subfolder5 in subfolders4:
                if 'benchmark' in subfolder5:
                    benchmark_subfolders.append(subfolder5)
    return benchmark_subfolders

folder_list = find_benchmark_subfolders('results')

# create a dataframe with the benchmark results
df_pepGMrun = pd.DataFrame()
df_graph_contruction = pd.DataFrame()
for folder in folder_list:
    #get all items in folder
    folder_items = [f.path for f in os.scandir(folder)]

    for item in folder_items:
        if os.path.isfile(item):
            df_temp = pd.read_csv(item,sep = '\t')
            df_graph = df_graph_contruction.append(df_temp)
        if os.path.isdir(item):
            files = [f.path for f in os.scandir(item)]
            for file in files:
                df_temp = pd.read_csv(file, sep='\t')
                df_pepGMrun = df_pepGMrun.append(df_temp)

#get min, max and average of all benchmark elements
#for the Peptonizer runs, per column
# Get the minimum, maximum, and average for each column
min_values = df_pepGMrun.min()
max_values = df_pepGMrun.max()
df_pepGMrun['h:m:s'] = pd.to_timedelta(df_pepGMrun['h:m:s'])
average_values = df_pepGMrun.mean()


# Create a new DataFrame to store the results
summary_df = pd.DataFrame({
    'Column': df_pepGMrun.columns,
    'Min': min_values.values,
    'Max': max_values.values,
    'Average': average_values.values
})

#svae results
summary_df.to_csv('results/summary_benchmark_pepGMrun.csv', sep='\t', index=False)

# Get the minimum, maximum, and average for each column of the graph construction
min_values = df_graph.min()
max_values = df_graph.max()
df_graph['h:m:s'] = pd.to_timedelta(df_graph['h:m:s'])
average_values = df_graph.mean()

# Create a new DataFrame to store the results
summary_df = pd.DataFrame({
    'Column': df_graph.columns,
    'Min': min_values.values,
    'Max': max_values.values,
    'Average': average_values.values
})

#svae results
summary_df.to_csv('results/summary_benchmark_graph_construction.csv', sep='\t', index=False)









#for the graph contruction runs, per column




