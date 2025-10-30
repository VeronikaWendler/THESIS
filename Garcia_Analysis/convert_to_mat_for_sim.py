# Veronika Wendler
# converts the csv file with all participants into and additional .mat file
import pandas as pd
import scipy.io as sio

# Define the file paths
csv_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
mat_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/Garcia_Eye_for_Simulation.mat"

# Load the CSV file into a Pandas DataFrame
df = pd.read_csv(csv_file, sep=',') 

# Convert DataFrame columns to a dictionary of lists (MATLAB-compatible table structure)
mat_table = {col: df[col].values for col in df.columns}

# Save to .mat file as a MATLAB table
sio.savemat(mat_file, {"TB": {"Garcia": mat_table}})
print(f"MAT file saved successfully at: {mat_file}")
