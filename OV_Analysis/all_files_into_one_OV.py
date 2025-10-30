# Veronika Wendler
# OV combined participants + .mat file

import os
import pandas as pd
import numpy as np
import scipy.io

# Directory containing the Participant_1, Participant_2, ..., Participant_99 folders
data_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data"

# Participant IDs to process: 1 through 27
# Also there is no EXP3 and EXP4 for participant 14 and issues with EXP4 for participant 20
participant_ids = [i for i in list(range(1, 28)) if i not in [99]]      #14,20

# Columns to retain from each file
columns_to_keep = ['ev1',
                   'ev2',
                   'p1',
                   'p2',
                   'cho',
                   'corr',
                   'trial', 
                   'cond',
                   'out',
                   'cfout',
                   'phase',
                   'sess',
                   'op1',
                   'op2',
                   'rew', 
                   'chose_right',
                   'rtime']

# data types
column_dtypes = {
    'SubID': 'float64',
    'phase': 'object',
    'p1': 'float64',
    'p2': 'float64',
    'rtime': 'float64',
    'out': 'float64',
    'cfout': 'float64',
    'cho': 'float64',
    'corr': 'float64',
    'trial': 'float64',
    'cond': 'float64',
    'chose_right': 'int64',
    'rew': 'float64',
    'sess': 'float64',
    'op1': 'object',
    'op2': 'object',
    'ev1': 'float64',
    'ev2': 'float64',
}

# preparing a list to store all data frames
all_dfs = []

for participant_id in participant_ids:
    print(f"Processing Participant {participant_id}...")
    participant_folder = os.path.join(data_dir, f"sub-{participant_id:02d}/beh")
    participant_dfs = []

    # checking EXP1, EXP2 (_re.csv), EXP3, and EXP4 files in the folder
    for exp_num in ['EXP1', 'EXP2', 'EXP3', 'EXP4']:
        exp_file = None
        for f in os.listdir(participant_folder):
            if exp_num in f and f.endswith(".csv") and ("_info" not in f):
                if exp_num == 'EXP2' and not f.endswith("_re.csv"): #  we need the _re files as they fit Garcia's code, re = reorganized (E on the left, S on the right)
                    continue  
                exp_file = os.path.join(participant_folder, f)
                break

        if exp_file is None:
            print(f"No {exp_num} CSV found for participant {participant_id}, skipping {exp_num}...")
            continue

        df_exp = pd.read_csv(exp_file)
        df_exp["SubID"] = participant_id
        df_exp = df_exp[["SubID"] + columns_to_keep]
        # get the right data types
        for col, dtype in column_dtypes.items():
            if col in df_exp.columns:
                df_exp[col] = df_exp[col].astype(dtype, errors='ignore')
                
        participant_dfs.append(df_exp)

    if participant_dfs:
        participant_combined = pd.concat(participant_dfs, ignore_index=True)
        all_dfs.append(participant_combined)

df_combined = pd.concat(all_dfs, ignore_index=True)
df_combined['index'] = range(1, len(df_combined) + 1)
df_combined.rename(columns={'SubID': 'sub_id'}, inplace=True)
df_combined['catch_trial'] = 0
df_combined['reversed'] = 0
df_combined['sess'] = 0

desired_order = [
    'index',
    'sub_id', 
    'phase',
    'p1',
    'p2',
    'rtime',
    'out', 
    'cfout',
    'cho',
    'corr',
    'trial',
    'cond', 
    'chose_right',
    'rew', 
    'sess', 
    'op1',
    'op2',
    'ev1',
    'ev2',
    'catch_trial',
    'reversed'
]

df_combined = df_combined[desired_order]
output_csv = os.path.join(data_dir, "OVParticipants_EXP_1_2_3_4.csv")
df_combined.to_csv(output_csv, index=False)
print(f"CSV saved at: {output_csv}")

#.MAT file creation (I know this is probably useless but I wanted to create it so I can really fit all of Garcia's code, to see if I can replicate it)
# .mat file processing
df_combined['elic'] = df_combined['phase'].map({
    'LE': -1,
    'ES': 0,
    'EE': 0,
    'SP': 2
})
df_combined['op1_mat'] = df_combined['op1'].map({
    'S': 0,
    'E': 1
})
df_combined['op2_mat'] = df_combined['op2'].map({
    'S': 0,
    'E': 1
}).fillna(-1)
df_combined['cont1'] = df_combined['p1'].apply(
    lambda x: int(x * 10) if pd.notna(x) and x not in [1, -1] else (-1 if pd.isna(x) else x)
)
df_combined['cont2'] = df_combined['p2'].apply(
    lambda x: int(x * 10) if pd.notna(x) and x not in [1, -1] else (-1 if pd.isna(x) else x)
)

df_combined['dist'] = -1
df_combined['plot'] = -1
df_combined['prolific_id'] = np.nan
df_combined['dbtime'] = np.nan

mat_columns = {
    1: 'sub_id',
    2: 'prolific_id',
    3: 'elic',
    4: 'p1',
    5: 'p2',
    6: 'rtime',
    7: 'out',
    8: 'cfout',
    9: 'cho',
    10: 'corr',
    12: 'trial',
    13: 'cond',
    14: 'cont1',
    15: 'cont2',
    19: 'rew',
    20: 'sess',
    21: 'op1_mat',
    22: 'op2_mat',
    23: 'ev1',
    24: 'ev2',
    25: 'catch_trial',
    28: 'dist',
    29: 'plot',
    30: 'dbtime'
}

mat_df = pd.DataFrame(index=df_combined.index, columns=range(1, 31))
for col_idx, col_name in mat_columns.items():
    if col_name in df_combined:
        mat_df[col_idx] = df_combined[col_name]
    else:
        mat_df[col_idx] = np.nan  
# dictionary for .mat file saving        
mat_data = {'data': mat_df.to_numpy()}
output_mat = os.path.join(data_dir, "OVParticipants_EXP_1_2_3_4.mat")
scipy.io.savemat(output_mat, mat_data)
print(f"MAT file saved at: {output_mat}")

# Variable	Description
# index	Index identifier for each row
# sub_id	Subject identifier, identifying the participant
# phase	Phase of the experiment (LE, ES, EE, SP)
# p1	p(win) option 1
# p2	p(win) option 2
# rtime	Response time
# out	Outcome of the choice (selected option)
# cfout	Counterfactual outcome (unchosen option outcome)
# cho	Choice made by the participant (1 or 2)
# corr	Correctness of the choice (1 if expected value of chosen option >= expected value of unchosen option, 0 else)
# trial	Trial number
# cond	Condition of the experiment (3=60/40, 2=70/30, 1=80/20, 0=90/10)
# chose_right	Whether the participant chose the rightmost option on screen
# rew	Total reward received (cumulated)
# sess	Session identifier (When sess in (-1, -2), it means it was a training session, when sess is 0 or 1, it means first or second session)
# op1	Option 1 presented to the participant (filename or identifier)
# op2	Option 2 presented to the participant ((filename or identifier)
# ev1	Expected value of option 1
# ev2	Expected value of option 2
# catch_trial	Indicates if it's a catch trial (testing attention)
# reversed  /nothing is reversed




