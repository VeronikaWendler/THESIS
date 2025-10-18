# Veronika Wendler
# code creates the important, beahavioural conditions (columns) which can be used for drift diffusion modelling or other behviourl analyses

import pandas as pd 
import numpy as np 


# data
data = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/Pilot_Analysis/data/data_sets/PilotParticipants_EXP_1_2_4.csv")
    
# participants to exclude
exclude_part = {}
data['sub_id'] = data['sub_id'].astype(int)  

# Initialize OV, OV_category, and absolute value difference columns as NaN
data['OV_num'] = np.nan
data['OV_num_2'] = np.nan
data['OVcate'] = np.nan
data['OVcate_2'] = np.nan
data['AbsValueDiff'] = np.nan
data['AbsValueDiff_2'] = np.nan
data['Abscate'] = np.nan
data['Abscate_2'] = np.nan
data['VD'] = np.nan
data['VD_2'] = np.nan
data['RLdiff'] = np.nan
data['feedback'] = np.nan  
data['split_by'] = np.nan  
data['q_init'] = np.nan  
data['cond'] = pd.to_numeric(data['cond'], errors='coerce')
data['out'] = pd.to_numeric(data['out'], errors='coerce')  
data['stim_chosen'] = np.nan
 
    
data.loc[data['phase'] == 'LE', 'feedback'] = data.loc[data['phase'] == 'LE', 'out'].replace({-1: 0, 1: 1})
data['feedback'] = data['feedback'].astype(float) 

data.loc[data['phase'] == 'LE', 'split_by'] = data.loc[data['phase'] == 'LE', 'cond'].astype(int)
data.loc[data['phase'] == 'LE', 'q_init'] = 0.5
data['q_init'] = data['q_init'].astype(float)  
    
data['stim_chosen'] = np.where(data['phase'] == 'ES', 
                               data['chose_right'].apply(lambda x: 'E' if x == 0 else 'S'), 
                               np.nan)

# function to calculate OV and Absolute Value Difference for a specific phase
def calculate_ov_and_abs_value_diff(data, exclude_part, phase):
    unique_participants = data['sub_id'].unique()
    for participant in unique_participants:
        if participant in exclude_part:
            continue

        participant_data = data[(data['sub_id'] == participant) & (data['phase'] == phase)].copy()

        if participant_data.empty:
            continue

        # get OV and round to 1 decimal place
        participant_data['OV_num'] = (participant_data['p1'] + participant_data['p2']).round(1)
        participant_data['AbsValueDiff'] = (participant_data['p1'] - participant_data['p2']).abs().round(1)
        participant_data['RLdiff'] = (participant_data['p1'] - participant_data['p2']).round(1)

#################### Quartiles ################################################################################
        # thresholds for OV (25th and 75th percentiles) and round
        T1_OV, T2_OV = participant_data['OV_num'].quantile([0.25, 0.75]).round(1)

        def categorize_ov(ov):
            if ov <= T1_OV:
                return 'low'
            elif ov <= T2_OV:
                return 'medium'
            else:
                return 'high'

        participant_data['OVcate'] = participant_data['OV_num'].apply(categorize_ov)

        T1_AbsDiff, T2_AbsDiff = participant_data['AbsValueDiff'].quantile([0.25, 0.75]).round(1)

        def categorize_abs_diff(abs_diff):
            if abs_diff <= T1_AbsDiff:
                return 'low'
            elif abs_diff <= T2_AbsDiff:
                return 'medium'
            else:
                return 'high'

        participant_data['Abscate'] = participant_data['AbsValueDiff'].apply(categorize_abs_diff)
        category_mapping = {'low': 1, 'medium': 2, 'high': 3}
        participant_data['OV'] = participant_data['OVcate'].map(category_mapping)
        participant_data['VD'] = participant_data['Abscate'].map(category_mapping)
            
        data.loc[participant_data.index, ['OV', 'VD', 'OVcate', 'OV_num', 'AbsValueDiff', 'Abscate', 'RLdiff']] = participant_data[['OV', 'VD','OVcate','OV_num', 'AbsValueDiff', 'Abscate', 'RLdiff']]
    
    return data
            
phases = ['ES', 'LE', 'EE']
for phase in phases:
    data = calculate_ov_and_abs_value_diff(data, exclude_part, phase)
        

######################## TERTILES #################################################################################
def calculate_ov_and_abs_value_diff_tertiles(data, exclude_part, phase):
    unique_participants = data['sub_id'].unique()
    for participant in unique_participants:
        if participant in exclude_part:
            continue

        participant_data = data[(data['sub_id'] == participant) & (data['phase'] == phase)].copy()

        if participant_data.empty:
            continue

        participant_data['OV_num_2'] = (participant_data['p1'] + participant_data['p2']).round(1)
        participant_data['AbsValueDiff_2'] = (participant_data['p1'] - participant_data['p2']).abs().round(1)

        T1_OV_tertile, T2_OV_tertile = participant_data['OV_num_2'].quantile([0.33, 0.66]).round(1)

        def categorize_ov_tertile(ov):
            if ov <= T1_OV_tertile:
                return 'low'
            elif ov <= T2_OV_tertile:
                return 'medium'
            else:
                return 'high'

        participant_data['OVcate_2'] = participant_data['OV_num_2'].apply(categorize_ov_tertile)

        T1_AbsDiff_tertile, T2_AbsDiff_tertile = participant_data['AbsValueDiff_2'].quantile([0.33, 0.66]).round(1)

        def categorize_abs_diff_tertile(abs_diff):
            if abs_diff <= T1_AbsDiff_tertile:
                return 'low'
            elif abs_diff <= T2_AbsDiff_tertile:
                return 'medium'
            else:
                return 'high'

        participant_data['Abscate_2'] = participant_data['AbsValueDiff_2'].apply(categorize_abs_diff_tertile)

        category_mapping2 = {'low': 1, 'medium': 2, 'high': 3}
        participant_data['OV_2'] = participant_data['OVcate_2'].map(category_mapping2)
        participant_data['VD_2'] = participant_data['Abscate_2'].map(category_mapping2)

        data.loc[participant_data.index, ['OV_2', 'VD_2', 'OVcate_2', 'OV_num_2', 'AbsValueDiff_2', 'Abscate_2']] = participant_data[['OV_2', 'VD_2','OVcate_2','OV_num_2', 'AbsValueDiff_2', 'Abscate_2']]
    
    return data

phases = ['ES', 'LE']
for phase in phases:
    data = calculate_ov_and_abs_value_diff_tertiles(data, exclude_part, phase)
        
    
output_path = "D:/Aberdeen_Uni_June24/cap/THESIS/Pilot_Analysis/data/data_sets/PilotParticipants_ddm_OV_Abs.csv"
data.to_csv(output_path, index=False)

print(data[data['phase'].isin(phases)].head(20))
