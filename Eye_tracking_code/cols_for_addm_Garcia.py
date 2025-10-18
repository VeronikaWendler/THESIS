# Veronika Wendler
# 22.01.25
# preparing columns important for the aDDM model in hdddm
# columns include PropDwell_opt, PropDwell_sub, ev_opt, ev_sub, OVcate, Abscate 


# libraries
import pandas as pd
import numpy as np


version = 2 # 1 runs first code, 2 runs second code, I don't know of some other idea but 1 needs to be proeccessed first


if version == 1:
    data = pd.read_csv('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_EXP_1_2_3_4.csv')  
    
    # modify the raw probabilities to have them on a scale, similar to CCT and Ian, i.e. 0.3 -> 30 
    data[['p1', 'p2']] = data[['p1', 'p2']] * 100
    data['ev_opt'] = data[['ev1', 'ev2']].max(axis=1)  # Optimal E(X)
    data['ev_sub'] = data[['ev1', 'ev2']].min(axis=1)  # Suboptimal E(X)

# ------RESPONSE eyetracking data (from start of stimuli presentation until response )------------------------------------------------------------------------------------

    data['DwellTime_opt'] = np.where(data['ev1'] == data['ev_opt'], data['DwellLeft'], data['DwellRight'])
    data['DwellTime_sub'] = np.where(data['ev1'] == data['ev_sub'], data['DwellLeft'], data['DwellRight'])
    # PropDwell
    data['PropDwell_opt'] = data['DwellTime_opt'] / data['DwellTotal']
    data['PropDwell_sub'] = data['DwellTime_sub'] / data['DwellTotal']
    #missing or zero dwell data
    data.loc[data['DwellTotal'] == 0, ['PropDwell_opt', 'PropDwell_sub']] = np.nan  
    
    # AttentionW & InattentionW computation
    # Compute V_corr (higher value) and V_sub (lower value)
    data['V_corr'] = data[['p1', 'p2']].max(axis=1)  
    data['V_sub'] = data[['p1', 'p2']].min(axis=1) 
    data['AttentionW'] = (data['PropDwell_opt'] * data['V_corr']) - (data['PropDwell_sub'] * data['V_sub'])
    data['InattentionW'] = (data['PropDwell_sub'] * data['V_corr']) - (data['PropDwell_opt'] * data['V_sub'])
    data['AttentionW'] = data['AttentionW'].round(3)
    data['InattentionW'] = data['InattentionW'].round(3)

    # better option based on the higher expected value for RESPONSE eyetracking data
    data['BetterOption'] = np.where(data['ev1'] > data['ev2'], 1, 2)  # 1 for left, 2 for right
    data['FixLocFirstCorr'] = (data['FirstFixLoc'] == data['BetterOption']).astype(int)

    #FixLocLastCorr: 1 if the final fixation was on the better option, 0 otherwise
    data['FixLocLastCorr'] = (data['FinalFixLoc'] == data['BetterOption']).astype(int)
    data['DwellTimeAdvantage'] = data['DwellRight'] - data['DwellLeft']
    data['PropDwell_Right'] = data['DwellRight'] / data['DwellTotal']
    data['PropDwell_Left'] = data['DwellLeft'] / data['DwellTotal']

    data.loc[data['DwellTotal'] == 0, ['PropDwell_Right', 'PropDwell_Left']] = np.nan  
    data['DwellPropAdvantage'] = data['PropDwell_Right'] - data['PropDwell_Left']

    # right-left
    data['FixationAdvantage'] = data['RightFixNR'] - data['LeftFixNR']
    
    data['Fix_opt'] = np.where(data['ev1'] == data['ev_opt'], data['LeftFixNR'], data['RightFixNR'])
    data['Fix_sub'] = np.where(data['ev1'] == data['ev_sub'], data['LeftFixNR'], data['RightFixNR'])
    # better - worse
    data['FixationAdvantage_Corr'] = data['Fix_opt'] - data['Fix_sub']
    
    
    data['DwelltimeAdvantageCorrect'] = data['DwellTime_opt'] - data['DwellTime_sub']
    data['DwellPropAdvantageCorrect'] = data['PropDwell_opt'] - data['PropDwell_sub']


    
#-------FEEDBACK eyetracking columns -----------------------------------------------------------------------------------
    
    data['DwellTime_opt_Feed'] = np.where(data['ev1'] == data['ev_opt'], data['DwellLeft_Feed'], data['DwellRight_Feed'])
    data['DwellTime_sub_Feed'] = np.where(data['ev1'] == data['ev_sub'], data['DwellLeft_Feed'], data['DwellRight_Feed'])
    # PropDwell
    data['PropDwell_opt_Feed'] = data['DwellTime_opt_Feed'] / data['DwellTotal_Feed']
    data['PropDwell_sub_Feed'] = data['DwellTime_sub_Feed'] / data['DwellTotal_Feed']
    #missing or zero dwell data
    data.loc[data['DwellTotal_Feed'] == 0, ['PropDwell_opt_Feed', 'PropDwell_sub_Feed']] = np.nan  
    
    data['AttentionW_Feed'] = (data['PropDwell_opt_Feed'] * data['V_corr']) - (data['PropDwell_sub_Feed'] * data['V_sub'])
    #  Compute InattentionW
    data['InattentionW_Feed'] = (data['PropDwell_sub_Feed'] * data['V_corr']) - (data['PropDwell_opt_Feed'] * data['V_sub'])
    data['AttentionW_Feed'] = data['AttentionW_Feed'].round(3)
    data['InattentionW_Feed'] = data['InattentionW_Feed'].round(3)
    data['BetterOption_Feed'] = np.where(data['ev1'] > data['ev2'], 1, 2)  # 1 for left, 2 for right
    data['FixLocFirstCorr_Feed'] = (data['FirstFixLoc_Feed'] == data['BetterOption_Feed']).astype(int)
    #FixLocLastCorr: 1 if the final fixation was on the better option, 0 otherwise
    data['FixLocLastCorr_Feed'] = (data['FinalFixLoc_Feed'] == data['BetterOption_Feed']).astype(int)
    data['DwellTimeAdvantage_Feed'] = data['DwellRight_Feed'] - data['DwellLeft_Feed']


#-------ALLFIX eyetracking data (all fixations in the trial; no restrictions)-----------------------------------------------    
    
    # optimal option for FEEDBACK eyetracking data
    data['DwellTime_opt_allfix'] = np.where(data['ev1'] == data['ev_opt'], data['DwellLeft_allfix'], data['DwellRight_allfix'])
    data['DwellTime_sub_allfix'] = np.where(data['ev1'] == data['ev_sub'], data['DwellLeft_allfix'], data['DwellRight_allfix'])
    # PropDwell
    data['PropDwell_opt_allfix'] = data['DwellTime_opt_allfix'] / data['DwellTotal_allfix']
    data['PropDwell_sub_allfix'] = data['DwellTime_sub_allfix'] / data['DwellTotal_allfix']
    data.loc[data['DwellTotal_allfix'] == 0, ['PropDwell_opt_allfix', 'PropDwell_sub_allfix']] = np.nan  
    data['AttentionW_allfix'] = (data['PropDwell_opt_allfix'] * data['V_corr']) - (data['PropDwell_sub_allfix'] * data['V_sub'])
    data['InattentionW_allfix'] = (data['PropDwell_sub_allfix'] * data['V_corr']) - (data['PropDwell_opt_allfix'] * data['V_sub'])
    data['AttentionW_allfix'] = data['AttentionW_allfix'].round(3)
    data['InattentionW_allfix'] = data['InattentionW_allfix'].round(3)
    data['BetterOption_allfix'] = np.where(data['ev1'] > data['ev2'], 1, 2)  # 1 for left, 2 for right
    data['FixLocFirstCorr_allfix'] = (data['FirstFixLoc_allfix'] == data['BetterOption_allfix']).astype(int)
    #FixLocLastCorr: 1 if the final fixation was on the better option, 0 otherwise
    data['FixLocLastCorr_allfix'] = (data['FinalFixLoc_allfix'] == data['BetterOption_allfix']).astype(int)
    data['DwellTimeAdvantage_allfix'] = data['DwellRight_allfix'] - data['DwellLeft_allfix']

    
    print(data[['p1', 'p2', 'V_corr', 'V_sub', 'PropDwell_Left', 'PropDwell_Right',
                'DwellTime_opt', 'DwellTime_sub','PropDwell_opt', 'PropDwell_sub', 'AttentionW', 'InattentionW', 'BetterOption','FixLocFirstCorr', 'FixLocLastCorr',
                'DwellTime_opt_Feed', 'DwellTime_sub_Feed','PropDwell_opt_Feed', 'PropDwell_sub_Feed', 'AttentionW_Feed', 'InattentionW_Feed', 'BetterOption_Feed','FixLocFirstCorr_Feed', 'FixLocLastCorr_Feed',
                'DwellTime_opt_allfix', 'DwellTime_sub_allfix','PropDwell_opt_allfix', 'PropDwell_sub_allfix', 'AttentionW_allfix', 'InattentionW_allfix', 'BetterOption_allfix','FixLocFirstCorr_allfix', 'FixLocLastCorr_allfix',
                'DwellTimeAdvantage','DwellPropAdvantage', 'DwellTimeAdvantage_Feed', 'DwellTimeAdvantage_allfix', 'FixationAdvantage', 'DwelltimeAdvantageCorrect',
                'Fix_opt', 'Fix_sub','FixationAdvantage_Corr']])

    data.to_csv('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_EXP_1_2_3_4_CCT.csv', index=False)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
elif version == 2:
    data = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_EXP_1_2_3_4_CCT.csv")
    # participants to exclude
    exclude_part = {4, 5, 6, 14}
    data['sub_id'] = data['sub_id'].astype(int)  
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

    # 'split_by' column is based on 'cond', converting to integers only where 'phase' == 'LE'
    data.loc[data['phase'] == 'LE', 'split_by'] = data.loc[data['phase'] == 'LE', 'cond'].astype(int)

    # 'q_init' column is where 'phase' == 'LE'
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

            # thresholds for absolute value difference (25th and 75th percentiles) and round
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

            # get thresholds (33rd and 66th percentiles) for OV
            T1_OV_tertile, T2_OV_tertile = participant_data['OV_num_2'].quantile([0.33, 0.66]).round(1)

            # terciles
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

    phases = ['ES', 'LE', 'EE']
    for phase in phases:
        data = calculate_ov_and_abs_value_diff_tertiles(data, exclude_part, phase)
        
    
    output_path = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
    data.to_csv(output_path, index=False)

    print(data[data['phase'].isin(phases)].head(20))
    
