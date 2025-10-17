# Fixation statistics for simulation (similar to Dr. Chih-Chung Ting's work)
# libraries
import pandas as pd
import numpy as np
import os

output_dir = r'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Fixation_Stats_for_Simulation_ES'
os.makedirs(output_dir, exist_ok=True)

# preprocess data
def load_and_preprocess_data(file_path, phase_filter, exclude_subjects):
    data = pd.read_csv(file_path)
    data = data[data["phase"] == phase_filter]
    data = data[~data['sub_id'].isin(exclude_subjects)]
    numeric_cols = ['p1', 'p2', 'FirstFixLoc', 'FirstFixDur', 'FinalFixDur', 'eachMiddleFixDur']
    for col in numeric_cols:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    ov_order = ['low', 'medium', 'high']
    data['OVcate_2'] = pd.Categorical(data['OVcate_2'], categories=ov_order, ordered=True)
    data = data.dropna(subset=['p1', 'p2', 'FirstFixLoc', 'FirstFixDur', 'FinalFixDur', 'eachMiddleFixDur', 'OVcate_2'])
    return data

# fixation probabilities and durations with log-normal transformation
def compute_fixation_statistics(data):
    results = []
    for ov_level in data['OVcate_2'].cat.categories:
        subset = data[data['OVcate_2'] == ov_level]
        vr_greater_subset = subset[subset['p2'] > subset['p1']]
        vl_greater_subset = subset[subset['p2'] < subset['p1']]
        
        p_fix_right_vr_greater = (vr_greater_subset['FirstFixLoc'] == 2).mean() if not vr_greater_subset.empty else np.nan
        p_fix_right_vl_greater = (vl_greater_subset['FirstFixLoc'] == 2).mean() if not vl_greater_subset.empty else np.nan
        
        sd_fix_right_vr_greater = (vr_greater_subset['FirstFixLoc'] == 2).std() if not vr_greater_subset.empty else np.nan
        sd_fix_right_vl_greater = (vl_greater_subset['FirstFixLoc'] == 2).std() if not vl_greater_subset.empty else np.nan
        
        if not np.isnan(p_fix_right_vr_greater) and not np.isnan(p_fix_right_vl_greater):
            total_prob = p_fix_right_vr_greater + p_fix_right_vl_greater
            p_fix_right_vr_greater /= total_prob
            p_fix_right_vl_greater /= total_prob
        
        mean_first_fix = subset['FirstFixDur'].mean()
        sd_first_fix = subset['FirstFixDur'].std()
        
        mean_remaining_fix = (subset['FinalFixDur'] + subset['eachMiddleFixDur']).mean()
        sd_remaining_fix = (subset['FinalFixDur'] + subset['eachMiddleFixDur']).std()
        
        # Log-normal transformation
        mu_first_fix = np.log(mean_first_fix)
        sigma_first_fix = np.sqrt(np.log(1 + (sd_first_fix / mean_first_fix) ** 2))
        
        mu_remaining_fix = np.log(mean_remaining_fix)
        sigma_remaining_fix = np.sqrt(np.log(1 + (sd_remaining_fix / mean_remaining_fix) ** 2))
        
        results.append({
            'OVcate_2': ov_level,
            'P_FixRight_VrGreater': p_fix_right_vr_greater,
            'SD_FixRight_VrGreater': sd_fix_right_vr_greater,
            'P_FixRight_VlGreater': p_fix_right_vl_greater,
            'SD_FixRight_VlGreater': sd_fix_right_vl_greater,
            'Mean_FirstFixDur': mean_first_fix,
            'SD_FirstFixDur': sd_first_fix,
            'Mu_FirstFixDur': mu_first_fix,
            'Sigma_FirstFixDur': sigma_first_fix,
            'Mean_RemainingFixDur': mean_remaining_fix,
            'SD_RemainingFixDur': sd_remaining_fix,
            'Mu_RemainingFixDur': mu_remaining_fix,
            'Sigma_RemainingFixDur': sigma_remaining_fix
        })
    fixation_stats_df = pd.DataFrame(results)
    return fixation_stats_df

# Main
file_path = r'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv'
exclude_subjects = [1, 4, 5, 6, 14, 99] 
phase_filter = 'ES'

data = load_and_preprocess_data(file_path, phase_filter, exclude_subjects)
fixation_stats = compute_fixation_statistics(data)
output_file = os.path.join(output_dir, 'Fixation_Statistics_LogNorm_mean.csv')
fixation_stats.to_csv(output_file, index=False)
print(f"Saved fixation statistics (with log-normal) to: {output_file}")
print("\n🔹 Fixation Statistics by OV Level:")
print(fixation_stats)
