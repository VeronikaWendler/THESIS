# # Veronika Wendler
# # 10.12.24

# #Automated conversion program
# # this program runs through the .mat files for EXP1 for the OV  experiment, 
# # participant 1-26 and creates panda frames and csv files of the varaibles stored in .mat 

# #Conversion code .mat -> pandas data frame -> csv


# importing libraries
import scipy.io
import pandas as pd
import os
import glob
import re
import numpy as np
from toolbox import clean_column, parse_nested_list, get_OV_condition

# Directory containing all .mat files
mat_files_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Veronika_Thesis_Exp2_OV/Data/Data_28_11_24/Data"
output_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/Data_OV_CSV"
os.makedirs(output_dir, exist_ok=True)

mat_files = sorted(glob.glob(os.path.join(mat_files_dir, "Exp2_participant*_dataEXP3_*.mat")))

for file_path in mat_files:
    participant_match = re.search(r"participant(\d+)", file_path)
    if participant_match:
        participant_id = participant_match.group(1)  
    else:
        raise ValueError(f"ID not in path: {file_path}")
    
    participant_folder = os.path.join(output_dir, f"Participant_{int(participant_id)}")
    os.makedirs(participant_folder, exist_ok=True)
    mat_data = scipy.io.loadmat(file_path)
    
    #'scoreEXP3', 'elapsed_timeEXP3', 'left_imagesEXP3', 'right_imagesEXP3', 'left_probabilityEXP3', 'right_probabilityEXP3', 'selected_probabilityEXP3', 'selected_imageEXP3', 'selected_sideEXP3', 'trialRunsEXP3', 'trialIcon1EXP3', 'trialIcon2EXP3', 'outcomeArrayEXP3', 'leftOutcomesEXP3', 'rightOutcomesEXP3', 'clicksArrayEXP3', 'utility_ArrayEXP3', 'scoreArrayEXP3')

    sub_id = [int(participant_id)] * len(mat_data['elapsed_timeEXP3'])
    elapsed_timeEXP3_vals = mat_data['elapsed_timeEXP3'].flatten()
    left_imagesEXP3_vals = mat_data['left_imagesEXP3'].flatten()
    left_probabilityEXP3_vals = mat_data['left_probabilityEXP3'].flatten()
    right_imagesEXP3_vals = mat_data['right_imagesEXP3'].flatten()
    right_probabilityEXP3_vals = mat_data['right_probabilityEXP3'].flatten()
    selected_probabilityEXP3_vals = mat_data['selected_probabilityEXP3'].flatten()
    selected_imageEXP3_vals = mat_data["selected_imageEXP3"].flatten()
    selected_sideEXP3_vals = mat_data["selected_sideEXP3"].flatten()
    trialRunsEXP3_vals = mat_data["trialRunsEXP3"].flatten()
    trialIcon1EXP3_vals = mat_data["trialIcon1EXP3"].flatten()
    trialIcon2EXP3_vals = mat_data["trialIcon2EXP3"].flatten()
    outcomeArrayEXP3_vals = mat_data['outcomeArrayEXP3'].flatten()
    leftOutcomesEXP3_vals = mat_data['leftOutcomesEXP3'].flatten()
    rightOutcomesEXP3_vals = mat_data['rightOutcomesEXP3'].flatten()
    clicksArrayEXP3_vals = mat_data['clicksArrayEXP3'].flatten()
    utility_ArrayEXP3_vals = mat_data['utility_ArrayEXP3'].flatten()
    scoreArrayEXP3_vals = mat_data['scoreArrayEXP3'].flatten()
    scoreEXP3_vals = mat_data["scoreEXP3"].flatten()
    if len(scoreEXP3_vals) == 1:
        scoreEXP3_vals = [scoreEXP3_vals[0]] * len(clicksArrayEXP3_vals)
    
    # cleaning
    left_imagesEXP3_vals = clean_column(left_imagesEXP3_vals)
    left_probabilityEXP3_vals = clean_column(left_probabilityEXP3_vals)
    right_imagesEXP3_vals = clean_column(right_imagesEXP3_vals)
    right_probabilityEXP3_vals = clean_column(right_probabilityEXP3_vals)
    selected_probabilityEXP3_vals = clean_column(selected_probabilityEXP3_vals)
    selected_imageEXP3_vals = clean_column(selected_imageEXP3_vals)
    selected_sideEXP3_vals = clean_column(selected_sideEXP3_vals)
    trialRunsEXP3_vals = clean_column(trialRunsEXP3_vals)
    trialIcon1EXP3_vals = clean_column(trialIcon1EXP3_vals)
    trialIcon2EXP3_vals = clean_column(trialIcon2EXP3_vals)
        
    #trial-level DataFrame
    df = pd.DataFrame({
        "SubID": sub_id,
        "elapsed_timeEXP3": elapsed_timeEXP3_vals,
        "left_imagesEXP3": left_imagesEXP3_vals,
        "left_probabilityEXP3": left_probabilityEXP3_vals,
        "right_imagesEXP3": right_imagesEXP3_vals,
        "right_probabilityEXP3": right_probabilityEXP3_vals,
        "selected_probabilityEXP3": selected_probabilityEXP3_vals,
        "selected_imageEXP3": selected_imageEXP3_vals,
        "selected_sideEXP3": selected_sideEXP3_vals,
        "trialRunsEXP3": trialRunsEXP3_vals,
        "trialIcon1EXP3": trialIcon1EXP3_vals,
        "trialIcon2EXP3": trialIcon2EXP3_vals,
        "outcomeArrayEXP3": outcomeArrayEXP3_vals ,
        "leftOutcomesEXP3": leftOutcomesEXP3_vals,
        "rightOutcomesEXP3": rightOutcomesEXP3_vals,
        "clicksArrayEXP3": clicksArrayEXP3_vals,
        "utility_ArrayEXP3": utility_ArrayEXP3_vals,
        "scoreArrayEXP3": scoreArrayEXP3_vals,
        "scoreEXP3": scoreEXP3_vals,

    })
    
    df["ev1"] = (2.0 * df["left_probabilityEXP3"] - 1.0).round(1)
    df["ev2"] = (2.0 * df["right_probabilityEXP3"] - 1.0).round(1)
    
    df["p1"] = df["left_probabilityEXP3"].round(1)
    df["p2"] = df["right_probabilityEXP3"].round(1)
    
    df["cho"] = df["clicksArrayEXP3"].apply(lambda x: 1 if x == 1 else 2)

    #1 if (chosen EV >= unchosen EV), else 0
    # ---------------------------------------------------------------
    
    df["corr"] = np.where(
        (
            ((df["cho"] == 1) & (df["ev1"] >= df["ev2"]))  # left chosen & ev1 >= ev2
            | 
            ((df["cho"] == 2) & (df["ev2"] >= df["ev1"]))  # right chosen & ev2 >= ev1
        ),
        1,  
        0   
    )
   
    #trial col
    df["trial"] = np.arange(1, len(df) + 1)
    # cond col
    df["cond"] = -1
    df["out"] = df["utility_ArrayEXP3"]
    df["cfout"] = -1
    
    # phase col
    df["phase"] = "EE"

    # sess col
    df["sess"] = 0
    df['op1'] = 'E'
    df['op2'] = 'E'

    # rew column, copy of scoreEXP
    df["rew"] = df["scoreArrayEXP3"]
    # chose_right column: 1 if clicksArray1EXP == 0, else 0
    #    Because in clicksArray1EXP: 0 = right, 1 = left
    df["chose_right"] = df["clicksArrayEXP3"].apply(lambda x: 1 if x == 0 else 0)
    df["rtime"] = df["elapsed_timeEXP3"]
    csv_file_trials = os.path.join(participant_folder, f"EXP3_OV_participant_{int(participant_id)}.csv")
    df.to_csv(csv_file_trials, index=False)
    print(f"participant {participant_id} files saved in {participant_folder}")