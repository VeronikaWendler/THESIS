# Veronika Wendler
# 08.01.25
# this program runs through the .mat files for EXP1 (experiential-symbolic (ES) phase) for the garcia replication, 
# participant 1-26 and creates panda frames and csv files of the variables stored in .mat 
# includes the creation of columns that match Garcia's code

import scipy.io
import pandas as pd
import os
import glob
import re
import numpy as np
# custom made functions
import toolbox
from toolbox import clean_column, rearrange_dataframe

# Directory containing all .mat files
mat_files_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Veronika_Thesis_Exp1_Garcia/Data/Data_28_11_24/Data/"
output_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Data_Garcia_CSV"
os.makedirs(output_dir, exist_ok=True)

mat_files = sorted(glob.glob(os.path.join(mat_files_dir, "Exp2_participant*_dataEXP2_*.mat")))

for file_path in mat_files:
    participant_match = re.search(r"participant(\d+)", file_path)
    if participant_match:
        participant_id = participant_match.group(1)  # Extract numeric ID
    else:
        raise ValueError(f"ID not in path: {file_path}")
    
    participant_folder = os.path.join(output_dir, f"Participant_{int(participant_id)}")
    os.makedirs(participant_folder, exist_ok=True)
    mat_data = scipy.io.loadmat(file_path)

    #'scoreEXP2', 'elapsed_timeEXP2', 'left_imagesEXP2', 'right_imagesEXP2', 'left_probabilityEXP2', 'right_probabilityEXP2', 'selected_probabilityEXP2', 'selected_imageEXP2', 'selected_sideEXP2', 'trialRunsEXP2', 'trialIconsEXP2', 'trialPiechartProbsEXP2', 'outcomeArrayEXP2', 'leftOutcomesEXP2', 'rightOutcomesEXP2', 'clicksArrayEXP2', 'utility_ArrayEXP2', 'scoreArrayEXP2'
    sub_id = [int(participant_id)] * len(mat_data['elapsed_timeEXP2'])
    elapsed_timeEXP2_vals = mat_data['elapsed_timeEXP2'].flatten()
    left_imagesEXP2_vals = mat_data['left_imagesEXP2'].flatten()
    left_probabilityEXP2_vals = mat_data['left_probabilityEXP2'].flatten()
    right_imagesEXP2_vals = mat_data['right_imagesEXP2'].flatten()
    right_probabilityEXP2_vals = mat_data['right_probabilityEXP2'].flatten()
    selected_probabilityEXP2_vals = mat_data['selected_probabilityEXP2'].flatten()
    selected_imageEXP2_vals = mat_data["selected_imageEXP2"].flatten()
    selected_sideEXP2_vals = mat_data["selected_sideEXP2"].flatten()
    trialRunsEXP2_vals = mat_data["trialRunsEXP2"].flatten()
    trialIconsEXP2_vals = mat_data["trialIconsEXP2"].flatten()
    trialPiechartProbsEXP2_vals = mat_data["trialPiechartProbsEXP2"].flatten()
    outcomeArrayEXP2_vals = mat_data['outcomeArrayEXP2'].flatten()
    leftOutcomesEXP2_vals = mat_data['leftOutcomesEXP2'].flatten()
    rightOutcomesEXP2_vals = mat_data['rightOutcomesEXP2'].flatten()
    clicksArrayEXP2_vals = mat_data['clicksArrayEXP2'].flatten()
    utility_ArrayEXP2_vals = mat_data['utility_ArrayEXP2'].flatten()
    scoreArrayEXP2_vals = mat_data['scoreArrayEXP2'].flatten()
    scoreEXP2_vals = mat_data["scoreEXP2"].flatten()
    if len(scoreEXP2_vals) == 1:
        scoreEXP2_vals = [scoreEXP2_vals[0]] * len(clicksArrayEXP2_vals)
        
    # cleaning all columns
    left_imagesEXP2_vals = clean_column(left_imagesEXP2_vals)
    left_probabilityEXP2_vals = clean_column(left_probabilityEXP2_vals)
    right_imagesEXP2_vals = clean_column(right_imagesEXP2_vals)
    right_probabilityEXP2_vals = clean_column(right_probabilityEXP2_vals)
    selected_probabilityEXP2_vals = clean_column(selected_probabilityEXP2_vals)
    selected_imageEXP2_vals = clean_column(selected_imageEXP2_vals)
    selected_sideEXP2_vals = clean_column(selected_sideEXP2_vals)
    trialRunsEXP2_vals = clean_column(trialRunsEXP2_vals)
    trialIconsEXP2_vals = clean_column(trialIconsEXP2_vals)
    trialPiechartProbsEXP2_vals = clean_column(trialPiechartProbsEXP2_vals)
        
    df = pd.DataFrame({
        "SubID": sub_id,
        "elapsed_timeEXP2": elapsed_timeEXP2_vals,
        "left_imagesEXP2": left_imagesEXP2_vals,
        "left_probabilityEXP2": left_probabilityEXP2_vals,
        "right_imagesEXP2": right_imagesEXP2_vals,
        "right_probabilityEXP2": right_probabilityEXP2_vals,
        "selected_probabilityEXP2": selected_probabilityEXP2_vals,
        "selected_imageEXP2": selected_imageEXP2_vals,
        "selected_sideEXP2": selected_sideEXP2_vals,
        "trialRunsEXP2": trialRunsEXP2_vals,
        "trialIconsEXP2": trialIconsEXP2_vals,
        "trialPiechartProbsEXP2": trialPiechartProbsEXP2_vals,
        "outcomeArrayEXP2": outcomeArrayEXP2_vals ,
        "leftOutcomesEXP2": leftOutcomesEXP2_vals,
        "rightOutcomesEXP2": rightOutcomesEXP2_vals,
        "clicksArrayEXP2": clicksArrayEXP2_vals,
        "utility_ArrayEXP2": utility_ArrayEXP2_vals,
        "scoreArrayEXP2": scoreArrayEXP2_vals,
        "scoreEXP2": scoreEXP2_vals,
    })
    
    # creating columns Garcia has
    # ---------------------------------------------------------------
    # ev1 and ev2:
    # ev(1) = p1 * (+1) + (1 - p1) * (-1) => 2*p1 - 1
    # ev(2) = p2 * (+1) + (1 - p2) * (-1) => 2*p2 - 1
    # ---------------------------------------------------------------
    
    df["ev1"] = (2.0 * df["left_probabilityEXP2"] - 1.0).round(1)
    df["ev2"] = (2.0 * df["right_probabilityEXP2"] - 1.0).round(1)
    
    df["p1"] = df["left_probabilityEXP2"].round(1)
    df["p2"] = df["right_probabilityEXP2"].round(1)
    
    df["cho"] = df["clicksArrayEXP2"].apply(lambda x: 1 if x == 1 else 2)

    #1 if (chosen EV >= unchosen EV), else 0
    # ---------------------------------------------------------------
    df["corr"] = np.where(
        (
            ((df["cho"] == 1) & (df["ev1"] >= df["ev2"]))  # left chosen & ev1 >= ev2
            | 
            ((df["cho"] == 2) & (df["ev2"] >= df["ev1"]))  # right chosen & ev2 >= ev1
        ),
        1,  # True
        0   # False
    )
   
    #trial col
    df["trial"] = np.arange(1, len(df) + 1)
    # cond col
    df["cond"] = -1
    df["out"] = df["utility_ArrayEXP2"]
    df["cfout"] = -1
    # phase column, constant "LE" for each row
    df["phase"] = "ES"
    # sess column: constant 1 for each row
    df["sess"] = 0
    # adding `op1` and `op2` based on image type
    df["op1"] = df["left_imagesEXP2"].apply(lambda x: 'S' if 'Pie' in x else 'E')
    df["op2"] = df["right_imagesEXP2"].apply(lambda x: 'S' if 'Pie' in x else 'E')
    # rew column, copy of scoreEXP
    df["rew"] = df["scoreArrayEXP2"]
    # chose_right column: 1 if clicksArray1EXP == 0, else 0
    #   Because in clicksArray1EXP: 0 = right, 1 = left
    df["chose_right"] = df["clicksArrayEXP2"].apply(lambda x: 1 if x == 0 else 0)
    df["rtime"] = df["elapsed_timeEXP2"]
    df = rearrange_dataframe(df)
    csv_file_trials = os.path.join(participant_folder, f"EXP2_garcia_participant_{int(participant_id)}_re.csv")
    df.to_csv(csv_file_trials, index=False)
    print(f"Processed participant {participant_id}, files saved in {participant_folder}")
