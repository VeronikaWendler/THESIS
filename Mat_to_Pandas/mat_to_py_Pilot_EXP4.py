# Veronika Wendler
# 10.12.24


#Conversion code .mat -> pandas data frame -> csv for the Slider Phase (EXP4) of Garcia's direct replication
# includes the creation of columns that match Garcia's code

# importing libraries
import scipy.io
import pandas as pd
import os
import glob
import re
import numpy as np
import toolbox
from toolbox import clean_column, parse_nested_list, clean_numeric_column, clean_categorical_column

# Directory containing all .mat files
mat_files_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Pilot_Analysis/Data_Collection_Start_VW_Laptop_20_10_23"
info_files_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Pilot_Analysis/Data_Pilot_CSV"
output_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Pilot_Analysis/Data_Pilot_CSV"
os.makedirs(output_dir, exist_ok=True)

mat_files = sorted(glob.glob(os.path.join(mat_files_dir, "participant*_SliderResponsesEXP_*.mat")))

for file_path in mat_files:
    participant_match = re.search(r"participant(\d+)", file_path)
    if participant_match:
        participant_id = participant_match.group(1)  # Extract numeric ID
    else:
        raise ValueError(f"Participant ID not found in file path: {file_path}")
    
    participant_folder = os.path.join(output_dir, f"Participant_{int(participant_id)}")
    os.makedirs(participant_folder, exist_ok=True)
    mat_data = scipy.io.loadmat(file_path)
    sub_id = [int(participant_id)] * len(mat_data['confidenceLevelsArrayEXP'])
    confidenceLevelsArrayEXP_vals = mat_data['confidenceLevelsArrayEXP'].flatten()
    selectedImageNamesArrayEXP_vals = mat_data['selectedImageNamesArrayEXP'].flatten()
    sliderResponsesArrayEXP_vals = mat_data['sliderResponsesArrayEXP'].flatten()

    # cleaning columns
    selectedImageNamesArrayEXP_vals = clean_categorical_column(selectedImageNamesArrayEXP_vals)
    sliderResponsesArrayEXP_vals = clean_numeric_column(sliderResponsesArrayEXP_vals)
   
    # info file
    info_file_path = os.path.join(info_files_dir, f"Participant_{int(participant_id)}", f"EXP1_Pilot_participant_{int(participant_id)}_info.csv")
    info_df = pd.read_csv(info_file_path)

    # parsing the `selectedImagesFromEXP` and `selectedProbabilitiesFromEXP`
    selected_images = parse_nested_list(info_df["selectedImagesFromEXP"].iloc[0])
    selected_probabilities = parse_nested_list(info_df["selectedProbabilitiesFromEXP"].iloc[0])

    # dictionaries for the mapping
    image_prob_map = {}
    for images, probs in zip(selected_images, selected_probabilities):
        for img, prob in zip(images, probs):
            if isinstance(img, list):
                img = img[0] 
            image_prob_map[img] = prob

    # 'p1'
    p1_values = []
    for image_name in selectedImageNamesArrayEXP_vals:
        if "Pie" in image_name:  
            match = re.search(r"(\d)$", image_name)   # get last digit -> it is the probability of the pie chart, the first digit shows the pie chart type
            if match:
                p1 = int(match.group(1)) / 10.0  # Convert last digit to probability
            else:
                p1 = np.nan 
        else:  
            p1 = image_prob_map.get(image_name, np.nan) 
        p1_values.append(p1)

    df = pd.DataFrame({
        "SubID": sub_id,
        "confidenceLevelsArrayEXP": confidenceLevelsArrayEXP_vals,
        "selectedImageNamesArrayEXP": selectedImageNamesArrayEXP_vals,
        "sliderResponsesArrayEXP": sliderResponsesArrayEXP_vals,
        "p1": p1_values, 
    })
    
    df['phase'] = "SP"
    df["p2"] = np.nan
    df["rtime"] = 0
    df["out"] = 1
    df["cfout"] = -1
    df["cho"] = df["sliderResponsesArrayEXP"] / 100
    df["corr"] = np.where(df["cho"] == df["p1"], 1, 0)
    df["trial"] = np.arange(1, len(df) + 1)
    df["cond"] = np.nan
    df["chose_right"] = 0
    df["rew"] = 0
    df["sess"] = 0
    df["op1"] = df["selectedImageNamesArrayEXP"].apply(lambda x: "S" if "Pie" in x else "E")
    df["op2"] = ""
    df["ev1"] = (2.0 * df["p1"] - 1.0).round(1)
    df["ev2"] = np.nan
    df["catch_trial"] = 0
    df["reversed"] = 0
    csv_file_trials = os.path.join(participant_folder, f"EXP4_Pilot_participant_{int(participant_id)}.csv")
    df.to_csv(csv_file_trials, index=False)
    print(f"Participant {participant_id} files saved in {participant_folder}")

