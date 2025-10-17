# Veronika Wendler
# 08.01.25
#Conversion code .mat -> pandas data frame -> csv for the Slider phase (Slider) of Garcia's direct replication
# includes the creation of columns that match Garcia's code

#importing libraries
import os
import pandas as pd
import numpy as np
import re
import scipy.io
from toolbox import get_OV_condition

# Directory containing all .mat files for OV
mat_files_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/Data_OV_CSV"
output_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/Data_OV_CSV" 
os.makedirs(output_dir, exist_ok= True)

participant_accuracies = []

for participant_folder in os.listdir(mat_files_dir):
    participant_path = os.path.join(mat_files_dir, participant_folder)
    if os.path.isdir(participant_path) and participant_folder.startswith("Participant_"):
        participant_id = re.search(r"Participant_(\d+)", participant_folder).group(1)
        csv_file_path = os.path.join(participant_path, f"EXP1_OV_participant_{participant_id}.csv")
        if os.path.exists(csv_file_path):
            df = pd.read_csv(csv_file_path)

            # new columns
            df["ev1"] = (2.0 * df["left_ProbabilityEXP"] - 1.0).round(1)
            df["ev2"] = (2.0 * df["right_ProbabilityEXP"] - 1.0).round(1)
            df["p1"] = df["left_ProbabilityEXP"]
            df["p2"] = df["right_ProbabilityEXP"]
            df["cho"] = df["clicksArray1EXP"].apply(lambda x: 1 if x == 1 else 2)
            df["corr"] = np.where(
                (
                    ((df["cho"] == 1) & (df["ev1"] >= df["ev2"]))
                    | ((df["cho"] == 2) & (df["ev2"] >= df["ev1"]))
                ),
                1,
                0,
            )
            df["trial"] = np.arange(1, len(df) + 1)
            df["cond"] = df.apply(lambda row: get_OV_condition(row["p1"], row["p2"]), axis=1)
            df["out"] = np.where(df["cho"] == 1, df["left_OutArrayEXP"], df["right_OutArrayEXP"])
            df["cfout"] = np.where(df["cho"] == 1, df["right_OutArrayEXP"], df["left_OutArrayEXP"])
            df["phase"] = "LE"
            df["sess"] = 0
            df["op1"] = "E"
            df["op2"] = "E"
            df["rew"] = df["scoreArrayEXP"]
            df["chose_right"] = df["clicksArray1EXP"].apply(lambda x: 1 if x == 0 else 0)
            df["rtime"] = df["elapsed_timeEXP"]
            accuracy = df["corr"].mean()
            participant_accuracies.append({"participant_id": int(participant_id), "accuracy": accuracy})

            updated_csv_file_path = os.path.join(participant_path, f"EXP1_OV_participant_{participant_id}.csv")
            df.to_csv(updated_csv_file_path, index=False)
            print(f"Updated CSV saved for participant {participant_id} at {updated_csv_file_path}")

# accuracy DataFrame to calculate median
accuracy_df = pd.DataFrame(participant_accuracies)
median_accuracy = accuracy_df["accuracy"].median()

# adding learner classification
for participant_folder in os.listdir(mat_files_dir):
    participant_path = os.path.join(mat_files_dir, participant_folder)
    if os.path.isdir(participant_path) and participant_folder.startswith("Participant_"):
        participant_id_match = re.search(r"Participant_(\d+)", participant_folder)
        if not participant_id_match:
            continue

        participant_id = participant_id_match.group(1)
        csv_file_trials = os.path.join(participant_path, f"EXP1_OV_participant_{participant_id}.csv")

        if os.path.isfile(csv_file_trials):
            df = pd.read_csv(csv_file_trials)

            participant_accuracy = accuracy_df[accuracy_df["participant_id"] == int(participant_id)]["accuracy"].values[0]
            learner_type = "good" if participant_accuracy >= median_accuracy else "poor"

            df["learners"] = learner_type
            final_csv_path = os.path.join(participant_path, f"EXP1_OV_participant_{participant_id}_with_learners.csv")
            df.to_csv(final_csv_path, index=False)
            print(f"CSV with learner saved for participant {participant_id} at {final_csv_path}")

