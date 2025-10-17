# Veronika Wendler
#12.12.24
# Automated conversion program
# this program runs through the .mat files for EXP1 (reinforcement learning (LE) phase) for the garcia replication, 
# participant 1-26 and creates panda frames and csv files of the variables stored in .mat 
# includes the creation of columns that match Garcia's code

# libraries
import scipy.io
import pandas as pd
import os
import glob
import re
import numpy as np
from toolbox import clean_categorical_column, get_condition_garcia

# Directory containing all .mat files for the Garcia replication
mat_files_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Veronika_Thesis_Exp1_Garcia/Data/Data_28_11_24/Data/"
output_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Data_Garcia_CSV"
os.makedirs(output_dir, exist_ok=True)

# list of all EXP1 .mat files for participants 1 to 26
mat_files = sorted(glob.glob(os.path.join(mat_files_dir, "Exp2_participant*_dataEXP1_*.mat")))


for file_path in mat_files:
    participant_match = re.search(r"participant(\d+)", file_path)
    if participant_match:
        participant_id = participant_match.group(1)  # Extract numeric ID
    else:
        raise ValueError(f"ID not in path: {file_path}")
 
    participant_folder = os.path.join(output_dir, f"Participant_{int(participant_id)}") # creates the participant folder
    os.makedirs(participant_folder, exist_ok=True)
    mat_data = scipy.io.loadmat(file_path)

    # creating the columns
    sub_id = [int(participant_id)] * len(mat_data['current_subblogEXP'])
    clicksArray1EXP_vals = mat_data['clicksArray1EXP'].flatten()
    current_chosenImEXP_vals = mat_data['current_chosenImEXP'].flatten()
    current_leftchosenEXP_vals = mat_data['current_leftchosenEXP'].flatten()
    current_subblogEXP_vals = mat_data['current_subblogEXP'].flatten()
    current_unchosenImEXP_vals = mat_data['current_unchosenImEXP'].flatten()
    elapsedTimeEXP_vals = mat_data['elapsedTimeEXP'].flatten()
    left_OutArrayEXP_vals = mat_data["left_OutArrayEXP"].flatten()
    right_OutArrayEXP_vals = mat_data["right_OutArrayEXP"].flatten()
    outcomeArrayEXP_vals = mat_data["outcomeArrayEXP"].flatten()
    scoreArrayEXP_vals = mat_data["scoreArrayEXP"].flatten()
    utility_ArrayEXP1_vals = mat_data["utility_ArrayEXP1"].flatten()
    leftimEXP_vals = mat_data['leftimEXP'].flatten()
    rightimEXP_vals = mat_data['rightimEXP'].flatten()
    scoreEXP_vals = mat_data["scoreEXP"].flatten()
    if len(scoreEXP_vals) == 1:
        scoreEXP_vals = [scoreEXP_vals[0]] * len(clicksArray1EXP_vals)

    #image probability pairing
    selected_images = mat_data['selectedImagesFromEXP']
    selected_probabilities = mat_data['selectedProbabilitiesFromEXP']

    image_prob_map = {}
    for row in range(selected_images.shape[0]):
        for col in range(selected_images.shape[1]):
            image = selected_images[row, col][0]
            probability = selected_probabilities[row, col]
            image_prob_map[image] = probability

    # Map probabilities to left and right images
    left_probabilities = [image_prob_map[img[0]] for img in leftimEXP_vals]
    right_probabilities = [image_prob_map[img[0]] for img in rightimEXP_vals]

    current_chosenImEXP_vals = clean_categorical_column(current_chosenImEXP_vals)
    current_unchosenImEXP_vals = clean_categorical_column(current_unchosenImEXP_vals)
    leftimEXP_vals = clean_categorical_column(leftimEXP_vals)
    rightimEXP_vals = clean_categorical_column(rightimEXP_vals)

    df = pd.DataFrame({
        "SubID": sub_id,
        "clicksArray1EXP": clicksArray1EXP_vals,
        "current_chosenImEXP": current_chosenImEXP_vals,
        "current_leftchosenEXP": current_leftchosenEXP_vals,
        "current_subblogEXP": current_subblogEXP_vals,
        "current_unchosenImEXP": current_unchosenImEXP_vals,
        "elapsed_timeEXP": elapsedTimeEXP_vals,
        "left_OutArrayEXP": left_OutArrayEXP_vals,
        "right_OutArrayEXP": right_OutArrayEXP_vals,
        "outcomeArrayEXP": outcomeArrayEXP_vals,
        "scoreArrayEXP": scoreArrayEXP_vals,
        "scoreEXP": scoreEXP_vals,
        "utility_ArrayEXP1": utility_ArrayEXP1_vals,
        "leftimEXP": leftimEXP_vals,
        "left_ProbabilityEXP": left_probabilities,
        "rightimEXP": rightimEXP_vals,
        "right_ProbabilityEXP": right_probabilities,
    })

    # creating columns Garcia has
    # ---------------------------------------------------------------
    # ev1 and ev2:
    # EV(1) = p1 * (+1) + (1 - p1) * (-1) => 2*p1 - 1
    # EV(2) = p2 * (+1) + (1 - p2) * (-1) => 2*p2 - 1
    # ---------------------------------------------------------------
    
    df["ev1"] = (2.0 * df["left_ProbabilityEXP"] - 1.0).round(1)
    df["ev2"] = (2.0 * df["right_ProbabilityEXP"] - 1.0).round(1)
    
    df["p1"] = df["left_ProbabilityEXP"].round(1)
    df["p2"] = df["right_ProbabilityEXP"].round(1)
    
    df["cho"] = df["clicksArray1EXP"].apply(lambda x: 1 if x == 1 else 2)

    #1 if (chosen EV >= unchosen EV), else 0
    # ---------------------------------------------------------------
    
    df["corr"] = np.where(
        (
            ((df["cho"] == 1) & (df["ev1"] >= df["ev2"]))  # left chosen & ev1 >= ev2
            | 
            ((df["cho"] == 2) & (df["ev2"] >= df["ev1"]))  # right chosen & ev2 >= ev1
        ),
        1,  # true
        0   # false
    )
   
    #trial col
    df["trial"] = np.arange(1, len(df) + 1)
    # cond col
    df["cond"] = df.apply(lambda row: get_condition_garcia(row["p1"], row["p2"]), axis=1)
    df["out"] = np.where(df["cho"] == 1,  # If choice=1 => left chosen
                         df["left_OutArrayEXP"],
                         df["right_OutArrayEXP"])
    df["cfout"] = np.where(df["cho"] == 1,  # If choice=1 => right unchosen
                           df["right_OutArrayEXP"],
                           df["left_OutArrayEXP"])
    # phase col
    df["phase"] = "LE"
    # sess col
    df["sess"] = 0
    df['op1'] = 'E'
    df['op2'] = 'E'
    # rew col: copy of scoreEXP
    df["rew"] = df["scoreArrayEXP"]

    # chose_right column: 1 if clicksArray1EXP == 0 else 0
    #    Because in clicksArray1EXP: 0 = right, 1 = left
    df["chose_right"] = df["clicksArray1EXP"].apply(lambda x: 1 if x == 0 else 0)
    
    df["rtime"] = df["elapsed_timeEXP"]
    
    df_participant = pd.DataFrame({
        "selectedImagesFromEXP": [selected_images.tolist()],
        "selectedProbabilitiesFromEXP": [selected_probabilities.tolist()],
    })

    csv_file_trials = os.path.join(participant_folder, f"EXP1_garcia_participant_{int(participant_id)}.csv")
    csv_file_participant = os.path.join(participant_folder, f"EXP1_garcia_participant_{int(participant_id)}_info.csv")
    df.to_csv(csv_file_trials, index=False)
    df_participant.to_csv(csv_file_participant, index=False)
    print(f"participant {participant_id} files saved in {participant_folder}")

#getting the learners column
participant_accuracies = []
for file_path in mat_files:
    participant_match = re.search(r"participant(\d+)", file_path)
    if participant_match:
        participant_id = participant_match.group(1)
    else:
        raise ValueError(f"ID not in path: {file_path}")
    mat_data = scipy.io.loadmat(file_path)
    participant_folder = os.path.join(output_dir, f"Participant_{int(participant_id)}")
    csv_file_trials = os.path.join(participant_folder, f"EXP1_garcia_participant_{int(participant_id)}.csv")
    df = pd.read_csv(csv_file_trials)
    overall_accuracy = df["corr"].mean()
    participant_accuracies.append({"participant_id": int(participant_id), "accuracy": overall_accuracy})

# get median accuracy
accuracy_df = pd.DataFrame(participant_accuracies)
median_accuracy = accuracy_df["accuracy"].median()
print(f"Median Accuracy: {median_accuracy}")
accuracy_df["learners"] = accuracy_df["accuracy"].apply(lambda x: "good" if x >= median_accuracy else "poor")
good_learners = accuracy_df[accuracy_df["learners"] == "good"]
poor_learners = accuracy_df[accuracy_df["learners"] == "poor"]
print("\nGood Learners:")
print(good_learners)
print("\nPoor Learners:")
print(poor_learners)

for _, row in accuracy_df.iterrows():
    participant_id = row["participant_id"]
    learner_type = row["learners"]
    participant_folder = os.path.join(output_dir, f"Participant_{int(participant_id)}")
    csv_file_trials = os.path.join(participant_folder, f"EXP1_garcia_participant_{int(participant_id)}.csv")
    df = pd.read_csv(csv_file_trials)
    df["learners"] = learner_type
    csv_file_with_learners = os.path.join(participant_folder, f"EXP1_garcia_participant_{int(participant_id)}_with_learners.csv")
    df.to_csv(csv_file_with_learners, index=False)
    print(f"file with learners saved for participant {participant_id} at {csv_file_with_learners}")