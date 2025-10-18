# Veronika Wendler
# 21.01.25
# code that selects and attaches eyetracking data to the behavioural Garcia file:GarciaParticipants_EXP_1_2_3_4.csv  that contains all participants

import pandas as pd

# File paths
eye_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data"
garcia_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_EXP_1_2_3_4.csv"
output_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_EXP_1_2_3_4.csv"

# participants to exclude based on insufficient data
exclude_part = {4, 5, 6, 14}

participant_indices = {
    1: (1, 660),
    2: (712, 1371),
    3: (1423, 2082),
    4: (2134, 2793),
    5: (2845, 3504),
    7: (3556, 4215),
    8: (4267, 4926),
    9: (4978, 5637),
    10: (5689, 6348),
    11: (6400, 7059),
    12: (7111, 7770),
    13: (7822, 8481),
    14: (8533, 9192),
    15: (9244, 9903),
    16: (9955, 10614),
    17: (10666, 11325),
    18: (11377, 12036),
    19: (12088, 12747),
    20: (12799, 13458),
    21: (13510, 14169),
    22: (14221, 14880),
    23: (14932, 15591),
    24: (15643, 16302),
    25: (16354, 17013),
    26: (17065, 17724),
    99: (17776, 18435)
}


# Load the Garcia dataset
garcia_data = pd.read_csv(garcia_file)

# Initialize eye-tracking columns in the original dataset
eye_tracking_columns = [
    "TrialID", 
    "Task",
    "DwellLeft_allfix",
    "DwellRight_allfix",
    "DwellTotal_allfix", 
    "Nfix_allfix", 
    "FirstFixLoc_allfix",
    "FirstFixDur_allfix",
    "FinalFixLoc_allfix",
    "FinalFixDur_allfix", 
    "MiddleFixDur_allfix",
    "eachMiddleFixDur_allfix",
    "GazeSwitch_allfix",
    "LeftFixNR_allfix",
    "RightFixNR_allfix",
    "GazeDiff_allfix", 
    "last_roi_allfix",
    "Messages",
    "Phase",
    "Fixations_allfix"
]

for col in eye_tracking_columns:
    if col not in garcia_data.columns:
        garcia_data[col] = None

# Process participants
participants = list(range(1, 27))  # Includes participants 1 to 26
participants.append(99)

for participant in participants:
    if participant in exclude_part:
        continue

    eye_file = f"{eye_dir}/sub-{str(participant).zfill(2)}/eyetracking/sub-{str(participant).zfill(2)}_eye_allfix_phase.csv"
    try:
        eye_data = pd.read_csv(eye_file)
    except FileNotFoundError:
        print(f"Eye-tracking file for participant {participant} not found. Skipping...")
        continue

    # Filter eye-tracking data for TrialID >= 49 (everything up to there is training data)
    eye_data_filtered = eye_data[eye_data["TrialID"] >= 49].reset_index(drop=True)

    if participant in participant_indices:
        start_idx, end_idx = participant_indices[participant]
        garcia_sub_indices = (
            (garcia_data["sub_id"] == float(participant)) &
            (garcia_data["index"] >= start_idx) &
            (garcia_data["index"] <= end_idx)
        )

        garcia_sub = garcia_data.loc[garcia_sub_indices]

        if len(eye_data_filtered) != len(garcia_sub):
            raise ValueError(f"Row count mismatch: Eye-tracking ({len(eye_data_filtered)}) and Garcia ({len(garcia_sub)}).")

        for col in eye_tracking_columns:
            if col not in eye_data_filtered.columns:
                eye_data_filtered[col] = None
            garcia_data.loc[garcia_sub_indices, col] = eye_data_filtered[col].values

        print(f"Participant {participant}: Eye-tracking data successfully merged.")
        print(garcia_data.loc[garcia_sub_indices].head())

garcia_data.to_csv(output_file, index=False)
print(f"Updated Garcia data saved to {output_file}")
garcia_data.update(garcia_sub)
garcia_data.to_csv(output_file, index=False)
print(f"Updated Garcia data saved to {output_file}")
