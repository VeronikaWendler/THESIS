# Veronika Wendler
# 21.01.25

# this code additionally appends eye tracking data from the feedback phase during the learning phase
# (that is, fixations were filtered for the duration of the feedback, e.g. gaining +1 or loosing -1 points)


import pandas as pd

# File paths
eye_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data"
garcia_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_EXP_1_2_3_4.csv"
output_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_EXP_1_2_3_4.csv"

# participants to exclude based on insufficient data
exclude_part = {4, 5, 6, 14}

# dictionary for participant indices, 
participant_indices = {
    1: (1, 160),
    2: (712, 871),
    3: (1423, 1582),
    4: (2134, 2293),
    5: (2845, 3004),
    7: (3556, 3715),
    8: (4267, 4426),
    9: (4978, 5137),
    10: (5689, 5848),
    11: (6400, 6559),
    12: (7111, 7270),
    13: (7822, 7981),
    14: (8533, 8692),
    15: (9244, 9403),
    16: (9955, 10114),
    17: (10666, 10825),
    18: (11377, 11536),
    19: (12088, 12247),
    20: (12799, 12958),
    21: (13510, 13669),
    22: (14221, 14380),
    23: (14932, 15091),
    24: (15643, 15802),
    25: (16354, 16513),
    26: (17065, 17224),
    99: (17776, 17935)
}

garcia_data = pd.read_csv(garcia_file)

eye_tracking_columns = [
    "DwellLeft_Feed",
    "DwellRight_Feed",
    "DwellTotal_Feed", 
    "Nfix_Feed", 
    "FirstFixLoc_Feed",
    "FirstFixDur_Feed",
    "FinalFixLoc_Feed",
    "FinalFixDur_Feed", 
    "MiddleFixDur_Feed",
    "eachMiddleFixDur_Feed",
    "GazeSwitch_Feed",
    "LeftFixNR_Feed",
    "RightFixNR_Feed",
    "GazeDiff_Feed", 
    "last_roi_Feed",
    "Messages_Feed",
    "Phase_Feed",
    "Fixations_Feed"                            
]

for col in eye_tracking_columns:
    if col not in garcia_data.columns:
        garcia_data[col] = None

participants = list(range(1, 27))  # Includes participants 1 to 26
participants.append(99)

for participant in participants:
    if participant in exclude_part:
        continue

    eye_file = f"{eye_dir}/sub-{str(participant).zfill(2)}/eyetracking/sub-{str(participant).zfill(2)}_eye_feedback_phase.csv"
    try:
        eye_data = pd.read_csv(eye_file)
    except FileNotFoundError:
        print(f"Eye-tracking file for participant {participant} not found. Skipping...")
        continue

    # Filter eye-tracking data for TrialID >= 49
    eye_data_filtered = eye_data[(eye_data["TrialID"] >= 49) & (eye_data["TrialID"] <= 208)].reset_index(drop=True)

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