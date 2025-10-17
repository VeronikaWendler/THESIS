# Veronika Wendler
# 21.01.25
# code that selects and attaches eyetracking data to the behavioural Garcia file:GarciaParticipants_EXP_1_2_3_4.csv  that contains all participants

import pandas as pd

# File paths
eye_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data"
ov_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_EXP_1_2_3_4.csv"
output_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_Eye_Response_EXP_1_2_3_4.csv"

# participants to exclude based on insufficient data
# participant 20 & 14 are already excluded

exclude_part = {26, 2, 6, 9, 20, 14, 18}            # participant 14 & 20 & 26 did not have a .edf file,  - exclude
                                                    # participant 2 & 9 has insufficient eye-tracking data (almsot non exisiting)  - exclude
                                                    # participant 6 only has data from phases LE to ES -exclude as well
                                                    # participant 18's eye-tracker stopped just before the end

# 20 particiapnts for the eye-tracking analysis

participant_indices = {
    1: (1, 660),                          
    2: (712, 1371),     # eye-tracker stopped at the start (not enough data)- exclude this participatn from eyetracking analysis                  
    3: (1423, 2082),
    4: (2134, 2793),                    
    5: (2845, 3504),   
    6: (3556, 4215),     #only has data from limited phases            
    7: (4267, 4926),
    8: (4978, 5637),
    9: (5689, 6348),     # no edf file
    10: (6400, 7059),
    11: (7111, 7770),
    12: (7822, 8481),
    13: (8533, 9192),
    14: (9244, 9763),   # exclude this one from analysis not enough eye-tracking data, missing phases
    15: (9764, 10423),
    16: (10475, 11134),
    17: (11186, 11845),
    18: (11897, 12556),  # eye-tracker stopped at the end
    19: (12608, 13267),
    20: (13319, 13978),   # no edf file
    21: (13979, 14638), 
    22: (14690, 15349),
    23: (15401, 16060), 
    24: (16112, 16771),
    25: (16823, 17482),
    26: (17534, 18193),   # 26 did not have a .edf file
    27: (18245, 18904)
}


# Load the Garcia dataset
ov_data = pd.read_csv(ov_file)

# Initialize eye-tracking columns in the original dataset
eye_tracking_columns = [
    "TrialID", 
    "Task",
    "DwellLeft",
    "DwellRight",
    "DwellTotal", 
    "Nfix", 
    "FirstFixLoc",
    "FirstFixDur",
    "FinalFixLoc",
    "FinalFixDur", 
    "MiddleFixDur",
    "eachMiddleFixDur",
    "GazeSwitch",
    "LeftFixNR",
    "RightFixNR",
    "GazeDiff", 
    "last_roi",
    "Messages",
    "Phase",
    "Fixations",
    "MiddleDominantLoc",
    "MiddleDominantDur", 
]

for col in eye_tracking_columns:
    if col not in ov_data.columns:
        ov_data[col] = None

# Process participants
participants = list(range(1, 28))  # Includes participants 1 to 27
#participants.append(99)

for participant in participants:
    if participant in exclude_part:
        continue

    # Load the eye-tracking file for the current participant
    eye_file = f"{eye_dir}/sub-{str(participant).zfill(2)}/eyetracking/sub-{str(participant).zfill(2)}_eye_response_phase.csv"
    try:
        eye_data = pd.read_csv(eye_file)
    except FileNotFoundError:
        print(f"Eye-tracking file for participant {participant} not found. Skipping...")
        continue

    # Filter eye-tracking data for TrialID >= 49
    eye_data_filtered = eye_data[eye_data["TrialID"] >= 49].reset_index(drop=True)

    # Check if participant has corresponding indices in the OV data
    if participant in participant_indices:
        start_idx, end_idx = participant_indices[participant]
        ov_sub_indices = (
            (ov_data["sub_id"] == float(participant)) &
            (ov_data["index"] >= start_idx) &
            (ov_data["index"] <= end_idx)
        )

        ov_sub = ov_data.loc[ov_sub_indices]

        # Ensure matching row counts
        if len(eye_data_filtered) != len(ov_sub):
            raise ValueError(f"Row count mismatch: Eye-tracking ({len(eye_data_filtered)}) and OV ({len(ov_sub)}).")

        # Assign eye-tracking data to the relevant rows in garcia_data
        for col in eye_tracking_columns:
            if col not in eye_data_filtered.columns:
                eye_data_filtered[col] = None
            ov_data.loc[ov_sub_indices, col] = eye_data_filtered[col].values

        # Debugging: Check assigned columns
        print(f"Participant {participant}: Eye-tracking data successfully merged.")
        print(ov_data.loc[ov_sub_indices].head())

# Save the updated dataset without modifying the original structure
ov_data.to_csv(output_file, index=False)
print(f"Updated OV data saved to {output_file}")


ov_data.update(ov_sub)

# Save the updated dataset
ov_data.to_csv(output_file, index=False)
print(f"Updated OV data saved to {output_file}")

# # Step 6: Resolve column duplication
# # # Keep `_x` columns and rename them to their original names
# garcia_df = garcia_df.rename(columns=lambda x: x[:-2] if x.endswith('_x') else x)

# # Drop all `_y` columns
# garcia_df = garcia_df[[col for col in garcia_df.columns if not col.endswith('_y')]]

# # Step 7: Save the updated Garcia dataset
# garcia_df.to_csv(output_file, index=False)
# print(f"Updated Garcia data saved to {output_file}")


