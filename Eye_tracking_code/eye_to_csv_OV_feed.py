# Veronika Wendler
# 21.01.25

# this code additionally appends eye tracking data from the feedback phase during the learning phase
# (that is, fixations were filtered for the duration of the feedback, e.g. gaining +1 or loosing -1 points)


import pandas as pd

# File paths
eye_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data"
ov_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_Eye_Response_EXP_1_2_3_4.csv"
output_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_Eye_Response_Feed_EXP_1_2_3_4.csv"

exclude_part = {26, 2, 9, 20, 14}            # participant 14 & 20 & 26 did not have a .edf file,  - exclude
                                                    # participant 2 & 9 has insufficient eye-tracking data (almsot non exisiting)  - exclude
                                                    # participant 6 only has data from phases LE to ES -exclude as well
                                                    # participant 18's eye-tracker stopped just before the end

# 20 particiapnts for the eye-tracking analysis

participant_indices = {
    1: (1, 160),                          
    2: (712, 871),     # eye-tracker stopped at the start (not enough data)- exclude this participatn from eyetracking analysis                  
    3: (1423, 1582),
    4: (2134, 2293),                    
    5: (2845, 3004),   
    6: (3556, 3715),     # can be used since has full data from LE phase       
    7: (4267, 4426),
    8: (4978, 5137),
    9: (5689, 5848),     # no edf file
    10: (6400, 6559),
    11: (7111, 7270),
    12: (7822, 7981),
    13: (8533, 8692),
    14: (9244, 9403),   # exclude this one from analysis no edf file
    15: (9764, 9923),
    16: (10475, 10634),
    17: (11186, 11345),
    18: (11897, 12056),  # can be used since LE has eye-tracking data
    19: (12608, 12767),
    20: (13319, 13478),   # no edf file
    21: (13979, 14138), 
    22: (14690, 14849),
    23: (15401, 15560), 
    24: (16112, 16271),
    25: (16823, 16982),  
    26: (17534, 17693),   # 26 did not have a .edf file
    27: (18245, 18404)
}


# Load the Garcia dataset
ov_data = pd.read_csv(ov_file)

# Initialize eye-tracking columns in the original dataset
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
    "eachMiddleFixDur",
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
    if col not in ov_data.columns:
        ov_data[col] = None

participants = list(range(1, 27))  
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
        ov_sub_indices = (
            (ov_data["sub_id"] == float(participant)) &
            (ov_data["index"] >= start_idx) &
            (ov_data["index"] <= end_idx)
        )

        ov_sub = ov_data.loc[ov_sub_indices]

        if len(eye_data_filtered) != len(ov_sub):
            raise ValueError(f"Row count mismatch: Eye-tracking ({len(eye_data_filtered)}) and OV ({len(ov_sub)}).")

        for col in eye_tracking_columns:
            if col not in eye_data_filtered.columns:
                eye_data_filtered[col] = None
            ov_data.loc[ov_sub_indices, col] = eye_data_filtered[col].values

        print(f"Participant {participant}: Eye-tracking data successfully merged.")
        print(ov_data.loc[ov_sub_indices].head())

ov_data.to_csv(output_file, index=False)
print(f"Updated OV data saved to {output_file}")