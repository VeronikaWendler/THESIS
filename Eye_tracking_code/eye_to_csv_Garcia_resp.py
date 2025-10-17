# Veronika Wendler
# 21.01.25
# code that selects and attaches eyetracking data to the behavioural Garcia file:GarciaParticipants_EXP_1_2_3_4.csv  that contains all participants

import pandas as pd

# File paths
eye_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data"
garcia_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_EXP_1_2_3_4.csv"
output_file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_EXP_1_2_3_4.csv"

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

#insufficient data for participant 6, 
# for participant 1, index <=660
# for participant 2, 712 <= index <= 1371
# for participant 3, 1423 <= index <= 2082
# for participant 4, 2134 <= index <= 2793
# for participant 5, 2845 <= index <= 3504
# for participant 7, 3556 <= index <= 4215
# for participant 8, 4267 <= index <= 4926
# for participant 9, 4978 <= index <= 5637
# for participant 10, 5689 <= index <= 6348
# for participant 11, 6400 <= index <= 7059
# for participant 12, 7111 <= index <= 7770
# for participant 13, 7822 <= index <= 8481
# for participant 14, 8533 <= index <= 9192
# for participant 15, 9244 <= index <= 9903
# for participant 16, 9955 <= index <= 10614
# for participant 17, 10666 <= index <= 11325
# for participant 18, 11377 <= index <= 12036
# for participant 19, 12088 <= index <= 12747
# for participant 20, 12799 <= index <= 13458
# for participant 21, 13510 <= index <= 14169
# for participant 22, 14221 <= index <= 14880
# for participant 23, 14932 <= index <= 15591
# for participant 24, 15643 <= index <= 16302
# for participant 25, 16354 <= index <= 17013
# for participant 26, 17065 <= index <= 17724
# for participant 99, 17776 <= index <= 18435

# Load the Garcia dataset
garcia_data = pd.read_csv(garcia_file)

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
    if col not in garcia_data.columns:
        garcia_data[col] = None

# Process participants
participants = list(range(1, 27))  # Includes participants 1 to 26
participants.append(99)

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

    # Check if participant has corresponding indices in the Garcia data
    if participant in participant_indices:
        start_idx, end_idx = participant_indices[participant]
        garcia_sub_indices = (
            (garcia_data["sub_id"] == float(participant)) &
            (garcia_data["index"] >= start_idx) &
            (garcia_data["index"] <= end_idx)
        )

        garcia_sub = garcia_data.loc[garcia_sub_indices]

        # Ensure matching row counts
        if len(eye_data_filtered) != len(garcia_sub):
            raise ValueError(f"Row count mismatch: Eye-tracking ({len(eye_data_filtered)}) and Garcia ({len(garcia_sub)}).")

        # Assign eye-tracking data to the relevant rows in garcia_data
        for col in eye_tracking_columns:
            if col not in eye_data_filtered.columns:
                eye_data_filtered[col] = None
            garcia_data.loc[garcia_sub_indices, col] = eye_data_filtered[col].values

        # Debugging: Check assigned columns
        print(f"Participant {participant}: Eye-tracking data successfully merged.")
        print(garcia_data.loc[garcia_sub_indices].head())

# Save the updated dataset without modifying the original structure
garcia_data.to_csv(output_file, index=False)
print(f"Updated Garcia data saved to {output_file}")

# for participant in participants:
#     if participant in exclude_part:
#         continue

#     # Load the eye-tracking file for the current participant
#     eye_file = f"{eye_dir}/sub-{str(participant).zfill(2)}/eyetracking/sub-{str(participant).zfill(2)}_eye_response_phase.csv"
#     try:
#         eye_data = pd.read_csv(eye_file)
#     except FileNotFoundError:
#         print(f"Eye-tracking file for participant {participant} not found. Skipping...")
#         continue

#     # Filter eye-tracking data for TrialID >= 49
#     eye_data_filtered = eye_data[eye_data["TrialID"] >= 49].reset_index(drop=True)

#     print(f"Participant {participant}: eye_data_filtered rows = {len(eye_data_filtered)}")

#     if participant in participant_indices:
#         start_idx, end_idx = participant_indices[participant]
#         garcia_sub = garcia_data[
#             (garcia_data["sub_id"] == float(participant)) &
#             (garcia_data["index"] >= start_idx) &
#             (garcia_data["index"] <= end_idx)
#         ].reset_index(drop=True)

#         print(f"Participant {participant}: garcia_sub rows = {len(garcia_sub)}")

#         # Reset indices to align rows
#         eye_data_filtered = eye_data_filtered.reset_index(drop=True)
#         garcia_sub = garcia_sub.reset_index(drop=True)

        
#         if len(eye_data_filtered) != len(garcia_sub):
#             raise ValueError(f"Row count mismatch: Eye-tracking ({len(eye_data_filtered)}) and Garcia ({len(garcia_sub)}).")

#         # adding eye-tracking columns to Garcia dataset
#         eye_tracking_columns = ["TrialID", 
#                                 "Task", 
#                                 "DwellLeft",
#                                 "DwellRight",
#                                 "DwellTotal", 
#                                 "Nfix", 
#                                 "FirstFixLoc",
#                                 "FirstFixDur",
#                                 "FinalFixLoc",
#                                 "FinalFixDur", 
#                                 "MiddleFixDur",
#                                 "GazeSwitch",
#                                 "LeftFixNR",
#                                 "RightFixNR", 
#                                 "GazeDiff", 
#                                 "last_roi", 
#                                 "Messages",
#                                 "Fixations",
#                                 "Phase"]

#         # Add new columns to the filtered Garcia rows
#         for col in eye_tracking_columns:
#             garcia_sub[col] = eye_data_filtered[col].values
#         # Assign columns dynamically
        
#         # Verify all required columns exist
#         missing_columns = [col for col in eye_tracking_columns if col not in eye_data_filtered.columns]
#         if missing_columns:
#             raise ValueError(f"Missing columns in eye-tracking data: {missing_columns}")

#         # # Assign eye-tracking columns dynamically
#         # garcia_sub = garcia_sub.assign(**{col: eye_data_filtered[col].values for col in eye_tracking_columns})

#         # Debugging: Check if data is populated
#         print(f"Participant {participant}: Assigned eye-tracking data to garcia_sub")
#         print(garcia_sub.head())
        
#Merge updated filtered rows back into the original Garcia dataset
#garcia_df = garcia_data.merge(garcia_sub, how="left", on="index")
# Update `garcia_data` with the modified rows from `garcia_sub`
garcia_data.update(garcia_sub)

# Save the updated dataset
garcia_data.to_csv(output_file, index=False)
print(f"Updated Garcia data saved to {output_file}")

# # Step 6: Resolve column duplication
# # # Keep `_x` columns and rename them to their original names
# garcia_df = garcia_df.rename(columns=lambda x: x[:-2] if x.endswith('_x') else x)

# # Drop all `_y` columns
# garcia_df = garcia_df[[col for col in garcia_df.columns if not col.endswith('_y')]]

# # Step 7: Save the updated Garcia dataset
# garcia_df.to_csv(output_file, index=False)
# print(f"Updated Garcia data saved to {output_file}")



       
        
        
        
# # Merge updated filtered rows back into the original Garcia dataset
# for col in eye_tracking_columns:
#     if col not in garcia_data.columns:
#         garcia_data[col] = pd.NA  # Add missing columns
        
# garcia_data = garcia_data.sort_values(by=['sub_id', 'phase', 'index']).reset_index(drop=True)


# # Update only the rows that match the indices in garcia_sub
# garcia_data.update(garcia_sub)

# # Save the updated dataset
# garcia_data.to_csv(output_file, index=False)
# print(f"Updated Garcia data saved to {output_file}")

        




# ######################################################################################################
#insufficient data for participant 6, 
# for participant 1, index <=660
# for participant 2, 712 <= index <= 1371
# for participant 3, 1423 <= index <= 2082
# for participant 4, 2134 <= index <= 2793
# for participant 5, 2845 <= index <= 3504
# for participant 7, 3556 <= index <= 4215
# for participant 8, 4267 <= index <= 4926
# for participant 9, 4978 <= index <= 5637
# for participant 10, 5689 <= index <= 6348
# for participant 11, 6400 <= index <= 7059
# for participant 12, 7111 <= index <= 7770
# for participant 13, 7822 <= index <= 8481
# for participant 14, 8533 <= index <= 9192
# for participant 15, 9244 <= index <= 9903
# for participant 16, 9955 <= index <= 10614
# for participant 17, 10666 <= index <= 11325
# for participant 18, 11377 <= index <= 12036
# for participant 19, 12088 <= index <= 12747
# for participant 20, 12799 <= index <= 13458
# for participant 21, 13510 <= index <= 14169
# for participant 22, 14221 <= index <= 14880
# for participant 23, 14932 <= index <= 15591
# for participant 24, 15643 <= index <= 16302
# for participant 25, 16354 <= index <= 17013
# for participant 26, 17065 <= index <= 17724
# for participant 99, 17776 <= index <= 18435

# +  659, 52




# import pandas as pd

# # File paths
# eye_tracking_dir = r"D:\Aberdeen_Uni_June24\cap\THESIS\Garcia_Analysis\data"
# garcia_file = r"D:\Aberdeen_Uni_June24\cap\THESIS\Garcia_Analysis\data\GarciaParticipants_EXP_1_2_3_4.csv"
# output_file = r"D:\Aberdeen_Uni_June24\cap\THESIS\Garcia_Analysis\data\Updated_2_GarciaParticipants_EXP_1_2_3_4.csv"

# # Participant-specific details
# participants_to_exclude = {4, 5, 6, 14}
# num_trials_per_participant = 660
# trial_offset = 52
# start_trial_id = 49

# # Load the datasets
# garcia_df = pd.read_csv(garcia_file)

# # Initialize an empty DataFrame to store updated Garcia data
# updated_garcia_df = pd.DataFrame()

# # Loop through participants
# for participant in range(1, 27 + 1):
#     if participant in participants_to_exclude:
#         continue

#     # Calculate the index range for the current participant
#     start_index = (participant - 1) * (num_trials_per_participant + trial_offset)
#     end_index = start_index + num_trials_per_participant

#     # Load the eye-tracking file for the current participant
#     eye_tracking_file = f"{eye_tracking_dir}\sub-{str(participant).zfill(2)}\eyetracking\sub-{str(participant).zfill(2)}_eye_response_phase.csv"
#     try:
#         eye_tracking_df = pd.read_csv(eye_tracking_file)
#     except FileNotFoundError:
#         print(f"Eye-tracking file for participant {participant} not found. Skipping...")
#         continue

#     # Step 1: Filter eye-tracking data for TrialID >= 49
#     eye_tracking_filtered = eye_tracking_df[eye_tracking_df["TrialID"] >= start_trial_id].reset_index(drop=True)

#     # Step 2: Filter Garcia data for the current participant's index range
#     garcia_filtered = garcia_df[(garcia_df["index"] >= start_index) & (garcia_df["index"] < end_index)].reset_index(drop=True)

#     # Step 3: Ensure the row count matches between datasets
#     if len(eye_tracking_filtered) != len(garcia_filtered):
#         raise ValueError(
#             f"Row count mismatch for participant {participant}: "
#             f"Eye-tracking ({len(eye_tracking_filtered)}) and Garcia ({len(garcia_filtered)})."
#         )

#     # Step 4: Add eye-tracking columns to Garcia dataset
#     eye_tracking_columns = [
#         "TrialID", 
#         "Task", 
#         "DwellLeft",
#         "DwellRight", 
#         "DwellTotal",
#         "Nfix",
#         "FirstFixLoc",
#         "FirstFixDur", 
#         "FinalFixLoc",
#         "FinalFixDur", 
#         "MiddleFixDur",
#         "GazeSwitch",
#         "LeftFixNR",
#         "RightFixNR", 
#         "GazeDiff",
#         "last_roi",
#         "Messages",
#         "Fixations",
#         "Phase"
#     ]

#     for col in eye_tracking_columns:
#         garcia_filtered[col] = eye_tracking_filtered[col].values

#     # Append the updated rows to the overall DataFrame
#     updated_garcia_df = pd.concat([updated_garcia_df, garcia_filtered], ignore_index=True)

# # Step 5: Merge updated filtered rows back into the original Garcia dataset
# garcia_df = garcia_df.merge(updated_garcia_df, how="left", on="index")

# # Step 6: Resolve column duplication
# # Keep `_x` columns and rename them to their original names
# garcia_df = garcia_df.rename(columns=lambda x: x[:-2] if x.endswith('_x') else x)

# # Drop all `_y` columns
# garcia_df = garcia_df[[col for col in garcia_df.columns if not col.endswith('_y')]]

# # Step 7: Save the updated Garcia dataset
# garcia_df.to_csv(output_file, index=False)
# print(f"Updated Garcia data saved to {output_file}")
