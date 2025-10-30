# Veronika Wendler
# this program takes the behavioural csv and attaches the eye-tracking metrics

import os
import pandas as pd

# paths
bids_root = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/"
garcia_csv_path = os.path.join(bids_root, "GarciaParticipants_EXP_1_2_3_4.csv")
updated_csv_path = os.path.join(bids_root, "GarciaParticipants_EXP_1_2_3_4_with_eyetracking.csv")

garcia_df = pd.read_csv(garcia_csv_path)

# Standardize sub_id in Garcia DataFrame
garcia_df["sub_id"] = garcia_df["sub_id"].astype(str).apply(lambda x: f"sub-{int(float(x)):02d}")

# list for collecting eye-tracking data
eye_tracking_data = []

for root, dirs, files in os.walk(bids_root):
    for file in files:
        if file.endswith("eye_response_phase.csv"):
            # get ID
            participant_id = file.split("_")[0]  # e.g., "sub-01"
            
            # load eye-tracking data
            eye_tracking_path = os.path.join(root, file)
            eye_df = pd.read_csv(eye_tracking_path)
            
            # standardize sub_id in eye-tracking DataFrame
            eye_df["sub_id"] = participant_id
            
            # map trials (trials 49 to 209 correspond to 1-160 in Garcia)
            eye_df_filtered = eye_df[(eye_df["TrialID"] >= 49) & (eye_df["TrialID"] <= 209)].copy()
            eye_df_filtered["trial"] = eye_df_filtered["TrialID"] - 48
            
            # keep columns
            cols_to_keep = ["sub_id",
                            "trial",
                            "DwellLeft",
                            "DwellRight",
                            "DwellTotal",
                            "Nfix", 
                            "FirstFixLoc",
                            "FinalFixLoc",
                            "FirstFixDur",
                            "FinalFixDur",
                            "MiddleFixDur", 
                            "GazeSwitch",
                            "LeftFixNR",
                            "RightFixNR",
                            "GazeDiff",
                            "last_roi",
                            "Fixations"]
            eye_tracking_data.append(eye_df_filtered[cols_to_keep])

# combine all
combined_eye_df = pd.concat(eye_tracking_data, ignore_index=True)

# trial should match the first 160 trials in Garcia (because this is the lenght of LE phase)
combined_eye_df = combined_eye_df[combined_eye_df["trial"].isin(range(1, 161))]

# merge with Garcia DataFrame on left
merged_df = pd.merge(
    garcia_df,
    combined_eye_df,
    how="left",
    on=["sub_id", "trial"]
)

merged_df.to_csv(updated_csv_path, index=False)
print(f"Updated data saved to {updated_csv_path}")
