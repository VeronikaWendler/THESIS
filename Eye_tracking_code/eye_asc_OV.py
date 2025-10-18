import os
import re
import pandas as pd

mode = 2 # 1 = filters fixations form start of trial until response 
          # 2 = filters fixations during the feedback phase of a trial (e.g. participants see outcome; only applicable in LE phase)
          # 3 = relaxed filetring, counts all fixations in the trial (most of the time I am not using this since I am explicitly interested in participatns fixations from the start of the trial until they made a response)

#directory
bids_root = 'D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/'

# screen parameters
image_size = [200, 200]
image_distance = 300
screen_width = 1919
screen_height = 1079
screen_center = (screen_width / 2, screen_height / 2)
middle_margin = 50
middle_left_bound = screen_center[0] - middle_margin
middle_right_bound = screen_center[0] + middle_margin

def assign_phase(messages):
    if not messages:
        return "Unknown"
    messages = messages.strip().lower()
    if "training 0" in messages and "task 1" in messages:
        return "LE"
    elif "training 0" in messages and "task 2" in messages:
        return "ES"
    elif "training 0" in messages and "task 3" in messages:
        return "EE"
    else:
        return "T"  


# Loop through participants' directories
for root, dirs, files in os.walk(bids_root):
    for file in files:
        if file.endswith('.asc'):
            participant_id = os.path.basename(os.path.dirname(root))
            file_path = os.path.join(root, file)
            print(f"Processing file: {file_path} for participant: {participant_id}")

            trial_data = []
            current_trial_id = None
            current_task_message = None
            fixation_events = []
            fixation_start_time = None
            trial_start_time = None
            trial_end_time = None

            # process .asc file
            with open(file_path, 'r') as asc:
                if mode == 1:
                    print("Running Code NR 1 (Trial to Response Period)...")

                    for line in asc:
                        # Detects trial ID
                        trial_match = re.search(r'MSG\s+(\d+)\s+TRIALID\s+(\d+)', line)
                        if trial_match:
                            # save trial data
                            if current_trial_id is not None:
                                
                                if len(fixation_events) > 2:
                                    mids = fixation_events[1:-1]              
                                    left_dur  = sum(f["duration"] for f in mids if f["roi"] == 1)
                                    right_dur = sum(f["duration"] for f in mids if f["roi"] == 2)
                                    if   left_dur  > right_dur:
                                        mid_dom_loc, mid_dom_dur = 1, left_dur
                                    elif right_dur > left_dur:
                                        mid_dom_loc, mid_dom_dur = 2, right_dur
                                    else:                                     
                                        mid_dom_loc, mid_dom_dur = 0, left_dur
                                else:
                                    mid_dom_loc, mid_dom_dur = None, 0
                                
                                
                                trial_data.append({
                                    "TrialID": current_trial_id,
                                    "Task": current_task_message,
                                    "DwellLeft": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1),
                                    "DwellRight": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2),
                                    "DwellTotal": sum(fix["duration"] for fix in fixation_events),
                                    "Nfix": len(fixation_events),
                                    "FirstFixLoc": fixation_events[0]["roi"] if fixation_events else None,
                                    "FirstFixDur": fixation_events[0]["duration"] if fixation_events else None,
                                    "FinalFixLoc": fixation_events[-1]["roi"] if fixation_events else None,
                                    "FinalFixDur": fixation_events[-1]["duration"] if fixation_events else None,
                                    "MiddleFixDur": max(
                                        0,
                                        sum(fix["duration"] for fix in fixation_events) -
                                        (fixation_events[0]["duration"] + fixation_events[-1]["duration"])
                                        if len(fixation_events) > 1 else 0
                                    ),
                                    "eachMiddleFixDur": (
                                        max(
                                            0,
                                            sum(fix["duration"] for fix in fixation_events) -
                                            (fixation_events[0]["duration"] + fixation_events[-1]["duration"])) / max(1, len(fixation_events) - 2)
                                        if len(fixation_events) > 2 else 0
                                    ),
                                    "MiddleDominantLoc": mid_dom_loc,      
                                    "MiddleDominantDur": mid_dom_dur,     
                                    "GazeSwitch": sum(
                                        1 for i in range(1, len(fixation_events))
                                        if fixation_events[i]["roi"] != fixation_events[i - 1]["roi"]
                                    ),
                                    "LeftFixNR": sum(1 for fix in fixation_events if fix["roi"] == 1),
                                    "RightFixNR": sum(1 for fix in fixation_events if fix["roi"] == 2),
                                    "GazeDiff": abs(
                                        sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1) -
                                        sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2)
                                    ),
                                    "last_roi": fixation_events[-1]["roi"] if fixation_events else None,
                                    "Messages": current_task_message,
                                    "Fixations": fixation_events,

                                })

                            # new trial
                            current_trial_id = int(trial_match.group(2))
                            trial_start_time = int(trial_match.group(1))
                            trial_end_time = None
                            fixation_events = []
                            current_task_message = None
                            continue

                        # Detect task message
                        msg_match = re.search(r'MSG\s+\d+\s+!V TRIAL_VAR\s+(.*)', line)
                        if msg_match:
                            current_task_message = msg_match.group(1).strip()
                            continue
                        
                        # Detect response message (end of trial)
                        response_match = re.search(r'MSG\s+(\d+)\s+(ResponseLeft|ResponseRight)', line)
                        if response_match:
                            trial_end_time = int(response_match.group(1))
                            continue

                        # Detect fixation start (SFIX)
                        sfix_match = re.search(r'SFIX\s+\w+\s+(\d+)', line)
                        if sfix_match:
                            fixation_start_time = int(sfix_match.group(1))
                            continue

                        # Detect fixation end (EFIX)
                        efix_match = re.search(r'EFIX\s+\w+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                        if efix_match and current_trial_id and trial_start_time:
                            if fixation_start_time is None:
                                continue
                            
                            fixation_end_time = int(efix_match.group(2))
                            fixation_duration = int(efix_match.group(3))
                            avg_x = float(efix_match.group(4))
                            
                             # process fixations before or overlapping with the response time (offers leeway for leftover EFIX)
                            if trial_end_time is not None and fixation_start_time > trial_end_time:
                                continue

                            roi = None
                            # ROI
                            if avg_x < middle_left_bound:
                                roi = 1
                            elif avg_x > middle_right_bound:
                                roi = 2
                                
                            # Debugging
                            print(f"Detected Fixation: avg_x={avg_x}, middle_left_bound={middle_left_bound}, middle_right_bound={middle_right_bound}, ROI={roi}")

                            if roi is not None:
                                fixation_events.append({
                                    "start_time": fixation_start_time,
                                    "end_time": fixation_end_time,
                                    "duration": fixation_duration,
                                    "roi": roi,
                                    "avg_x": avg_x
                                })

                            fixation_start_time = None

                    # save trial data
                    if current_trial_id is not None:
                        trial_data.append({
                        "TrialID": current_trial_id,
                        "Task": current_task_message,
                        "DwellLeft": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1),
                        "DwellRight": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2),
                        "DwellTotal": sum(fix["duration"] for fix in fixation_events),
                        "Nfix": len(fixation_events),
                        "FirstFixLoc": fixation_events[0]["roi"] if fixation_events else None,
                        "FirstFixDur": fixation_events[0]["duration"] if fixation_events else None,
                        "FinalFixLoc": fixation_events[-1]["roi"] if fixation_events else None,
                        "FinalFixDur": fixation_events[-1]["duration"] if fixation_events else None,
                        "MiddleFixDur": max(
                            0,
                            sum(fix["duration"] for fix in fixation_events) -
                            (fixation_events[0]["duration"] + fixation_events[-1]["duration"])
                            if len(fixation_events) > 1 else 0
                        ),
                        "eachMiddleFixDur": (
                                        max(
                                            0,
                                            sum(fix["duration"] for fix in fixation_events) -
                                            (fixation_events[0]["duration"] + fixation_events[-1]["duration"])) / max(1, len(fixation_events) - 2)
                                        if len(fixation_events) > 2 else 0
                                    ),
                        "MiddleDominantLoc": mid_dom_loc,      
                        "MiddleDominantDur": mid_dom_dur,    
                        
                        "GazeSwitch": sum(
                            1 for i in range(1, len(fixation_events))
                            if fixation_events[i]["roi"] != fixation_events[i - 1]["roi"]
                        ),
                        "LeftFixNR": sum(1 for fix in fixation_events if fix["roi"] == 1),
                        "RightFixNR": sum(1 for fix in fixation_events if fix["roi"] == 2),
                        "GazeDiff": abs(
                            sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1) -
                            sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2)
                        ),
                        "last_roi": fixation_events[-1]["roi"] if fixation_events else None,
                        "Fixations": fixation_events,
                        "Messages": current_task_message,

                    })
#------------------------------------------------------------------------------------------------------------------------------------------------                    
#------------------------------------------------------------------------------------------------------------------------------------------------
                elif mode == 2:
                    print("Running Code NR 2 (Feedback Phase)...")
                    trial_data = {}
                    fixation_events = []  # different approach; collect fixations for the entire file and filter them for feedback
                    feedback_start_time = None
                    feedback_end_time = None
                    current_trial_id = None
                    current_task_message = None
                    fixation_start_time = None

                    # valid trial ranges - matchin gwith LE phase
                    valid_trial_ranges = list(range(1, 25)) + list(range(49, 210))

                    with open(file_path, 'r') as asc:
                        for line in asc:
                            trial_match = re.search(r'MSG\s+(\d+)\s+TRIALID\s+(\d+)', line)
                            if trial_match:
                                current_trial_id = int(trial_match.group(2))
                                trial_start_time = int(trial_match.group(1))

                                if current_trial_id not in valid_trial_ranges:
                                    current_trial_id = None
                                    continue

                                # Initializing data structure for that trial
                                trial_data[current_trial_id] = {
                                    "DwellLeft_Feed": 0,
                                    "DwellRight_Feed": 0,
                                    "DwellTotal_Feed": 0,
                                    "Nfix_Feed": 0,
                                    "FirstFixLoc_Feed": None,
                                    "FinalFixLoc_Feed": None,
                                    "FirstFixDur_Feed": 0,
                                    "FinalFixDur_Feed": 0,
                                    "MiddleFixDur_Feed": 0,
                                    "eachMiddleFixDur_Feed": 0,
                                    "GazeSwitch_Feed": 0,
                                    "LeftFixNR_Feed": 0,
                                    "RightFixNR_Feed": 0,
                                    "GazeDiff_Feed": 0,
                                    "Fixations_Feed": [],
                                    "last_roi_Feed": None,
                                    "Messages": None,
                                    "Task": None,
                                }
                                # Reset 
                                feedback_start_time = None
                                feedback_end_time = None
                                fixation_events = []
                                continue

                            msg_match = re.search(r'MSG\s+\d+\s+!V TRIAL_VAR\s+(.*)', line)
                            if msg_match and current_trial_id is not None:
                                current_task_message = msg_match.group(1).strip()
                                trial_data[current_trial_id]["Messages"] = current_task_message
                                trial_data[current_trial_id]["Task"] = assign_phase(current_task_message)
                                continue
                            
                            # filter for FeedbackOn
                            feedback_on_match = re.search(r'MSG\s+(\d+)\s+FeedbackOn', line)
                            if feedback_on_match and current_trial_id:
                                feedback_start_time = int(feedback_on_match.group(1))
                                continue
                           
                            # filter for FeedbackOff
                            feedback_off_match = re.search(r'MSG\s+(\d+)\s+FeedbackOff', line)
                            if feedback_off_match and current_trial_id:
                                feedback_end_time = int(feedback_off_match.group(1))

                                if feedback_start_time and fixation_events:
                                    trial_info = trial_data[current_trial_id]

                                    for fix in fixation_events:
                                        if (fix["start_time"] >= feedback_start_time and
                                            fix["end_time"] <= feedback_end_time):

                                            trial_info["Nfix_Feed"] += 1
                                            trial_info["Fixations_Feed"].append(fix)
                                            roi = fix["roi"]
                                            duration = fix["duration"]

                                            if roi == 1:
                                                trial_info["DwellLeft_Feed"] += duration
                                                trial_info["LeftFixNR_Feed"] += 1
                                            elif roi == 2:
                                                trial_info["DwellRight_Feed"] += duration
                                                trial_info["RightFixNR_Feed"] += 1

                                            # First fix
                                            if trial_info["FirstFixLoc_Feed"] is None:
                                                trial_info["FirstFixLoc_Feed"] = roi
                                                trial_info["FirstFixDur_Feed"] = duration
                                            #  update "final" to the last one processed
                                            trial_info["FinalFixLoc_Feed"] = roi
                                            trial_info["FinalFixDur_Feed"] = duration

                                            # Gaze Switch (compare to last ROI)
                                            if (trial_info["last_roi_Feed"] is not None 
                                                and trial_info["last_roi_Feed"] != roi):
                                                trial_info["GazeSwitch_Feed"] += 1

                                            # Update last ROI
                                            trial_info["last_roi_Feed"] = roi

                                    # Summarize dwell totals
                                    trial_info["DwellTotal_Feed"] = (trial_info["DwellLeft_Feed"]
                                                                 + trial_info["DwellRight_Feed"])
                                    first_fix_dur = trial_info["FirstFixDur_Feed"] or 0
                                    final_fix_dur = trial_info["FinalFixDur_Feed"] or 0
                                    trial_info["MiddleFixDur_Feed"] = max(
                                        0,
                                        trial_info["DwellTotal_Feed"] - (first_fix_dur + final_fix_dur)
                                    )
                                    trial_info["eachMiddleFixDur_Feed"] = (
                                        trial_info["MiddleFixDur_Feed"] / max(1, trial_info["Nfix_Feed"] - 2)
                                        if trial_info["Nfix_Feed"] > 2 else 0
                                        )
                                    trial_info["GazeDiff_Feed"] = abs(trial_info["DwellLeft_Feed"]
                                                                  - trial_info["DwellRight_Feed"])

                                # Reset 
                                fixation_events = []
                                feedback_start_time = None
                                feedback_end_time = None
                                continue

                            # Detect fixation start (SFIX)
                            sfix_match = re.search(r'SFIX\s+\w+\s+(\d+)', line)
                            if sfix_match and feedback_start_time:
                                fixation_start_time = int(sfix_match.group(1))
                                continue

                            # Detect fixation end (EFIX)
                            efix_match = re.search(r'EFIX\s+\w+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                            if efix_match and current_trial_id and feedback_start_time:
                                if fixation_start_time is None:
                                    continue

                                fixation_end_time = int(efix_match.group(2))
                                fixation_duration = int(efix_match.group(3))
                                avg_x = float(efix_match.group(4))
                                avg_y = float(efix_match.group(5))

                                # If we already have feedback_end_time, skip fixations that go beyond it
                                if feedback_end_time and fixation_end_time > feedback_end_time:
                                    continue

                                # ROI
                                roi = None
                                if avg_x < middle_left_bound:
                                    roi = 1  # Left ROI
                                elif avg_x > middle_right_bound:
                                    roi = 2  # Right ROI

                                if roi is not None:
                                    # Collect the fixation; we will filter it by time window on FeedbackOff
                                    fixation_events.append({
                                        "start_time": fixation_start_time,
                                        "end_time": fixation_end_time,
                                        "duration": fixation_duration,
                                        "roi": roi,
                                        "avg_x": avg_x,
                                        "avg_y": avg_y
                                    })

                                fixation_start_time = None  # Reset
                                
                                
                    # Convert trial_data to a list of dictionaries
                    trial_data_list = []
                    for trial_id, data in trial_data.items():
                        trial_data_list.append({"TrialID": trial_id, **data})

                    trial_df = pd.DataFrame(trial_data_list)

                    # # Ensure required columns exist
                    # if "Messages" not in trial_df.columns:
                    #     trial_df["Messages"] = ""
                    # if "Task" not in trial_df.columns:
                    #     trial_df["Task"] = "Unknown"

#-------------------------------------------------------------------------------------------------------------------------------------------------------*32                            
                if mode == 3:
                    print("Running Code NR 3 (all fixations without restrictions per trial)...")

                    for line in asc:
                        # Detect trial ID
                        trial_match = re.search(r'MSG\s+(\d+)\s+TRIALID\s+(\d+)', line)
                        if trial_match:
                            # Save trial data
                            if current_trial_id is not None:
                                trial_data.append({
                                    "TrialID": current_trial_id,
                                    "Task": current_task_message,
                                    "DwellLeft_allfix": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1),
                                    "DwellRight_allfix": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2),
                                    "DwellTotal_allfix": sum(fix["duration"] for fix in fixation_events),
                                    "Nfix_allfix": len(fixation_events),
                                    "FirstFixLoc_allfix": fixation_events[0]["roi"] if fixation_events else None,
                                    "FirstFixDur_allfix": fixation_events[0]["duration"] if fixation_events else None,
                                    "FinalFixLoc_allfix": fixation_events[-1]["roi"] if fixation_events else None,
                                    "FinalFixDur_allfix": fixation_events[-1]["duration"] if fixation_events else None,
                                    "MiddleFixDur_allfix": max(
                                        0,
                                        sum(fix["duration"] for fix in fixation_events) -
                                        (fixation_events[0]["duration"] + fixation_events[-1]["duration"])
                                        if len(fixation_events) > 1 else 0
                                    ),
                                    "eachMiddleFixDur_allfix": (
                                        max(
                                            0,
                                            sum(fix["duration"] for fix in fixation_events) -
                                            (fixation_events[0]["duration"] + fixation_events[-1]["duration"])) / max(1, len(fixation_events) - 2)
                                        if len(fixation_events) > 2 else 0
                                        ),
                                    "GazeSwitch_allfix": sum(
                                        1 for i in range(1, len(fixation_events))
                                        if fixation_events[i]["roi"] != fixation_events[i - 1]["roi"]
                                    ),
                                    "LeftFixNR_allfix": sum(1 for fix in fixation_events if fix["roi"] == 1),
                                    "RightFixNR_allfix": sum(1 for fix in fixation_events if fix["roi"] == 2),
                                    "GazeDiff_allfix": abs(
                                        sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1) -
                                        sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2)
                                    ),
                                    "last_roi_allfix": fixation_events[-1]["roi"] if fixation_events else None,
                                    "Messages": current_task_message,
                                    "Fixations_allfix": fixation_events,

                                })

                            # Start a new trial
                            current_trial_id = int(trial_match.group(2))
                            fixation_events = []
                            current_task_message = None
                            continue

                        # Detect task message
                        msg_match = re.search(r'MSG\s+\d+\s+!V TRIAL_VAR\s+(.*)', line)
                        if msg_match:
                            current_task_message = msg_match.group(1).strip()
                            continue

                        # Detect fixation start (SFIX)
                        sfix_match = re.search(r'SFIX\s+\w+\s+(\d+)', line)
                        if sfix_match:
                            fixation_start_time = int(sfix_match.group(1))
                            continue

                        # Detect fixation end (EFIX)
                        efix_match = re.search(r'EFIX\s+\w+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
                        if efix_match and fixation_start_time is not None:
                            fixation_end_time = int(efix_match.group(2))
                            fixation_duration = int(efix_match.group(3))
                            avg_x = float(efix_match.group(4))
                            roi = None

                            # Determine ROI
                            if avg_x < middle_left_bound:
                                roi = 1
                            elif avg_x > middle_right_bound:
                                roi = 2

                            if roi is not None:
                                fixation_events.append({
                                    "start_time": fixation_start_time,
                                    "end_time": fixation_end_time,
                                    "duration": fixation_duration,
                                    "roi": roi,
                                })

                            fixation_start_time = None

                    # Save the last trial's data
                    if current_trial_id is not None:
                        trial_data.append({
                        "TrialID": current_trial_id,
                        "Task": current_task_message,
                        "DwellLeft_allfix": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1),
                        "DwellRight_allfix": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2),
                        "DwellTotal_allfix": sum(fix["duration"] for fix in fixation_events),
                        "Nfix_allfix": len(fixation_events),
                        "FirstFixLoc_allfix": fixation_events[0]["roi"] if fixation_events else None,
                        "FirstFixDur_allfix": fixation_events[0]["duration"] if fixation_events else None,
                        "FinalFixLoc_allfix": fixation_events[-1]["roi"] if fixation_events else None,
                        "FinalFixDur_allfix": fixation_events[-1]["duration"] if fixation_events else None,
                        "MiddleFixDur_allfix": max(
                            0,
                            sum(fix["duration"] for fix in fixation_events) -
                            (fixation_events[0]["duration"] + fixation_events[-1]["duration"])
                            if len(fixation_events) > 1 else 0
                        ),
                        "eachMiddleFixDur_allfix": (
                                        max(
                                            0,
                                            sum(fix["duration"] for fix in fixation_events) -
                                            (fixation_events[0]["duration"] + fixation_events[-1]["duration"])) / max(1, len(fixation_events) - 2)
                                        if len(fixation_events) > 2 else 0
                                        ),
                        "GazeSwitch_allfix": sum(
                            1 for i in range(1, len(fixation_events))
                            if fixation_events[i]["roi"] != fixation_events[i - 1]["roi"]
                        ),
                        "LeftFixNR_allfix": sum(1 for fix in fixation_events if fix["roi"] == 1),
                        "RightFixNR_allfix": sum(1 for fix in fixation_events if fix["roi"] == 2),
                        "GazeDiff_allfix": abs(
                            sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1) -
                            sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2)
                        ),
                        "last_roi_allfix": fixation_events[-1]["roi"] if fixation_events else None,
                        "Fixations_allfix": fixation_events,
                        "Messages": current_task_message,

                        })


#------------------------------------------------------------------------------------------------------------
            # trial_data_list = []
            # for trial_id, data in trial_data.items():
            #     trial_data_list.append({"TrialID": trial_id, **data})  # Add TrialID to each row

            # trial_df = pd.DataFrame(trial_data_list)

            # if "Messages" not in trial_df.columns:
            #     trial_df["Messages"] = ""
            # if "Task" not in trial_df.columns:
            #     trial_df["Task"] = "Unknown"
            
            
            
            trial_df = pd.DataFrame(trial_data)
            trial_df["Messages"] = trial_df["Messages"].fillna("")
            trial_df["Phase"] = trial_df["Messages"].apply(assign_phase)

            output_file = os.path.join(root,
                           f'{participant_id}_{"eye_response_phase" if mode == 1 else "eye_feedback_phase" if mode == 2 else "eye_allfix_phase"}.csv')
            trial_df.to_csv(output_file, index=False)
            print(f"Saved processed file {output_file}")



























# import os
# import re
# import pandas as pd

# mode = 1  # Mode 1 for trial-to-response period processing

# # Define base BIDS directory
# bids_root = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/'

# # Screen parameters
# image_size = [200, 200]
# image_distance = 300
# screen_width = 1919
# screen_height = 1079
# screen_center = (screen_width / 2, screen_height / 2)
# middle_margin = 50
# middle_left_bound = screen_center[0] - middle_margin
# middle_right_bound = screen_center[0] + middle_margin

# def assign_phase(messages):
#     # Handle None or empty strings
#     if not messages:
#         return "Unknown"
#     messages = messages.strip().lower()
#     if "training 0" in messages and "task 1" in messages:
#         return "LE"
#     elif "training 0" in messages and "task 2" in messages:
#         return "ES"
#     elif "training 0" in messages and "task 3" in messages:
#         return "EE"
#     else:
#         return "T"  # Default to training if no match


# # Loop through participants' directories
# for root, dirs, files in os.walk(bids_root):
#     for file in files:
#         if file.endswith('.asc'):
#             participant_id = os.path.basename(os.path.dirname(root))
#             file_path = os.path.join(root, file)
#             print(f"Processing file: {file_path} for participant: {participant_id}")

#             trial_data = []
#             current_trial_id = None
#             current_task_message = None
#             fixation_events = []
#             trial_start_time = None

#             # Open and process the .asc file
#             with open(file_path, 'r') as asc:
#                 for line in asc:
#                     # Detect trial ID and task message
#                     trial_match = re.search(r'MSG\s+(\d+)\s+TRIALID\s+(\d+)', line)
#                     if trial_match:
#                         # Save previous trial data if applicable
#                         if current_trial_id is not None:
#                             trial_data.append({
#                                 "TrialID": current_trial_id,
#                                 "Task": current_task_message,
#                                 "Fixations": fixation_events,
#                                 "DwellLeft": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1),
#                                 "DwellRight": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2),
#                                 "DwellTotal": sum(fix["duration"] for fix in fixation_events),
#                                 "Nfix": len(fixation_events),
#                                 "Messages": current_task_message,
#                             })

#                         # Start a new trial
#                         current_trial_id = int(trial_match.group(2))
#                         fixation_events = []
#                         current_task_message = None
#                         continue

#                     # Detect task message
#                     msg_match = re.search(r'MSG\s+\d+\s+!V TRIAL_VAR\s+(.*)', line)
#                     if msg_match:
#                         current_task_message = msg_match.group(1).strip()
#                         continue

#                     # Detect fixation start (SFIX)
#                     sfix_match = re.search(r'SFIX\s+\w+\s+(\d+)', line)
#                     if sfix_match:
#                         fixation_start_time = int(sfix_match.group(1))
#                         continue

#                     # Detect fixation end (EFIX)
#                     efix_match = re.search(r'EFIX\s+\w+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
#                     if efix_match and fixation_start_time is not None:
#                         fixation_end_time = int(efix_match.group(2))
#                         fixation_duration = int(efix_match.group(3))
#                         avg_x = float(efix_match.group(4))
#                         roi = None

#                         # Determine ROI
#                         if avg_x < middle_left_bound:
#                             roi = 1
#                         elif avg_x > middle_right_bound:
#                             roi = 2

#                         if roi is not None:
#                             fixation_events.append({
#                                 "start_time": fixation_start_time,
#                                 "end_time": fixation_end_time,
#                                 "duration": fixation_duration,
#                                 "roi": roi,
#                             })

#                         fixation_start_time = None

#                 # Save the last trial data
#                 if current_trial_id is not None:
#                     trial_data.append({
#                         "TrialID": current_trial_id,
#                         "Task": current_task_message,
#                         "Fixations": fixation_events,
#                         "DwellLeft": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 1),
#                         "DwellRight": sum(fix["duration"] for fix in fixation_events if fix["roi"] == 2),
#                         "DwellTotal": sum(fix["duration"] for fix in fixation_events),
#                         "Nfix": len(fixation_events),
#                         "Messages": current_task_message,
#                     })

#             # Convert trial data to a DataFrame
#             # Convert trial data to a DataFrame
#             trial_df = pd.DataFrame(trial_data)
#             trial_df["Messages"] = trial_df["Messages"].fillna("")  # Ensure no None values in Messages
#             trial_df["Phase"] = trial_df["Messages"].apply(assign_phase)

#             # Save the processed data
#             output_file = os.path.join(
#             root, f'{participant_id}_eye_response_phase.csv'
#             )
#             trial_df.to_csv(output_file, index=False)
#             print(f"Saved processed file {output_file}")


























# import os
# import re
# import pandas as pd

# # Define base BIDS directory
# bids_root = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/'

# # Helper function to assign phase based on messages
# def assign_phase(messages):
#     if not messages:
#         return "Unknown"
#     messages = messages.strip().lower()
#     if "training 0" in messages and "task 1" in messages:
#         return "LE"
#     elif "training 0" in messages and "task 2" in messages:
#         return "ES"
#     elif "training 0" in messages and "task 3" in messages:
#         return "EE"
#     else:
#         return "T"  # Default to training if no match

# # Initialize the processed trials list
# processed_trials = []

# # Loop through participants' directories
# for root, dirs, files in os.walk(bids_root):
#     for file in files:
#         if file.endswith(".asc"):
#             participant_id = os.path.basename(os.path.dirname(root))
#             file_path = os.path.join(root, file)
#             print(f"Processing file: {file_path} for participant: {participant_id}")

#             current_trial_id = None
#             trial_messages = {}

#             # Open and process the .asc file
#             with open(file_path, "r") as asc:
#                 for line in asc:
#                     # Detect trial IDs
#                     trial_match = re.search(r"MSG\s+(\d+)\s+TRIALID\s+(\d+)", line)
#                     if trial_match:
#                         current_trial_id = int(trial_match.group(2))
#                         if current_trial_id not in trial_messages:
#                             trial_messages[current_trial_id] = []
#                         continue

#                     # Collect trial messages
#                     msg_match = re.search(r"MSG\s+\d+\s+!V TRIAL_VAR\s+(.*)", line)
#                     if msg_match and current_trial_id is not None:
#                         trial_messages[current_trial_id].append(msg_match.group(1).strip())

#             # Process trial messages into rows
#             for trial_id, messages_list in trial_messages.items():
#                 for message in messages_list:
#                     # Split the message into separate entries if it contains multiple messages
#                     split_messages = message.split("Training 0")
#                     for sub_message in split_messages:
#                         if sub_message.strip():  # Ignore empty parts
#                             full_message = "Training 0 " + sub_message.strip() if "Training 0" not in sub_message else sub_message.strip()
#                             processed_trials.append({
#                                 "TrialID": trial_id,
#                                 "Messages": full_message.strip(),
#                                 "Phase": assign_phase(full_message.strip())
#                             })

# # Create a DataFrame from the processed trial data
# trial_df = pd.DataFrame(processed_trials)

# # Save the processed data
# output_file = os.path.join(
#     root, f"{participant_id}_eye_response_phase.csv"
# )
# trial_df.to_csv(output_file, index=False)
# print(f"Saved processed file {output_file}")





# import os
# import re
# import pandas as pd

# # Define base BIDS directory
# bids_root = 'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/'

# # Helper function to assign phase based on messages
# def assign_phase(messages):
#     if not messages:
#         return "Unknown"
#     messages = messages.strip().lower()
#     if "training 0" in messages and "task 1" in messages:
#         return "LE"
#     elif "training 0" in messages and "task 2" in messages:
#         return "ES"
#     elif "training 0" in messages and "task 3" in messages:
#         return "EE"
#     else:
#         return "T"  # Default to training if no match

# # Initialize the processed trials list
# processed_trials = []

# # Loop through participants' directories
# for root, dirs, files in os.walk(bids_root):
#     for file in files:
#         if file.endswith(".asc"):
#             # Get participant ID based on the parent folder of the file
#             participant_id = os.path.basename(os.path.dirname(root))
            
#             # Construct participant's directory path
#             participant_folder = os.path.dirname(os.path.dirname(root))
            
#             file_path = os.path.join(root, file)
#             print(f"Processing file: {file_path} for participant: {participant_id}")

#             current_trial_id = None
#             trial_messages = {}

#             # Open and process the .asc file
#             with open(file_path, "r") as asc:
#                 for line in asc:
#                     # Detect trial IDs
#                     trial_match = re.search(r"MSG\s+(\d+)\s+TRIALID\s+(\d+)", line)
#                     if trial_match:
#                         current_trial_id = int(trial_match.group(2))
#                         if current_trial_id not in trial_messages:
#                             trial_messages[current_trial_id] = []
#                         continue

#                     # Collect trial messages
#                     msg_match = re.search(r"MSG\s+\d+\s+!V TRIAL_VAR\s+(.*)", line)
#                     if msg_match and current_trial_id is not None:
#                         trial_messages[current_trial_id].append(msg_match.group(1).strip())

#             # Process trial messages into rows
#             for trial_id, messages_list in trial_messages.items():
#                 for message in messages_list:
#                     # Split the message into separate entries if it contains multiple messages
#                     split_messages = message.split("Training 0")
#                     for sub_message in split_messages:
#                         if sub_message.strip():  # Ignore empty parts
#                             full_message = "Training 0 " + sub_message.strip() if "Training 0" not in sub_message else sub_message.strip()
#                             processed_trials.append({
#                                 "TrialID": trial_id,
#                                 "Messages": full_message.strip(),
#                                 "Phase": assign_phase(full_message.strip())
#                             })

#             # Create a DataFrame from the processed trial data
#             trial_df = pd.DataFrame(processed_trials)

#             # Save the processed data to the participant's folder
#             output_file = os.path.join(participant_folder, f"{participant_id}_eye_response_phase.csv")
#             trial_df.to_csv(output_file, index=False)
#             print(f"Saved processed file {output_file}")
