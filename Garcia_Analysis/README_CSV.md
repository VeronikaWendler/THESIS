# Readme file for understanding the CSV file corresponding to **Study 1**  
`GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv`

This folder contains raw behavioural data from Study 1.

This means the CSV file above contains already preprocessed column data – ready to use!  
If you have questions regarding the preprocessing please let me know, at this stage I still need to clean some code until I can put it up on GitHub.

---

## Columns in the DataFrame

### Core behavioural variables

| Column | Explanation |
|------|------------|
| `index` | row index |
| `sub_id` | participant identifier |
| `phase` | LE = Learning phase, ES = Experiential-Symbolic phase, EE = Experiential-Experiential phase, SP = Stated Prob. Phase = no useful data in this frame |
| `p1` | value of the option corresponding to op1 |
| `p2` | value of the option corresponding to op2 |
| `rtime` | RT |
| `out` | outcome of the chosen option |
| `cfout` | outcome of the unchosen option – counterfactual outcome |
| `cho` | side of option chosen – 1 for left/p1 – 2 for right/p2 |
| `corr` | 1 if correct, 0 if choice was incorrect == Accuracy |
| `trial` | trial identifier – resets with every phase |
| `cond` | probability condition classifier only used in LE phase – 3 for hardest (0.6/0.4), 2 for (0.7/0.3), 1 for (0.8/0.2), 0 for (0.9/0.1) |
| `chose_right` | 0 if left chosen, 1 if right chosen |
| `rew` | cumulative reward |
| `sess` | not used |
| `op1` | E option |
| `op2` | E option or S option (**IMPORTANT** – for analysis purposes, data from ES-phase of the experiment was rearranged so that the S-option is always on the right in the dataframe (not in the real experiment), all other columns account for this change) |
| `ev1` | expected value of the option corresponding to op1 |
| `ev2` | expected value of the option corresponding to op2 |
| `catch_trial` | not used |
| `reversed` | not used |

---

## Response fixation data

These are fixations measured from the start of the trial until the participant presses a response  
(if EFIX – end of fixation identifier – falls slightly after the response was made and SFIX (start of fixation identifier) is still within the valid range it still counts as a valid fixation).

| Column | Explanation |
|------|------------|
| `TrialID` | corresponds to trial identifier from the eye-tracking .asc file – starts at 49 – everything beforehand is training |
| `Task` | task message from the eye-tracking file – Task 1 == LE, 2 == ES, 3 == EE |
| `DwellLeft` | Dwell time on the left option (fixation time on left opt) |
| `DwellRight` | Dwell time on the right option (fixation time on right opt) |
| `DwellTotal` | Total dwell time (DwellRight + DwellLeft) |
| `Nfix` | number of fixations |
| `FirstFixLoc` | Location of first fixation (1 = left, 2 = right) |
| `FirstFixDur` | Duration of the first fixation |
| `FinalFixLoc` | Location of final fixation (1 = left, 2 = right) |
| `FinalFixDur` | Duration of final fixation |
| `MiddleFixDur` | Duration of middle fixation (total) |
| `eachMiddleFixDur` | Duration of an average middle fixation |
| `GazeSwitch` | Number of gaze switches |
| `LeftFixNR` | Number of fixations to the left |
| `RightFixNR` | Number of fixations to the right |
| `GazeDiff` | Absolute difference in gaze time between the options |
| `last_roi` | in which ROI was the last fixation on trial? Bit redundant – Left = 1, Right = 2 |
| `Messages` | not used |
| `Phase` | not used |

---

## Feedback fixation data

Columns for feedback eye-tracking measurements are the same as for response measurements, except for the fact that these fixations were measured from the start of **Feedback onset** (participant sees the outcome) until the **Feedback offset** (feedback disappears).

Same columns, just with a `_Feed` ending.

---

## Allfix data

I did not use these as they contain absolutely all fixations within a trial period – even those from **Feedback offset until just before the next trial** (in these periods the screen is white and I am not interested in this).

Simply contains too much noise but does not differ much from the `Response` fixations in ES and EE phases.  
Same columns, just with a `_allfix` ending.

---

## Columns useful for creating the aDDM calculations

| Column | Explanation |
|------|------------|
| `ev_opt` | expected value of the optimal option |
| `ev_sub` | expected value of the suboptimal option |
| `DwellTime_opt` | Dwell time (that is total fixation duration) on the optimal option |
| `DwellTime_sub` | Dwell time on the suboptimal option |
| `PropDwell_opt` | proportion of fixation time on optimal option |
| `PropDwell_sub` | proportion of fixation time on suboptimal option |
| `V_corr` | value of the optimal option (either p1 or p2) |
| `V_sub` | value of the suboptimal option (either p1 or p2) |
| `AttentionW` | Attentional weight parameter of the aDDM – calculated as `(PropDwell_opt × V_corr) − (PropDwell_sub × V_sub)` |
| `InattentionW` | Inattentional weight parameter of the aDDM – calculated as `(PropDwell_sub × V_corr) − (PropDwell_opt × V_sub)` |
| `BetterOption` | np.where(data['ev1'] > data['ev2'], 1, 2) – 1 for left, 2 for right |
| `FixLocFirstCorr` | 1 if the first fixation was to the correct option, 0 otherwise |
| `FixLocLastCorr` | 1 if the last fixation was to the correct option, 0 otherwise |
| `DwellTimeAdvantage` | dwell time on the right option − dwell time on the left option |
| `PropDwell_Right` | proportion of dwell time on the right option = left / total |
| `PropDwell_Left` | proportion of dwell time on the left option = right / total |
| `FixationAdvantage` | fixation count right − fixation count left |
| `DwellPropAdvantageCorrect` | proportion of dwell time on correct option − proportion of dwell time on incorrect option |

---

## Overall value (OV) and value-difference variables

| Column | Explanation |
|------|------------|
| `OV_num` | raw overall value (p1 + p2) |
| `OV_num_2` | same as OV_num – bit redundant |
| `OVcate` | categorical overall value identifier (low, medium, high) – according to quantiles – used pandas `.qcut(q=3)` – would not suggest using this binning |
| `OVcate_2` | categorical overall value identifier (low, medium, high) – according to tertile binning – used pandas `.qcut(q=4)` – more equal |
| `AbsValueDiff` | absolute value difference (raw, trial-by-trial) |
| `AbsValueDiff_2` | absolute value difference – redundant |
| `Abscate` | categorical value difference – would not use since quantiles |
| `Abscate_2` | categorical value difference – would use since tertile binning |
| `VD` | numeric value difference – 1 = low, 2 = medium, 3 = high |
| `VD_2` | numeric value difference – 1 = low, 2 = medium, 3 = high – would use |
| `RLdiff` | p1 − p2 rounded to 1 dec. |
| `feedback` | Feedback received on every trial (1 = +1 point, 0 = −1 point) |
| `split_by` | same as `cond`, just integer |
| `q_init` | initial q-value |
| `stim_chosen` | categorical – stimulus chosen, E or S |
| `OV` | numeric OV – would not use since quantile binning |
| `OV_2` | numeric OV – would use |

---

## OLD

| Variable | Description |
|---------|-------------|
| `index` | Index identifier for each row |
| `sub_id` | Subject identifier, identifying the participant |
| `phase` | Phase of the experiment (LE, ES, EE, SP) |
| `p1` | p(win) option 1 |
| `p2` | p(win) option 2 |
| `rtime` | Response time |
| `out` | Outcome of the choice (selected option) |
| `cfout` | Counterfactual outcome (unchosen option outcome) |
| `cho` | Choice made by the participant (1 or 2) |
| `corr` | Correctness of the choice (1 if expected value of chosen option ≥ expected value of unchosen option, 0 else) |
| `trial` | Trial number |
| `cond` | Condition of the experiment (3 = 60/40, 2 = 70/30, 1 = 80/20, 0 = 90/10) |
| `chose_right` | Whether the participant chose the rightmost option on screen |
| `rew` | Total reward received (cumulated) |
| `sess` | Session identifier (when `sess` is −1 or −2 it was a training session; when `sess` is 0 or 1 it was the first or second session) |
| `op1` | Option 1 presented to the participant (filename or identifier) |
| `op2` | Option 2 presented to the participant (filename or identifier) |
| `ev1` | Expected value of option 1 |
| `ev2` | Expected value of option 2 |
| `catch_trial` | Indicates if it’s a catch trial (testing attention) |
| `reversed` | nothing is reversed |

