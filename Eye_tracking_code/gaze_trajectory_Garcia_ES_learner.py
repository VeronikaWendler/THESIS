# Veronika Wendler
# 19.06.2025
# program that creates gaze trajectories during the learning phase (similar to plots in Checchi et al.,2025) to understand dynamic gaze allocation from stimuöi preseontation to outcome fixation

#libraries
import pandas as pd
import numpy as np 
import ast
import matplotlib.pyplot as plt
import itertools
from matplotlib.lines import Line2D
import os
from collections import defaultdict


# big file that has everything
df_raw = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv")

# figures dir
fig_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Figures/gaze_trajectory"
if not fig_dir:
    os.mkdir(fig_dir)

BIN_MS    = 25                       # bin width in ms
ROI   = {1: 'left', 2: 'right'}      # ROI

raw = df_raw.query("phase == 'ES'").copy()   # keep only Learning‐phase trials
rt_glob_mean = raw.rtime.mean() * 1000          
rt_glob_sem  = raw.rtime.std(ddof=1) / np.sqrt(len(raw)) * 1000
raw['corr'] = raw['corr'].astype(int)


# Define your groups
#poor_learners = [2,3,10,12,13,15,16,23,24,25]
#good_learners = [7,8,9,11,17,18,19,20,21,22,26]
#exclude_sub = [1, 4, 5, 6, 14, 99]

# Remove any excluded IDs from good learners
#good_learners_filtered = [i for i in poor_learners if i not in exclude_sub]

# Now filter the DataFrame to keep only these participants
#raw = raw[raw['sub_id'].isin(good_learners_filtered)]


# helper for response phase -----------------------------------------------------------------------

def trial_bins(row):
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    hi_side = 'left' if row.p1 > row.p2 else 'right'
    t0      = fixes[0]['start_time']
    bins    = {}

    for f in fixes:
        side   = ROI.get(f['roi'],'none')
        is_hi  = side == hi_side
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0,0]          # [high, low]
            bins[b][0 if is_hi else 1] += 1

    out = pd.DataFrame.from_dict(bins, orient='index', columns=['hi','lo'])
    out['sub'] = row.sub_id
    out['trial'] = row.TrialID
    out['bin'] = out.index
    out['hi_ratio'] = out.hi / (out.hi+out.lo)
    return out[['sub','trial','bin','hi_ratio']]


# helper for ES phase -----------------------------------------------------------------------
def trial_bins_ES(row):
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side = 'left'
    t0 = fixes[0]['start_time']
    bins = {}

    for f in fixes:
        side  = ROI.get(f['roi'], 'none')
        is_E = side == E_side
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]
            bins[b][0 if is_E else 1] += 1

    out = pd.DataFrame.from_dict(bins, orient='index', columns=['E','S'])
    out['sub'] = row.sub_id
    out['trial'] = row.TrialID
    out['bin'] = out.index
    out['ES_ratio'] = out.E / (out.E+out.S)
    return out[['sub','trial','bin','ES_ratio']]


def trial_bins_EvS_byChoice(row):
    """
    Return per-25 ms bins:
        • E_cnt, S_cnt
        • identify trial group: 'E-chosen' (cho==1) or 'S-chosen' (cho==2)
    """
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side, S_side = 'left', 'right'
    t0, bins = fixes[0]['start_time'], {}

    for f in fixes:
        side = ROI.get(f['roi'], 'none')            # 'left' | 'right' | 'none'
        idx  = 0 if side == E_side else 1           # 0 = E, 1 = S
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]                    # [E_cnt, S_cnt]
            bins[b][idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['E_cnt', 'S_cnt'])
    out['sub']      = row.sub_id
    out['trial']    = row.TrialID
    out['bin']      = out.index
    out['cho_grp']  = 'E-chosen' if row.cho == 1 else 'S-chosen'
    out['ratio_E']  = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S']  = 1 - out['ratio_E']
    return out[['sub','trial','bin','cho_grp','ratio_E','ratio_S']]




def trial_bins_EvS_byAcc(row):
    """
    For every 25-ms bin during ES trials count E- and S-fixations.
    Attach a trial label:
        • 'Correct'    (row.corr == 1)
        • 'Incorrect'  (row.corr == 0)
    """
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side, S_side = 'left', 'right'
    t0, bins = fixes[0]['start_time'], {}

    for f in fixes:
        side = ROI.get(f['roi'], 'none')          # 'left' / 'right' / 'none'
        idx  = 0 if side == E_side else 1         # 0 = E, 1 = S
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]                  # [E_cnt, S_cnt]
            bins[b][idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['E_cnt', 'S_cnt'])
    out['sub']   = row.sub_id
    out['trial'] = row.TrialID
    out['bin']   = out.index

    # ── FIX here ───────────────────────
    out['acc_grp'] = 'Correct' if row['corr'] == 1 else 'Incorrect'
    # ───────────────────────────────────

    out['ratio_E'] = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S'] = 1 - out['ratio_E']
    return out[['sub','trial','bin','acc_grp','ratio_E','ratio_S']]



########################################################################
# === RESPONSE-LOCKED, FIXATION-BY-FIXATION TABLE =====================#
########################################################################
def seq_resp_locked(row):
    try:
        seq = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not seq:
        return None

    t_resp = int(row.rtime * 1000)   # click (ms from stim onset)
    E_side = 'left'

    # ‼️ keep every fixation whose **start** is before the click
    pre_click = [f for f in seq if f['start_time'] < t_resp]
    if not pre_click:
        return None

    # --- build rows, starting with the LAST pre-click fixation -----------
    rows = []
    for i_back, f in enumerate(reversed(pre_click)):   # 0,-1,-2...
        side = ROI.get(f['roi'], 'none')
        rows.append(dict(
            sub          = row.sub_id,
            trial        = row.TrialID,
            rel_fix_idx  = -i_back,                 # 0 = last fixation
            is_E         = int(side == E_side),
            is_S         = int(side != E_side)
        ))
    return pd.DataFrame(rows)

# ── build the master table ────────────────────────────────────────────

fix_parts = raw.apply(seq_resp_locked, axis=1).dropna().tolist()

if not fix_parts:
    raise RuntimeError("No pre-click fixations found – check time units!")

fix_tbl = pd.concat(fix_parts, ignore_index=True)

# how many fixations do we keep?
mean_fix = fix_tbl.groupby(['sub','trial']).size().mean()   # float
N_back   = int(np.ceil(mean_fix))                           # round UP

fix_tbl = fix_tbl[fix_tbl.rel_fix_idx >= -N_back]          

#######################################################################
# === AGGREGATE & SEM ================================================#
#######################################################################
fix_agg = (fix_tbl
           .groupby('rel_fix_idx')
           .agg(E_m   = ('is_E','mean'),
                E_sem = ('is_E', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                S_m   = ('is_S','mean'),
                S_sem = ('is_S', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
           .reset_index()
           .sort_values('rel_fix_idx'))  # ascending: -N … 0



# bin table per trial -------------------------------------------------------------------------------------
tbl = raw.apply(trial_bins, axis=1).dropna()
tbl = pd.concat(tbl.tolist(), ignore_index=True)

# mask bins after each individual RT 
rts = (raw[['sub_id','TrialID','rtime']]
       .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)

tbl = tbl.merge(rts, on=['sub','trial'])
tbl = tbl[tbl.bin <= tbl.rt_bin]           # keep up to the click

# time in ms relative to onset
tbl['time_ms'] = tbl['bin'] * BIN_MS

# subject-average trajectory 
sub_traj = (tbl.groupby(['sub','time_ms']).hi_ratio.mean().reset_index(name='hi_ratio'))

# grand mean & SEM across participants
grand = (sub_traj.groupby('time_ms')
                   .agg(mean=('hi_ratio','mean'),
                        sem =('hi_ratio', lambda x: x.std()/np.sqrt(len(x))))
                   .reset_index())
grand['lo_mean'] = 1 - grand['mean']   # complementary curve
grand['lo_sem' ] = grand['sem']        # same SEM

# global mean RT & its SEM 
rt_mean = raw.rtime.mean()*1000            # ms
rt_sem  = raw.rtime.std(ddof=1)/np.sqrt(len(raw))*1000


# ─── build table for ES phase ─────────────────────────────────────
ES_tbl = raw.apply(trial_bins_ES, axis=1).dropna()
ES_tbl = pd.concat(ES_tbl.tolist(), ignore_index=True)

rts_ES = (raw[['sub_id','TrialID','rtime']]
       .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
rts_ES['rtES_bin'] = (rts_ES.rt*1000 // BIN_MS).astype(int)

ES_tbl = ES_tbl.merge(rts_ES, on=['sub','trial'])
ES_tbl = ES_tbl[ES_tbl.bin <= ES_tbl.rtES_bin]           # keep up to the click

# time in ms relative to onset
ES_tbl['timeES_ms'] = ES_tbl['bin'] * BIN_MS

# subject-average trajectory 
sub_trajES = (ES_tbl.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index(name='ES_ratio'))

# grand mean & SEM across participants
grandES = (sub_trajES.groupby('timeES_ms')
                   .agg(mean=('ES_ratio','mean'),
                        sem =('ES_ratio', lambda x: x.std()/np.sqrt(len(x))))
                   .reset_index())
grandES['S_mean'] = 1 - grandES['mean']   # complementary curve
grandES['S_sem' ] = grandES['sem']        # same SEM

# global mean RT & its SEM 
rtES_mean = raw.rtime.mean()*1000            # ms
rtES_sem  = raw.rtime.std(ddof=1)/np.sqrt(len(raw))*1000

#----------------------------------------------------------------------------
# chocie E vs S

t_cap = int(rt_glob_mean + 100)      # ≈ 1 × mean RT + 100 ms
choice_tbl = raw.apply(trial_bins_EvS_byChoice, axis=1).dropna()
choice_tbl = pd.concat(choice_tbl.tolist(), ignore_index=True)
choice_tbl['time_ms'] = choice_tbl['bin'] * BIN_MS
choice_tbl = choice_tbl[choice_tbl.time_ms <= rt_glob_mean + 100]  # trim

sub_choice = (choice_tbl
              .groupby(['sub','cho_grp','time_ms'])
              .agg(E_mean=('ratio_E','mean'),
                   S_mean=('ratio_S','mean'))
              .reset_index())

grp_choice = (sub_choice
              .groupby(['cho_grp','time_ms'])
              .agg(E_m   = ('E_mean','mean'),
                   E_sem = ('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                   S_m   = ('S_mean','mean'),
                   S_sem = ('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
              .reset_index())

#-----------------------------------------------------------------------
# ----------  accuracy-based E vs S  ---------- #
acc_tbl = raw.apply(trial_bins_EvS_byAcc, axis=1).dropna()
acc_tbl = pd.concat(acc_tbl.tolist(), ignore_index=True)
acc_tbl['time_ms'] = acc_tbl['bin'] * BIN_MS
acc_tbl = acc_tbl[acc_tbl.time_ms <= rt_glob_mean + 100]   # trim to mean RT + 100 ms

sub_acc = (acc_tbl
           .groupby(['sub','acc_grp','time_ms'])
           .agg(E_mean=('ratio_E','mean'),
                S_mean=('ratio_S','mean'))
           .reset_index())

grp_acc = (sub_acc
           .groupby(['acc_grp','time_ms'])
           .agg(E_m   = ('E_mean','mean'),
                E_sem = ('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                S_m   = ('S_mean','mean'),
                S_sem = ('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
           .reset_index())

#-----------------
MAX_FIX = 10           # cut once <70 % of trials still have data

def extract_seq(row, max_fix=MAX_FIX):
    try: seq = ast.literal_eval(row.Fixations) or []
    except Exception: return None
    if not seq: return None

    E_side = 'left'
    rec = []
    for i, f in enumerate(seq[:max_fix], 1):
        side = ROI.get(f['roi'],'none')
        rec.append(dict(sub=row.sub_id, trial=row.TrialID,
                        fix_idx=i,
                        is_E = int(side==E_side),
                        is_S = int(side!=E_side),
                        dur   = f['end_time']-f['start_time']))
    return pd.DataFrame(rec)

fix_tbl = raw.apply(extract_seq, axis=1).dropna()
fix_tbl = pd.concat(fix_tbl.tolist(), ignore_index=True)

# keep fixation positions present in ≥70 % of trials
valid = fix_tbl.groupby('fix_idx').size()
fix_tbl = fix_tbl[fix_tbl.fix_idx.isin(valid[valid>=0.7*raw.TrialID.nunique()].index)]

fix_agg = (fix_tbl
           .groupby('fix_idx')
           .agg(E_m=('is_E','mean'),
                E_sem=('is_E', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                S_m=('is_S','mean'),
                S_sem=('is_S', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                dur_m=('dur','mean'),
                dur_sem=('dur', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
           .reset_index())

# -----------------------------------------------------------------------
# helper: small, tidy legend kwargs
legend_kw = dict(loc='upper right', frameon=False,
                 fontsize=8,  title_fontsize=9,
                 borderpad=0.2, labelspacing=0.25, handlelength=2)

line_handles = [
    Line2D([], [], color='darkorchid',  lw=2, label='higher'),
    Line2D([], [], color='darkorchid', ls='--', lw=2, label='lower')
]

line_handles2 = [
    Line2D([], [], color='royalblue',  lw=2, label='E'),
    Line2D([], [], color='darkorchid', ls='--', lw=2, label='S')
]

# ────────────────────────────────────────────────────────────────────────
# CHOICE  (stimulus-locked) – its own figure
# ────────────────────────────────────────────────────────────────────────
fig_choice, ax_c = plt.subplots(figsize=(4.5, 3))

# RT band  (mean ± SEM)
ax_c.axvspan(rt_mean-rt_sem, rt_mean+rt_sem,
             color='lightgrey', alpha=.5, zorder=0)
# RT mean line on top of the band
ax_c.axvline(rt_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_c.plot(grand.time_ms, grand['mean'],
          color='darkorchid', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['mean']-grand['sem'],
                  grand['mean']+grand['sem'],
                  color='darkorchid', alpha=.25, label='_nolegend_')

# lower-EV
ax_c.plot(grand.time_ms, grand['lo_mean'],
          color='darkorchid', ls='--', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['lo_mean']-grand['lo_sem'],
                  grand['lo_mean']+grand['lo_sem'],
                  color='darkorchid', alpha=.20, label='_nolegend_')

ax_c.axhline(.5, ls='--', color='k')
ax_c.axvline(0,  ls='--', color='k')
ax_c.set(xlim=(0, rt_mean+100), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Choice (stimulus-locked)')

ax_c.legend(handles=line_handles, title='Gaze on', **legend_kw)
fig_choice.tight_layout()
plot_path= os.path.join(fig_dir, "gaze_choice_ES_EV.png")
fig_choice.savefig(plot_path, dpi=300, bbox_inches='tight')
fig_choice.show()

# ────────────────────────────────────────────────────────────────────────
# FEEDBACK  (outcome-locked) – separate figure
# ────────────────────────────────────────────────────────────────────────

fig_choiceES, ax_ES = plt.subplots(figsize=(4.5, 3))

# RT band  (mean ± SEM)
ax_ES.axvspan(rtES_mean-rtES_sem, rtES_mean+rtES_sem,
             color='lightgrey', alpha=.5, zorder=0)
# RT mean line on top of the band
ax_ES.axvline(rtES_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_ES.plot(grandES.timeES_ms, grandES['mean'],
          color='royalblue', lw=2)
ax_ES.fill_between(grandES.timeES_ms,
                  grandES['mean']-grandES['sem'],
                  grandES['mean']+grandES['sem'],
                  color='royalblue', alpha=.25, label='_nolegend_')

# lower-EV
ax_ES.plot(grandES.timeES_ms, grandES['S_mean'],
          color='darkorchid', ls='--', lw=2)
ax_ES.fill_between(grandES.timeES_ms,
                  grandES['S_mean']-grandES['S_sem'],
                  grandES['S_mean']+grandES['S_sem'],
                  color='darkorchid', alpha=.20, label='_nolegend_')

ax_ES.axhline(.5, ls='--', color='k')
ax_ES.axvline(0,  ls='--', color='k')
ax_ES.set(xlim=(0, rtES_mean+100), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Choice E vs S (stimulus-locked)')

ax_ES.legend(handles=line_handles2, title='Gaze on', **legend_kw)
fig_choiceES.tight_layout()
plot_path2= os.path.join(fig_dir, "gaze_choice_ES_ES.png")
fig_choiceES.savefig(plot_path2, dpi=300, bbox_inches='tight')
fig_choiceES.show()

#-----------------------------------------------------------------------------------------------------------------
# ── colours / styles
fig, ax = plt.subplots(figsize=(4.8,3.3))

# RT cue
ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
           color='lightgrey', alpha=.35, zorder=0)
ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=1)

colour_grp = {'E-chosen':'royalblue', 'S-chosen':'darkorange'}

for grp, df in grp_choice.groupby('cho_grp'):
    col = colour_grp[grp]

    # solid  = gaze to E
    ax.plot(df.time_ms, df.E_m,  color=col, lw=2,
            label=f'{grp} – gaze E')
    ax.fill_between(df.time_ms,
                    df.E_m-df.E_sem, df.E_m+df.E_sem,
                    color=col, alpha=.20)

    # dashed = gaze to S
    ax.plot(df.time_ms, df.S_m,  color=col, ls='--', lw=2,
            label=f'{grp} – gaze S')
    ax.fill_between(df.time_ms,
                    df.S_m-df.S_sem, df.S_m+df.S_sem,
                    color=col, alpha=.20)

ax.axhline(.5, ls='--', color='k')
ax.axvline(0,  ls='--', color='k')
ax.set(xlim=(0, rt_glob_mean+100), ylim=(0,1),
       xlabel='Time (ms)', ylabel='Fixation ratio',
       title='ES phase - gaze to E vs S, split by chosen item')
ax.legend(frameon=False, fontsize=7, ncol=1,
          bbox_to_anchor=(1.02,1), loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, 'ES_EvS_byChoice.png'),
            dpi=300, bbox_inches='tight')
plt.show()

#--------------------------------------

fig, ax = plt.subplots(figsize=(4.8,3.3))

# grey RT window + mean
ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
           color='lightgrey', alpha=.35, zorder=0)
ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=1)

colour_grp = {'Correct':'seagreen', 'Incorrect':'crimson'}

for grp, df in grp_acc.groupby('acc_grp'):
    col = colour_grp[grp]

    # solid  = gaze to E
    ax.plot(df.time_ms, df.E_m,  color=col, lw=2,
            label=f'{grp} – gaze E')
    ax.fill_between(df.time_ms,
                    df.E_m-df.E_sem, df.E_m+df.E_sem,
                    color=col, alpha=.20)

    # dashed = gaze to S
    ax.plot(df.time_ms, df.S_m,  color=col, ls='--', lw=2,
            label=f'{grp} – gaze S')
    ax.fill_between(df.time_ms,
                    df.S_m-df.S_sem, df.S_m+df.S_sem,
                    color=col, alpha=.20)

ax.axhline(.5, ls='--', color='k')
ax.axvline(0,  ls='--', color='k')
ax.set(xlim=(0, rt_glob_mean+100), ylim=(0,1),
       xlabel='Time (ms)', ylabel='Fixation ratio',
       title='ES phase - gaze to E vs S, split by accuracy')
ax.legend(frameon=False, fontsize=7, ncol=1,
          bbox_to_anchor=(1.02,1), loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, 'ES_EvS_byAccuracy.png'),
            dpi=300, bbox_inches='tight')
plt.show()


#-----------------------------------------------------------------------------------------
########################################################################
# === PLOT ============================================================#
########################################################################
fig, ax = plt.subplots(figsize=(5.2,3.4))

# solid = gaze E  dashed = gaze S
ax.plot(fix_agg.rel_fix_idx, fix_agg.E_m,   lw=2, c='royalblue')
ax.fill_between(fix_agg.rel_fix_idx,
                fix_agg.E_m-fix_agg.E_sem, fix_agg.E_m+fix_agg.E_sem,
                color='royalblue', alpha=.25, label='gaze E')

ax.plot(fix_agg.rel_fix_idx, fix_agg.S_m,   lw=2, ls='--', c='darkorange')
ax.fill_between(fix_agg.rel_fix_idx,
                fix_agg.S_m-fix_agg.S_sem, fix_agg.S_m+fix_agg.S_sem,
                color='darkorange', alpha=.25, label='gaze S')

# cosmetics -----------------------------------------------------------------
ax.axvline(0, ls='-', c='k', lw=.8)                    # response moment
ax.axhline(.5, ls=':', c='k', lw=.8)                   # 50 %
ax.set(xlim=(-N_back, 0), ylim=(0,1),
       xlabel='Fixations before response\n(0 = last fixation before click)',
       ylabel='Proportion of fixations',
       title=f'Response-locked gaze trajectory\n(average ES fixations ≈ {mean_fix:.1f})')
ax.set_xticks(range(-N_back,1))
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, 'RespLocked_FixationTrajectory.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ──────────────────────────────────────────────────────────────────────────
#  EXTRA: plots *per learning block*  (cond = 0‒3)
# ──────────────────────────────────────────────────────────────────────────

# ---------------------------------------------------------------
OV_level = {1: "low",
            2: "medium",
            3: "high"}  

block_cols = ['pink', 'violet', 'indigo']

for OV_val, col in zip(range(1, 4), block_cols):
    label = OV_level[OV_val]
    blk_raw = raw[raw['OV_2'] == OV_val] 
    if blk_raw.empty:
        print(f'Block {OV_val} is empty - skipped.')
        continue

    # ===========================================================
    # ---------- RE-RUN the pipeline on *blk_raw*  ---------------
    # ===========================================================
    # 1)  choice-phase bins  ------------------------------------
    tbl = blk_raw.apply(trial_bins, axis=1).dropna()
    tbl = pd.concat(tbl.tolist(), ignore_index=True)

    rts = (blk_raw[['sub_id','TrialID','rtime']]
           .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
    rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)
    tbl = tbl.merge(rts, on=['sub','trial'])
    tbl = tbl[tbl.bin <= tbl.rt_bin]
    tbl['time_ms'] = tbl['bin'] * BIN_MS

    sub_traj = (tbl.groupby(['sub','time_ms']).hi_ratio
                    .mean().reset_index(name='hi_ratio'))
    grand = (sub_traj.groupby('time_ms')
                     .agg(mean=('hi_ratio','mean'),
                          sem =('hi_ratio', lambda x: x.std()/np.sqrt(len(x))))
                     .reset_index())
    grand['lo_mean'] = 1 - grand['mean']
    grand['lo_sem']  = grand['sem']

    rt_mean = blk_raw.rtime.mean()*1000
    rt_sem  = blk_raw.rtime.std(ddof=1)/np.sqrt(len(blk_raw))*1000



    # for E vs S ######################################################################
    ES_tbl = blk_raw.apply(trial_bins_ES, axis=1).dropna()
    ES_tbl = pd.concat(ES_tbl.tolist(), ignore_index=True)

    rts_ES = (raw[['sub_id','TrialID','rtime']]
              .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
    rts_ES['rtES_bin'] = (rts_ES.rt*1000 // BIN_MS).astype(int)

    ES_tbl = ES_tbl.merge(rts_ES, on=['sub','trial'])
    ES_tbl = ES_tbl[ES_tbl.bin <= ES_tbl.rtES_bin]           # keep up to the click

    # time in ms relative to onset
    ES_tbl['timeES_ms'] = ES_tbl['bin'] * BIN_MS

    # subject-average trajectory 
    sub_trajES = (ES_tbl.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index(name='ES_ratio'))

    grandES = (sub_trajES.groupby('timeES_ms')
                     .agg(mean=('ES_ratio','mean'),
                          sem =('ES_ratio', lambda x: x.std()/np.sqrt(len(x))))
                     .reset_index())
    grandES['S_mean'] = 1 - grandES['mean']
    grandES['S_sem']  = grandES['sem']

    rtES_mean = blk_raw.rtime.mean()*1000
    rtES_sem  = blk_raw.rtime.std(ddof=1)/np.sqrt(len(blk_raw))*1000

    # ===========================================================
    # ----------------------  PLOTS  -----------------------------
    # ===========================================================

    # ----- figure 1: choice -----------------------------------
    fig_c, ax_c = plt.subplots(figsize=(4.3,3))
    ax_c.axvspan(rt_mean-rt_sem, rt_mean+rt_sem,
                 color='lightgrey', alpha=.5, zorder=0)
    ax_c.axvline(rt_mean, color='grey', lw=1.2)

    ax_c.plot(grand.time_ms, grand['mean'],  color=col, lw=2)
    ax_c.fill_between(grand.time_ms,
                      grand['mean']-grand['sem'],
                      grand['mean']+grand['sem'],
                      color=col, alpha=.25)
    ax_c.plot(grand.time_ms, grand['lo_mean'], color=col,
              ls='--', lw=2)
    ax_c.fill_between(grand.time_ms,
                      grand['lo_mean']-grand['lo_sem'],
                      grand['lo_mean']+grand['lo_sem'],
                      color=col, alpha=.25)

    ax_c.axhline(.5, ls='--', color='k')
    ax_c.set(xlim=(0, rt_mean+100), ylim=(0,1),
             xlabel='Time (ms)', ylabel='Fixation ratio',
             title=f'OV-level {label}: Choice')

    # ----- figure 2: ES  E VS S  -------------------
    fig_f, ax_f = plt.subplots(figsize=(4.3,3))
    ax_f.axvspan(rtES_mean-rtES_sem, rtES_mean+rtES_sem,
                 color='lightgrey', alpha=.5, zorder=0)
    ax_f.axvline(rtES_mean, color='grey', lw=1.2)    
    ax_f.plot(grandES.timeES_ms, grandES['mean'],  color=col, lw=2)
    
    ax_f.fill_between(grandES.timeES_ms,
                      grandES['mean']-grandES['sem'],
                      grandES['mean']+grandES['sem'],
                      color=col, alpha=.25)
    ax_f.plot(grandES.timeES_ms, grandES['S_mean'], color=col,
              ls='--', lw=2)
    ax_f.fill_between(grandES.timeES_ms,
                      grandES['S_mean']-grandES['S_sem'],
                      grandES['S_mean']+grandES['S_sem'],
                      color=col, alpha=.25)

    ax_f.axhline(.5, ls='--', color='k')
    ax_f.set(xlim=(0, rtES_mean+100), ylim=(0,1),
             xlabel='Time (ms)', ylabel='Fixation ratio',
             title=f'OV-level {label}: Choice')
    
    # tidy legends (same style for every plot)
    # ── tidy legends (same style for every plot) ─────────────────────────────
    solid_hi  = Line2D([], [], color=col, lw=2,      label='higher EV')
    dashed_lo = Line2D([], [], color=col, ls='--', lw=2, label='lower EV')

    solid_E = Line2D([], [], color=col, lw=2,      label='E')
    dash_S   = Line2D([], [], color=col, ls='--', lw=2, label='S')

    ax_c.legend(handles=[solid_hi, dashed_lo], frameon=False, fontsize=7)
    ax_f.legend(handles=[solid_E, dash_S], frameon=False, fontsize=7)

    fig_c.savefig(os.path.join(fig_dir, f'OV{label}_ES_EV.png'),dpi=300, bbox_inches='tight')
    fig_f.savefig(os.path.join(fig_dir, f'OV{label}_ES_ES.png'),dpi=300, bbox_inches='tight')




# 1) -------------------------------------------------------------
all_curves = defaultdict(list)          # keys → 'choice' | 'feed' | 'outcf'


for OV_val, col in zip(range(1, 4), block_cols):      # 1, 2, 3  (= low-medium-high)
    label   = OV_level[OV_val]
    blk_raw = raw[raw['OV_2'] == OV_val]              # ← FIX 1
    if blk_raw.empty:
        continue
    # ---------- a) choice: hi vs lo EV -------------------------
    tbl   = blk_raw.apply(trial_bins, axis=1).dropna()
    tbl   = pd.concat(tbl.tolist(), ignore_index=True)
    rts   = blk_raw[['sub_id','TrialID','rtime']].rename(
              columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'})
    rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)
    tbl   = tbl.merge(rts, on=['sub','trial'])
    tbl   = tbl[tbl.bin <= tbl.rt_bin]
    tbl['time_ms'] = tbl['bin']*BIN_MS
    sub   = tbl.groupby(['sub','time_ms']).hi_ratio.mean().reset_index()
    g     = sub.groupby('time_ms').hi_ratio.agg(['mean','sem']).reset_index()
    g['lo_mean'] = 1-g['mean'];  g['lo_sem'] = g['sem']
    all_curves['choice'].append((label, col, g, 'higher EV', 'lower EV'))

    # ---------- b) feedback: hi vs lo EV -----------------------
    
    tbl_ES   = blk_raw.apply(trial_bins_ES, axis=1).dropna()
    tbl_ES   = pd.concat(tbl_ES.tolist(), ignore_index=True)

    rts_ES   = blk_raw[['sub_id','TrialID','rtime']].rename(
              columns={'sub_id':'sub',
                       'TrialID':'trial',
                       'rtime':'rt'})
    rts_ES['rtES_bin'] = (rts_ES.rt * 1000 // BIN_MS).astype(int)

    tbl_ES = tbl_ES.merge(rts_ES, on=['sub','trial'])
    tbl_ES = tbl_ES[tbl_ES.bin <= tbl_ES.rtES_bin]   # ← now this works
    tbl_ES['timeES_ms'] = tbl_ES['bin'] * BIN_MS
    sub_ES   = tbl_ES.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index()
    g_ES     = sub_ES.groupby('timeES_ms').ES_ratio.agg(['mean','sem']).reset_index()
    g_ES['S_mean'] = 1-g_ES['mean'];  g_ES['S_sem'] = g_ES['sem']
    g_ES.rename(columns={'timeES_ms': 'time_ms'}, inplace=True)   # ← unify x-axis name
    all_curves['ES'].append((label, col, g_ES, 'E', 'S'))   


# ---------------------------------------------------------------
# 2)  helper to draw one overlay figure
# ---------------------------------------------------------------
def overlay_plot(key, title, ylab, fname):
    """
    key        : 'choice', 'ES'
    title      : figure title
    ylab       : y-axis label
    fname      : file name (saved in fig_dir)
    """
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    
    ax.axvspan(rt_glob_mean-rt_glob_sem,
               rt_glob_mean+rt_glob_sem,
               color='lightgrey', alpha=.35, zorder=10)
    ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=11)

    for label, col, df, solid_lbl, dash_lbl in all_curves[key]:
        # solid line (+ SEM band)
        ax.plot(df.time_ms, df['mean'],
                color=col, lw=2, label=f'{label} – {solid_lbl}')
        ax.fill_between(df.time_ms,
                        df['mean']-df['sem'],
                        df['mean']+df['sem'],
                        color=col, alpha=.20)

        # dashed line (+ SEM band)
        second   = 'lo_mean' if key == 'choice' else 'S_mean'
        second_s = 'lo_sem'  if key == 'choice' else 'S_sem'
        ax.plot(df.time_ms, df[second],
                color=col, ls='--', lw=2, label=f'{label} – {dash_lbl}')
        ax.fill_between(df.time_ms,
                        df[second]-df[second_s],
                        df[second]+df[second_s],
                        color=col, alpha=.20)

    ax.axhline(.5, ls='--', color='k')
    ax.axvline(0 , ls='--', color='k')
    ax.set(xlim=(0, 1300), ylim=(0,1),
           xlabel='Time (ms)', ylabel=ylab, title=title)

    # put legend outside (to the right)
    ax.legend(frameon=False, fontsize=7, ncol=1,
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, fname),
                dpi=300, bbox_inches='tight')
    plt.show()

# ---------------------------------------------------------------
# 3)  make the three overlay figures
# ---------------------------------------------------------------
overlay_plot('choice',
             'Choice - gaze on higher vs lower EV (OV levels)',
             'Fixation ratio',
             'ALL_OV_levels_ES_EV.png')

overlay_plot('ES',
             'Choice - gaze on E vs S (OV levels)',
             'Fixation ratio',
             'ALL_OV_levels_ES_ES.png')

