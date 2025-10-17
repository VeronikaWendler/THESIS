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



BIN_MS    = 25                       # bin width in ms
ROI   = {1: 'left', 2: 'right'}      # ROI

raw = df_raw.query("phase == 'LE'").copy()   # keep only Learning‐phase trials

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


# helper for FEEDBACK phase -----------------------------------------------------------------------
def trial_bins_feed(row):
    try:
        fixes = ast.literal_eval(row.Fixations_Feed) or []
    except Exception:
        return None
    if not fixes:
        return None

    hi_side = 'left' if row.p1 > row.p2 else 'right'
    t0 = fixes[0]['start_time']
    bins = {}

    for f in fixes:
        side  = ROI.get(f['roi'], 'none')
        is_hi = side == hi_side
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]
            bins[b][0 if is_hi else 1] += 1

    out = pd.DataFrame.from_dict(bins, orient='index', columns=['hi','lo'])
    out['sub'] = row.sub_id
    out['trial'] = row.TrialID
    out['bin'] = out.index
    out['hi_ratio'] = out.hi / (out.hi+out.lo)
    return out[['sub','trial','bin','hi_ratio']]


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


# ─── build table for FEEDBACK phase ─────────────────────────────────────
feed_tbl = raw.apply(trial_bins_feed, axis=1).dropna()
feed_tbl = pd.concat(feed_tbl.tolist(), ignore_index=True)
feed_tbl['time_ms'] = feed_tbl['bin'] * BIN_MS

feed_sub = (feed_tbl.groupby(['sub','time_ms'])
                      .hi_ratio.mean()
                      .reset_index(name='hi_ratio'))

feed_grand = (feed_sub.groupby('time_ms')
                        .agg(mean=('hi_ratio','mean'),
                             sem =('hi_ratio',
                                   lambda x: x.std()/np.sqrt(len(x))))
                        .reset_index())
feed_grand['lo_mean'] = 1 - feed_grand['mean']
feed_grand['lo_sem']  = feed_grand['sem']

# decide how far to plot the feedback axis (e.g. 700 ms or the 95 % percentile)
feed_end = min(700, feed_grand.time_ms.quantile(.95))

# -----------------------------------------------------------------------
# helper: small, tidy legend kwargs
legend_kw = dict(loc='upper right', frameon=False,
                 fontsize=8,  title_fontsize=9,
                 borderpad=0.2, labelspacing=0.25, handlelength=2)

line_handles = [
    Line2D([], [], color='royalblue',  lw=2, label='higher'),
    Line2D([], [], color='royalblue', ls='--', lw=2, label='lower')
]

# ────────────────────────────────────────────────────────────────────────
# ❶ CHOICE  (stimulus-locked) – its own figure
# ────────────────────────────────────────────────────────────────────────
fig_choice, ax_c = plt.subplots(figsize=(4.5, 3))

# RT band  (mean ± SEM)
ax_c.axvspan(rt_mean-rt_sem, rt_mean+rt_sem,
             color='lightgrey', alpha=.5, zorder=0)
# RT mean line on top of the band
ax_c.axvline(rt_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_c.plot(grand.time_ms, grand['mean'],
          color='royalblue', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['mean']-grand['sem'],
                  grand['mean']+grand['sem'],
                  color='royalblue', alpha=.25, label='_nolegend_')

# lower-EV
ax_c.plot(grand.time_ms, grand['lo_mean'],
          color='royalblue', ls='--', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['lo_mean']-grand['lo_sem'],
                  grand['lo_mean']+grand['lo_sem'],
                  color='royalblue', alpha=.20, label='_nolegend_')

ax_c.axhline(.5, ls='--', color='k')
ax_c.axvline(0,  ls='--', color='k')
ax_c.set(xlim=(0, rt_mean+100), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Choice (stimulus-locked)')

ax_c.legend(handles=line_handles, title='Gaze on', **legend_kw)
fig_choice.tight_layout()
plot_path= os.path.join(fig_dir, "gaze_choice.png")
fig_choice.savefig(plot_path, dpi=300, bbox_inches='tight')
fig_choice.show()

# ────────────────────────────────────────────────────────────────────────
# ❷ FEEDBACK  (outcome-locked) – separate figure
# ────────────────────────────────────────────────────────────────────────
fg = feed_grand.query('time_ms <= @feed_end')

fig_feed, ax_f = plt.subplots(figsize=(4.5, 3))

ax_f.plot(fg.time_ms, fg['mean'],
          color='royalblue', lw=2)
ax_f.fill_between(fg.time_ms,
                  fg['mean']-fg['sem'],
                  fg['mean']+fg['sem'],
                  color='royalblue', alpha=.25, label='_nolegend_')

ax_f.plot(fg.time_ms, fg['lo_mean'],
          color='royalblue', ls='--', lw=2)
ax_f.fill_between(fg.time_ms,
                  fg['lo_mean']-fg['lo_sem'],
                  fg['lo_mean']+fg['lo_sem'],
                  color='royalblue', alpha=.20, label='_nolegend_')

ax_f.axhline(.5, ls='--', color='k')
ax_f.axvline(0,  ls='--', color='k')
ax_f.set(xlim=(0, feed_end), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Feedback (outcome-locked)')

ax_f.legend(handles=line_handles, title='Gaze on', **legend_kw)
fig_feed.tight_layout()
plot_path2= os.path.join(fig_dir, "gaze_feed.png")
fig_feed.savefig(plot_path2, dpi=300, bbox_inches='tight')
fig_feed.show()   


def trial_bins_outcf(row):
    """
    During feedback, count 50-ms bins that fall on the *obtained* outcome
    (‘out’) vs. the *counter-factual* outcome (‘cfout’).
    Which screen side shows which number is given by `chose_right`.
    """
    try:
        fixes = ast.literal_eval(row.Fixations_Feed) or []
    except Exception:
        return None
    if not fixes:                      # no gaze in feedback
        return None

    out_side = 'right' if row.chose_right == 1 else 'left'
    cf_side  = 'left'  if row.chose_right == 1 else 'right'

    t0, bins = fixes[0]['start_time'], {}        # first sample = 0 ms

    for f in fixes:
        side = ROI.get(f['roi'], 'none')
        idx  = 0 if side == out_side else 1       # 0 = out, 1 = cf
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]
            bins[b][idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['out_cnt', 'cf_cnt'])
    out['sub']        = row.sub_id
    out['trial']      = row.TrialID
    out['bin']        = out.index
    out['out_ratio']  = out.out_cnt / (out.out_cnt + out.cf_cnt)
    return out[['sub', 'trial', 'bin', 'out_ratio']]


# ---------- build table & average ----------------------------------------
outcf_tbl = raw.apply(trial_bins_outcf, axis=1).dropna()
outcf_tbl = pd.concat(outcf_tbl.tolist(), ignore_index=True)
outcf_tbl['time_ms'] = outcf_tbl['bin'] * BIN_MS

outcf_sub = (outcf_tbl.groupby(['sub', 'time_ms'])
                       .out_ratio.mean()
                       .reset_index(name='out_ratio'))

outcf_grand = (outcf_sub.groupby('time_ms')
                          .agg(mean=('out_ratio','mean'),
                               sem =('out_ratio',
                                     lambda x: x.std()/np.sqrt(len(x))))
                          .reset_index())
outcf_grand['cf_mean'] = 1 - outcf_grand['mean']
outcf_grand['cf_sem']  = outcf_grand['sem']

# time window (same rule as earlier)
outcf_end = min(700, outcf_grand.time_ms.quantile(.95))


# ---------- plot ----------------------------------------------------------
fig_outcf, ax_o = plt.subplots(figsize=(4.5, 3))

ax_o.plot(outcf_grand.time_ms, outcf_grand['mean'],
          color='royalblue', lw=2, label='obtained (out)')
ax_o.fill_between(outcf_grand.time_ms,
                  outcf_grand['mean']-outcf_grand['sem'],
                  outcf_grand['mean']+outcf_grand['sem'],
                  color='royalblue', alpha=.25)

ax_o.plot(outcf_grand.time_ms, outcf_grand['cf_mean'],
          color='royalblue', ls='--', lw=2, label='counter-factual (cfout)')
ax_o.fill_between(outcf_grand.time_ms,
                  outcf_grand['cf_mean']-outcf_grand['cf_sem'],
                  outcf_grand['cf_mean']+outcf_grand['cf_sem'],
                  color='royalblue', alpha=.25)

ax_o.axhline(.5, ls='--', color='k')
ax_o.axvline(0,  ls='--', color='k')
ax_o.set(xlim=(0, outcf_end), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Feedback – obtained vs counter-factual')

ax_o.legend(loc='upper right', frameon=False, fontsize=8, title_fontsize=9)
fig_outcf.tight_layout()
plot_path3= os.path.join(fig_dir, "gaze_feed_out_cfout.png")
fig_outcf.savefig(plot_path3, dpi=300, bbox_inches='tight')
fig_outcf.show()    # shows the third figure


#-----------------------------------------------------------------------------------------------------------------

# plots per LE block

# ──────────────────────────────────────────────────────────────────────────
#  EXTRA: plots *per learning block*  (cond = 0‒3)
# ──────────────────────────────────────────────────────────────────────────

# ---------------------------------------------------------------
difficulty = {0: "90_10",
              1: "80_20",
              2: "70_30",
              3: "60_40"}  

block_cols = ['darkred', 'darkorange', 'seagreen', 'mediumpurple']

for blk_idx, (cond_val, col) in enumerate(zip(range(4), block_cols), start=1):
    blk_raw = raw[raw['cond'] == cond_val].copy()
    if blk_raw.empty:
        print(f'Block {cond_val} is empty – skipped.')
        continue

    label = difficulty[cond_val]
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

    # 2)  feedback (higher vs lower EV) -------------------------
    feed_tbl = blk_raw.apply(trial_bins_feed, axis=1).dropna()
    feed_tbl = pd.concat(feed_tbl.tolist(), ignore_index=True)
    feed_tbl['time_ms'] = feed_tbl['bin'] * BIN_MS
    feed_sub = (feed_tbl.groupby(['sub','time_ms']).hi_ratio
                         .mean().reset_index(name='hi_ratio'))
    feed_grand = (feed_sub.groupby('time_ms')
                            .agg(mean=('hi_ratio','mean'),
                                 sem =('hi_ratio', lambda x: x.std()/np.sqrt(len(x))))
                            .reset_index())
    feed_grand['lo_mean'] = 1 - feed_grand['mean']
    feed_grand['lo_sem']  = feed_grand['sem']
    feed_end = min(700, feed_grand.time_ms.quantile(.95))

    # 3)  feedback (obtained vs counter-factual) ----------------
    outcf_tbl = blk_raw.apply(trial_bins_outcf, axis=1).dropna()
    outcf_tbl = pd.concat(outcf_tbl.tolist(), ignore_index=True)
    outcf_tbl['time_ms'] = outcf_tbl['bin'] * BIN_MS
    outcf_sub = (outcf_tbl.groupby(['sub','time_ms']).out_ratio
                           .mean().reset_index(name='out_ratio'))
    outcf_grand = (outcf_sub.groupby('time_ms')
                              .agg(mean=('out_ratio','mean'),
                                   sem =('out_ratio', lambda x: x.std()/np.sqrt(len(x))))
                              .reset_index())
    outcf_grand['cf_mean'] = 1 - outcf_grand['mean']
    outcf_grand['cf_sem']  = outcf_grand['sem']
    outcf_end = min(700, outcf_grand.time_ms.quantile(.95))

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
             title=f'Block {label}: Choice')

    # ----- figure 2: feedback (hi vs lo EV) -------------------
    fig_f, ax_f = plt.subplots(figsize=(4.3,3))
    fg = feed_grand.query('time_ms <= @feed_end')
    ax_f.plot(fg.time_ms, fg['mean'],  color=col, lw=2)
    ax_f.fill_between(fg.time_ms,
                      fg['mean']-fg['sem'],
                      fg['mean']+fg['sem'],
                      color=col, alpha=.25)
    ax_f.plot(fg.time_ms, fg['lo_mean'], color=col,
              ls='--', lw=2)
    ax_f.fill_between(fg.time_ms,
                      fg['lo_mean']-fg['lo_sem'],
                      fg['lo_mean']+fg['lo_sem'],
                      color=col, alpha=.25)
    ax_f.axhline(.5, ls='--', color='k')
    ax_f.set(xlim=(0, feed_end), ylim=(0,1),
             xlabel='Time (ms)', ylabel='Fixation ratio',
             title=f'Block {label}: Feedback (hi/lo EV)')

    # ----- figure 3: feedback (out vs cf) ---------------------
    fig_o, ax_o = plt.subplots(figsize=(4.3,3))
    ax_o.plot(outcf_grand.time_ms, outcf_grand['mean'],
              color=col, lw=2)
    ax_o.fill_between(outcf_grand.time_ms,
                      outcf_grand['mean']-outcf_grand['sem'],
                      outcf_grand['mean']+outcf_grand['sem'],
                      color=col, alpha=.25)
    ax_o.plot(outcf_grand.time_ms, outcf_grand['cf_mean'],
              color=col, ls='--', lw=2)
    ax_o.fill_between(outcf_grand.time_ms,
                      outcf_grand['cf_mean']-outcf_grand['cf_sem'],
                      outcf_grand['cf_mean']+outcf_grand['cf_sem'],
                      color=col, alpha=.25)
    ax_o.axhline(.5, ls='--', color='k')
    ax_o.set(xlim=(0, outcf_end), ylim=(0,1),
             xlabel='Time (ms)', ylabel='Fixation ratio',
             title=f'Block {label}: Feedback (out/cf)')

    # tidy legends (same style for every plot)
    # ── tidy legends (same style for every plot) ─────────────────────────────
    solid_hi  = Line2D([], [], color=col, lw=2,      label='higher EV')
    dashed_lo = Line2D([], [], color=col, ls='--', lw=2, label='lower EV')

    solid_out = Line2D([], [], color=col, lw=2,      label='obtained (out)')
    dash_cf   = Line2D([], [], color=col, ls='--', lw=2, label='counter-factual')

    ax_c.legend(handles=[solid_hi, dashed_lo], frameon=False, fontsize=7)
    ax_f.legend(handles=[solid_hi, dashed_lo], frameon=False, fontsize=7)
    ax_o.legend(handles=[solid_out, dash_cf ], frameon=False, fontsize=7)

    fig_c.savefig(os.path.join(fig_dir, f'block{label}_choice.png'),dpi=300, bbox_inches='tight')
    fig_f.savefig(os.path.join(fig_dir, f'block{label}_feedbackEV.png'),dpi=300, bbox_inches='tight')
    fig_o.savefig(os.path.join(fig_dir, f'block{label}_feedbackOutCF.png'),dpi=300, bbox_inches='tight')




# 1) -------------------------------------------------------------
all_curves = defaultdict(list)          # keys → 'choice' | 'feed' | 'outcf'

for cond_val, col in zip(range(4), block_cols):
    label = difficulty[cond_val]        # 90_10 … 60_40
    blk_raw = raw[raw['cond'] == cond_val]
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
    ftbl  = blk_raw.apply(trial_bins_feed, axis=1).dropna()
    ftbl  = pd.concat(ftbl.tolist(), ignore_index=True)
    ftbl['time_ms'] = ftbl['bin']*BIN_MS
    fsub  = ftbl.groupby(['sub','time_ms']).hi_ratio.mean().reset_index()
    fg    = fsub.groupby('time_ms').hi_ratio.agg(['mean','sem']).reset_index()
    fg['lo_mean'] = 1-fg['mean'];  fg['lo_sem']=fg['sem']
    all_curves['feed'].append((label, col, fg, 'higher EV', 'lower EV'))

    # ---------- c) feedback: obtained vs counter-factual -------
    otbl  = blk_raw.apply(trial_bins_outcf, axis=1).dropna()
    otbl  = pd.concat(otbl.tolist(), ignore_index=True)
    otbl['time_ms'] = otbl['bin']*BIN_MS
    osub  = otbl.groupby(['sub','time_ms']).out_ratio.mean().reset_index()
    og    = osub.groupby('time_ms').out_ratio.agg(['mean','sem']).reset_index()
    og['cf_mean'] = 1-og['mean']; og['cf_sem']=og['sem']
    all_curves['outcf'].append((label, col, og,
                                'obtained (out)', 'counter-factual'))

# ---------------------------------------------------------------
# 2)  helper to draw one overlay figure
# ---------------------------------------------------------------
def overlay_plot(key, title, ylab, fname):
    """
    key        : 'choice', 'feed', or 'outcf'
    title      : figure title
    ylab       : y-axis label
    fname      : file name (saved in fig_dir)
    """
    fig, ax = plt.subplots(figsize=(6.2, 3.8))

    for label, col, df, solid_lbl, dash_lbl in all_curves[key]:
        # solid line (+ SEM band)
        ax.plot(df.time_ms, df['mean'],
                color=col, lw=2, label=f'{label} – {solid_lbl}')
        ax.fill_between(df.time_ms,
                        df['mean']-df['sem'],
                        df['mean']+df['sem'],
                        color=col, alpha=.20)

        # dashed line (+ SEM band)
        second   = 'lo_mean' if key in ('choice','feed') else 'cf_mean'
        second_s = 'lo_sem'  if key in ('choice','feed') else 'cf_sem'
        ax.plot(df.time_ms, df[second],
                color=col, ls='--', lw=2, label=f'{label} – {dash_lbl}')
        ax.fill_between(df.time_ms,
                        df[second]-df[second_s],
                        df[second]+df[second_s],
                        color=col, alpha=.20)

    ax.axhline(.5, ls='--', color='k')
    ax.axvline(0 , ls='--', color='k')
    ax.set(xlim=(0, 700), ylim=(0,1),
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
             'Choice – gaze on higher vs lower EV (all blocks)',
             'Fixation ratio',
             'ALL_blocks_choice.png')

overlay_plot('feed',
             'Feedback – gaze on higher vs lower EV (all blocks)',
             'Fixation ratio',
             'ALL_blocks_feedbackEV.png')

overlay_plot('outcf',
             'Feedback – gaze on obtained vs counter-factual (all blocks)',
             'Fixation ratio',
             'ALL_blocks_feedbackOutCF.png')