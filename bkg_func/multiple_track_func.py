'''set of functions to run the analysis on multiple tracks (XLM file) at once '''

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from bkg_func.conf_ratio_func import ComputeConfinementRatioScore, FindTransitionPoints, DenoiseTransitionPoints
from bkg_func.conf_ratio_func import ComputeSubSegmentStats


def AnalyzeAllTracks(all_tracks, frame_rate, tracks = -1, w = 10, p_thres = 1500, t_thres = 0.5):
    '''incorporates functions from single track analysis into a for loop to read all tracks
    tracks: number of tracks to analyze. default = -1 (all tracks)
    w: window size;
    p_thres: value above which a subsegment is considered confined (dimensionless)
    t_thres: minimal time (sec) a confined event needs to last to be considered '''
    
    #print('rolling window: {0:.4} sec'.format(window_size * frame_rate))
    
    table = all_tracks.copy()
    
    # truncate number of tracks to analyze; default is -1 to analyze eveything
    
    if tracks == -1: 
        table = table

    else:
        trunc = table[table['TRACK_ID'] == (tracks-1)].index[-1]
        table = table.iloc[:trunc]
    
    # run all function for each track
    
    all_tracks_stats = []
    
    for i, track in table.groupby('TRACK_ID'):
        # check if trajectory is just noise (i.e empty tracks)
            
        noise = (track[~track.duplicated(['POSITION_X', 'POSITION_Y'])].shape[0] < 10)
        
        if noise == False:
        
            print('track {}/{} - {} steps'.format(i+1, len(table.TRACK_ID.unique()), track.shape[0]), end = '\r')

            try: 
                track_score = ComputeConfinementRatioScore(track, conf_ratio_thres = p_thres, rol_window = w, show_progress = False)
                track_score_denoised = DenoiseTransitionPoints(track_score , FindTransitionPoints(track_score), frame_rate, t_thres)
                track_score_denoised_stats = ComputeSubSegmentStats(track_score_denoised, frame_rate)

                all_tracks_stats.append(track_score_denoised_stats)
            
            except:
                print('track '+ str(i) + ' is weird: skipping')
                continue
    
    # add two extra columns to the final table: track_id and subtrajectory_id
    for i, track in enumerate(all_tracks_stats):
        track.insert(0, 'subtraj_id', np.arange(track.shape[0]))
        track.insert(0, 'track_id', i)
    
    return all_tracks_stats

def PlotOutputStuff(all_tracks_stats):
    
    import warnings
    warnings.filterwarnings("ignore") 

    mode_col = 'mode'
    lifetimes_col = 'lifetime_s'
    area_col = 'area_nm2'
    colors = ['crimson', 'steelblue', 'seagreen']
    
    # extract lifetimes information
    lifetime_free_events = all_tracks_stats[all_tracks_stats[mode_col] == 'free_diff'][lifetimes_col].tolist()
    lifetime_confined_events = all_tracks_stats[all_tracks_stats[mode_col] == 'confined'][lifetimes_col].tolist()
    lifetimes_df = pd.DataFrame([lifetime_confined_events, lifetime_free_events], index = ['confined','free']).T
    
    # extract confinement area information
    area_free_events = all_tracks_stats[all_tracks_stats[mode_col] == 'free_diff'][area_col].tolist()
    area_confined_events = all_tracks_stats[all_tracks_stats[mode_col] == 'confined'][area_col].tolist()
    areas_df = pd.DataFrame([area_confined_events, area_free_events], index = ['confined','free']).T
    
    # extract number of confined events per track - to plot distribution
    conf_events_per_tracks = []
    for i, track in all_tracks_stats.groupby('track_id'):
        conf_events_per_tracks.append(track[track[mode_col] == 'confined'].shape[0])
    
    # plot stuff
    fig, ax = plt.subplots(1, 3, figsize = (8, 2.5), dpi = 150)
    fig.subplots_adjust(wspace = 0.4)
    
    # lifetimes 
    sns.barplot(data = lifetimes_df, ax = ax[0], 
                 palette = colors, edgecolor = 'black', 
                 linewidth = 1, errwidth=1, capsize=0.1);
    
    # confinement area distribution
    sns.histplot(areas_df.confined.dropna().tolist(), ax = ax[1], shrink = 0.8, 
                 edgecolor = 'black', linewidth = 0.8); 
    
    #sns.kdeplot(data = conf_events_per_tracks);
    sns.countplot(conf_events_per_tracks, ax = ax[2], #shrink = 0.8, 
                  #bins = np.arange(np.max(conf_events_per_tracks) + 1),
                  color = colors[2], edgecolor = 'black', linewidth = 0.8);
    
    ax[0].set_title('time under confinement \n vs. free diffusion', fontsize = 8)
    ax[0].set_ylabel('subdiffusion lifetime (s)', fontsize = 9); #ax.set_xlabel('type of motion', fontsize = 11);
    ax[0].tick_params(direction = 'out', bottom = False, labelsize=9)
    
    ax[1].set_title('confinement area \n for all sub-trajectories', fontsize = 8)
    ax[1].set_ylabel('confinement area (nm2)', fontsize = 9); 
    ax[1].tick_params(direction = 'out', bottom = False, labelsize=9)
    
    ax[2].set_title('number of confined \n events per track', fontsize = 8)
    ax[2].set_ylabel('number of tracks', fontsize = 9); 
    ax[2].set_xlabel('number of confined events', fontsize = 8); 
    ax[2].tick_params(direction = 'out', bottom = False, labelsize=9)
    #ax[2].set_xlim([-1, np.max(conf_events_per_tracks) + 1]);
    
    
def TrimEnds(all_tracks_stats, start = True, end = True):
    '''discards confined motions identified at the end or start of each trajectory
    if start == True: discards confinement motions identified at the begninning of each trajectory
    if end == True discards confinement motions identified at the end of each trajectory'''
    
    all_tracks_stats_filtered = []
    
    for i, traj in enumerate(all_tracks_stats):
        
        traj_filtered = traj.copy()

        if start == True:
            if traj.iloc[0]['mode'] == 'confined':
                traj_filtered = traj_filtered.drop(0)
            else: 
                traj_filtered = traj_filtered
            
        if end == True:
            if traj.iloc[-1]['mode'] == 'confined':
                traj_filtered = traj_filtered.drop(traj.index[-1])

            else: 
                traj_filtered = traj_filtered

        all_tracks_stats_filtered.append(traj_filtered.reset_index(drop = True))
        
    return all_tracks_stats_filtered