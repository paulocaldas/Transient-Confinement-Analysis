'''set of functions to run the analysis on multiple tracks (XLM file) at once '''

import pandas as pd
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
    
    return all_tracks_stats

def ExtractConfinementLifetimes(all_tracks_stats):
    '''takes the output from the function AnalyzeAllTracks and returns
    a dataframe with the lifetimes of confinement vs. free diffusion'''
    
    tracks_lifetimes_confined = []
    tracks_lifetimes_free = []

    for track in all_tracks_stats:

        confined_lifetimes = track[track['mode'] == 'confined']['lifetime (s)'].values
        free_lifetimes = track[track['mode'] == 'free_diff']['lifetime (s)'].values

        tracks_lifetimes_confined.append(confined_lifetimes)
        tracks_lifetimes_free.append(free_lifetimes)

    # unpack lists
    tracks_lifetimes_confined = [i for j in tracks_lifetimes_confined for i in j]
    tracks_lifetimes_free = [i for j in tracks_lifetimes_free for i in j]
    
    # final dataframe
    lifetimes_df = pd.DataFrame([tracks_lifetimes_confined, tracks_lifetimes_free],
                                 index = ['confined','free']).T
    
    return lifetimes_df

def ExtractConfinementCounts(all_tracks_stats):
    '''takes the output from the function AnalyzeAllTracks and returns
    a dataframe with the number and % of confiinement events'''
    
    tracks_w_confinement = 0
    tracks_wo_confinement = 0
    
    for track in all_tracks_stats:

        if 'confined' in track['mode'].values:
            tracks_w_confinement += 1
        
        if 'confined' not in track['mode'].values: 
            tracks_wo_confinement += 1
    
    print('tracks w confinement: {}/{} ~ {:.3}% '.format(tracks_w_confinement, len(all_tracks_stats), 
                                              (tracks_w_confinement/len(all_tracks_stats)*100)))
    return tracks_w_confinement

def ExtractConfinementCountsII(all_tracks_stats):
    '''returns a distrbuiton of confined motion '''
    
    conf_dist = {'one':0, 'two':0, '>three':0}
    
    tracks_wo_confinement = 0
    
    for track in all_tracks_stats:

        if 'confined' in track['mode'].values:
            
            if track['mode'].value_counts()['confined'] == 1:
                conf_dist['one'] += 1
            
            if track['mode'].value_counts()['confined'] == 2:
                conf_dist['two'] += 1
                
            if track['mode'].value_counts()['confined'] >= 3:
                conf_dist['>three'] += 1
        
        if 'confined' not in track['mode'].values: 
            tracks_wo_confinement += 1
        
    return conf_dist

def ExtractConfinementArea(all_tracks_stats):
    '''takes the output from the function AnalyzeAllTracks and returns
    a dataframe with the cage area of each subsegment analyzed'''
    
    tracks_area_confined = []
    tracks_area_free = []

    for track in all_tracks_stats:

        confined_areas = track[track['mode'] == 'confined']['cage_area (nm2)'].values
        free_areas = track[track['mode'] == 'free_diff']['cage_area (nm2)'].values

        tracks_area_confined.append(confined_areas)
        tracks_area_free.append(free_areas)

    # unpack lists
    tracks_area_confined = [i for j in tracks_area_confined for i in j]
    tracks_area_free = [i for j in tracks_area_free for i in j]
    
    # make a dataframe
    
    tracks_area = pd.DataFrame([tracks_area_confined, tracks_area_free], index = ['confined','free']).T
    
    return tracks_area

def PlotLifetimes(lifetimes, colors = ['steelblue', 'crimson']):
    '''takes the output from the ExtractConfinementLifetimes to 
    plot confinement vs. free diffusion lifetimes'''
    
    fig, ax1 = plt.subplots(figsize = (2,2.5), dpi = 150)
    
    sns.barplot(data = lifetimes,
                ax = ax1, palette = colors, edgecolor = 'black', 
                linewidth = 1, errwidth=1, capsize=0.1)
    
    # fancy stuff
    ax1.set_ylabel('subdiffusion time (s)', fontsize = 9); #ax.set_xlabel('type of motion', fontsize = 11);
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.tick_params(direction = 'out', top=False, right = False, bottom = False, labelsize=8)
    
    #ax.set_ylim([0, DefineYmax(ax)])

def PlotConfinementCounts(conf_counts, colors = ['peachpuff', 'salmon', 'orange']):
    
    data_to_plot = pd.DataFrame.from_dict(conf_counts, orient = 'index').T
    
    fig, ax = plt.subplots(figsize = (2,1.5), dpi = 150)
    
    sns.barplot(data = data_to_plot, ax = ax, palette = colors, edgecolor = 'black', 
                    linewidth = 1, errwidth=1, capsize=0.1,)
    
    # fancy stuff
    ax.set_ylabel('number of tracks', fontsize = 9); ax.set_xlabel('number of confined events', fontsize = 8);
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.tick_params(direction = 'out', top=False, right = False, bottom = False, labelsize=8)
    
   # ax.set_ylim([0, DefineYmax(ax)])
    
    return conf_counts

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