'''single track - main function '''

#import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
#from matplotlib.ticker import FormatStrFormatter

from bkg_func.utils import read_trackmate_xml_tracks, FilterTracks, SingleTrack
from bkg_func.conf_ratio_func import ComputeConfinementRatioScore, FindTransitionPoints, DenoiseTransitionPoints
from bkg_func.conf_ratio_func import PlotConfinedRegions, ComputeSubSegmentStats, ShowStats
from bkg_func.multiple_track_func import AnalyzeAllTracks, ExtractConfinementLifetimes, PlotLifetimes, TrimEnds


def ReadTracks(file, minlen):
    '''reads file into a dataframe, extracts frame rate and number of tracks'''
    table, rate, n_tracks = read_trackmate_xml_tracks(file)
    
    # apply filter for minimum lenght
    table = FilterTracks(table, minlen)
    n_tracks = table.TRACK_ID.unique().max() 
    print('{:} tracks (filter = {} frames)'.format(int(n_tracks), minlen))    

    return table, rate, n_tracks

def TrajectoryClassification(all_tracks, track, thres, w, frame_rate, t_thresh):
        
    track = SingleTrack(all_tracks, track);
    track_score = ComputeConfinementRatioScore(track, conf_ratio_thres = thres, rol_window = w)
    track_score_denoised = DenoiseTransitionPoints(track_score , FindTransitionPoints(track_score), frame_rate, t_thresh)
    PlotConfinedRegions(track_score_denoised, track_score_denoised)
        
    print('time threshold: {} frames ~ {}s'.format(np.ceil(t_thresh / frame_rate), t_thresh))
    
    track_score_denoised_stats = ComputeSubSegmentStats(track_score_denoised, frame_rate)
    ShowStats(track_score_denoised_stats)
    
    return track_score_denoised, track_score_denoised_stats
	
def TrajectoryClassificationAllTracks(file, min_track_len = 30, tracks = 30, window = 5, p_thres = 1500, t_thres = 0.5, trim_trajectory_ends = False):
    
    all_tracks, frame_rate, n_tracks = ReadTracks(file, minlen = min_track_len)
    
    print('confinement threshold: {} '.format(p_thres))
    print('time threshold: {} frames ~ {}s'.format(np.ceil(t_thres / frame_rate), t_thres))
    
    all_tracks_stats = AnalyzeAllTracks(all_tracks, frame_rate, tracks = tracks, w = window, p_thres = p_thres, t_thres = t_thres)
    
    # discard all confined motions identified at the begnning and end of each track
    if trim_trajectory_ends == True:
        all_tracks_stats = TrimEnds(all_tracks_stats)
        
    lifetimes = ExtractConfinementLifetimes(all_tracks_stats)
    PlotLifetimes(lifetimes)
    
    # save lifetimes table in the output folder
    lifetimes.to_csv('outputs/Lifetimes_ConfinedMotion{}.csv'.format(file.split('/')[-1][:-4]), index = False)
    
    return lifetimes