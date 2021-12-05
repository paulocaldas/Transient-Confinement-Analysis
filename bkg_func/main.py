'''single track - main function '''

import pandas as pd
import numpy as np
import os, glob, sys
#import matplotlib.pyplot as plt
#import seaborn as sns
#from matplotlib.ticker import FormatStrFormatter

from bkg_func.utils import read_trackmate_xml_tracks, FilterTracks, SingleTrack
from bkg_func.conf_ratio_func import ComputeConfinementRatioScore, FindTransitionPoints, DenoiseTransitionPoints
from bkg_func.conf_ratio_func import ComputeSubSegmentStats, PlotConfinedRegions, ShowStats
from bkg_func.multiple_track_func import AnalyzeAllTracks, PlotOutputStuff, TrimEnds


def ReadTracks(file, minlen):
    '''reads file into a dataframe, extracts frame rate and number of tracks'''
    table, rate, n_tracks = read_trackmate_xml_tracks(file)
    
    # apply filter for minimum lenght
    table = FilterTracks(table, minlen)
    n_tracks = table.TRACK_ID.unique().max() 
    print('{} total tracks; filtering shorter than {} frames ({:.3} sec)'.format(int(n_tracks), minlen, minlen*rate))    

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
	
def TrajectoryClassificationAllTracks(file, min_track_len = 30, tracks = 30, window = 5,
                                      p_thres = 1500, t_thres = 0.5, trim_trajectory_ends = False,
                                      output_folder = None):
    
    '''save all info from each track in an output table
    by default output_folder is the same as the file directory'''
    
    all_tracks, frame_rate, n_tracks = ReadTracks(file, minlen = min_track_len)
    
    print('p_thresh = {}nm-2; t_thresh = {}s ({} frames)'.format(p_thres, t_thres,
                                                              np.ceil(t_thres / frame_rate)))
    
    all_tracks_stats = AnalyzeAllTracks(all_tracks, frame_rate, tracks = tracks, w = window, p_thres = p_thres, t_thres = t_thres)
    
    # discard all confined motions identified at the begnning and end of each track
    if trim_trajectory_ends == True:
        all_tracks_stats = TrimEnds(all_tracks_stats)
        
    # merge all tracks into a master table and save output
    all_tracks_stats = pd.concat(all_tracks_stats, ignore_index = True)
    
    if output_folder == None:
        filename = file.rsplit('.', maxsplit=1)[0] + '_TransientConfinementClassification.tab'
        all_tracks_stats.to_csv(filename, index = False, sep = '\t')
        
    else:
        get_dir = output_folder
        filename = os.path.basename(file).rsplit('.', maxsplit=1)[0] + '_TransientConfinementClassification.tab'
        all_tracks_stats.to_csv(get_dir + '/' + filename, index = False, sep = '\t')
    
    # plot some of the info for good looks
    PlotOutputStuff(all_tracks_stats)
    
    return all_tracks_stats

def TrajectoryClassificationAllFiles(folder, min_track_len = 30, tracks = 30, window = 5,
                                     p_thres = 1500, t_thres = 0.5, trim_trajectory_ends = False,
                                     output_folder = None):
    
    '''save output table for all xml files in folder
    by default output_folder is the same as the file directory'''
    
    #print('files to analyze : {}'.format(glob.glob(folder + '/*.xml')))
    
    for file in glob.glob(folder + '/*.xml'):
        
        print('Analyze All Tracks for {}'.format(os.path.basename(file)))
        
        with HiddenPrints():
            
            all_tracks_stats = TrajectoryClassificationAllTracks(file, min_track_len = min_track_len, tracks = tracks, 
                            window = window, p_thres = p_thres, t_thres = t_thres, 
                            trim_trajectory_ends = trim_trajectory_ends,
                            output_folder = output_folder)
    print('done')
    
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout    

    
    
    
    