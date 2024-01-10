'''set of functions to run the confinement ratio analysis on a single track 
from an XLM file with trajectories '''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.spatial import ConvexHull
from bkg_func.utils import ComputeDirectionality

# FUNCTIONS TO COMPUTE THE CONFINEMENT RATIO ALONG A GIVEN TRACK

def ComputeSquareDisplacement(track_table, coord = ['POSITION_X','POSITION_Y']):
    '''estimates the sum of square distances along step size of a given track '''
    
    steps = []
    
    for frame in range(1,len(track_table)):
        step_dist = np.linalg.norm(track_table[coord].iloc[frame] - (track_table[coord].iloc[frame-1]))
        steps.append(step_dist**2) #append squared distances between time points
    
    # sum all squared distances between time points
    square_displacement = np.sum(steps)
    
    return square_displacement

def FitElipse(track_table, X='POSITION_X', Y = 'POSITION_Y'):
    '''Fits and elipse around all data points
    returns elipse coordinates and area of the convex hull surroding data points'''
    
    # get convex hull for data points
    hull = ConvexHull(track_table[[X,Y]]) 
    
    # use hull vertices to get area around the data points
    x = track_table.reset_index(drop = True).iloc[hull.vertices].POSITION_X
    y = track_table.reset_index(drop = True).iloc[hull.vertices].POSITION_Y
    
    # smooth area around the data points with an elipse
    
    xmean, ymean = x.mean(), y.mean()
    x -= xmean
    y -= ymean
    U, S, V = np.linalg.svd(np.stack((x, y)))

    tt = np.linspace(0, 2*np.pi, 1000)
    circle = np.stack((np.cos(tt), np.sin(tt)))    # unit circle
    transform = np.sqrt(2/len(y)) * U.dot(np.diag(S))   # transformation matrix
    fit = transform.dot(circle) + np.array([[xmean], [ymean]])
    
    return fit, hull.volume

def ConfinementRatio(track_table, X='POSITION_X', Y = 'POSITION_Y'):
    ''' relates the sum of the square displacements between successive
    time points by the total surface area of the track '''
    
    squared_distances = ComputeSquareDisplacement(track_table)
    cage, cage_area = FitElipse(track_table, X = X, Y = Y) 
    
    segment_pc = squared_distances/(cage_area**2)
    
    return segment_pc
	
# FUNCTIONS TO COMPUTE AND VISUALIZE CONFINEMENT RATIO ALONG A GIVEN TRACK

def TuneConfinementThreshold(track, frame_rate, conf_thres = 1500, windows = [5, 8, 12]):
    '''computes the confinement ratio using different window sizes (windows)
    and returns a list with conf_ratio per position for each window size'''
    
    print('rolling windows: {} frames ~ {} sec'.format(windows, list(np.array(windows) * frame_rate)))
    
    conf_ratios = []  # save all confinement ratios and respective rolling window size

    for i, window_size in enumerate(windows, start = 1):
    
        print('rolling window {0} : {1:.4} sec'.format(i, window_size * frame_rate), end ="\r")
    
        for window in np.arange(0, track.shape[0] - window_size, 1):
       
            X = track.POSITION_X.iloc[window: window + window_size + 1]
            Y = track.POSITION_Y.iloc[window: window + window_size + 1]

            segment = pd.DataFrame([X,Y]).T # segment to analyze

            # compute confinement ratio
            conf_ratio = ConfinementRatio(segment)
            conf_ratios.append([conf_ratio, window_size])
    
    return conf_ratios

def PlotParameterTunnig(param, frame_rate, thres = 800, ylim = None, log = False,):
    '''takes the output of TuneConfinementThreshold to plot the results'''
    
    rol_windows = np.unique([j for i,j in param]) # reads each par coef_ratio, window_size
    df = pd.DataFrame(param, columns=['conf_ratio','window_size']) # converts into a pretty table
    
    fig, ax = plt.subplots(figsize = (8, 3), dpi = 120)
    # plot conf_ratio vs. position along the track, for each window_size
    for w in rol_windows:
        x = np.arange(0, df[df['window_size'] == w].shape[0]) * frame_rate
        y = df[df['window_size'] == w]['conf_ratio']
        
        plt.plot(x,y, lw= 1.5, label = w, color = plt.cm.tab20b(w))
    
    plt.legend(title= "wind_size", fontsize=7, title_fontsize=7, loc = 0, frameon = False)
    ax.set_xlabel('start position (sec)'); 
    ax.set_ylabel('confinement ratio ($\mathregular{nm^{-2}}$)');
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'));
    
    # threhsold for confinement ratio
    plt.hlines(y = thres, xmin = np.min(x), xmax = np.max(x), alpha = 0.2, ls = '--', color = 'crimson')
    
    # fancy stuff
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.tick_params(direction = 'inout', top=False, right = False, bottom = False, labelsize=8)
    
    # set ylim 
    
    if ylim == None: 
        ax.set_ylim([0, np.max(y) + thres])
    else:
        ax.set_ylim([0, ylim])
        #'tracks with confined regions: {}/{}'.format(tracks_w_confinement, len(all_tracks_stats))
        
    if log == True: # set y-axis to log scale 
        plt.yscale('symlog')
	

def ComputeConfinementRatioScore(track, conf_ratio_thres = 1500, rol_window = 5, coords = ['POSITION_X','POSITION_Y'], show_progress = True):
    '''computes de confinement ratio for each subsegment with lenght = rol_window (using step = 1) 
    
    track: input table track; columns should contain coords
    conf_ratio_thres: subsegments with a conf_ratio above this value are considered confined
    rol_window: lenght of each subsegment
    coords: column names of each X,Y coordinate pair
    '''
    conf_ratio_col = 'conf_ratio'
    
    # create a list that will contain the confinement ratio
    # for each position along the track
    conf_ratio_matrix = []
    
    #print('rolling window: {0:.4} sec'.format(window_size * frame_interval))

    for window in np.arange(0, track.shape[0] - rol_window, 1):
        
        if show_progress == True:
            print('rolling {0}/{1}'.format(window, track.shape[0]), end ="\r")
        
        # create subsegments using rol_window size
        X = track.POSITION_X.iloc[window: window + rol_window + 1]
        Y = track.POSITION_Y.iloc[window: window + rol_window + 1]
        segment = pd.DataFrame([X,Y]).T 

        # compute confinement ratio for each segment
        conf_ratio = ConfinementRatio(segment)

        # add conf_ratio to a list; keep only the first row for the main table
        conf_ratio_segment = segment.copy()
        conf_ratio_segment[conf_ratio_col] = conf_ratio
        conf_ratio_segment_first = conf_ratio_segment.iloc[:1] #conf_ratio_segment[~conf_ratio_segment.duplicated(conf_ratio_segment[conf_ratio_col])]
        conf_ratio_matrix.append(conf_ratio_segment_first)
    
        # Combine all matrix in one; colpase NaN values
        #conf_ratio_matrix_merged = reduce(lambda X,Y: pd.merge(X,Y, on = coords, how = 'outer'), conf_ratio_matrix)
        #conf_ratio_matrix_reduced = pd.concat([conf_ratio_matrix_merged.iloc[:,:2], conf_ratio_matrix_merged.iloc[:,2:].mean(axis = 1)], axis = 1)
        #conf_ratio_matrix_reduced.columns = coords + list([conf_ratio_col])

        # add column with a mode label for each data point: free diffusion or confined

    conf_ratio_matrix = pd.concat(conf_ratio_matrix)
    conf_ratio_matrix = conf_ratio_matrix.reset_index(drop = True)
    
    # add boolean a column: is the conf_ratio above the treshold?
    isitconfined = conf_ratio_matrix[conf_ratio_col].map(lambda x : True if (x > conf_ratio_thres) else False)
    conf_ratio_matrix['confined'] = isitconfined
           
    return conf_ratio_matrix

def FindTransitionPoints(track_score, col = 'confined'):
    ''''finds transition points along column col'''
    
    score_matrix = track_score.copy()

    # first and last datapoint are always added to the list
    transition_points = [0, score_matrix.shape[0]]
    
    # for each transition between modes, consider a transition point
    for i, row in score_matrix.iloc[:-1].iterrows():
        if score_matrix[col].iloc[i+1] != score_matrix[col].iloc[i]:
            transition_points.append(i+1)

    transition_points = np.unique(transition_points)

    return transition_points

def DenoiseTransitionPoints(track_score, transition_points,  frame_rate, t_thres = 0.5,  col = 'confined',):
    '''finds transition points along column 'col', discarding short subsegments (shorter than t_thres in sec)
    uses orginal transition_points as input. returns a new track_score with only significant transitions'''

    track_score = track_score.copy()
    min_len = np.ceil(t_thres / frame_rate)
    
    confined_tracks_index = []
    
    # rationale
    # if a confined subsegment is longer than min_len, we keep those coordinates 
    # the 'confined' column is True for this coordinates and False anywhere else
    
    for i, pos in enumerate(transition_points[:-1]):

        sub_segment = track_score.iloc[transition_points[i]:transition_points[i+1]]
        
        if len(sub_segment) >= min_len and sub_segment[col].iloc[-1] == True:
            confined_tracks_index.append(sub_segment.index.values)
        
    # combine multiple segment indexes into a list
    confined_tracks_index = [i for j in confined_tracks_index for i in j]

    # re-name column col for true confined spots
    track_score.loc[confined_tracks_index, col] = True
    track_score.loc[~track_score.index.isin(confined_tracks_index), col] = False # False for everything else
    
    return track_score

def PlotConfinedRegions(track, track_score, col = 'confined', colors = ['red','steelblue']):
    '''takes original track table and track_score table to plot trajectory color code'''
    
    fig, ax = plt.subplots(figsize = (4,3), dpi = 150)
    
    # plot all points from raw track in gray 
    plt.plot(track.POSITION_X, track.POSITION_Y, '-', lw = 1, markersize = 3,
             markeredgecolor = 'black', color = colors[1], alpha = 0.4)
    
    # color in red the confined data points from the score track table
    
    for i, row in track_score.iterrows():
            
        if row[col] == True: 
            color = colors[0]
            size = 3.5
        
        else: 
            color = colors[1]
            size = 3
    
        plt.plot(row.POSITION_X, row.POSITION_Y, 'o', lw = 1.2, markersize = size, 
                 markeredgecolor = 'black', markeredgewidth = 0.2, color = color, alpha = 0.5)
        
    # add starting point and end point
    kwargs = {'markeredgecolor': 'black', 'markersize':3, 'markeredgewidth': 0.2}
    plt.plot(track.POSITION_X.iloc[0], track.POSITION_Y.iloc[0], 'o', markerfacecolor = 'darkgreen', **kwargs)
    plt.plot(track.POSITION_X.iloc[-1], track.POSITION_Y.iloc[-1], 'o', markerfacecolor = 'red', **kwargs)
    
    # format plot
    ax.set_xlabel('X ($\mu m$)'); ax.set_ylabel('Y ($\mu m$)'); ax.axis('off');

def ComputeSubSegmentStats(track_score_denoised, frame_rate, col = 'confined'):
    '''computes the lifetime, cage area and other parameters of each subsegment in track_score_denoised (in seconds)
    using the transition points in the column col. 
    
    returns a data frame with all parameters for each segment along the track
    
    params:
    track_score_denoised: trajectory table containing a column col
    col: column, boolean, whether a given spot is in confined motion or not
    '''
    
    params = []
    transition_points = FindTransitionPoints(track_score_denoised)
    
    for i, pos in enumerate(transition_points[:-1]):

        sub_segment = track_score_denoised.iloc[transition_points[i]:transition_points[i+1]]
        
        # compute confinement ratio of each segment
        if sub_segment.shape[0] > 5: 
            conf_ratio = ConfinementRatio(sub_segment)
        else:
            conf_ratio = np.nan
        
        # compute lifetime of each subsegment
        lifetime = sub_segment.shape[0] * frame_rate

        # compute cage area of each subsegment - only if segment is larger than 3 steps 
        # free difussion segments can have less than 5 frames sometimes ...
        
        if sub_segment.shape[0] > 5: 
            cage_area = FitElipse(sub_segment)[1] * 1000 # in nanometer
        else: 
            cage_area = np.nan
            
        # compute total distance of each subsegment
        
        if sub_segment.shape[0] > 5:
            dist = ComputeDirectionality(sub_segment)[0] * 1000 # in nanometer
        else: 
            dist = np.nan
            
        # add parameters to list
        
        if sub_segment[col].iloc[0] == True:
            params.append(['confined', conf_ratio, lifetime, round(cage_area,0), round(dist,0), sub_segment.shape[0]])
        else:
            params.append(['free_diff', conf_ratio, lifetime, round(cage_area,0), round(dist,0), sub_segment.shape[0]])      
            
    return pd.DataFrame(params, columns=['mode', 'p_coeff_um-2', 'lifetime_s', 'area_nm2', 'distance_nm', 'steps'])

def ShowStats(stats_table):
    
    if any(stats_table['mode'] == 'confined'):
        
        n_conf = stats_table['mode'].value_counts().confined
        tau_conf = stats_table[stats_table['mode'] == 'confined']['lifetime_s'].mean()
        area_conf = stats_table[stats_table['mode'] == 'confined']['area_nm2'].mean()
        
        print('{} confined events; avg lifetime = {:.3}s; avg cage_area = {:.3} nm2'.format(n_conf, tau_conf, area_conf))
        
    else:
        print('no confined events') 
    return 