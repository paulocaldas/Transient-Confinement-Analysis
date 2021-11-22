'''support functions'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_trackmate_xml_tracks(xml_file):
    """import tracks from trackmate xml track file and converts into a user-friendly DataFrame """
    from xml.etree import cElementTree as ET
    tracks = ET.parse(xml_file)
    frame_interval = float(tracks.getroot().attrib["frameInterval"])
    n_tracks = float(tracks.getroot().attrib["nTracks"])
    
    attributes = []
    for ti, track in enumerate(tracks.iterfind('particle')):
        for spots in track.iterfind('detection'):
            attributes.append([ti, int(spots.attrib.get('t')),
                                   float(spots.attrib.get('x')),
                                   float(spots.attrib.get('y'))])

    track_table = pd.DataFrame(attributes, columns=['TRACK_ID','FRAME','POSITION_X','POSITION_Y'])
    return track_table, frame_interval, n_tracks

def FilterTracks(traj_table, minlen):
    '''filter tracks shorter than minlen (in frames)'''
    
    traj_table = traj_table.groupby('TRACK_ID').filter(lambda track: track.TRACK_ID.count() > minlen)
    
    # create a dictionary to attribute a new number to each track ID
    new_ids = {}
    for i, value in enumerate(traj_table.TRACK_ID.unique()):
        new_ids[value] = i
    
    # replace existent IDs with new ID's using the dicitonary above
    traj_table['TRACK_ID'].replace(new_ids, inplace = True)
    
    # reset index (just in case)
    traj_table.reset_index(drop = True, inplace = True)
    
    return traj_table

def SingleTrack(data, track):
    '''select one track from all_tracks table'''
    track = data[data['TRACK_ID'] == track]
    print('track with {} steps'.format(len(track)))
    
    return track

def plotSingleTrack(track, color = 'peachpuff', start = 0 , end = -1):
   '''plot a single track, from start to end position'''
   
   fig, ax = plt.subplots(figsize = (4,3), dpi = 120)
    
   # plot window to add sub segments
   if end == -1: end = track.shape[0]
   track = track.reset_index(drop=True)
   track = track.loc[start:end]
    
   plt.plot(track.POSITION_X, track.POSITION_Y, '-o', lw = 2, color = color,
             markersize = 6, markeredgecolor = 'black', markeredgewidth = 0.8); plt.axis('off');
     
def color_pick(p, thres = 1000, colors = ['steelblue', 'crimson']):
    if p >= thres:
        color = colors[1]
    elif p < thres:
        color = colors[0]
    return color

def ComputeDirectionality(track_table, coord = ['POSITION_X','POSITION_Y']):
    '''estimates distance, displacement and directionality for a given track '''

    steps = []
    
    for frame in range(1,len(track_table)):
        step_dist = np.linalg.norm(track_table[coord].iloc[frame] - (track_table[coord].iloc[frame-1]))
        steps.append(step_dist)

    displacement = np.linalg.norm(track_table[coord].iloc[0] - track_table[coord].iloc[-1])
    distance = np.sum(steps)
    directionality = displacement/np.sum(steps)
    
    return distance, displacement, directionality

def DefineYmax(ax):
    '''defines ylim max - this function is irrelevant'''
    
    if np.ceil(max(ax.get_ylim()))%2 == 0:
        ymax = np.ceil(max(ax.get_ylim()))
    else: 
        ymax = np.ceil(max(ax.get_ylim())) + 1
    
    return ymax