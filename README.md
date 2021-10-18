# Transient-Confinement-Analysis
### Identify when a given molecule switches between free diffusion and confined motion

**Tune-Parameters-Single-Track-Confinement-Classification Notebook** <br>
Inspect Tracks Individually and Tune Parameters

**Multi-Track-Confinement-Classification Notebook** <br> 
Run Analysis for all tracks in a XML file: returns the lifetime of confined segments vs. free diffusion segments

#### Overview

Here, we use the confinement ratio (p) to identify periods of confined motion along a free diffusion trajectory (when they exist).
p is defined as the **length of the trajectory in a short time window and the surface area that it occupies.**
This gives an estimate of the degree of free movement that a molecule displays in a period independently of its global diffusivity.
A rolling window with size w is used to compute p for each subsegment. Each time point will have a characteristic p, based on the behavior of the 
following w time points. Periods of confinement are identified by setting a threshold (p_thres) corresponding to a certain confinement area size, 
since p scales with the size of the confinement area. Then it is possible to calculate the frequency and duration of confinement periods and to localize them in space
This approach is adapted from Renner et al. (2017) and implemented here as an easy-to-use python script.

The use of a threshold value of p (p_thresh) and a minimal duration above this threshold (t_thresh) can suppress the detection of apparent nonrandom behaviors without excluding the detection 
of real confinement. These parameters depend on the acquisition frequency and the characteristic time of confinement. 
Too large windows will not detect properly the confinement period, while the statistical uncertainty increases in shorter windows. 
Thus, to accurately detect confinement periods, the window size should be adjusted accordingly to the acquisition rate and other experimental features.

### Parameters

The script runs on a simple jupyther notebook using a XML file (contains trajectories and coordinates - typical TrackMate's output file Fiji)

- window (w): is the length of the time window (in frames)
- confinement ratio threshold (p_thres): value above which a subsegment is considered confined (dimensionless)
- confinement time threshold (t_thres): minimum time (in sec) a confined event needs to last to be considered (converted in frames along the script)

*optional*
- min_track_len: discards trajectories shorter than this value (in frames) when reading the XML file
- trim_trajectory_ends: if True, discards confined events identified at the end or start of each trajectory

### References
Analysis is based on Renner et al. (2017) approach: https://doi.org/10.1016/j.bpj.2017.09.017 <br>
*"A Simple and Powerful Analysis of Lateral Subdiffusion Using Single Particle Tracking"*