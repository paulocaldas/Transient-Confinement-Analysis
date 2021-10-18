# Transient-Confinement-Analysis
### Identify when a given molecule switches between free diffusion and confined motion

**Tune-Parameters-Single-Track-Confinement-Classification Notebook** <br>
Inspect Tracks Individually and Tune Parameters

**Multi-Track-Confinement-Classification Notebook** <br> 
Run Analysis for all tracks in a XML file: returns the lifetime of confined segments vs. free diffusion segments

### Overview

Here, we use the confinement ratio (p) to identify periods of confined motion along a free diffusion trajectory (if they exist).
p is defined as the **length of the trajectory in a short time window and the surface area that it occupies.**
This gives an estimate of the degree of free movement that a molecule displays in a period independently of its global diffusivity.
This approach is adapted from Renner et al. (2017) and implemented here as an easy-to-use python script.

### Parameters

p= ∑_i^(i+n-1)▒((x_(i+1)- x_i )^2- (y_(i+1)- y_i )^2)/〖S_i〗^2 



[write some documentation here]
