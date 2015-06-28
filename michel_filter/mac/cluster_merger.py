import numpy as np
import ROOT as rr
from scipy.stats.stats import pearsonr


# Check endpoints, if endoints of two clusters are ``near" each other
# and they have good covariance in their two respective windows, 
# then merge the two clusters

# no covariance first...

def check_boundaries(clusters):
    
    yes = True
    poop = 0
    while yes:
        poop = 0
        final_objects = []
        for c in clusters:
            for k in clusters:
                if(c != k and
                   c in clusters and k in clusters):
                    if(c.touching(k)) :
                        final_objects.append(c + k)
                        clusters.remove(c)
                        clusters.remove(k)
                        poop += 1

        clusters += final_objects
        if(poop == 0): 
            yes = False

    
    return clusters

