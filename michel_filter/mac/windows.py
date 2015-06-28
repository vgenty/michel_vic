import scipy.stats
from   scipy.misc import comb
import numpy as np
import sys

def windowed_means(data,window_size,p_down,p_up) :
    if(window_size % 2 == 0):
        print "Give an odd window size dear"
        sys.exit()

    w = window_size+2
    w = (w - 1)/2
    num = len(data)
    
    mean_window = []
    means       = []

    for i in xrange(1,num+1):
        if  (i < w) :
            means = np.asarray(data[0:2*(i%w) - 1])
        elif( i > num - w + 1) :
            means = np.asarray(data[num - 2*((num - i)%w)-1:num])
        else :
            means = np.asarray(data[i - w: i + w - 1])

        lower_limit = np.percentile(means, p_down)
        upper_limit = np.percentile(means, p_up)
        mean_window.append(scipy.stats.tmean(means,
                                             limits=(lower_limit,upper_limit),
                                             inclusive=(True, True)))

    return mean_window

def coeff(k,N) :
    m = (N - 3.0)/2.0
    return 1.0/np.power(2,2*m+1) * (comb(2*m,m-k+1) - comb(2*m,m-k-1))

def smooth_derive(f,x,N) :
    M   = int((N - 1.0)/2.0)
    tot = 0.0
    
    for k in xrange(M) :
        tot += coeff(k+1,N) * (f[k+M] - f[M - 1 - k])/(x[k + M] - x[M - 1 - k]) * 2 * (k+1)
        
    return tot
