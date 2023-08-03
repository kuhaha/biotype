from BaselineRemoval import BaselineRemoval
from scipy import signal

def find_peaks(y, cwt=False, smooth=False, baseline=False):
    window, deg = 15, 2     
    z = y
    if smooth:
        z = signal.savgol_filter(y, window, deg, deriv=0)

    if baseline:
        brm = BaselineRemoval(z)
        z = brm.ZhangFit(lambda_=400,repitition=15, porder=1)

        if cwt == True:
        peaks = signal.find_peaks_cwt(z, [20])  
    else:
        dist, prom = 100, 600
        peaks,_ = signal.find_peaks(z, distance=dist, prominence=prom)
    
    return peaks, z