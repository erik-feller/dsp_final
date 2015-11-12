###############################################################################
#
# butterworth.py
# Erik Feller
# 2015-11-11
# a function to create a butterworth filter for the specified stop and pass 
# bands with lowest order possible. 
#
###############################################################################
import math

def butterworth_build(cutoff_freq, stopband_freq):
    cutoff_ang = math.pi*cutoff_freq
    stopban_ang = math.pi*stopband_freq

