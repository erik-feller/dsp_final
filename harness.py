###############################################################################
# 
# This harness should utilize all of the various filters. That are created from
# the design steps. 
#
###############################################################################

#Setup PyQt for using with the matplotlib
from PyQt4 import QtGui, QtCore
import matplotlib
matplotlib.use("Qt4Agg") #force matplotlib to use Qt4
import matplotlib.pyplot as plt

import scipy.signal as sig
import scipy.io.wavfile as wav
import math
import numpy as np

#Set variables for the stopband frequencies in case they might change
leftstop = 1600
rightstop = 1700
leftpass = 1400
rightpass = 1900
nyquist = 5512.5

#Define values for the tolerances allowed
pass_tol = 0.5
stop_val = 100 #given in negative dB

#Now change the values into ratios so that SciPy functions can accept them
stopband_norm = [(leftstop/nyquist),(rightstop/nyquist)]
passband_norm = [(leftpass/nyquist),(rightpass/nyquist)]
print(passband_norm)
print(stopband_norm)

#import the dirty audio
audio_input = wav.read("rawaudio/noisy.wav")

#apply some filters to the audio
##############################################################################
# Butterworth Filter
##############################################################################
test = "butterworth"
N, Wn = sig.buttord(passband_norm,stopband_norm, 3, stop_val, analog=False)
n, d = sig.butter(N, Wn, btype='bandstop',analog=False, output='ba') 
w, h = sig.freqz(n, d)
plt.plot(w/max(w)*nyquist, 20*np.log10(abs(h)))
plt.show()


#Define a function to create an output file with order and other information.
def output(order):
    return 0



