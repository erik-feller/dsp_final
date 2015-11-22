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

import scipy.signal as sig
import scipy.io.wavfile as wav
import math

#import the dirty audio
audio_input = wav.read("rawaudio/noisy.wav")

#apply some filters to the audio






###############################################################################
#
# Define a set of functions to design each of the required filters. 
#
###############################################################################

#Butterworth 
def butter_filt(stopband_freq, passband_freq):
    stopband_rads = stopband_freq/(math.pi*2)
    passband_rads = passband_freq/(math.pi*2)
    N, Wn = sig.buttord(
