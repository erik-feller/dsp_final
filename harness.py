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

#Set variables for the stopband frequencies in case they might change
leftstop = 1600
rightstop = 1700
leftpass = 1400
rightpass = 1900
nyquist = 5512.5

#Now change the values into ratios so that SciPy functions can accept them
stopband_ratio = [(leftstop/nyquist),(rightstop/nyquist)]
passband_ratio = [(leftpass/nyquist),(rightpass/nyquist)]

#import the dirty audio
audio_input = wav.read("rawaudio/noisy.wav")

#apply some filters to the audio
##############################################################################
# Butterworth Filter
##############################################################################



