#############################################################################
# 
# This harness should utilize all of the various filters. That are created from
# the design steps. 
#
###############################################################################

#Check OS version and set newline character accordingly
from sys import platform
if platform == 'linux' or platform == 'linux2':
    newline = '\n'
if platform == 'darwin':
    newline = '\n'
elif platform == 'win32':
    newline = '\r\n'

#Setup PyQt for using with the matplotlib
from PyQt4 import QtGui, QtCore
import matplotlib
matplotlib.use("Qt4Agg") #force matplotlib to use Qt4
import matplotlib.pyplot as plt

from grpdelay import group_delay
import scipy as sp
import scipy.signal as sig
import scipy.io.wavfile as wav
import math
import numpy as np

#define a visual display switch
visual = True

#Define a variable to tweak computations
computations = 512

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

#import the dirty audio
audio_input = wav.read("rawaudio/noisy.wav")

#Output formatting junk
divide_line = "############################################################################" + newline

#Define a function to create output graphs as well as other information.
def output(order, n, d, z, p, k):
    out_file = open('outputs/output.txt', 'a') #open file for output information
    out_file.write(divide_line + test + newline)
    out_file.write("The order of the filter is: " + str(2*N) + newline)
    #now work on the graphing
    #start with magnitude and phase graphing
    plt.figure(1)
    plt.clf()
    ax = plt.axes()
    axes = [ax, ax.twinx()]
    w, h = sig.freqz(n, d, computations)
    phase = sig.savgol_filter(np.unwrap(np.angle(h)),3,1)
    #h = sig.savgol_filter(h, 11, 2) #Smooth that data
    axes[0].plot(w/max(w)*nyquist, 20*np.log10(abs(h)),color='Blue') #plotphase
    axes[0].set_ylim(-160, 5)
    axes[0].set_ylabel('Magnitude in dB')
    axes[0].grid(True,linestyle='--')
    axes[1].plot(w/max(w)*nyquist, phase, color='Green') #plot phase
    axes[1].set_ylabel('Phase')
    plt.xlim(1300, 2000)
    plt.xlabel('Frequency in Hz')
    plt.title('Frequency Response')
    plt.savefig('outputs/' + test + '_magn_phase.png') #write out image
    #Plot zero-pole
    plt.figure(2)
    plt.clf()
    plt.plot(np.real(z),np.imag(z),'or') #plot the zeroes
    plt.plot(np.real(p),np.imag(p),'xb') #plot the poles
    t = np.linspace(0,np.pi*2,100) #generate a time interval for unit circle 
    plt.plot(np.cos(t),np.sin(t),'b--')
    plt.xlim(-1.1, 1.1) #Set the x axis
    plt.ylim(-1.1, 1.1) #Set the y axis
    plt.gca().grid(True, linestyle='--')
    plt.savefig('outputs/' + test + '_pz.png')
    #Plot group delay
    plt.figure(3)
    plt.clf()
    angle_h = np.unwrap(np.angle(h))
    #angle_h = np.angle(n)-np.angle(d)
    h_diff = (angle_h[0:len(h)-2]-angle_h[2:len(h)])/2
    grpdelay = -np.diff(np.unwrap(np.angle(h)))/np.diff(w) #calculating the group delay using linear derivative approx
    grpdelay = sig.savgol_filter(grpdelay, 29, 5) #Smooth that data
    plt.plot((w/max(w)*nyquist)[0:len(grpdelay)],grpdelay)
    plt.gca().grid(True, linestyle='--')
    plt.savefig('outputs/' + test + '_grpdelay.png')
    #Plot impulse response
    plt.figure(4)
    plt.clf()
    t,(im,) = sig.dimpulse([n,d,1])
    plt.stem(t,im)
    plt.xlim(-1, 100)
    plt.gca().grid(True, linestyle='--')
    plt.savefig('outputs/' + test + '_impulse.png')
    #Plot the graphs
    if(visual):
        plt.show()
    return 0



#apply some filters to the audio
##############################################################################
# Butterworth Filter
##############################################################################
test = "butterworth"
N, Wn = sig.buttord(passband_norm,stopband_norm, .5, stop_val, analog=False)
n, d = sig.butter(N, Wn, btype='bandstop',analog=False, output='ba') 
z, p, k = sig.butter(N, Wn, btype='bandstop',analog=False, output='zpk') 
output(N, n, d, z, p, k)
##############################################################################
# Chebyshev 1
##############################################################################
test = "chebyshev1"
N, Wn = sig.cheb1ord(passband_norm,stopband_norm,.5, stop_val, analog=False)
n, d = sig.cheby1(N, Wn, btype='bandstop', analog=False, output='ba')
z, p, k = sig.cheby1(N, Wn, btype='bandstop', analog=False, output='zpk')
output(N, n, d, z, p, k)
