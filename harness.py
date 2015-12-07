###############################################################################
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

import scipy.signal as sig
import scipy.io.wavfile as wav
import math
import numpy as np

#define a visual display switch
visual = True

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
    ax = plt.axes()
    axes = [ax, ax.twinx()]
    w, h = sig.freqz(n, d)
    axes[0].plot(w/max(w)*nyquist, 20*np.log10(abs(h)),color='Blue') #plotphase
    axes[0].set_ylim(-160, 5)
    axes[0].set_ylabel('Magnitude in dB')
    axes[1].plot(w/max(w)*nyquist,np.unwrap(np.arctan2(np.imag(h),np.real(h))), color='Green') #plot phase
    axes[1].set_ylabel('Phase')
    plt.xlim(1300, 2000)
    plt.xlabel('Frequency in Hz')
    plt.title('Frequency Response')
    plt.savefig('outputs/' + test + '_magn_phase.png') #write out image
    #Plot zero-pole
    plt.figure(2)
    plt.plot(np.real(z),np.imag(z),'or') #plot the zeroes
    plt.plot(np.real(p),np.imag(p),'xb') #plot the poles
    t = np.linspace(0,np.pi*2,100) #generate a time interval for unit circle 
    plt.plot(np.cos(t),np.sin(t),'b--')
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,1.1)
    plt.savefig('outputs/' + test + '_pz.png')
    #Plot group delay
    plt.figure(3)
    grpdelay = -np.diff(np.unwrap(np.angle(h)))/np.diff(w) #calculating the group delay using linear derivative approx
    #print(len(-np.diff(np.unwrap(np.angle(h)))))
    #print(len(np.diff(w)))
    plt.plot(w[0:len(grpdelay)]/max(w)*nyquist,grpdelay)
    #Plot impulse response
    plt.figure(4)
    t,(im,) = sig.dimpulse([n,d,1])
    plt.stem(t,im)
    plt.xlim(-1, 100)
    if(visual):
        plt.show()
    return 0



#apply some filters to the audio
##############################################################################
# Butterworth Filter
##############################################################################
test = "butterworth"
N, Wn = sig.buttord(passband_norm,stopband_norm, 3, stop_val, analog=False)
n, d = sig.butter(N, Wn, btype='bandstop',analog=False, output='ba') 
z, p, k = sig.butter(N, Wn, btype='bandstop',analog=False, output='zpk') 
output(N+1, n, d, z, p, k)
