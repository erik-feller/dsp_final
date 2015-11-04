"""
Copyright (c) 2008 Christopher Terman

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import numpy,wx
import blocks,wxplot
from filters import compute_window
from sinks import sample_sink

# helper function looks at range of values and returns a new range,
# an engineering prefix and the scale factor
def eng_notation(range):
    x = max(abs(range[0]),abs(range[1]))   # find largest value
    if x <= 1e-12:
        scale = 1e15
        units = "f"
    elif x <= 1e-9:
        scale = 1e12
        units = "p"
    elif x <= 1e-6:
        scale = 1e9
        units = "n"
    elif x <= 1e-3:
        scale = 1e6
        units = "u"
    elif x <= 1:
        scale = 1e3
        units = "m"
    elif x >= 1e9:
        scale = 1e-9
        units = "G"
    elif x >= 1e6:
        scale = 1e-6
        units = "M"
    elif x >= 1e3:
        scale = 1e-3
        units = "k"
    else:
        scale = 1
        units = ""
    return ((range[0]*scale,range[1]*scale),units,scale)

################################################################################
# plot_frame: base class for all plots
################################################################################

# Override GetPlotData and your method return a list wxplot.PolyXXX objects to 
# be plotted; data is in self.buffer
class plot_frame(wx.Frame):
    def __init__(self,buffer,title='plot_frame',position=None,size=wx.DefaultSize,
                 plot_title='',xlabel=None,ylabel=None,
                 xaxis=None,yaxis=None,step=None,
                 peak_hold = False, average = False,
                 zoom=True,grid=True):
        wx.Frame.__init__(self,None,-1,title=title,size=size,pos=position)

        # remember various parameters
        self.plot_title = plot_title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.step = step

        # where we'll keep the arriving data
        self.buffer = buffer

        wx.EVT_CLOSE(self,self.OnClose)   # turn off data delivery when closed

        self.average = average
        self.peak_hold = peak_hold
        self.peak_vals = None

        # set up plot window
        self.plot = wxplot.PlotCanvas(self)
        self.plot.SetEnableGrid(grid)
        self.plot.SetEnableZoom(zoom)

        vsizer = wx.BoxSizer(wx.VERTICAL)
        vsizer.Add(self.plot,1,wx.EXPAND)
        self.SetSizer(vsizer)
        self.Layout()
        if position is None:
            self.SetPosition(blocks.current_testbench.get_position(self.GetSize()))
        self.Show(True)

    def OnClose(self,event):
        # tell whoever is sending us events to stop!
        pass

    # return list of plot.PolyXXX objects, override as appropriate
    def GetPlotData(self):
        return []

    def ProcessSample(self,sampler):
        sampler.get_sample(self.buffer)
        pdata = self.GetPlotData()
        graphics = wxplot.PlotGraphics(pdata,title=self.plot_title,
                                       xLabel=self.xlabel,yLabel=self.ylabel)
        self.plot.Draw(graphics,xAxis=self.xaxis,yAxis=self.yaxis,step=self.step)

################################################################################
# plot_fft: show fft of sampled data
################################################################################

# customized plot_frame for ffts
class fft_frame(plot_frame):
    def __init__(self,buffer,xrange,**params):
        plot_frame.__init__(self,buffer,**params)

        self.fft_size = len(buffer)
        self.points = numpy.zeros((self.fft_size, 2),numpy.float64)
        xstep = (xrange[1] - xrange[0])/float(self.fft_size-1)
        self.points[:,0] = numpy.arange(self.fft_size)*xstep + xrange[0]
        self.window = compute_window('blackman-harris',self.fft_size)
        power = 0.0
        for tap in self.window: power += tap*tap
        self.power_adjust = -10*numpy.log10(self.fft_size) - 10*numpy.log10(power/self.fft_size)

    def GetPlotData(self):
        numpy.multiply(self.buffer,self.window,self.buffer)
        fft = numpy.fft.fft(self.buffer)
        # convert to power spectrum in db
        numpy.absolute(fft,fft)
        numpy.log10(fft,fft)
        numpy.multiply(fft,20,fft)

        l = self.fft_size/2

        # adjust power for number of bins and windowing
        numpy.add(fft,self.power_adjust,fft)
        # ??? add 3dB to all bins except 0 and midpoint
        #numpy.add(fft[1:l],3,fft[1:l])
        #numpy.add(fft[l+1:],3,fft[l+1:])

        # fftshift without consing...
        self.points[0:l,1] = fft[l:]
        self.points[l:,1] = fft[0:l]
        del fft

        return [wxplot.PolyLine(self.points,colour='BLUE'),]

def plot_fft(input,
             title = 'FFT',
             size = wx.DefaultSize,
             position = None,
             plot_rate = 5,
             sample_size = 512,
             xaxis = None,
             yaxis = (-160,0),
             peak_hold = False,
             average = False
             ):
    # set up buffer to hold sample points
    buffer = numpy.empty(sample_size,numpy.float64)

    # set up xaxis
    nyquist = input.sample_rate/2
    x_range,x_prefix,x_scale = eng_notation((-nyquist,nyquist))
    if xaxis is None: xaxis = x_range
    else: xaxis = (xaxis[0]*x_scale,xaxis[1]*x_scale)

    # make a plot window for user to see results
    frame = fft_frame(buffer,x_range,
                      title=title, size=size, position=position,
                      plot_title = title,
                      xaxis=xaxis, xlabel=x_prefix+'Hz',
                      yaxis=yaxis, ylabel='dB',
                      peak_hold=peak_hold, average=average)

    # set up sampler block to deliver data samples to plot window
    return sample_sink(input,sample_size=sample_size,interval=(1.0/plot_rate),wxclient=frame)

################################################################################
# plot_data: oscilloscope-like display of plot data
################################################################################

# customized plot_frame for scope-like displays
class scope_frame(plot_frame):
    def __init__(self,buffer,xrange,**params):
        plot_frame.__init__(self,buffer,**params)

        self.sample_size = len(buffer)
        self.points = numpy.zeros((self.sample_size, 2),numpy.float64)
        xstep = (xrange[1] - xrange[0])/float(self.sample_size-1)
        self.points[:,0] = numpy.arange(self.sample_size)*xstep + xrange[0]

    def GetPlotData(self):
        self.points[:,1] = self.buffer
        return [wxplot.PolyLine(self.points,colour='BLUE'),]

def plot_data(input,
              title = 'Scope',
              size = wx.DefaultSize,
              position = None,
              interval = .1,
              sample_size = 1000,
              xaxis = None,
              yaxis = (-1.0,1.0),
              peak_hold = False,
              average = False
              ):
    # set up buffer to hold sample points
    buffer = numpy.empty(sample_size,numpy.float64)

    # set up xaxis
    x_range,x_prefix,x_scale = eng_notation((0,float(sample_size - 1)/input.sample_rate))
    if xaxis is None: xaxis = x_range
    else: xaxis = (xaxis[0]*x_scale,xaxis[1]*x_scale)

    # make a plot window for user to see results
    frame = scope_frame(buffer,x_range,
                        title=title, size=size, position=position,
                        plot_title = title,
                        xaxis=xaxis, xlabel=x_prefix+'s',
                        yaxis=yaxis, ylabel='data',
                        peak_hold=peak_hold, average=average)

    # set up sampler block to deliver data samples to plot window
    return sample_sink(input,sample_size=sample_size,interval=interval,wxclient=frame)



__all__ = [
    'plot_fft',
    'plot_data',
    ]
