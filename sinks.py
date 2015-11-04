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

from blocks import *
import numpy

################################################################################
# sink: base class for all sink blocks
#  sink blocks are special since when they've all reached end-of-stream
#  and their _run() methods report that to the testbench, the simulation
#  will exit.
################################################################################

import blocks
class sink(dsp_block):
    def _reset(self):
        blocks.current_testbench.add_sink(self)
        dsp_block._reset(self)

################################################################################
# null_sink: consume inputs at approximately real-time
################################################################################

import time
class null_sink(sink):
    def _initialize(self):
        self._seconds_per_buffer = self._buffer_size/self.sample_rate
        self._wait_until = time.time() + self._seconds_per_buffer

    def _buffers_ready(self,inputs,outputs):
        now = time.time()
        if now >= self._wait_until:
            self._wait_until = now + self._seconds_per_buffer
            return True
        else:
            return False

################################################################################
# sample_sink: gather a sample of size N at a specified interval
################################################################################

# the sample is stored in an internal buffer and is updated atomically. It can
# be copied into a user-supplied numpy array by calling get_sample(array); the
# code worries about the appropriate synchronization.  User can supply a wx
# object which will get EVT_SAMPLE events when a new sample is ready.

import Queue,threading,wx
class sample_sink(sink):
    _default_parameters = {
        'sample_size': None,    # number of data items in a sample
        'interval': None,       # how long (in seconds) between start of samples
        'wxclient': None,       # (optional) who receives notifications
        'queue': None           # (optional) Queue for samples
        }
    _expected_number_of_inputs = (1,1)

    def _initialize(self):
        assert isinstance(self.sample_size,int) and self.sample_size > 0,\
               "%s: sample_size must be a positive integer" % self._name
        if self.interval:
            self._nskip = max(0,int(self.interval * self.sample_rate) - self.sample_size)
        else:
            # continuous sampling
            self._nskip = 0

        # buffer where we collect a sample
        self._temp = numpy.empty(self.sample_size,dtype=self._inputs[0].dtype)
        # buffer to hold sample for user to copy
        self._sample = numpy.empty(self.sample_size,dtype=self._inputs[0].dtype)

    # make a copy of the current sample for the user.  Grab the lock so that
    # we synchronize properly with this block running in the background
    def get_sample(self,dest):
        self._mutex.acquire()
        dest[:] = self._sample
        self._mutex.release()
        return dest

    def _reset(self):
        self._skip = self._nskip  # start by skipping
        sink._reset(self)

    def _buffers_ready(self,inputs,outputs):
        # handle the case where we take multiple samples from a block....
        ioffset = 0
        while (True):
            #print "    %s: skip=%d, record=%d" % (self._name,self._skip,self._record)
            if self._skip >= self._buffer_size - ioffset:
                # we're done with this block
                self._skip -= self._buffer_size - ioffset
                return True

            # finishing skipping and start sampling
            ioffset += self._skip
            self._skip = 0
            self._record = self.sample_size
            self._soffset = 0

            # recording
            count = min(self._buffer_size - ioffset,self._record)
            self._temp[self._soffset:self._soffset+count] = inputs[0][ioffset:ioffset+count]
            self._record -= count
            if self._record == 0:
                #print "    %s: done recording" % self._name
                ioffset += count
                # finished taking the sample, save it away and start skipping.
                # we're running locked, so we can just do this without worrying
                # about atomicity
                self._sample[:] = self._temp
                if self.queue:
                    try:
                        self.queue.put(self._sample.copy(),False)
                    except Queue.Full:   # don't care if sample gets dropped
                        pass
                if self.wxclient:
                    # send client a notification that there's data to look at.
                    # they can call get_sample(ndarray) to retrieve the samples
                    #print "    %s: notifying wxclient" % self._name
                    wx.CallAfter(self.wxclient.ProcessSample,self)
                self._skip = self._nskip
            else:
                self._soffset += count
                return True

################################################################################
# audio_sink
#       play samples using pyaudio.  Expects float values between -1.0 and +1.0
################################################################################

audio = None
try:
    import ossaudiodev
    audio = 'oss'
except ImportError:
    try:
        import pyaudio
        audio = 'pyaudio'
    except ImportError:
        pass

class audio_sink(sink):
    _default_parameters = {
        'gain': 1.0,
        }
    _expected_number_of_inputs = (1,2)

    def _initialize(self):
        self._stream = None
        if audio == 'oss':
            pass
        elif audio == 'pyaudio':
            self._pyaudio = pyaudio.PyAudio()
            devinfo = self._pyaudio.get_default_output_device_info()
            assert self._ninputs <= devinfo['maxOutputChannels'],\
                   "%s: maximum number of output channels is %s" % \
                   (self._name,devinfo['maxOutputChannels'])
            assert self._pyaudio.is_format_supported(self.sample_rate,
                                                     output_device = devinfo['index'],
                                                     output_channels = self._ninputs,
                                                     output_format = pyaudio.paInt16),\
                     "%s: %g is not a supported sample rate" % (s._name,self.sample_rate)
        else:
            assert False,"Sorry, no audio device library could be found"

        # buffer for scaling of input samples
        self._temp1 = numpy.empty(self._buffer_size,numpy.float32)
        # buffer for assembling int16 output stream
        self._temp2 = numpy.empty(self._buffer_size*self._ninputs,numpy.int16)

    def set_gain(self,gain):
        self._mutex.acquire()
        self.gain = gain
        self._mutex.release()

    def _reset(self):
        if audio == 'oss':
            self._stream = ossaudiodev.open('w')
            (fmt,chans,rate) = self._stream.setparameters(ossaudiodev.AFMT_S16_LE,self._ninputs,self.sample_rate)
            assert rate == self.sample_rate,\
                   "%s: %g is not a supported sample rate" % (s._name,self.sample_rate)
        elif audio == 'pyaudio':
            if self._stream: self._stream.close()
            self._stream = self._pyaudio.open(int(self.sample_rate),self._ninputs,pyaudio.paInt16,
                                              output=True)
        self._offset = -1   # >= 0 if we're emptying internal buffer to audio output
        sink._reset(self)

    def _can_send(self):
        if audio == 'oss':
            return self._stream.obuffree()
        elif audio == 'pyaudio':
            return self._stream.get_write_available()
        else:
            # doesn't matter, we're throwing the samples away!
            return 32768

    def _audio_write(self,data,nframes):
        if audio == 'oss':
            self._stream.write(data)
        elif audio == 'pyaudio':
            self._stream.write(data,num_frames = nframes)
        else:
            # dump samples into the bit bucket
            pass

    def _buffers_ready(self,inputs,outputs):
        if self._offset < 0:
            # fill internal buffer with interleaved samples
            for i in xrange(self._ninputs):
                numpy.multiply(inputs[i],self.gain*32767.0,self._temp1)
                self._temp2[i::self._ninputs] = self._temp1
            self._offset = 0

        # see how many frames we can write without blocking
        nframes = (len(self._temp2) - self._offset)/self._ninputs
        can_send = min(self._can_send(),nframes)
        if can_send > 0:
            end = self._offset + can_send*self._ninputs
            self._audio_write(self._temp2[self._offset:end].tostring(),can_send)
            # update offset into internal buffer
            self._offset = end

        if self._offset == len(self._temp2):
            self._offset = -1    # we need another input buffer
            return True
        else:
            return False  # still working on this input buffer

################################################################################
# usrp_sink
#       send samples to USRP.  Expects float values between -1.0 and +1.0
#       accepts one (I only) or two (I and Q) input streams
################################################################################

try:
    import usrp

    class usrp_sink(sink):
        _default_parameters = {
            'gain': 1.0,
            'board': None,          # None: use board 0, string: match serial number
            'dbid':  usrp.USRP_DBID_LF_TX,   # which daughterboard to use
            'pga_gain_db': 0,       # gain of programmable gain amplifier (-20..0)
            'frequency': 0,         # set freq for digital up converter (-44e6..44e6)
            }
        _expected_number_of_inputs = (1,2)

        NBOARDS = 4   # how many boards to look for

        def _initialize(self):
            for board in xrange(self.NBOARDS):
                if not usrp.usrp_find_device(board): continue
                self._tx = usrp.usrp_tx(500,     # interp rate updated below
                                        nchan = 1,
                                        board = board)
                if self.board is None or self.board.lower() == self._tx.serial_number.lower():
                    self.board = self._tx.serial_number
                    break
                del self._tx
                self._tx = None
            else:
                assert False,"%s: can't find USRP board %s" % (self._name,self.board)

            # set correct interpolation rate
            interp = self._tx.dac_rate/self.sample_rate
            assert interp <= 512,\
                   "%s: sorry, sample rate must be at least %g" % (self._name,self._tx.dac_rate/512)
            self._tx.set_interp_rate(interp)
            self.interp_rate = self._tx.interp_rate

            # send data to correct daughterboard
            if self._tx.daughterboard_id(0) == self.dbid:
                self._chan = 0
                self._tx.set_mux(0x0098)  # send chan 0 to DACs 0 and 1
            elif self._tx.daughterboard_id(1) == self.dbid:
                self._chan = 1
                self._tx.set_mux(0x9800)  # send chan 0 to DACs 2 and 3
            else:
                assert False,"%s: can't find daughterboard with id 0x%x" % (self._name,self.dbid)

            # set PGA gain
            assert self.pga_gain_db >= self._tx.pga_min and self.pga_gain_db <= self._tx.pga_max,\
                   "%s: pga_gain_db must be in the range %g to %g" % \
                   (self._name,self._tx.pga_min,self._tx.pga_max)
            # DACs 0 and 1 share a gain setting, as do DACs 2 and 3
            self._tx.set_pga(self._chan << 1,self.pga_gain_db)

            # set frequency of DUC
            assert self._tx.set_tx_freq(self._chan,self.frequency),\
                   "%s: set_tx_freq failed with value %g" % (self._name,self.frequency)

            # buffer for scaling of output samples
            self._temp1 = numpy.empty(self._buffer_size,numpy.float32)
            # buffer for assembling int16 output stream (interleaved I & Q samples)
            # start with all zeros in case Q channel isn't connected
            self._temp2 = numpy.zeros(self._buffer_size*2,numpy.int16)

        def set_gain(self,gain):
            self._mutex.acquire()
            self.gain = gain
            self._mutex.release()

        def set_frequency(self,freq):
            self._mutex.acquire()
            self.frequency = freq
            assert self._tx.set_tx_freq(self._chan,freq),\
                   "%s: set_tx_freq failed with value %g" % (self._name,freq)
            self._mutex.release()

        def _reset(self):
            self._tx.stop()
            assert self._tx.start(),"%s: usrp failed to start" % self._name
            self._offset = -1   # >= 0 if we're emptying internal buffer to audio output
            sink._reset(self)

        def _buffers_ready(self,inputs,outputs):
            if self._offset < 0:
                # fill internal buffer with interleaved samples
                for i in xrange(self._ninputs):
                    numpy.multiply(inputs[i],self.gain*32767.0,self._temp1)
                    self._temp2[i::2] = self._temp1
                self._offset = 0

            # try sending what we have
            data = self._temp2[self._offset:].tostring()
            bytes_sent,underrun = self._tx.write(data)
            self._offset += bytes_sent/2

            if self._offset == len(self._temp2):
                self._offset = -1    # call us again when there are more inputs
                return True
            else:
                return False  # more to send, call us again soon

except ImportError:
    class usrp_sink(sink):
        def __init__(self,*inputs,**parameters):
            raise NotImplementedError,"sorry, no usrp library is available"

################################################################################
# wavfile_sink
#   expects inputs to be in the range (-1,+1) after scale_factor is applied
################################################################################

import wave
class wavfile_sink(sink):
    _default_parameters = {
        'filename': None,
        'sample_width': 2,    # sample width in bytes (1, 2, or 4)
        'scale_factor': 1.0,  # each sample value multiplied by this
        }
    _expected_number_of_inputs = (1,16)
    _output_ports = (())

    def _initialize(self):
        assert self.filename,"%s: filename must be specified" % self._name
        if self.sample_width == 1:
            dtype = numpy.uint8
            self._scale = 127.0 * self.scale_factor
            self._offset = 127.0
        elif self.sample_width == 2:
            dtype = numpy.int16
            self._scale = 32767.0 * self.scale_factor
            self._offset = 0
        elif self.sample_width == 4:
            dtype = numpy.int32
            self._scale = 2147483647.0 * self.scale_factor
            self._offset = 0
        else:
            assert False,"%s: sample_width must be 1, 2, or 4 bytes" % self._name
        # array to hold scaled data for 1 channel
        self._temp = numpy.empty(self._buffer_size,dtype=numpy.float64)
        # array to hold frame data all channels
        self._data = numpy.empty(self._ninputs * self._buffer_size,dtype=dtype)

    def _reset(self):
        self._wav = wave.open(self.filename,'wb')
        self._wav.setnchannels(self._ninputs)
        self._wav.setsampwidth(self.sample_width)
        self._wav.setframerate(int(self.sample_rate))
        sink._reset(self)

    def _buffers_ready(self,inputs,outputs):
        for i in xrange(self._ninputs):
            # apply appropriate scale and offset
            numpy.multiply(inputs[i],self._scale,self._temp)
            if self._offset != 0: numpy.add(self._temp,self._offset,self._temp)
            # interleave channel samples in the output array
            self._data[i::self._ninputs] = self._temp[:]
        # send frames to wav file
        self._wav.writeframes(self._data.tostring())
        return True

################################################################################
# print_sink
#       just calls print on samples in the input buffer
################################################################################

class print_sink(sink):
    _default_parameters = {}
    _expected_number_of_inputs = (1,1)

    def _buffers_ready(self,inputs,outputs):
        print self
        print inputs[0]
        return True

################################################################################
# rate_sink
#       computes the rate at which samples are delivered, printing result
#       at specified interval
################################################################################

import time
class rate_sink(sink):
    _default_parameters = {
        'interval': 10.0,
        'label': None,
        }
    _expected_number_of_inputs = (1,1)

    def _reset(self):
        self.count = 0
        self.then = time.time()
        sink._reset(self)

    def _buffers_ready(self,inputs,outputs):
        self.count += self._buffer_size
        now = time.time()
        interval = now - self.then
        if interval >= self.interval:
            # print out stats every interval seconds
            self._report(self.label or self,self.count,interval)
            self.then = now
            self.count = 0
        return True

    # override for a different printout...
    def _report(self,label,count,interval):
        print "%s: %3.0f samples/sec" % (label,count/interval)
        

__all__ = [
    'null_sink',
    'sample_sink',
    'audio_sink',
    'usrp_sink',
    'wavfile_sink',
    'print_sink',
    'rate_sink',
    ]
