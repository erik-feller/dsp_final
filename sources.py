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
import numpy,os,wave

################################################################################
# sin_source: sinusoid on port 0, quadrature (+90 degrees) on port 1
################################################################################

class sin_source(dsp_block):
    _default_parameters = {
        'frequency': None,
        'sample_rate': None,
        'phase_offset': 0.0,
        'amplitude': 1.0,
        'offset': 0.0
        }
    _expected_number_of_inputs = (0,0)
    _output_ports = (numpy.float64,numpy.float64)

    def _initialize(self):
        assert isinstance(self.frequency,(int,float)),\
               "%s: frequency must be specified" % self._name
        assert self.frequency < self.sample_rate/2.0,\
               "%s: frequency (%g) must be <= half of sample rate (%g)" % \
               (self._name,self.frequency,self.sample_rate)
        self.set_frequency(self.frequency)

    def _reset(self):
        self._phase = self.phase_offset
        dsp_block._reset(self)

    def set_frequency(self,frequency):
        assert frequency <= self.sample_rate/2,\
               "%s attempt to set frequency (%g) that's greater than the Nyquist rate (%g)" % \
               (self._name,frequency,self.sample_rate/2)
        self._mutex.acquire()
        self.frequency = frequency
        # radians/sec / (samples/sec) = radians/sample
        phase_step = 2*numpy.pi*self.frequency/self.sample_rate
        self._phase_inc = self._buffer_size * phase_step
        self._phase_angles = numpy.arange(0,self._phase_inc,phase_step)
        self._mutex.release()

    def _buffers_ready(self,inputs,outputs):
        # compute sin of next batch of phase angles
        numpy.add(self._phase_angles,self._phase,outputs[0])
        numpy.sin(outputs[0],outputs[0])
        if self.amplitude != 1.0:
            numpy.multiply(outputs[0],self.amplitude,outputs[0])
        if self.offset != 0.0:
            numpy.add(outputs[0],self.offset,outputs[0])

        # output cos if anyone is listening
        if self._output_connected[1]:
            numpy.add(self._phase_angles,self._phase,outputs[1])
            numpy.cos(outputs[1],outputs[1])
            if self.amplitude != 1.0:
                numpy.multiply(outputs[1],self.amplitude,outputs[1])
            if self.offset != 0.0:
                numpy.add(outputs[1],self.amplitude,outputs[1])

        # update phase accumulator
        self._phase = numpy.fmod(self._phase + self._phase_inc,2*numpy.pi)

        # finally, propagate the new output values and report we've used inputs
        self._propagate_outputs()
        return True

################################################################################
# wavefile_source
#       produce stream of float64 values from samples in WAV file
#       left channel on port 0, right channel on port 1
#       mono signals appear on both outputs
################################################################################

class wavfile_source(dsp_block):
    _default_parameters = {
        'filename': None,
        'amplitude': 1.0,
        'repeat': True,
        }
    _expected_number_of_inputs = (0,0)
    _output_ports = (())   # will be adjusted in initialize

    def _initialize(self):
        assert self.filename,"%s: filename must be specified" % self._name
        assert os.path.exists(self.filename),\
               "%s: file %s doesn't exist" % (self._name,self.filename)
        self._wav = wave.open(self.filename,'rb')
        self._total_frames = self._wav.getnframes()
        self._nchan = self._wav.getnchannels()
        assert self._total_frames > 0,\
               "%s: %s doesn't have any audio data!" % (self._name,self.filename)
        self.sample_rate = self._wav.getframerate()
        sample_width = self._wav.getsampwidth()

        # set up an output port for each channel
        self._setup_outputs(output_ports=[numpy.float64,] * self._nchan)

        # see http://ccrma.stanford.edu/courses/422/projects/WaveFormat/
        if sample_width == 1:
            # data is unsigned bytes, 0 to 255
            self._dtype = numpy.uint8
            self._scale = self.amplitude / 127.0
            self._offset = -1.0
        elif sample_width == 2:
            # data is signed 2's complement 16-bit samples (little-endian byte order)
            self._dtype = numpy.int16
            self._scale = self.amplitude / 32767.0
            self._offset = 0.0
        elif sample_width == 4:
            # data is signed 2's complement 16-bit samples (little-endian byte order)
            self._dtype = numpy.int32
            self._scale = self.amplitude / 2147483647.0
            self._offset = 0.0
        else:
            assert False,"%s: unrecognized sample width %d" % (self._name,sample_width)

    def _reset(self):
        self._wav.rewind()
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        count = 0
        while count < self._buffer_size:
            # read data into numpy array
            # we receive data in frames, each frame has nchannels*samplewidth bytes
            want = (self._buffer_size - count)
            frames = self._wav.readframes(want)
            audio = numpy.frombuffer(frames,dtype=self._dtype)
            got = len(audio)
            if got < want:
                if self.repeat:
                    self._wav.rewind()
                else:
                    # fill remainder of output buffer with zeros
                    temp = numpy.zeros(want,dtype=self._dtype)
                    temp[0:got] = audio[:]
                    audio = temp
                    got = want
                    # indicate we've reached the end-of-stream
                    self._eos = True

            for i in xrange(self._nchan):
                outputs[i][count:count+got] = audio[i::self._nchan]

            count += got

        # scale data appropriately
        for i in xrange(self._nchan):
            numpy.multiply(outputs[i],self._scale,outputs[i])
            if self._offset != 0: numpy.add(outputs[i],self._offset,outputs[i])

        self._propagate_outputs()
        return True

# like wavfile_source except you can specify the desired sample rate (must be
# correspond to an integer decimation or interpolation factor of the WAV file
# sample rate).  Additionally you can band-limit the output (must be <= the
# lesser of the WAV nyquist rate and the sample nyquist rate).
import blocks
from filters import fir_low_pass,decimate,interpolate
def resampled_wavfile_source(filename=None,amplitude=1.0,sample_rate=None,bw=None,ntaps=100,repeat=True):
    # start by creating the wavfile_source block
    wav = wavfile_source(filename=filename,repeat=repeat)

    # use default sample rate if none specified
    if sample_rate is None:
        sample_rate = blocks.current_testbench.sample_rate

    # if the wavfile_source meets the specs, we're done!
    if wav.sample_rate == sample_rate and bw == wav.sample_rate/2:
        return wav

    if wav.sample_rate >= sample_rate:
        # decimation required, start with low-pass filter to avoid aliasing
        cutoff = sample_rate/2
        if bw: cutoff = min(bw,cutoff)

        fdecimation = float(wav.sample_rate)/sample_rate
        decimation = int(fdecimation)
        assert sample_rate * decimation == wav.sample_rate,\
               "%s: required decimation (%g) must be an integer" % (self._name,fdecimation)

        filter = fir_low_pass(wav,gain=amplitude,
                              cutoff_frequency=cutoff,ntaps=ntaps)

        # now add the decimator, if needed
        if decimation == 1: return filter
        else: return decimate(filter,decimation=decimation)
    else:
        # interpolation required
        cutoff = wav.sample_rate/2
        if bw: cutoff = min(bw,cutoff)

        finterpolation = float(sample_rate)/wav.sample_rate
        interpolation = int(finterpolation)
        assert wav.sample_rate * interpolation == sample_rate,\
               "%s: required interpolation (%g) must be an integer" % (self._name,finterpolation)

        interpolator = interpolate(wav,interpolation=interpolation)
        return fir_low_pass(interpolator,gain=interpolation*amplitude,
                            cutoff_frequency=cutoff,ntaps=ntaps)

################################################################################
# usrp_source
#       get samples from USRP.  Delivers float values between -1.0 and +1.0
#       has two output ports (I and Q)
################################################################################

try:
    import usrp

    class usrp_source(dsp_block):
        _default_parameters = {
            'gain': 1.0,
            'board': None,          # None: use board 0, string: match serial number
            'dbid':  usrp.USRP_DBID_LF_RX,   # which daughterboard to use
            'pga_gain_db': 0,       # gain of programmable gain amplifier (-20..0)
            'frequency': 0,         # set freq for digital down converter
            }
        _expected_number_of_inputs = (0,0)
        _output_ports = (numpy.float64,numpy.float64)

        NBOARDS = 4   # how many boards to look for

        def _initialize(self):
            for board in xrange(self.NBOARDS):
                if not usrp.usrp_find_device(board): continue
                self._rx = usrp.usrp_rx(250,     # decim rate updated below
                                        nchan = 1,
                                        board = board)
                if self.board is None or self.board.lower() == self._rx.serial_number.lower():
                    self.board = self._rx.serial_number
                    break
                del self._rx
                self._rx = None
            else:
                assert False,"%s: can't find USRP board %s" % (self._name,self.board)

            # set correct decim rate
            decim = self._rx.adc_rate/self.sample_rate
            assert decim <= 256,\
                   "%s: sorry, sample rate must be at least %g" % (self._name,self._rx.adc_rate/256)
            self._rx.set_decim_rate(decim)
            self.decim_rate = self._rx.decim_rate

            # send data to correct daughterboard
            if self._rx.daughterboard_id(0) == self.dbid:
                self._chan = 0
                self._rx.set_mux(0x10101010)  # use first pair of ADCs
            elif self._rx.daughterboard_id(1) == self.dbid:
                self._chan = 1
                self._rx.set_mux(0x32323232)  # use second pait of ADCs
            else:
                assert False,"%s: can't find daughterboard with id 0x%x" % (self._name,self.dbid)

            # set PGA gain
            assert self.pga_gain_db >= self._rx.pga_min and self.pga_gain_db <= self._rx.pga_max,\
                   "%s: pga_gain_db must be in the range %g to %g" % \
                   (self._name,self._rx.pga_min,self._rx.pga_max)
            self._rx.set_pga(2*self._chan,self.pga_gain_db)    # set it for both ADCs
            self._rx.set_pga(2*self._chan+1,self.pga_gain_db)

            # set frequency of DDC (we're only using the first DDC)
            assert self._rx.set_rx_freq(0,self.frequency),\
                   "%s: set_rx_freq failed with value %g" % (self._name,self.frequency)

        def set_gain(self,gain):
            self._mutex.acquire()
            self.gain = gain
            self._mutex.release()

        def set_frequency(self,freq):
            self._mutex.acquire()
            self.frequency = freq
            assert self._rx.set_rx_freq(0,freq),\
                   "%s: set_rx_freq failed with value %g" % (self._name,freq)
            self._mutex.release()

        def _reset(self):
            self._count = 0   # where we are in filling output buffers
            self._rx.stop()
            assert self._rx.start(),"%s: usrp failed to start" % self._name
            dsp_block._reset(self)

        def _buffers_ready(self,inputs,outputs):
            want = self._buffer_size - self._count
            # each sample is 2 bytes of I and 2 bytes of Q
            buffer,overrun = self._rx.read(want * 4)
            data = numpy.frombuffer(buffer,dtype=numpy.int16)
            got = len(data)/2   # contains I and Q interleaved

            if got > 0:
                end = self._count + got
                outputs[0][self._count:end] = data[0::2]  # copy I
                outputs[1][self._count:end] = data[1::2]  # copy Q
                self._count += got
                
                if self._count == self._buffer_size:
                    # we've filled our buffers, send them downstream
                    numpy.multiply(outputs[0],self.gain/32768.0,outputs[0])
                    numpy.multiply(outputs[1],self.gain/32768.0,outputs[1])
                    self._propagate_outputs()
                    self._count = 0

            return True

except ImportError:
    class usrp_source(dsp_block):
        def __init__(self,*inputs,**parameters):
            raise NotImplementedError,"sorry, no usrp library is available"

################################################################################
# noise_source
################################################################################

class noise_source(dsp_block):
    _default_parameters = {
        'sample_rate': None,
        'distribution': 'normal',   # 'normal', 'laplace', 'raleigh', 'uniform', 'triangular', 'impulse'
        'amplitude': 1.0,
        'loc': 0.0,
        'scale': 1.0,
        }
    _expected_number_of_inputs = (0,0)
    _output_ports = (numpy.float64,)

    def _initialize(self):
        self._temp = numpy.empty(self._buffer_size,dtype=bool)

    def _buffers_ready(self,inputs,outputs):
        if self.distribution in ('normal','gaussian'):
            noise = numpy.random.normal(self.loc,self.scale,size=self._buffer_size)
        elif self.distribution == 'uniform':
            noise = numpy.random.uniform(self.loc-self.scale,self.loc+self.scale,size=self._buffer_size)
        elif self.distribution == 'triangular':
            noise = numpy.random.triangular(self.loc-self.scale,self.loc,self.loc+self.scale,size=self._buffer_size)
        elif self.distribution == 'impulse':
            # from gr_random.cc in the gnuradio code
            # scale: 5 => scratchy, 8 => geiger
            noise = numpy.random.uniform(size=self._buffer_size)
            numpy.log(noise,noise)
            numpy.multiply(noise,-1.4142135623730951,noise)   # * -sqrt(2)
            numpy.less_equal(noise,self.scale,self._temp)
            noise[self._temp] = 0.0
        elif self.distribution == 'laplace':
            noise = numpy.random.laplace(self.loc,self.scale,size=self._buffer_size)
        elif self.distribution == 'raleigh':
            noise = numpy.random.raleigh(self.scale,size=self._buffer_size)
        else:
            assert False,"%s: unrecognized distribution %s" % (self._name,self.distribution)
        outputs[0][:] = noise
        del noise
        if self.amplitude != 1.0:
            numpy.multiply(outputs[0],self.amplitude,outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# array_source
################################################################################

import array

# create a data stream from an array.array 
class array_source(dsp_block):
    _default_parameters = {
        'data': None,    # array.array containing data
        'repeat': True,  # repeat data sequence, or just send 0's when done
        }
    _expected_number_of_inputs = (0,0)
    _output_ports = (())   # will be adjusted in _initialize

    def _initialize(self):
        assert isinstance(self.data,array.array),\
               "%s: data must be supplied in an array.array" % self._name
        out_type = {
            'c': numpy.uint8,
            'b': numpy.int8,
            'B': numpy.uint8,
            'u': numpy.uint16,
            'h': numpy.int16,
            'H': numpy.uint16,
            'i': numpy.int32,
            'I': numpy.uint32,
            'l': numpy.int64,
            'L': numpy.uint64,
            'f': numpy.float32,
            'd': numpy.float64,
            }.get(self.data.typecode,None)
        assert out_type is not None,\
               "%s: unrecognized array type %s" % (self._name,self.data.typecode)
        self._setup_outputs(output_ports=(out_type,))

        # for efficiency, make sure data array is long enough to generate a buffer's
        # worth of samples.
        if len(self.data) < self._buffer_size:
            temp = array.array(self.data.typecode,self.data)  # copy of data
            while len(self.data) < self._buffer_size:
                self.data.extend(temp)

        # an array of zeros to use when we run out of data
        if not self.repeat:
            self._zeros =  array.array(self.data.typecode)
            for k in xrange(1024): self._zeros.append(0)

    def _reset(self):
        self._doffset = 0    # start at beginning of data
        self._dlen = len(self.data)
        self._dzeros = False
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        # produce another output buffer...
        # start by restoring state into local vars for speed
        data = self.data
        dlen = self._dlen
        doffset = self._doffset

        index = 0   # next word to be filled in output buffer
        while index < self._buffer_size:
            count = min(self._buffer_size - index,dlen - doffset)
            end_in = doffset + count
            end_out = index + count
            outputs[0][index:end_out] = data[doffset:end_in]
            index = end_out
            doffset = end_in
            
            if doffset >= dlen:
                # we've reached the end of the array...
                if not self.repeat and not self._dzeros:
                    # fill remaining output buffer locations with zeros
                    data = self._zeros; self._data = data
                    dlen = len(data); self._dlen = dlen
                    self._dzeros = True
                    # indicate that we've reached the end of stream
                    self._eos = True
                doffset = 0

        # remember state for next time
        self._doffset = doffset

        self._propagate_outputs()
        return True

# build data stream from a file
def file_source(filename=None,repeat=True):
    f = open(filename,'rb')
    data = f.read()
    f.close()
    return array_source(data=array.array('B',data),repeat=repeat)

################################################################################
# queue_source
################################################################################

# create a data stream from a queue of sequence objects
class queue_source(dsp_block):
    _default_parameters = {
        'producer': None,  # instance of blocks.queue_block with _produce() method
        }
    _expected_number_of_inputs = (0,0)
    _output_ports = (())   # will be adjusted in _initialize

    def _initialize(self):
        self._setup_outputs(output_ports=(self.producer._dtype(),))
        # array to supply values when producer has none
        self._zeros = numpy.zeros(self._buffer_size,self.producer._dtype())

    def _reset(self):
        self._in = 0     # where we are in input sequence
        self._seq = None # current input sequence
        self._out = 0    # where we are in output buffer
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        while self._out < self._buffer_size:
            # if no input sequence, grab another from the queue
            if self._seq is None:
                self._seq = self.producer._get()
                # nothing in queue, finish off current buffer with zeros
                if self._seq is None:
                    self._seq = self._zeros
                    self._in = self._out
                    if self.producer._eos:
                        self._eos = True
                else:
                    self._in = 0

            # so how much of input sequence we can move to output buffer
            count = min(len(self._seq) - self._in,self._buffer_size - self._out)
            end_out = self._out + count
            end_in = self._in + count
            outputs[0][self._out:end_out] = self._seq[self._in:end_in]
            self._out = end_out
            self._in = end_in
            if self._in == len(self._seq): self._seq = None
        # we've filled up the output buffer, send it to fanouts...
        self._propagate_outputs()
        return True

__all__ = [
    'sin_source',
    'wavfile_source',
    'resampled_wavfile_source',
    'usrp_source',
    'noise_source',
    'array_source',
    'file_source',
    'queue_source',
    ]
