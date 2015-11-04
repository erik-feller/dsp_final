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
# convolve: basic convolution of stream with supplied coefficients
################################################################################

class convolve(dsp_block):
    _default_parameters = {
        'taps': None,
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = (numpy.float64,)

    def _initialize(self):
        assert self.taps,"%s: taps must be specified" % self._name
        self.ntaps = len(self.taps)

        # construct a numpy array for use during convolution
        self._taps = numpy.array(self.taps,dtype=numpy.float64)

        # to avoid problems as we move from block to block, set up an internal
        # buffer that holds the last _overlap samples from the previous block so
        # that the convolution operation will produce exactly _buffer_size valid outputs
        # when operating on the current input block
        self._overlap = self.ntaps - 1
        self._temp = numpy.empty(self._buffer_size+self._overlap,dtype=numpy.float64)

    def _reset(self):
        self._temp[0:self._overlap] = 0.0   # initialize our "memory"
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        # add incoming samples after the ones we saved from last buffer
        self._temp[self._overlap:] = inputs[0]

        # perform the convolution with the taps, send result to fanouts
        conv = numpy.convolve(self._temp,self._taps,mode='valid')
        outputs[0][:] = conv
        del conv
        self._propagate_outputs()

        # now remember the last _overlap samples for next time
        self._temp[0:self._overlap] = self._temp[self._buffer_size:]

        # mark inputs as consumed
        return True

################################################################################
# fir_filter: basic convolution with support for calculating the taps
################################################################################

class fir_filter(convolve):
    def _initialize(self):
        # set up taps
        if self.__dict__.has_key('taps') and self.taps:
            self.ntaps = len(self.taps)
        else:
            # self.ntaps might be set by a subclass
            if self.ntaps:
                self.ntaps = (self.ntaps & ~0x1) + 1   # make it odd
            else:
                self.ntaps = self._compute_ntaps()
            window = compute_window(self.window,self.ntaps)
            self._compute_tap_coefficients(window)
        convolve._initialize(self)

    # figure out filter taps
    def _compute_tap_coefficients(self):
        assert False,\
               "%s: _compute_tap_coefficients must be define if taps not specified" % self._name

    # compute number of taps given sample_rate and transition_width.
    # Stolen from the gnuradio firdes routines
    _width_factor = {
        'hamming': 3.3,
        'hann': 3.1,
        'blackman': 5.5,
        'rectangular': 2.0,
        }
    def _compute_ntaps(self):
        assert self.transition_width,\
               "%s: transition_width must be specified if ntaps is not" % self._name
        delta_f = float(self.transition_width)/self.sample_rate
        width_factor = self._width_factor.get(self.window,None)
        assert width_factor,\
               "%s: unrecognized window type %s" % (self._name,self.window)
        ntaps = int(width_factor/delta_f + 0.5)
        return (ntaps & ~0x1) + 1   # ensure it's odd

# compute specified window given number of taps
# formulae from Wikipedia
def compute_window(window,ntaps):
    order = float(ntaps - 1)
    if window == 'hamming':
        return [0.53836 - 0.46164*numpy.cos((2*numpy.pi*i)/order)
                for i in xrange(ntaps)]
    elif window == 'hann' or window == 'hanning':
        return [0.5 - 0.5*numpy.cos((2*numpy.pi*i)/order)
                for i in xrange(ntaps)]
    elif window == 'bartlett':
        return [1.0 - abs(2*i/order - 1)
                for i in xrange(ntaps)]
    elif window == 'blackman':
        alpha = .16
        return [(1-alpha)/2 - 0.50*numpy.cos((2*numpy.pi*i)/order)
                - (alpha/2)*numpy.cos((4*numpy.pi*i)/order)
                for i in xrange(ntaps)]
    elif window == 'nuttall':
        return [0.355768 - 0.487396*numpy.cos(2*numpy.pi*i/order)
                         + 0.144232*numpy.cos(4*numpy.pi*i/order)
                         - 0.012604*numpy.cos(6*numpy.py*i/order)
                for i in xrange(ntaps)]
    elif window == 'blackman-harris':
        return [0.35875 - 0.48829*numpy.cos(2*numpy.pi*i/order)
                        + 0.14128*numpy.cos(4*numpy.pi*i/order)
                        - 0.01168*numpy.cos(6*numpy.pi*i/order)
                for i in xrange(ntaps)]
    elif window == 'blackman-nuttall':
        return [0.3635819 - 0.4891775*numpy.cos(2*numpy.pi*i/order)
                          + 0.1365995*numpy.cos(4*numpy.pi*i/order)
                          - 0.0106411*numpy.cos(6*numpy.py*i/order)
                for i in xrange(ntaps)]
    elif window == 'flat top':
        return [1 - 1.93*numpy.cos(2*numpy.pi*i/order)
                  + 1.29*numpy.cos(4*numpy.pi*i/order)
                  - 0.388*numpy.cos(6*numpy.py*i/order)
                  + 0.032*numpy.cos(8*numpy.py*i/order)
                for i in xrange(ntaps)]
    elif window == 'rectangular' or window == 'dirichlet':
        return [1 for i in xrange(ntaps)]
    else:
        assert False,"compute_window: unrecognized window type %s" % window

################################################################################
# fir_low_pass: fir_filter with low-pass taps
################################################################################

class fir_low_pass(fir_filter):
    _default_parameters = {
        'gain': 1.0,
        'ntaps': None,
        'transition_width': None,
        'cutoff_frequency': None,
        'window': 'hamming',
        }

    # Stolen from the gnuradio firdes routines
    def _compute_tap_coefficients(self,window):
        assert self.cutoff_frequency,\
               "%s: cutoff_frequency must be specified" % self._name
        fc = float(self.cutoff_frequency) / self.sample_rate
        wc = 2 * numpy.pi * fc
        middle = (self.ntaps - 1)/2
        self.taps = [0] * self.ntaps
        fmax = 0  # for low pass, gain @ DC = 1.0
        for i in xrange(self.ntaps):
            if i == middle:
                coeff = (wc/numpy.pi) * window[i]
                fmax += coeff
            else:
                n = i - middle
                coeff = (numpy.sin(n*wc)/(n*numpy.pi)) * window[i]
                fmax += coeff
            self.taps[i] = coeff
        gain = self.gain / fmax
        for i in xrange(self.ntaps):
            self.taps[i] *= gain

################################################################################
# fir_high_pass: fir_filter with high-pass taps
################################################################################

class fir_high_pass(fir_filter):
    _default_properties = {
        'gain': 1.0,
        'ntaps': None,
        'transition_width': None,
        'cutoff_frequency': None,
        'window': 'hamming',
        }

    # Stolen from the gnuradio firdes routines
    def _compute_tap_coefficients(self,window):
        assert self.cutoff_frequency,\
               "%s: cutoff_frequency must be specified" % self._name
        fc = float(self.cutoff_frequency) / self.sample_rate
        wc = 2 * numpy.pi * fc
        middle = (self.ntaps - 1)/2
        self.taps = [0] * self.ntaps
        fmax = 0  # for high pass gain @ nyquist freq = 1.0
        for i in xrange(self.ntaps):
            if i == middle:
                coeff = (1.0 - wc/numpy.pi) * window[i]
                fmax += coeff
            else:
                n = i - middle
                coeff = (-numpy.sin(n*wc)/(n*numpy.pi)) * window[i]
                fmax += coeff * numpy.cos(n*numpy.pi)
            self.taps[i] = coeff
        gain = self.gain / fmax
        for i in xrange(self.ntaps):
            self.taps[i] *= gain

################################################################################
# fir_band_pass: fir_filter with band-pass taps
################################################################################

class fir_band_pass(fir_filter):
    _default_parameters = {
        'gain': 1.0,
        'ntaps': None,
        'transition_width': None,
        'low_cutoff_frequency': None,
        'high_cutoff_frequency': None,
        'window': 'hamming',
        }

    def _compute_tap_coefficients(self,window):
        assert self.low_cutoff_frequency,\
               "%s: low_cutoff_frequency must be specified" % self._name
        assert self.high_cutoff_frequency,\
               "%s: high_cutoff_frequency must be specified" % self._name
        fc_lo = float(self.low_cutoff_frequency) / self.sample_rate
        wc_lo = 2 * numpy.pi * fc_lo
        fc_hi = float(self.high_cutoff_frequency) / self.sample_rate
        wc_hi = 2 * numpy.pi * fc_hi
        middle = (self.ntaps - 1)/2
        self.taps = [0] * self.ntaps
        fmax = 0  # for band pass gain @ (fc_lo + fc_hi)/2 = 1.0
        # a band pass filter is simply the combination of
        #   a high-pass filter at fc_lo  in series with
        #   a low-pass filter at fc_hi
        # so convolve taps to get the effect of composition in series
        for i in xrange(self.ntaps):
            if i == middle:
                coeff = ((wc_hi - wc_lo)/numpy.pi) * window[i]
                fmax += coeff
            else:
                n = i - middle
                coeff = ((numpy.sin(n*wc_hi) - numpy.sin(n*wc_lo))/(n*numpy.pi)) * window[i]
                fmax += coeff * numpy.cos(n*(wc_lo + wc_hi)*0.5)
            self.taps[i] = coeff
        gain = self.gain / fmax
        for i in xrange(self.ntaps):
            self.taps[i] *= gain

################################################################################
# fir_band_reject: fir_filter with band-reject taps
################################################################################

class fir_band_reject(fir_filter):
    _default_parameters = {
        'gain': 1.0,
        'ntaps': None,
        'transition_width': None,
        'low_cutoff_frequency': None,
        'high_cutoff_frequency': None,
        'window': 'hamming',
        }

    def _compute_tap_coefficients(self,window):
        assert self.low_cutoff_frequency,\
               "%s: low_cutoff_frequency must be specified" % self._name
        assert self.high_cutoff_frequency,\
               "%s: high_cutoff_frequency must be specified" % self._name
        fc_lo = float(self.low_cutoff_frequency) / self.sample_rate
        wc_lo = 2 * numpy.pi * fc_lo
        fc_hi = float(self.high_cutoff_frequency) / self.sample_rate
        wc_hi = 2 * numpy.pi * fc_hi
        middle = (self.ntaps - 1)/2
        self.taps = [0] * self.ntaps
        fmax = 0  # for band reject gain @ DC = 1.0
        # a band reject filter is simply the combination of
        #   a low-pass filter at fc_lo   in series with a
        #   a high-pass filter at fc_hi
        # so convolve taps to get the effect of composition in series
        for i in xrange(self.ntaps):
            if i == middle:
                coeff = (1.0 - ((wc_hi - wc_lo)/numpy.pi)) * window[i]
                fmax += coeff
            else:
                n = i - middle
                coeff = ((numpy.sin(n*wc_lo) - numpy.sin(n*wc_hi))/(n*numpy.pi)) * window[i]
                fmax += coeff
            self.taps[i] = coeff
        gain = self.gain / fmax
        for i in xrange(self.ntaps):
            self.taps[i] *= gain

################################################################################
# decimate: select every Nth value
################################################################################

# This block is usually preceded by a low-pass filter with a cutoff at the nyquist
# rate of the output (to avoid aliasing problems).

class decimate(dsp_block):
    _default_parameters = {
        'decimation': None,
        'offset': 0,     # where to start collecting
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _initialize(self):
        assert isinstance(self.decimation,int) and self.decimation > 0,\
               "%s: decimation factor must be a positive integer" % self._name
        assert self.decimation < self._buffer_size,\
               "%s: Sorry, decimation factor must be less than block size" % self._name
        assert self.offset < self.decimation,\
               "%s: Sorry, offset must be less than decimation factor" % self._name
        self.sample_rate /= self.decimation

        # we'll need an array in which to collect output samples.  In processing
        # each input block we'll collect X samples we want to keep.  So depending
        # on X and _buffer_size we might need up to X extra locations (maybe that's
        # X-1 locations, but it's okay if we have extra).
        X = int(numpy.ceil(self._buffer_size / self.decimation))
        self._temp = numpy.empty(self._buffer_size + X,dtype=self._input_types[0])

    def _reset(self):
        # decimation factor may not be a multiple of the block size, so we need
        # to keep track of where we are in collecting samples
        self._index = self.offset
        # also keep track of where we are in filling temp buffer
        self._offset = 0
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        # we'll be slicing using [index::decimation], ie, collect samples from
        # index, index+decimation, index+2*decimation, ..., index+n*decimation
        # where n is the largest integer such that index+n*decimation <= buffer_size-1.
        # The number of samples collected is thus n+1.

        # compute number of samples we'll collect from this input buffer
        nsamples = int(float(self._buffer_size - 1 - self._index)/self.decimation) + 1
        
        # collect the samples
        self._temp[self._offset:self._offset+nsamples] = inputs[0][self._index::self.decimation]
        self._offset += nsamples
        self._index = (self._index + self._buffer_size) % self._buffer_size

        # if we've got enough samples, send them off!
        if self._offset >= self._buffer_size:
            outputs[0][:] = self._temp[:self._buffer_size]
            self._propagate_outputs()
            # now move what we didn't send down to the front of our internal buffer
            new_offset = self._offset - self._buffer_size
            if new_offset > 0:
                self._temp[0:new_offset] = self._temp[self._buffer_size:self._offset]
            self._offset = new_offset

        # finished with this input buffer!
        return True

################################################################################
# interpolate: insert N-1 zeros between each value
################################################################################

# This block is usually followed by a low-pass filter with a cutoff at the nyquist
# rate of the input (the fiter will fill in the zeroes with appropriately
# interpolated values).

class interpolate(dsp_block):
    _default_parameters = {
        'interpolation': None,
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _initialize(self):
        assert isinstance(self.interpolation,int) and self.interpolation > 0,\
               "%s: interpolation factor must be a positive integer" % self._name
        self.sample_rate *= self.interpolation
        self._nzeros = 0
        self._offset = 0

    def _reset(self):
        # remember how many extra zeros we need generate at beginning of
        # next output buffer
        self._nzeros = 0
        # also keep track of where we are in processing current input buffer
        self._offset = 0
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        # fill output buffer with zeros
        outputs[0][:] = 0

        # keep outputing buffers of 0 until next input sample is needed
        if self._nzeros >= self._buffer_size:
            self._nzeros -= self._buffer_size
            self._propagate_outputs()
            return False

        # we'll be slicing using [nzeros::interpolation], ie, placing samples at
        # nzeros, nzeros+interpolation, nzeros+2*interpolation, ..., nzeros+n*interpolation
        # where n is the largest integer such that nzeros+n*interpolation <= buffer_size-1.
        # The number of samples placed is thus n+1.
        nsamples = int(float(self._buffer_size - 1 - self._nzeros)/self.interpolation) + 1

        # move input samples to output buffer using a stride of interpolation
        outputs[0][self._nzeros::self.interpolation] = \
          inputs[0][self._offset:(self._offset+nsamples)]
        self._propagate_outputs()

        # figure out how many zeros in the next buffer(s)
        self._nzeros = self._nzeros + nsamples*self.interpolation - self._buffer_size

        # now update what we've consumed from input buffer
        self._offset += nsamples

        assert self._offset <= self._buffer_size,"%s: oops, we overconsumed!" % self._name
        if self._offset == self._buffer_size:
            self._offset = 0
            return True
        else:
            return False

################################################################################
# transmit: convert symbols to a stream of values
################################################################################

class transmit(dsp_block):
    _default_parameters = {
        'table': None,
        'sample_rate': None,    # rate at which to generate values
        'symbol_rate': None,    # how long each symbol lasts, in secs.
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = (numpy.float64,)

    def _initialize(self):
        assert isinstance(self.table,(tuple,list)),\
               "%s: lookup table must be a tuple or list" % self._name
        self._table = numpy.fromiter(self.table,dtype=numpy.float64)
        print 'hi there'
        assert self.symbol_rate > 0,\
               "%s: symbol_rate must be positive" % self._name
        assert self.symbol_rate <= self.sample_rate/2,\
               "%s: symbol_rate must be less than or equal to Nyquist rate" % self._name
        self._duration = int(float(self.sample_rate)/self.symbol_rate)

    def _reset(self):
        self._in = 0   # where we are in input buffer
        self._out = 0   # where we are in output buffer
        self._remain = self._duration     # how many samples of current value remain to be output
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        # keep looping until we consume input buffer
        while self._in < self._buffer_size:
            # see how many output samples we can generate
            count = min(self._buffer_size - self._out,self._remain)
            end = self._out + count
            outputs[0][self._out:end] = self._table[inputs[0][self._in]]
            if self._remain == count:
                # finished with this symbol, move on to next
                self._remain = self._duration
                self._in += 1
            else:
                # ran of out output before we finished, pick it up next time
                self._remain -= count
            self._out = end
            if self._out == self._buffer_size:
                # output this buffer, start a new one next time
                self._out = 0
                self._propagate_outputs()
                break
        # we're here because input was exhausted or output buffer got full
        if self._in == self._buffer_size:
            # input exhausted, now we need a new buffer
            self._in = 0
            return True
        else:
            return False

################################################################################
# slice: convert ranges of values into symbols
################################################################################

# inputs will divided into N+1 buckets according to the N levels supplied.
# the bucket index (0 for min value bucket) is used to look up symbol in table
class slice(dsp_block):
    _default_parameters = {
        'levels': None,     # slicing levels
        'table': None,      # symbols (min value bucket first)
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = (numpy.int32,)

    def _initialize(self):
        assert isinstance(self.levels,(tuple,list)),\
               "%s: levels must be a tuple or list" % self._name
        self._levels = [float(x) for x in self.levels]
        self._levels.sort()   # min to max
        assert isinstance(self.table,(tuple,list)),\
               "%s: lookup table must be a tuple or list" % self._name
        assert len(self.levels)+1 == len(self.table),\
               "%s: lookup table must have one more entry than levels" % self._name
        # _compare identifies elements <= current slicing level
        self._compare = numpy.empty(self._buffer_size,dtype=bool)
        # _mask identifies elements yet to be categorized
        self._mask = numpy.empty(self._buffer_size,dtype=bool)

    def _buffers_ready(self,inputs,outputs):
        self._mask[:] = True   # all elements need processing

        # work our way from min slicing level to max
        for index in xrange(len(self._levels)):
            # find elements <= current slicing level
            numpy.less_equal(inputs[0],self._levels[index],self._compare)
            # mask out those we've already categorized
            numpy.logical_and(self._mask,self._compare,self._compare)
            # set symbol value for outputs in this bucket
            outputs[0][self._compare] = self.table[index]
            # now mark the newly buckets values as processed
            self._mask[self._compare] = False
        # remaining elements are in last bucket
        outputs[0][self._mask] = self.table[-1]

        self._propagate_outputs()
        return True

################################################################################
# repack: repack symbols into different size symbols
################################################################################

class repack(dsp_block):
    _default_parameters = {
        'in_size': None,       # number of bits in each incoming sample
        'out_size': None,      # number of bits consumed for each output symbol
        }
    _expected_inputs = (1,1)
    _output_ports = (numpy.int32,)

    def _initialize(self):
        assert self.in_size >= 1 and self.in_size <= 32,\
               "%s: incoming data samples must contain between 1 and 32 bits each" % self._name
        assert self.out_size >= 1 and self.out_size <= 32,\
               "%s: outgoing data samples must contain between 1 and 32 bits each" % self._name

    def _reset(self):
        self._in = 0     # where we are in input buffer
        self._out = 0    # where we are in output buffer
        self._word = 0   # a one-word unpacking buffer
        self._bits = 0   # number of bits in our one-word buffer
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        iindex = self._in
        oindex = self._out
        word = self._word
        bits = self._bits

        ssize = self.out_size
        mask = (1 << ssize) - 1
        dinc = self.in_size
        len = self._buffer_size

        # repack data ssize bits at a time, a process that is slightly less
        # trivial when dinc isn't an integer multiple of ssize
        while oindex < len:
            while ssize > bits:
                if iindex == len: break
                # add bits to our one-word buffer until we have enough
                # for another output chunk
                word |= int(inputs[0][iindex]) << bits
                bits += dinc
                iindex += 1
            outputs[0][oindex] = word & mask
            oindex += 1
            word >>= ssize
            bits -= ssize

        self._word = word
        self._bits = bits
        if oindex == len:
            self._out = 0
            self._propagate_outputs()
        else:
            self._out = oindex
        if iindex == len:
            self._in = 0
            return True
        else:
            self._in = iindex
            return False

################################################################################
# lookup: convert value via table lookup
################################################################################

import operator
class lookup(dsp_block):
    _default_parameters = {
        'table': None,
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = (())   # will be adjusted in _initialize

    def _initialize(self):
        assert isinstance(self.table,(tuple,list)),\
               "%s: lookup table must be a tuple or list" % self._name

        # see what type of value we should produce
        all_values = reduce(operator.add,self.table)
        if isinstance(all_values,int) : dtype = numpy.int32
        else: dtype = numpy.float64
        self._setup_outputs(output_ports=(dtype,))
        self._table = numpy.fromiter(self.table,dtype=dtype)

    def _buffers_ready(self,inputs,outputs):
        # use each input value as an index into the table to lookup output value
        outputs[0][:] = self._table[inputs[0]]
        self._propagate_outputs()
        return True

################################################################################
# delay: delay output by specified number of samples
################################################################################

class delay(dsp_block):
    _default_parameters = {
        'delay': None,    # number of samples to delay
        }
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _initialize(self):
        assert isinstance(self.delay,int) and self.delay >= 0,\
               "%s: delay must be a non-negative integer" % self._name
        # how many buffers of 0's to output first
        self._zcount = int(self.delay / self.buffer_size)
        # how many samples to shift each input buffer
        self._delay = self.delay - self.buffer_size*self._zcount
        # where we remember samples from last input buffer
        self._buffer = numpy.empty(self._delay,dtype=self._input_types[0])

    def _reset(self):
        self._zcount = self.delay // self.buffer_size   # ooo! modern python....
        self._buffer[:] = 0
        dsp_block._reset(self)

    def _buffers_ready(self,inputs,outputs):
        if self._zcount > 0:
            # for large delays, we have to start with some number of
            # output buffers full of zeros
            self._zcount -= 1
            outputs[0][:] = 0
            self._propagate_outputs()
            # no inputs consumed yet
            return False

        # output data starts with values saved from last time
        outputs[0][0:self._delay] = self._buffer[:]
        # how much of current input buffer we'll consume now
        end = self._buffer_size - self._delay
        outputs[0][self._delay:] = inputs[0][0:end]
        self._buffer[:] = inputs[0][end:]

        self._propagate_outputs()
        return True

__all__=[
    'convolve',
    'fir_filter',
    'fir_low_pass',
    'fir_high_pass',
    'fir_band_pass',
    'fir_band_reject',
    'decimate',
    'interpolate',
    'transmit',
    'slice',
    'repack',
    'lookup',
    'delay',
    ]
