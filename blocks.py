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

import Queue,sys,threading
import numpy

current_testbench = None   # instance of testbench
instance_counts = {}       # maps class names -> count of instances

def make_instance_name(class_name):
    count = instance_counts.get(class_name,0) + 1
    instance_counts[class_name] = count
    return "%s%d" % (class_name,count)

################################################################################
# block: instances of this appear on the run queue
################################################################################

class block(object):
    def __init__(self):
        self._name = make_instance_name(self.__class__.__name__)
        self._run_queue = current_testbench.run_queue
        self._mutex = threading.Lock()  # for synchronizing access to this block
        self._queued = False

    def __repr__(self):
        return '<%s>' % self._name

    # should be called from within _run (ie, with muxtex acquired)
    def _enqueue(self):
        assert not self._queued,"%s: oops, already in the run queue!" % self._name
        self._queued = True
        self._run_queue.put(self)

    def _enqueue_if_necessary(self):
        self._mutex.acquire()
        if not self._queued:
            self._enqueue()
        self._mutex.release()

    # this method is called when block is taken off the run queue
    def _execute(self):
        self._mutex.acquire()
        #print "    %s: acquired lock" % self._name
        self._queued = False
        self._run()
        self._mutex.release()
        #print "    %s: released lock" % self._name

    # called once at beginning of simulation
    def _reset(self):
        raise NotImplementedError

    # override to do processing.  Call _enqueue if block should be
    # rescheduled for execution.
    def _run(self):
        raise NotImplementedError

################################################################################
# queue_block: block + queue, base class for python producers, consumers
################################################################################

class queue_block(block):
    def __init__(self,qsize = 10):
        block.__init__(self)
        self._qsize = qsize
        self._buffer_size = current_testbench.buffer_size

        # finally add this block to the test bench
        current_testbench.add_block(self)

    # called before starting a simulation run
    def _reset(self):
        self._eos = False
        self._queue = Queue.Queue(self._qsize)
        self._enqueue()

    def _run(self):
        # handles both production and consumption cases.
        # only one of the methods should be meaningful.
        if not self._queue.full():
            self._produce()
        if not self._queue.empty():
            self._consume()

    # producers should override to return data type for output stream
    def _dtype(self):
        raise NotImplemented

    # override this method to include your Python production code
    # which should deliver data by calling self._queue.put(data).
    # also override _dtype()...
    def _produce(self):
        pass

    # override this method to include your Python consumption code
    # which should consume data by calling self._queue.get().
    def _consume(self):
        pass

    # called by external agents to get something from the queue, return None if empty.
    # self should be an instance of a producer...
    def _get(self):
        try:
            result = self._queue.get(False)
        except Queue.Empty:
            result = None
        self._enqueue_if_necessary()  # wake up producer
        return result

    # called by external agents to put something into the queue, return True if succeeded
    # self should be an instance of a consumer...
    def _put(self,item):
        try:
            self._queue.put(item,False)
            result = True
        except Queue.Full:
            result = False
        self._enqueue_if_necessary()  # wake up consumer
        return result

################################################################################
# block: base class for all DSP blocks
################################################################################

class dsp_block(block):
    _default_parameters = {}  # dictionary of parameters and their default values
    _expected_number_of_inputs = (1,sys.maxint)   # range for number of inputs
    _output_ports = (())  # sequence of types for output buffers
                          # if an element is 'inherit' it means same as first input
    
    # output_list should be a list of numpy types, one for each output port
    # input_list should be a list of numpy types, one for each input port
    def __init__(self,*inputs,**parameters):
        block.__init__(self)
        self._setup_parameters(parameters)
        self._setup_inputs(inputs)
        self._setup_outputs()
        
        # keep track of who consumes our outputs
        self._nfanouts = 0    # total number of modules who consume our outputs
        self._fanouts = []    # list of blocks that consume our outputs

        self._input_wait = self._nfanins    # wait for input values to arrive
        self._output_wait = 0   # output buffers are empty

        # complete any block-specific initialization
        self._initialize()

        # finally add this block to the test bench
        current_testbench.add_block(self)

    # convert user-supplied and default parameters into attributes of this instance
    def _setup_parameters(self,parameters):
        self._buffer_size = current_testbench.buffer_size

       # set up our attributes by first copying default parameters
        for k,v in self._default_parameters.items():
            assert not self.__dict__.has_key(k),\
                   "%s: parameter name %s conflicts with existing attribute" % (self._name,k)
            self.__dict__[k] = v
        # now update with user-supplied properties
        for k,v in parameters.items():
            assert self.__dict__.has_key(k),\
                   "%s: unrecognized parameter %s" % (self._name, k)
            self.__dict__[k] = v

        # make sure we have sample_rate attribute
        if not self.__dict__.has_key('sample_rate'):
            self.sample_rate = None

    # process inputs to block
    def _setup_inputs(self,inputs):
        self._ninputs = len(inputs)
        assert self._ninputs >= self._expected_number_of_inputs[0] and \
               self._ninputs <= self._expected_number_of_inputs[1],\
               "%s: expected between %d and %d inputs, got %d" % \
               (self._name,self._expected_number_of_inputs[0],self._expected_number_of_inputs[1])

        self._inputs = []   # list of buffers/constant values
        self._input_types = []

        self._fanins = []   # list of blocks who supply us values
        self._nfanins = 0   # len(_fanins)

        self._input_sample_rate = None
        for i in xrange(self._ninputs):
            input = inputs[i]

            if isinstance(input,block_port):
                driver = input.block
                port = input.port
            elif isinstance(input,block):
                driver = input
                port = 0
            else:
                driver = None
                port = None

            if driver:
                # handle an input supplied by another block
                assert port < driver._noutputs,\
                       "%s: %s doesn't have a port %d" % (self._name,driver._name,port)
                buffer,buffer_type = driver._add_fanout(port,self)
                self._nfanins += 1
                self._fanins.append(driver)
                self._inputs.append(buffer)
                self._input_types.append(buffer_type)
                if self._input_sample_rate:
                    assert self._input_sample_rate==driver.sample_rate,\
                           "%s: sample rates of inputs must match" % self._name
                else:
                    self._input_sample_rate = driver.sample_rate
            else:
                # handle a constant input
                self._inputs.append(input)
                self._input_types.append(None)

        # if user hasn't specified a sample rate, inherit from the inputs or default
        # to whatever the testbench has been told
        if self.sample_rate is None:
            self.sample_rate = self._input_sample_rate or current_testbench.sample_rate

    # set up output buffers using supplied list of port types, or self._output_ports
    def _setup_outputs(self,output_ports=None):
        if output_ports is None: output_ports = self._output_ports
        self._noutputs = len(output_ports)
        self._outputs = []
        self._output_types = []
        self._output_connected = []
        for output_type in output_ports:
            if output_type == 'inherit':
                output_type = self._input_types[0]
            self._outputs.append(numpy.empty(current_testbench.buffer_size,output_type))
            self._output_types.append(output_type)
            self._output_connected.append(False)

    # allow block[i] to refer to the ith output port
    def __getitem__(self,index):
        assert index < self._noutputs,\
               "%s: port index (%d) too large, only have %d ports" % (self._name,index,self._noutputs)
        if index == 0: return self
        else: return block_port(self,index)

    # override for block-specific initialization (eg, allocate internal buffers)
    def _initialize(self):
        pass

    # add a new fanout block to specified output port
    # return buffer and its type for that output
    def _add_fanout(self,port_number,fanout):
        assert port_number < self._noutputs,\
               "%s: doesn't have output port %d" % (self._name,port_number)
        self._nfanouts += 1
        self._output_connected[port_number] = True
        self._fanouts.append(fanout)
        return (self._outputs[port_number],self._output_types[port_number])

    # called by fanin blocks when they've filled one of our input buffers
    def _fanin_ready(self):
        self._mutex.acquire()
        assert self._input_wait > 0,"%s oops, input_wait already 0!" % self._name
        self._input_wait -= 1
        self._enqueue_if_runnable()
        self._mutex.release()

    # called by fanout blocks when they've emptied one of our output buffers.
    # return True if we've reached end-of-stream
    def _fanout_ready(self):
        self._mutex.acquire()
        assert self._output_wait > 0,"%s oops, output_wait already 0!" % self._name
        self._output_wait -= 1
        self._enqueue_if_runnable()
        self._mutex.release()
        return self._eos

    # internal routine that adds this block to run queue if it's runnable.
    # should be called with mutex locked...
    def _enqueue_if_runnable(self):
        if not self._eos and self._input_wait == 0 and self._output_wait == 0:
            self._enqueue()
            #print "queueing %s" % self._name

    # called before starting a simulation run
    def _reset(self):
        self._eos = False       # haven't reached end of stream yet
        self._input_wait = self._nfanins    # wait for input values to arrive
        self._output_wait = 0   # output buffers are empty
        self._enqueue_if_runnable()

    # called by testbench when it's our chance to run
    def _run(self):
        assert self._input_wait==0 and self._output_wait==0,\
               "%s: oops, running but input_wait=%d, output_wait=%d" % \
               (self._name,self._input_wait,self._output_wait)

        eos = self._eos
        if self._buffers_ready(self._inputs,self._outputs):
            #print "    %s: notifying fanins" % self._name
            # if we consumed the inputs, tell our suppliers to send more
            self._input_wait = self._nfanins  # wait until all inputs have arrived
            for b in self._fanins:
                # remember if one or more of our inputs has reached end-of-stream
                eos |= b._fanout_ready()

        if not self._eos and eos:
            # we've just reached the end-of-stream.  This will be
            # propagated to our outputs when they call _fanout_ready()
            # to report finishing with the current output buffer.
            self._eos = True
            current_testbench.reached_eos(self)

        # handle case where block hasn't consumed inputs or produced outputs
        # but is able to continue running (ie, the block is just yielding).
        self._enqueue_if_runnable()

    # called by _run when inputs are full and outputs are empty.  Return True
    # if input values were consumed.  Call _propagate_outputs if new output
    # values were produced
    def _buffers_ready(self,inputs,outputs):
        raise NotImplementedError

    # called by _buffers_ready if there are new outputs to propagate
    def _propagate_outputs(self):
        #print "    %s: notifying fanouts" % self._name
        self._output_wait = self._nfanouts   # wait until all outputs have been consumed
        for b in self._fanouts: b._fanin_ready()

    # wow! signal arithmetic
    def __add__(self,other): return add(self,other)
    def __sub__(self,other): return subtract(self,other)
    def __mul__(self,other): return multiply(self,other)
    def __rmul__(self,other): return multiply(self,other)
    def __div__(self,other): return divide(self,other)
    def __truediv__(self,other): return true_divide(self,other)
    def __floordiv__(self,other): return floor_divide(self,other)
    def __mod__(self,other): return mod(self,other)
    def __pow__(self,other): return power(self,other)
    def __lshift__(self,other): return left_shift(self,other)
    def __rshift__(self,other): return right_shift(self,other)
    def __and__(self,other): return logical_and(self,other)
    def __or__(self,other): return logical_or(self,other)
    def __xor__(self,other): return logical_xor(self,other)

    def __radd__(self,other): return add(other,self)
    def __rsub__(self,other): return subtract(other,self)
    def __rmul__(self,other): return mul(other,self)
    def __rdiv__(self,other): return divide(other,self)
    def __rtruediv__(self,other): return true_divide(other,self)
    def __rfloordiv__(self,other): return floor_divide(other,self)
    def __rmod__(self,other): return mod(other,self)
    def __rpow__(self,other): return power(other,self)
    def __rlshift__(self,other): return left_shift(other,self)
    def __rrshift__(self,other): return right_shift(other,self)
    def __rand__(self,other): return logical_and(other,self)
    def __ror__(self,other): return logical_or(other,self)
    def __rxor__(self,other): return logical_xor(other,self)

    def __neg__(self): return negative(self)
    def __pos__(self): return self
    def __invert__(self): return logical_not(self)

    def __lt__(self,other): return less(self,other)
    def __le__(self,other): return less_equal(self,other)
    def __gt__(self,other): return greater(self,other)
    def __ge__(self,other): return greater_equal(self,other)
    def __eq__(self,other): return equal(self,other)
    def __ne__(self,other): return not_equal(self,other)

################################################################################
# block_port: a block object representing an output port
################################################################################

class block_port(dsp_block):
    def __init__(self,block,port):
        self.block = block
        self.port = port
        self.sample_rate = block.sample_rate

    # no one should be calling the usual methods...
    def _setup_parameters(self,parameters): raise NotImplementedError
    def _setup_inputs(self,inputs): raise NotImplementedError
    def _setup_outputs(self,output_ports=None): raise NotImplementedError
    def __getitem__(self,index): raise NotImplementedError
    def _initialize(self): raise NotImplementedError
    def _add_fanout(self,port_number,fanout): raise NotImplementedError
    def _fanin_ready(self): raise NotImplementedError
    def _fanout_ready(self): raise NotImplementedError
    def _enqueue_if_runnable(self): raise NotImplementedError
    def _reset(self): raise NotImplementedError
    def _run(self): raise NotImplementedError
    def _buffers_ready(self,inputs,outputs): raise NotImplementedError
    def _propagate_outputs(self): raise NotImplementedError

    # but we'll inherit the other __xxx__ to support arithmetic operations...

################################################################################
# add: out = in1 + in2 + ...
################################################################################

class add(dsp_block):
    _output_ports = ('inherit',)
    
    def _buffers_ready(self,inputs,outputs):
        if self._ninputs == 1:
            outputs[0][:] = inputs[0]
        else:
            numpy.add(inputs[0],inputs[1],outputs[0])
            for i in xrange(2,self._ninputs):
                numpy.add(outputs[0],inputs[i],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# subtract: out = in1 - in2 - ...
################################################################################

class subtract(dsp_block):
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        if self._ninputs == 1:
            outputs[0][:] = inputs[0]
        else:
            numpy.subtract(inputs[0],inputs[1],outputs[0])
            for i in xrange(2,self._ninputs):
                numpy.subtract(outputs[0],inputs[i],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# multiply: out = in1 * in2 * ...
################################################################################

class multiply(dsp_block):
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        if self._ninputs == 1:
            outputs[0][:] = inputs[0]
        else:
            numpy.multiply(inputs[0],inputs[1],outputs[0])
            for i in xrange(2,self._ninputs):
                numpy.multiply(outputs[0],inputs[i],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# divide: out = in1 / in2 / ...
################################################################################

class divide(dsp_block):
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        if self._ninputs == 1:
            outputs[0][:] = inputs[0]
        else:
            numpy.divide(inputs[0],inputs[1],outputs[0])
            for i in xrange(2,self._ninputs):
                numpy.divide(outputs[0],inputs[i],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# true_divide: out = in1 / in2 / ...
################################################################################

class true_divide(dsp_block):
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        if self._ninputs == 1:
            outputs[0][:] = inputs[0]
        else:
            numpy.true_divide(inputs[0],inputs[1],outputs[0])
            for i in xrange(2,self._ninputs):
                numpy.true_divide(outputs[0],inputs[i],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# floor_divide: out = floor(in1 // in2)
################################################################################

class floor_divide(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.floor_divide(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# negative: out = -in1
################################################################################

class negative(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.negative(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# power: out = in1**in2
################################################################################

class power(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.power(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# remainder: out = in1 - floor(in1/in2)*in2
################################################################################

class remainder(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.remainder(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# mod: like remainder
################################################################################

class mod(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.mod(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# fmod: like remainder but output has sign of in1
################################################################################

class fmod(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.fmod(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# absolute: out = absolute(in1)
################################################################################

class absolute(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.absolute(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# rint: out = int(in1)
################################################################################

class rint(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.rint(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# sign: out = sign(in1)
################################################################################

class sign(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.sign(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# conj: out = conj(in1)
################################################################################

class conj(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.conj(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# exp: out = e**in1
################################################################################

class exp(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.exp(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# log: out = natural log of in1
################################################################################

class log(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.log(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# expm1: out = e**in1 - 1
################################################################################

class expm1(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.expm1(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# lop1p: out = log(1+in)
################################################################################

class log1p(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.log1p(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# log10: out = log base 10 in1
################################################################################

class log10(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.log10(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# sqrt: out = sqrt(in1)
################################################################################

class sqrt(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.sqrt(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# square: out = square(in1)
################################################################################

class square(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.square(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# reciprocal: out = 1/in
################################################################################

class reciprocal(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.reciprocal(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# angle: out = angle(in1 + j*in2)
# counterclockwise angle from the positive real axis (range: 0 to pi, -pi to 0)
################################################################################

class angle(dsp_block):
    _expected_number_of_inputs = (1,2)
    _output_ports = (numpy.float64,)

    def _initialize(self):
        self._temp = numpy.empty(self._buffer_size,dtype=numpy.complex_)

    def _buffers_ready(self,inputs,outputs):
        if self._ninputs == 2:
            numpy.multiply(inputs[1],1j,self._temp)
        else:
            self._temp[:] = 0 + 0j
        numpy.add(inputs[0],self._temp,self._temp)
        angle = numpy.angle(self._temp)
        outputs[0][:] = angle
        del angle
        self._propagate_outputs()
        return True

################################################################################
# sin
################################################################################

class sin(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.sin(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# cos
################################################################################

class cos(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.cos(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# tan
################################################################################

class tan(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.tan(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arcsin
################################################################################

class arcsin(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arcsin(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arccos
################################################################################

class arccos(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arccos(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arctan
################################################################################

class arctan(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arctan(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arctan2
################################################################################

class arctan2(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arctan2(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# hypot
################################################################################

class hypot(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.hypot(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# sinh
################################################################################

class sinh(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.sinh(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# cosh
################################################################################

class cosh(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.cosh(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# tanh
################################################################################

class tanh(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.tanh(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arcsinh
################################################################################

class arcsinh(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arcsinh(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arccosh
################################################################################

class arccosh(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arccosh(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# arctanh
################################################################################

class arctanh(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.arctanh(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# bitwise_and
################################################################################

class bitwise_and(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.bitwise_and(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# bitwise_or
################################################################################

class bitwise_or(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.bitwise_or(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# bitwise_xor
################################################################################

class bitwise_xor(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.bitwise_xor(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# invert
################################################################################

class invert(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.invert(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# left_shift
################################################################################

class left_shift(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.left_shift(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# right_shift
################################################################################

class right_shift(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.right_shift(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# greater
################################################################################

class greater(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.greater(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# greater_equal
################################################################################

class greater_equal(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.greater_equal(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# less
################################################################################

class less(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.less(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# less_equal
################################################################################

class less_equal(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.less_equal(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# not_equal
################################################################################

class not_equal(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.not_equal(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# equal
################################################################################

class equal(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.equal(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# logical_and
################################################################################

class logical_and(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.logical_and(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# logical_or
################################################################################

class logical_or(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.logical_or(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# logical_xor
################################################################################

class logical_xor(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.logical_xor(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# logical_not
################################################################################

class logical_not(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.logical_not(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# maximum
################################################################################

class maximum(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.maximum(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# minimum
################################################################################

class minimum(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.minimum(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# signbit
################################################################################

class signbit(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = (numpy.bool,)

    def _buffers_ready(self,inputs,outputs):
        numpy.signbit(inputs[0],inputs[1],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# modf
################################################################################

class modf(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit','inherit')

    def _buffers_ready(self,inputs,outputs):
        numpy.modf(inputs[0],outputs[0],outputs[1])
        self._propagate_outputs()
        return True

################################################################################
# ldexp
################################################################################

class ldexp(dsp_block):
    _expected_number_of_inputs = (2,2)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        result = numpy.ldexp(inputs[0],inputs[1])
        outputs[0][:] = result
        del result
        self._propagate_outputs()
        return True

################################################################################
# frexp
################################################################################

class frexp(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit','inherit')

    def _buffers_ready(self,inputs,outputs):
        numpy.frexp(inputs[0],outputs[0],outputs[1])
        self._propagate_outputs()
        return True

################################################################################
# floor
################################################################################

class floor(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.floor(inputs[0],outputs[0])
        self._propagate_outputs()
        return True

################################################################################
# ceil
################################################################################

class ceil(dsp_block):
    _expected_number_of_inputs = (1,1)
    _output_ports = ('inherit',)

    def _buffers_ready(self,inputs,outputs):
        numpy.ceil(inputs[0],outputs[0])
        self._propagate_outputs()
        return True


__all__ = [
    'block','queue_block','dsp_block',
    'add','subtract','multiply','divide','true_divide','floor_divide',
    'negative','power','remainder','mod','fmod','absolute',
    'rint','sign','conj','exp','log','expm1','log1p','log10',
    'sqrt','square','reciprocal','angle',
    'sin','cos','tan','arcsin','arccos','arctan','arctan2','hypot',
    'sinh','cosh','tanh','arcsinh','arccosh','arctanh',
    'bitwise_and','bitwise_or','bitwise_xor','invert','left_shift','right_shift',
    'greater','greater_equal','less','less_equal','not_equal','equal',
    'logical_and','logical_or','logical_xor','logical_not',
    'maximum','minimum','signbit','modf','ldexp','frexp','floor','ceil'
    ]
