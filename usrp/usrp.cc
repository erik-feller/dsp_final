/**
 * usrp: Python wrapper for libusrp
 *
 * Copyright (c) 2008 Chris Terman
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "Python.h"
#include "structmember.h"

#include <usrp_standard.h>	// includes usrp_basic.h
#include <usrp_slots.h>
#include <usrp_dbid.h>
#include <usrp_prims.h>

// setter for read-only values
static int noset(void *self,PyObject *value, void *closure) {
  PyErr_SetString(PyExc_AttributeError,"usrp: read-only attribute");
  return NULL;
}

//////////////////////////////////////////////////
// usrp.rx: USRP receiver object  (built on usrp_standard_rx in libusrp)
//
// rx = usrp.rx(decim_rate,which_board=0,nchan=2,mux=-1,mode=FPGA_MODE_NORMAL,
//              fusb_block_size=0,fusb_nblocks=0,fpga_filename="",firmware_filename="")
//   defaults used:
//   mux = 0x32103210
//   firmware_filename = "/usr/local/share/usrp/rev4/std.ihx"
//   fpga_filename = "/usr/local/share/usrp/rev4/std_2rxhb_2tx.rbf"
//   fusb_block_size = 4096   # max is 16K
//   fusb_nblocks = 256       # FUSB_BUFFER_SIZE/FUSB_BLOCK_SIZE, BUFSIZE = 1M
//
// int = rx.start()   # start receiver, bool success flag
// (string,int) = rx.read([byte_bufsize])    # returns buffer,overrun flag
// int = rx.stop()  # stop receiver, bool success flag
//
// int = rx.daughterboard_id(which_dboard)  # returns USRP_DBID_xxx
// int = rx.set_decim_rate(rate)  # set decimator rate: even, 4..256, boolean success flag
// int = rx.decim_rate  # decimator rate
// int = rx.set_nchannels(nchannels)  # set number of active channels: 1, 2 or 4, boolean success flag
// int = rx.nchannels  # number of active channels
// int = rx.set_mux(value)   # set input mux configuration, boolean success flag
//   value is four pairs of 4-bit fields: Q3I3Q2I2Q1I1Q0I0
//   Ix,Qx select which physical channel is routed to which DDC
//   for Qx can specify 0xf for constant 0 value
//   all Qx must 0xf, or none of them may be 0xf
// int = rx.mux    # input mux selection
// int = rx.set_format(shift=0,width=16,want_q=1,bypass_halfband=0)  # returns new format word
//   at the moment only accepts:
//   shift=0,width=16,want_q=1,bypass_halfband=0  (16-bit data, 2 channels)
//   shift=8,width=8,want_q=1,bypass_halfband=0  (shifted 8-bit data, 2 channels)
// int = rx.format # rx data format
// int = rx.set_rx_freq(channel,freq)  # set rx_freq, boolean success flag
// double = rx.rx_freq(channel)  # get rx_freq, channel is 0..3
// int = rx.set_fpga_mode(mode)   # set fpga_mode, one of FPGA_MODE_xxx
// int = rx.set_ddc_phase(channel,phase)  # set digital down converter phase register
// int = rx.set_pga(which,gain_in_dB)  # set programmable amplifier: which=0..3
// double = rx.get_pga(which)  # get programmable amplifier
// double = rx.pga_min # min legal PGA gain in dB (0)
// double = rx.pga_max # max legal PGA gain in dB (20.0)
// double = rx.pga_db_per_step # hardware step size of PGA (20.0/20)
// int = rx.block_size  # fusb block size
// int = rx.adc_rate # ADC conversion rate
// string = rx.serial_number # serial number of USRP board
//////////////////////////////////////////////////

typedef struct {
  PyObject_HEAD
  usrp_standard_rx *rx;
} rx;

// properties defined in usrp_standard_rx

// int = rx.decim_rate  # decimator rate
static PyObject *rx_get_decim_rate(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->decim_rate());
}

// int = rx.nchannels  # number of active channels
static PyObject *rx_get_nchannels(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->nchannels());
}

// int = rx.mux    # input mux selection
static PyObject *rx_get_mux(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->mux());
}

// int = rx.format # rx data format
static PyObject *rx_get_format(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->format());
}

// double = rx.pga_min  # min legal PGA gain in dB
static PyObject *rx_get_pga_min(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("d",self->rx->pga_min());
}

// double = rx.pga_max  # max legal PGA gain in dB
static PyObject *rx_get_pga_max(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("d",self->rx->pga_max());
}

// double = rx.pga_db_per_step  # hardware PGA step size in dB
static PyObject *rx_get_pga_db_per_step(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("d",self->rx->pga_db_per_step());
}

// int = rx.block_size  # fusb block size
static PyObject *rx_get_block_size(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->block_size());
}

// int = rx.adc_rate  # ADC conversion rate
static PyObject *rx_get_adc_rate(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->adc_rate());
}

// string = rx.serial_number  # USRP board serial number
static PyObject *rx_get_serial_number(rx *self,void *closure) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("s",self->rx->serial_number().c_str());
}

static PyGetSetDef rx_getseters[] = {
  {"decim_rate",(getter)rx_get_decim_rate,(setter)noset,"decimator rate",NULL},
  {"nchannels",(getter)rx_get_nchannels,(setter)noset,"number of active channels",NULL},
  {"mux",(getter)rx_get_mux,(setter)noset,"input mux control",NULL},
  {"format",(getter)rx_get_format,(setter)noset,"rx data format",NULL},
  {"pga_min",(getter)rx_get_pga_min,(setter)noset,"min legal PGA gain in dB",NULL},
  {"pga_max",(getter)rx_get_pga_max,(setter)noset,"max legal PGA gain in dB",NULL},
  {"pga_db_per_step",(getter)rx_get_pga_db_per_step,(setter)noset,"PGA step size in dB",NULL},
  {"block_size",(getter)rx_get_block_size,(setter)noset,"fusb block size",NULL},
  {"adc_rate",(getter)rx_get_adc_rate,(setter)noset,"fusb block size",NULL},
  {"serial_number",(getter)rx_get_serial_number,(setter)noset,"USRP serial number",NULL},
  {NULL}  // sentinel
};

// rx = usrp.rx(decim_rate,which_board=0,nchan=2,mux=-1,mode=FPGA_MODE_NORMAL,
//              fusb_block_size=0,fusb_nblocks=0,fpga_filename="",firmware_filename="")
// defaults used:
//   mux = 0x32103210
//   firmware_filename = /usr/local/share/usrp/rev4/std.ihx
//   fpga_filename = /usr/local/share/usrp/rev4/std_2rxhb_2tx.rbf
//   fusb_block_size = 4096   # max is 16K
//   fusb_nblocks = 256       # FUSB_BUFFER_SIZE/FUSB_BLOCK_SIZE, BUFSIZE = 1M
static PyObject *rx_new(PyTypeObject *type,PyObject *args,PyObject *kwds) {
  rx *self;

  self = (rx *)type->tp_alloc(type,0);
  if (self != NULL) {
    self->rx = NULL;   // do the heavy lifting in init
  }
  return (PyObject *)self;
}

// initialize rx object
static int rx_init(rx *self,PyObject *args,PyObject *kwds) {
  unsigned int decim_rate;	// required, even in range 8..256
  int which_board = 0;		// optional, 0 is what you usually want
  int nchan = 2;		// 1, 2 or 4
  int mux = -1;			// use default
  int mode = usrp_standard_rx::FPGA_MODE_NORMAL;
  int fusb_block_size = 0;	// use default
  int fusb_nblocks = 0;		// use default
  char *fpga_filename = "";	// use default
  char *firmware_filename = ""; // use default

  static char *kwlist[] = {
    "decim_rate",
    "board",
    "nchan",
    "mux",
    "mode",
    "fusb_block_size",
    "fusb_nblocks",
    "fpga_filename",
    "firmware_filename",
    NULL};

  // decim_rate is required, rest of args are optional
  if (!PyArg_ParseTupleAndKeywords(args,kwds,"i|iiiiiiss",kwlist,
				   &decim_rate,
				   &which_board,
				   &nchan,
				   &mux,
				   &mode,
				   &fusb_block_size,
				   &fusb_nblocks,
				   &fpga_filename,
				   &firmware_filename))
    return -1;

  if ((decim_rate & 1) || decim_rate < 4 || decim_rate > 256) {
    PyErr_SetString(PyExc_ValueError,"decim_rate must be even and in the range 4..256");
    return -1;
  }

  if (self->rx != NULL) delete self->rx;
  self->rx = usrp_standard_rx::make(which_board,decim_rate,nchan,mux,mode,
				    fusb_block_size,fusb_nblocks,
				    fpga_filename,firmware_filename);

  if (self->rx == NULL) {
    PyErr_SetString(PyExc_IOError,"failed to create usrp_standard_rx input");
    return -1;
  }

  return 0;
}

// deallocate rx object
static void rx_dealloc(rx *self) {
  if (self->rx) delete self->rx;
  self->ob_type->tp_free((PyObject*)self);  // free the object
}

static PyMemberDef rx_members[] = {
  // members accesible from Python
  // {namestring,type,offsetof(rx,member),flags,docstring}
  // type: T_INT for integer, T_OBJECT_EX for PyObject *
  // flags: READONLY,READ_RESTRICTED,WRITE_RESTRICTED,RESTRICTED
  {NULL}  // sentinel
};

// methods defined in usrp_standard_rx

// int = rx.set_decim_rate(rate)  # set decimator rate: even, 8..256
static PyObject *rx_set_decim_rate(rx *self,PyObject *args) {
  int decim_rate;

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&decim_rate))
    return NULL;

  if ((decim_rate & 1) || decim_rate < 4 || decim_rate > 256) {
    PyErr_SetString(PyExc_ValueError,"decim_rate must be in the range 4..256");
    return NULL;
  }
    
  return Py_BuildValue("i",self->rx->set_decim_rate(decim_rate));
}

// int = rx.set_nchannels(nchannels)  # set number of active channels: 1, 2 or 4, boolean success flag
static PyObject *rx_set_nchannels(rx *self,PyObject *args) {
  int nchannels; 

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&nchannels))
    return NULL;
    
  if (nchannels != 1 && nchannels != 2 && nchannels != 4) {
    PyErr_SetString(PyExc_ValueError,"The nchannels attribute must be 1, 2 or 4");
    return NULL;
  }
    
  return Py_BuildValue("i",self->rx->set_nchannels(nchannels));
}

// int = rx.set_mux(value)   # set input mux configuration, boolean success flag
//   value is four pairs of 4-bit fields: Q3I3Q2I2Q1I1Q0I0
//   Ix,Qx select which physical channel is routed to which DDC
//   for Qx can specify 0xf for constant 0 value
//   all Qx must 0xf, or none of them may be 0xf
static PyObject *rx_set_mux(rx *self,PyObject *args) {
  int mux; 

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&mux))
    return NULL;
    
  // check I values
  for (int i = 0; i < 4; i += 1) {
    int v = (mux >> (i*8)) & 0xF;
    if (v > 3) {
      PyErr_SetString(PyExc_ValueError,"I Channel settings must be 0, 1, 2 or 3");
      return NULL;
    }
  }

  // check Q values
  int select_const_zero = 0;
  for (int q = 0; q < 4; q += 1) {
    int v = (mux >> (q*8 + 4)) & 0xF;
    if (v == 0xF) select_const_zero |= (1 << q);
    else if (v > 3) {
      PyErr_SetString(PyExc_ValueError,"Q Channel settings must be 0, 1, 2, 3, or 0xF");
      return NULL;
    }
  }
  if (select_const_zero != 0 && select_const_zero != 0xF) {
    PyErr_SetString(PyExc_ValueError,"Q Channel settings must all be 0xF or none of them 0xF");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->set_mux(mux));
}

// int = rx.set_rx_freq(channel,freq)  # set rx_freq, boolean success flag
static PyObject *rx_set_rx_freq(rx *self,PyObject *args) {
  int channel;
  double freq;

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"id",&channel,&freq))
    return NULL;

  if (channel < 0 || channel > 3) {
    PyErr_SetString(PyExc_ValueError,"channel must be 0, 1, 2 or 3");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->set_rx_freq(channel,freq));
}

// double = rx.get_rx_freq(channel)  # get rx_freq, channel is 0..3
static PyObject *rx_get_rx_freq(rx *self,PyObject *args) {
  int channel;

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&channel))
    return NULL;

  if (channel < 0 || channel > 3) {
    PyErr_SetString(PyExc_ValueError,"channel must be 0, 1, 2 or 3");
    return NULL;
  }

  return Py_BuildValue("d",self->rx->rx_freq(channel));
}

// int = rx.set_fpga_mode(mode)   # set fpga_mode, one of FPGA_MODE_xxx
static PyObject *rx_set_fpga_mode(rx *self,PyObject *args) {
  int mode;

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&mode))
    return NULL;

  return Py_BuildValue("i",self->rx->set_fpga_mode(mode));
}

// int = rx.set_ddc_phase(channel,phase)  # set digital down converter phase register
static PyObject *rx_set_ddc_phase(rx *self,PyObject *args) {
  int channel;
  int phase;

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"ii",&channel,&phase))
    return NULL;

  if (channel < 0 || channel > 3) {
    PyErr_SetString(PyExc_ValueError,"channel must be 0, 1, 2 or 3");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->set_ddc_phase(channel,phase));
}

// int = rx.set_format(shift=0,width=16,want_q=1,bypass_halfband=0)  # returns new format word
// at the moment only accepts:
//   shift=0,width=16,want_q=1,bypass_halfband=0  (16-bit data, 2 channels)
//   shift=8,width=8,want_q=1,bypass_halfband=0  (shifted 8-bit data, 2 channels)
static PyObject *rx_set_format(rx *self,PyObject *args,PyObject *kwds) {
  int shift = 0;
  int width = 16;
  int want_q = 1;
  int bypass_halfband = 0;

  static char *kwdlist[] = {"shift","width","want_q","bypass_halfband",NULL};

  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  if (!PyArg_ParseTupleAndKeywords(args,kwds,"|iiii",kwdlist,
				   &shift,&width,&want_q,&bypass_halfband))
    return NULL;

  // check for legal combinations
  unsigned int format = usrp_standard_rx::make_format(width,shift,want_q,bypass_halfband);
  if (format != 0x300 && format != 0x288) {
    PyErr_SetString(PyExc_ValueError,"format must be b=0,q=1,w=16,s=0 or b=0,q=1,w=8,s=8");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->set_format(format));
}

// bool = rx.start()   # start receiver, bool success flag
static PyObject *rx_start(rx *self) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->start());
}

// bool = rx.stop()  # stop receiver, bool success flag
static PyObject *rx_stop(rx *self) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->stop());
}

// methods defined in usrp_basic_rx

// (string,int) = rx.read([byte_bufsize])    # returns buffer,overrun flag
static PyObject *rx_read(rx *self,PyObject *args) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  int bufsize = 512;
  if (!PyArg_ParseTuple(args,"|i",&bufsize))
    return NULL;

  if (bufsize < 512 || (bufsize % 512) != 0) {
    PyErr_SetString(PyExc_ValueError,"bufsize must be a multiple of 512");
    return NULL;
  }

  char *buf = (char *)alloca(bufsize);
  bool overrun;
  int status;

  Py_BEGIN_ALLOW_THREADS
  status = self->rx->read((void *)buf,bufsize,&overrun);
  Py_END_ALLOW_THREADS

  if (status <= 0) return Py_BuildValue("(s,i)","",overrun);
  else return Py_BuildValue("(N,i)",PyString_FromStringAndSize(buf,status),overrun);
}

// int = rx.daughterboard_id(which_dboard)  # returns USRP_DBID_xxx
static PyObject *rx_daughterboard_id(rx *self,PyObject *args) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  int which_dboard;

  if (!PyArg_ParseTuple(args,"i",&which_dboard))
    return NULL;

  int dbid = self->rx->daughterboard_id(which_dboard);

  return Py_BuildValue("i",dbid);
}

// int = rx.set_pga(which,gain_in_dB)  # set programmable amplifier: which=0..3, bool success flag
static PyObject *rx_set_pga(rx *self,PyObject *args) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  int which;
  double gain;

  if (!PyArg_ParseTuple(args,"id",&which,&gain))
    return NULL;

  if (which < 0 || which > 3) {
    PyErr_SetString(PyExc_ValueError,"which PGA must be in range 0..3");
    return NULL;
  }

  if (gain < self->rx->pga_min() || gain > self->rx->pga_max()) {
    PyErr_SetString(PyExc_ValueError,"PGA gain outside of allowable range (0..20)");
    return NULL;
  }

  return Py_BuildValue("i",self->rx->set_pga(which,gain));
}

// double = rx.get_pga(which)  # get programmable amplifier
static PyObject *rx_get_pga(rx *self,PyObject *args) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  int which;

  if (!PyArg_ParseTuple(args,"i",&which))
    return NULL;

  if (which < 0 || which > 3) {
    PyErr_SetString(PyExc_ValueError,"which PGA must be in range 0..3");
    return NULL;
  }

  return Py_BuildValue("d",self->rx->pga(which));
}


// int = rx.set_dc_offset_cl_enable(bits,mask)  # enable/disable automatic DC offset control loop
//   bool success flag, enabled by default...
static PyObject *rx_set_dc_offset_cl_enable(rx *self,PyObject *args) {
  if (!self->rx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp receiver info available");
    return NULL;
  }

  int bits,mask;

  if (!PyArg_ParseTuple(args,"ii",&bits,&mask))
    return NULL;

  return Py_BuildValue("i",self->rx->set_dc_offset_cl_enable(bits,mask));
}

static PyMethodDef rx_methods[] = {
  //{namestring,(PyCFunction)func,METH_xxx,docstring}
  //METH_NOARGS, METH_VARARGS, METH_KEYWORDS
  {"set_decim_rate",(PyCFunction)rx_set_decim_rate,METH_VARARGS,"set decimator rate"},
  {"set_nchannels",(PyCFunction)rx_set_nchannels,METH_VARARGS,"set number of active channels"},
  {"set_mux",(PyCFunction)rx_set_mux,METH_VARARGS,"set input mux configuration"},
  {"rx_freq",(PyCFunction)rx_get_rx_freq,METH_VARARGS,"get channel's DDC center frequency"},
  {"set_rx_freq",(PyCFunction)rx_set_rx_freq,METH_VARARGS,"set channel's DDC center frequency"},
  {"set_fpga_mode",(PyCFunction)rx_set_fpga_mode,METH_VARARGS,"set fpga mode"},
  {"set_ddc_phase",(PyCFunction)rx_set_ddc_phase,METH_VARARGS,"set ddc phase register"},
  {"set_format",(PyCFunction)rx_set_format,METH_VARARGS|METH_KEYWORDS,"set rx data format"},
  {"start",(PyCFunction)rx_start,METH_NOARGS,"start receiver"},
  {"stop",(PyCFunction)rx_stop,METH_NOARGS,"stop receiver"},
  {"read",(PyCFunction)rx_read,METH_VARARGS,"read from usrp"},
  {"daughterboard_id",(PyCFunction)rx_daughterboard_id,METH_VARARGS,"get id of a daughterboard"},
  {"set_pga",(PyCFunction)rx_set_pga,METH_VARARGS,"set gain of PGA"},
  {"get_pga",(PyCFunction)rx_get_pga,METH_VARARGS,"get gain of PGA"},
  {"set_dc_offset_cl_enable",(PyCFunction)rx_set_dc_offset_cl_enable,METH_VARARGS,"enable/disable automatic DC offset removal control loop"},
  {NULL}  // sentinel
};

static PyTypeObject rx_type = {
    PyObject_HEAD_INIT(NULL)
    0,                  /*ob_size*/
    "usrp.usrp_rx",	/*tp_name*/
    sizeof(rx),		/*tp_basicsize*/
    0,                  /*tp_itemsize*/
    (destructor)rx_dealloc, /*tp_dealloc*/
    0,                  /*tp_print*/
    0,                  /*tp_getattr*/
    0,                  /*tp_setattr*/
    0,                  /*tp_compare*/
    0,                  /*tp_repr*/
    0,                  /*tp_as_number*/
    0,                  /*tp_as_sequence*/
    0,                  /*tp_as_mapping*/
    0,                  /*tp_hash */
    0,                  /*tp_call*/
    0,                  /*tp_str*/
    0,                  /*tp_getattro*/
    0,                  /*tp_setattro*/
    0,                  /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT, /*tp_flags*/
    "USRP receiver",    /* tp_doc */
    0,			/* tp_traverse */
    0,			/* tp_clear */
    0,			/* tp_richcompare */
    0,			/* tp_weaklistoffset */
    0,			/* tp_iter */
    0,			/* tp_iternext */
    rx_methods,		/* tp_methods */
    rx_members,		/* tp_members */
    rx_getseters,	/* tp_getset */
    0,			/* tp_base */
    0,			/* tp_dict */
    0,			/* tp_descr_get */
    0,			/* tp_descr_set */
    0,			/* tp_dictoffset */
    (initproc)rx_init,	/* tp_init */
    0,			/* tp_alloc */
    rx_new,		/* tp_new */
};

//////////////////////////////////////////////////
// usrp.tx: USRP transmitter object  (built on usrp_standard_tx in libusrp)
//
// tx = usrp.tx(interp_rate,which_board=0,nchan=2,mux=-1,
//              fusb_block_size=0,fusb_nblocks=0,fpga_filename="",firmware_filename="")
//   defaults used:
//   mux = (1 chan) 0x0098, (2 chan) 0xba98
//   firmware_filename = "/usr/local/share/usrp/rev4/std.ihx"
//   fpga_filename = "/usr/local/share/usrp/rev4/std_2rxhb_2tx.rbf"
//   fusb_block_size = 4096   # max is 16K
//   fusb_nblocks = 256       # FUSB_BUFFER_SIZE/FUSB_BLOCK_SIZE, BUFSIZE = 1M
//
// int = tx.start()   # start transmitter, bool success flag
// (int,int) = tx.write(string)    # returns number of bytes written,underrun flag
// int = tx.stop()  # stop transmitter, bool success flag
//
// int = tx.daughterboard_id(which_dboard)  # returns USRP_DBID_xxx
// int = tx.set_interp_rate(rate)  # set decimator rate: multiple of 4, 4..512, boolean success flag
// int = tx.interp_rate  # decimator rate
// int = tx.set_nchannels(nchannels)  # set number of active channels: 1 or 2, boolean success flag
// int = tx.nchannels  # number of active channels
// int = tx.set_mux(value)   # set output mux configuration, boolean success flag
// int = tx.mux    # output mux selection
// int = tx.set_tx_freq(channel,freq)  # set tx_freq [-44M..44M] channel [0..1], boolean success flag
// double = tx.tx_freq(channel)  # get tx_freq, channel is 0..1
// int = tx.set_pga(which,gain_in_dB)  # set programmable amplifier: which=0..3
// double = tx.set_pga(which)  # get programmable amplifier
// double = tx.pga_min # min legal PGA gain in dB (-20.0)
// double = tx.pga_max # max legal PGA gain in dB (0.0)
// double = tx.pga_db_per_step # hardware step size of PGA (20.0/255)
// int = tx.dac_rate # DAC sample rate
// string = tx.serial_number  # USRP board serial number
//////////////////////////////////////////////////

typedef struct {
  PyObject_HEAD
  usrp_standard_tx *tx;
} tx;

// properties defined in usrp_standard_tx

// int = tx.interp_rate  # interpolator rate
static PyObject *tx_get_interp_rate(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->interp_rate());
}

// int = tx.nchannels  # number of active channels
static PyObject *tx_get_nchannels(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->nchannels());
}

// int = tx.mux    # output mux selection
static PyObject *tx_get_mux(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->mux());
}

// double = tx.pga_min  # min legal PGA gain in dB
static PyObject *tx_get_pga_min(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("d",self->tx->pga_min());
}

// double = tx.pga_max  # max legal PGA gain in dB
static PyObject *tx_get_pga_max(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("d",self->tx->pga_max());
}

// double = tx.pga_db_per_step  # hardware PGA step size in dB
static PyObject *tx_get_pga_db_per_step(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("d",self->tx->pga_db_per_step());
}

// int = tx.block_size  # fusb block size
static PyObject *tx_get_block_size(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->block_size());
}

// int = tx.dac_rate  # sampling rate of D/A converter
static PyObject *tx_get_dac_rate(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->dac_rate());
}

// string = tx.serial_number  # USRP board serial number
static PyObject *tx_get_serial_number(tx *self,void *closure) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("s",self->tx->serial_number().c_str());
}

static PyGetSetDef tx_getseters[] = {
  {"interp_rate",(getter)tx_get_interp_rate,(setter)noset,"interpolator rate",NULL},
  {"nchannels",(getter)tx_get_nchannels,(setter)noset,"number of active channels",NULL},
  {"mux",(getter)tx_get_mux,(setter)noset,"output mux control",NULL},
  {"pga_min",(getter)tx_get_pga_min,(setter)noset,"min legal PGA gain in dB",NULL},
  {"pga_max",(getter)tx_get_pga_max,(setter)noset,"max legal PGA gain in dB",NULL},
  {"pga_db_per_step",(getter)tx_get_pga_db_per_step,(setter)noset,"PGA step size in dB",NULL},
  {"block_size",(getter)tx_get_block_size,(setter)noset,"fusb block size",NULL},
  {"dac_rate",(getter)tx_get_dac_rate,(setter)noset,"sampling rate of D/A converter",NULL},
  {"serial_number",(getter)tx_get_serial_number,(setter)noset,"USRP serial number",NULL},
  {NULL}  // sentinel
};

// tx = usrp.tx(interp_rate,which_board=0,nchan=2,mux=-1,
//              fusb_block_size=0,fusb_nblocks=0,fpga_filename="",firmware_filename="")
// defaults used:
//   mux = (1 chan) 0x0098 (2 chan) 0xba98
//   firmware_filename = /usr/local/share/usrp/rev4/std.ihx
//   fpga_filename = /usr/local/share/usrp/rev4/std_2rxhb_2tx.rbf
//   fusb_block_size = 4096   # max is 16K
//   fusb_nblocks = 256       # FUSB_BUFFER_SIZE/FUSB_BLOCK_SIZE, BUFSIZE = 1M
static PyObject *tx_new(PyTypeObject *type,PyObject *args,PyObject *kwds) {
  tx *self;

  self = (tx *)type->tp_alloc(type,0);
  if (self != NULL) {
    self->tx = NULL;   // do the heavy lifting in init
  }
  return (PyObject *)self;
}

// initialize tx object
static int tx_init(tx *self,PyObject *args,PyObject *kwds) {
  unsigned int interp_rate;	// required, even in range 8..256
  int which_board = 0;		// optional, 0 is what you usually want
  int nchan = 2;		// 1, 2 or 4
  int mux = -1;			// use default
  int fusb_block_size = 0;	// use default
  int fusb_nblocks = 0;		// use default
  char *fpga_filename = "";	// use default
  char *firmware_filename = ""; // use default

  static char *kwlist[] = {
    "interp_rate",
    "board",
    "nchan",
    "mux",
    "fusb_block_size",
    "fusb_nblocks",
    "fpga_filename",
    "firmware_filename",
    NULL};

  // interp_rate is required, rest of args are optional
  if (!PyArg_ParseTupleAndKeywords(args,kwds,"i|iiiiiss",kwlist,
				   &interp_rate,
				   &which_board,
				   &nchan,
				   &mux,
				   &fusb_block_size,
				   &fusb_nblocks,
				   &fpga_filename,
				   &firmware_filename))
    return -1;

  if ((interp_rate & 3) || interp_rate < 4 || interp_rate > 512) {
    PyErr_SetString(PyExc_ValueError,"intep_rate must be a multiple of 4 and in the range 4..512");
    return -1;
  }

  if (self->tx != NULL) delete self->tx;
  self->tx = usrp_standard_tx::make(which_board,interp_rate,nchan,mux,
				    fusb_block_size,fusb_nblocks,
				    fpga_filename,firmware_filename);

  if (self->tx == NULL) {
    PyErr_SetString(PyExc_IOError,"failed to create usrp_standard_tx output");
    return -1;
  }

  return 0;
}

// deallocate tx object
static void tx_dealloc(tx *self) {
  if (self->tx) delete self->tx;
  self->ob_type->tp_free((PyObject*)self);  // free the object
}

static PyMemberDef tx_members[] = {
  // members accesible from Python
  // {namestring,type,offsetof(tx,member),flags,docstring}
  // type: T_INT for integer, T_OBJECT_EX for PyObject *
  // flags: READONLY,READ_RESTRICTED,WRITE_RESTRICTED,RESTRICTED
  {NULL}  // sentinel
};

// methods defined in usrp_standard_tx

// int = tx.set_interp_rate(rate)  # set decimator rate: multiple of 4, 4..512
static PyObject *tx_set_interp_rate(tx *self,PyObject *args) {
  int interp_rate;

  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&interp_rate))
    return NULL;

  if ((interp_rate & 3) || interp_rate < 4 || interp_rate > 512) {
    PyErr_SetString(PyExc_ValueError,"interp_rate must be a multiple of 4 in the range 4..512");
    return NULL;
  }
    
  return Py_BuildValue("i",self->tx->set_interp_rate(interp_rate));
}

// int = tx.set_nchannels(nchannels)  # set number of active channels: 1 or 2, boolean success flag
static PyObject *tx_set_nchannels(tx *self,PyObject *args) {
  int nchannels; 

  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&nchannels))
    return NULL;
    
  if (nchannels != 1 && nchannels != 2) {
    PyErr_SetString(PyExc_ValueError,"The nchannels attribute must be 1 or 2");
    return NULL;
  }
    
  return Py_BuildValue("i",self->tx->set_nchannels(nchannels));
}

// int = tx.set_mux(value)   # set output mux configuration, boolean success flag
//   value is four 4-bit fields: DAC3 | DAC2 | DAC1 | DAC0
//   4-bit field is: ENNN
//   E is 1 if DAC is enable, NNN specifies which interpolator output connected to this DAC
//   NNN = 000   chan 0 I
//   NNN = 001   chan 0 Q
//   NNN = 010   chan 1 I
//   NNN = 011   chan 1 Q
static PyObject *tx_set_mux(tx *self,PyObject *args) {
  int mux; 

  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&mux))
    return NULL;
    
  // check NNN values
  for (int i = 0; i < 4; i += 1) {
    int v = (mux >> (i*4)) & 0x3;
    if (v > 3) {
      PyErr_SetString(PyExc_ValueError,"interpolator settings must be 0, 1, 2 or 3");
      return NULL;
    }
  }

  return Py_BuildValue("i",self->tx->set_mux(mux));
}

// int = tx.set_tx_freq(channel,freq)  # set tx_freq [-44M..44M] , boolean success flag
static PyObject *tx_set_tx_freq(tx *self,PyObject *args) {
  int channel;
  double freq;

  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"id",&channel,&freq))
    return NULL;

  if (channel < 0 || channel > 1) {
    PyErr_SetString(PyExc_ValueError,"channel must be 0 or 1");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->set_tx_freq(channel,freq));
}

// double = tx.get_tx_freq(channel)  # get tx_freq, channel is 0..3
static PyObject *tx_get_tx_freq(tx *self,PyObject *args) {
  int channel;

  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  if (!PyArg_ParseTuple(args,"i",&channel))
    return NULL;

  if (channel < 0 || channel > 2) {
    PyErr_SetString(PyExc_ValueError,"channel must be 0 or 1");
    return NULL;
  }

  return Py_BuildValue("d",self->tx->tx_freq(channel));
}

// bool = tx.start()   # start transmitter, bool success flag
static PyObject *tx_start(tx *self) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->start());
}

// bool = tx.stop()  # stop transmitter, bool success flag
static PyObject *tx_stop(tx *self) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->stop());
}

// methods defined in usrp_basic_tx

// (int,int) = tx.write(string)  # returns number of bytes sent, underrun flag
static PyObject *tx_write(tx *self,PyObject *args) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  char *buffer;
  int len;
  if (!PyArg_ParseTuple(args,"s#",&buffer,&len))
    return NULL;

  if ((len % 512) != 0) {
    PyErr_SetString(PyExc_ValueError,"transmit data must be a multiple of 512 bytes");
    return NULL;
  }

  bool underrun;
  int status;

  Py_BEGIN_ALLOW_THREADS
  status = self->tx->write(buffer,len,&underrun);
  Py_END_ALLOW_THREADS

  return Py_BuildValue("(i,i)",status,underrun);
}

// int = tx.daughterboard_id(which_dboard)  # returns USRP_DBID_xxx
static PyObject *tx_daughterboard_id(tx *self,PyObject *args) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  int which_dboard;

  if (!PyArg_ParseTuple(args,"i",&which_dboard))
    return NULL;

  int dbid = self->tx->daughterboard_id(which_dboard);

  return Py_BuildValue("i",dbid);
}

// int = tx.set_pga(which,gain_in_dB)  # set programmable amplifier: which=0..3, bool success flag
static PyObject *tx_set_pga(tx *self,PyObject *args) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  int which;
  double gain;

  if (!PyArg_ParseTuple(args,"id",&which,&gain))
    return NULL;

  if (which < 0 || which > 3) {
    PyErr_SetString(PyExc_ValueError,"which PGA must be in range 0..3");
    return NULL;
  }

  if (gain < self->tx->pga_min() || gain > self->tx->pga_max()) {
    PyErr_SetString(PyExc_ValueError,"PGA gain outside of allowable range (-20..0)");
    return NULL;
  }

  return Py_BuildValue("i",self->tx->set_pga(which,gain));
}

// double = tx.get_pga(which)  # get programmable gain amplifier gain
static PyObject *tx_get_pga(tx *self,PyObject *args) {
  if (!self->tx) {	// sanity check
    PyErr_SetString(PyExc_AttributeError,"No usrp transmitter info available");
    return NULL;
  }

  int which;

  if (!PyArg_ParseTuple(args,"i",&which))
    return NULL;

  if (which < 0 || which > 3) {
    PyErr_SetString(PyExc_ValueError,"which PGA must be in range 0..3");
    return NULL;
  }

  return Py_BuildValue("d",self->tx->pga(which));
}

static PyMethodDef tx_methods[] = {
  //{namestring,(PyCFunction)func,METH_xxx,docstring}
  //METH_NOARGS, METH_VARARGS, METH_KEYWORDS
  {"set_interp_rate",(PyCFunction)tx_set_interp_rate,METH_VARARGS,"set interpolator rate"},
  {"set_nchannels",(PyCFunction)tx_set_nchannels,METH_VARARGS,"set number of active channels"},
  {"set_mux",(PyCFunction)tx_set_mux,METH_VARARGS,"set output mux configuration"},
  {"tx_freq",(PyCFunction)tx_get_tx_freq,METH_VARARGS,"get channel's DUC center frequency"},
  {"set_tx_freq",(PyCFunction)tx_set_tx_freq,METH_VARARGS,"set channel's DUC center frequency"},
  {"start",(PyCFunction)tx_start,METH_NOARGS,"start transmitter"},
  {"stop",(PyCFunction)tx_stop,METH_NOARGS,"stop transmitter"},
  {"write",(PyCFunction)tx_write,METH_VARARGS,"write to usrp"},
  {"daughterboard_id",(PyCFunction)tx_daughterboard_id,METH_VARARGS,"get id of a daughterboard"},
  {"set_pga",(PyCFunction)tx_set_pga,METH_VARARGS,"set gain of PGA"},
  {"get_pga",(PyCFunction)tx_get_pga,METH_VARARGS,"get gain of PGA"},
  {NULL}  // sentinel
};

static PyTypeObject tx_type = {
    PyObject_HEAD_INIT(NULL)
    0,                  /*ob_size*/
    "usrp.usrp_tx",	/*tp_name*/
    sizeof(tx),		/*tp_basicsize*/
    0,                  /*tp_itemsize*/
    (destructor)tx_dealloc, /*tp_dealloc*/
    0,                  /*tp_print*/
    0,                  /*tp_getattr*/
    0,                  /*tp_setattr*/
    0,                  /*tp_compare*/
    0,                  /*tp_repr*/
    0,                  /*tp_as_number*/
    0,                  /*tp_as_sequence*/
    0,                  /*tp_as_mapping*/
    0,                  /*tp_hash */
    0,                  /*tp_call*/
    0,                  /*tp_str*/
    0,                  /*tp_getattro*/
    0,                  /*tp_setattro*/
    0,                  /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT, /*tp_flags*/
    "USRP transmitter", /* tp_doc */
    0,			/* tp_traverse */
    0,			/* tp_clear */
    0,			/* tp_richcompare */
    0,			/* tp_weaklistoffset */
    0,			/* tp_iter */
    0,			/* tp_iternext */
    tx_methods,		/* tp_methods */
    tx_members,		/* tp_members */
    tx_getseters,	/* tp_getset */
    0,			/* tp_base */
    0,			/* tp_dict */
    0,			/* tp_descr_get */
    0,			/* tp_descr_set */
    0,			/* tp_dictoffset */
    (initproc)tx_init,	/* tp_init */
    0,			/* tp_alloc */
    tx_new,		/* tp_new */
};

//////////////////////////////////////////////////
// module setup
//////////////////////////////////////////////////

static PyObject *find_device(PyObject *self,PyObject *args,PyObject *kwds) {
  int nth = 0;
  int fx2_ok_p = 0;

  static char *kwdlist[] = {"nth","fx2_ok_p",NULL};

  if (!PyArg_ParseTupleAndKeywords(args,kwds,"|ii",kwdlist,&nth,&fx2_ok_p))
    return NULL;

  return Py_BuildValue("i",usrp_find_device(nth,fx2_ok_p));
}

static PyMethodDef usrpmethods[] = {
  {"usrp_find_device",(PyCFunction)find_device,METH_VARARGS|METH_KEYWORDS,"locate nth USRP device"},
  {NULL,NULL,0,NULL}  // sentinel
};

PyMODINIT_FUNC initusrp(void) {
  PyObject *m;

  if (PyType_Ready(&rx_type) < 0) return;
  if (PyType_Ready(&tx_type) < 0) return;

  m = Py_InitModule("usrp",usrpmethods);
  if (m == NULL) return;

  Py_INCREF(&rx_type);
  PyModule_AddObject(m,"usrp_rx",(PyObject *)&rx_type);

  Py_INCREF(&tx_type);
  PyModule_AddObject(m,"usrp_tx",(PyObject *)&tx_type);

  // from usrp_slots.h
  PyModule_AddIntConstant(m,"SLOT_TX_A",SLOT_TX_A);
  PyModule_AddIntConstant(m,"SLOT_RX_A",SLOT_RX_A);
  PyModule_AddIntConstant(m,"SLOT_TX_B",SLOT_TX_B);
  PyModule_AddIntConstant(m,"SLOT_RX_B",SLOT_RX_B);

  // from usrp_standard.h
  PyModule_AddIntConstant(m,"FPGA_MODE_NORMAL",usrp_standard_rx::FPGA_MODE_NORMAL);
  PyModule_AddIntConstant(m,"FPGA_MODE_LOOPBACK",usrp_standard_rx::FPGA_MODE_LOOPBACK);
  PyModule_AddIntConstant(m,"FPGA_MODE_COUNTING",usrp_standard_rx::FPGA_MODE_COUNTING);
  PyModule_AddIntConstant(m,"FPGA_MODE_COUNTING_32BIT",usrp_standard_rx::FPGA_MODE_COUNTING_32BIT);

  PyModule_AddIntConstant(m,"CM_NEG_FDAC_OVER_4",usrp_standard_tx::CM_NEG_FDAC_OVER_4);
  PyModule_AddIntConstant(m,"CM_NEG_FDAC_OVER_8",usrp_standard_tx::CM_NEG_FDAC_OVER_8);
  PyModule_AddIntConstant(m,"CM_OFF",usrp_standard_tx::CM_OFF);
  PyModule_AddIntConstant(m,"CM_POS_FDAC_OVER_8",usrp_standard_tx::CM_POS_FDAC_OVER_8);
  PyModule_AddIntConstant(m,"CM_POS_FDAC_OVER_4",usrp_standard_tx::CM_POS_FDAC_OVER_4);

  // from usrp_dbid.h
  PyModule_AddIntConstant(m,"USRP_DBID_BASIC_TX",USRP_DBID_BASIC_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_BASIC_RX",USRP_DBID_BASIC_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_DBS_RX",USRP_DBID_DBS_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_TV_RX",USRP_DBID_TV_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_400_RX",USRP_DBID_FLEX_400_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_900_RX",USRP_DBID_FLEX_900_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1200_RX",USRP_DBID_FLEX_1200_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_2400_RX",USRP_DBID_FLEX_2400_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_400_TX",USRP_DBID_FLEX_400_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_900_TX",USRP_DBID_FLEX_900_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1200_TX",USRP_DBID_FLEX_1200_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_2400_TX",USRP_DBID_FLEX_2400_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_TV_RX_REV_2",USRP_DBID_TV_RX_REV_2);
  PyModule_AddIntConstant(m,"USRP_DBID_DBS_RX_REV_2_1",USRP_DBID_DBS_RX_REV_2_1);
  PyModule_AddIntConstant(m,"USRP_DBID_LF_TX",USRP_DBID_LF_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_LF_RX",USRP_DBID_LF_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_400_RX_MIMO_A",USRP_DBID_FLEX_400_RX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_900_RX_MIMO_A",USRP_DBID_FLEX_900_RX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1200_RX_MIMO_A",USRP_DBID_FLEX_1200_RX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_2400_RX_MIMO_A",USRP_DBID_FLEX_2400_RX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_400_TX_MIMO_A",USRP_DBID_FLEX_400_TX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_900_TX_MIMO_A",USRP_DBID_FLEX_900_TX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1200_TX_MIMO_A",USRP_DBID_FLEX_1200_TX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_2400_TX_MIMO_A",USRP_DBID_FLEX_2400_TX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_400_RX_MIMO_B",USRP_DBID_FLEX_400_RX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_900_RX_MIMO_B",USRP_DBID_FLEX_900_RX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1200_RX_MIMO_B",USRP_DBID_FLEX_1200_RX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_2400_RX_MIMO_B",USRP_DBID_FLEX_2400_RX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_400_TX_MIMO_B",USRP_DBID_FLEX_400_TX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_900_TX_MIMO_B",USRP_DBID_FLEX_900_TX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1200_TX_MIMO_B",USRP_DBID_FLEX_1200_TX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_2400_TX_MIMO_B",USRP_DBID_FLEX_2400_TX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1800_RX",USRP_DBID_FLEX_1800_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1800_TX",USRP_DBID_FLEX_1800_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1800_RX_MIMO_A",USRP_DBID_FLEX_1800_RX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1800_TX_MIMO_A",USRP_DBID_FLEX_1800_TX_MIMO_A);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1800_RX_MIMO_B",USRP_DBID_FLEX_1800_RX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_FLEX_1800_TX_MIMO_B",USRP_DBID_FLEX_1800_TX_MIMO_B);
  PyModule_AddIntConstant(m,"USRP_DBID_TV_RX_REV_3",USRP_DBID_TV_RX_REV_3);
  PyModule_AddIntConstant(m,"USRP_DBID_DTT754",USRP_DBID_DTT754);
  PyModule_AddIntConstant(m,"USRP_DBID_DTT768",USRP_DBID_DTT768);
  PyModule_AddIntConstant(m,"USRP_DBID_WBX_LO_TX",USRP_DBID_WBX_LO_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_WBX_LO_RX",USRP_DBID_WBX_LO_RX);
  PyModule_AddIntConstant(m,"USRP_DBID_EXPERIMENTAL_TX",USRP_DBID_EXPERIMENTAL_TX);
  PyModule_AddIntConstant(m,"USRP_DBID_EXPERIMENTAL_RX",USRP_DBID_EXPERIMENTAL_RX);
}
