# cjt (Nov. 08)

import blocks
import os,Queue,sys,threading,wx

# crossplatform way of determining number of cpus
# from http://cheeseshop.python.org/pypi/processing/0.34
def cpuCount():
    num = 0
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            pass
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            pass
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            pass

    if num >= 1:
        return num
    else:
        raise NotImplementedError

# The worker threads repeatedly take a block from the run queue and
# calls its _execute method.
class worker_process(threading.Thread):
    def __init__(self,queue,id):
        threading.Thread.__init__(self)
        self.queue = queue
        self.id = id

    def run(self):
        self.enabled = True   # might be set to false externally
        while (self.enabled):
            try:
                # set a .1 second timeout so we can check the self.enabled flag
                # and exit gracefully if simulation is finished
                b = self.queue.get(True,.1)
                #print "%s: running %s, qsize=%d" % (self.id,b._name,self.queue.qsize())
                b._execute()
            except Queue.Empty:
                pass

# testbench makes an instance of this frame which interleaves wx event
# processing with running the diagram blocks.  If nthreads > 1, additional
# background threads are started to also process blocks.
# NB. starting up too many background threads seems to starve the wx
# event processing and the GUI stops functioning.
# This frame also holds the various control widget created by the user.
class ControlFrame(wx.Frame):
    def __init__(self,sinks,queue,iterations=20,nthreads=cpuCount()):
        wx.Frame.__init__(self,None,title='Controls')
        self.sinks = sinks
        self.queue = queue
        self.iterations = iterations
        self.nthreads = nthreads
        self.workers = []
        self.enabled = False

        # a vertical stack of controls
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.ncontrols = 0

        self.Bind(wx.EVT_IDLE,self.OnIdle)

    def AddControl(self,label,control):
        self.ncontrols += 1
        self.sizer.Add(wx.StaticText(self,label=label),0,wx.LEFT|wx.TOP|wx.RIGHT,10)
        # controls expand horizontally but not vertically
        self.sizer.Add(control,0,wx.EXPAND|wx.LEFT|wx.RIGHT,10)

    def ShowIfControls(self,pos):
        if self.ncontrols > 0:
            self.sizer.Add((1,10))   # a bit of space at the bottom
            self.SetSizer(self.sizer)
            self.Fit()
            self.SetPosition(pos)
            self.Show(True)

    def StartProcessing(self):
        # we'll be one of the threads and interleave block processing
        # with wx event processing.  Create additional background
        # processing threads if user wants more than 1 thread.
        self.workers = [worker_process(self.queue,'thread %d'%i)
                        for i in xrange(1,self.nthreads)]

        # start up background threads and schedule ourselves on the event queue
        for worker in self.workers: worker.start()
        self.enabled = True

    def StopProcessing(self):
        # stop the worker threads
        for worker in self.workers: worker.enabled = False
        # and wait for them to exit
        for worker in self.workers: worker.join()
        self.workers = []
        self.enabled = False
       
        # for now, just exit when we're done processing
        sys.exit()

    def OnIdle(self,event):
        if self.enabled:
            # stop when all the sinks have reached end-of-stream
            if len(self.sinks) == 0:
                self.StopProcessing()
            
            # process several blocks and then return control to wx runtime
            for i in xrange(self.iterations):
                try:
                    b = self.queue.get(False)
                    #print "thread 0: running %s, qsize=%d" % (b._name,queue.qsize())
                    b._execute()
                except Queue.Empty:
                    break

        # ask for another IDLE event so we can run some more,
        # but only after GUI has a shot at processing events
        event.RequestMore()

# subclass this and supply a block_diagram() method that describes
# modules and their connections.  Make an instance of the subclass
# to instatiate the application and then call the .MainLooop() method
# of the instance.
class testbench(wx.App):
  def __init__(self,sample_rate=256000,buffer_size=None,nthreads=None,
               title="6.02 Testbench",nstatus=2):
      self.sample_rate = sample_rate
      self.buffer_size = buffer_size or 8192
      self.title = title
      self.nstatus = nstatus
      self.nthreads = nthreads or cpuCount()
      print 'using', self.nthreads, 'threads'
      wx.App.__init__(self,redirect=False)

  TOP = 50      # constants for positioning instrumentation windows
  RIGHT = 50
  GAP = 10

  def OnInit(self):
      blocks.current_testbench = self
      self.next_position = (self.RIGHT,self.TOP)   # default positioning of instruments
      self.maxw = 0

      self.run_queue = Queue.Queue(0)
      self.diagram_blocks = []
      self.sink_blocks = []

      # give controls a frame to live in.
      self.control_frame = ControlFrame(self.sink_blocks,self.run_queue,
                                        nthreads=self.nthreads)

      # let user add dsp blocks to the simulation
      self.block_diagram()

      assert len(self.diagram_blocks) > 0,"Oops, block_diagram() didn't define any blocks!"

      # show control frame if there are any controls.  Place it in its own column on the right.
      if self.next_position[1] != self.TOP:
          self.next_position = (self.next_position[0] + self.maxw + self.GAP,self.TOP)
          self.control_frame.ShowIfControls(self.next_position)

      # Reset diagram blocks to their starting state. Hopefully some blocks (eg, the
      # sources) will add themselves to the run queue...
      self.reset()
      self.control_frame.StartProcessing()

      return True

  # Reset diagram blocks to their starting state. Hopefully some blocks (eg, the
  # sources) will add themselves to the run queue...
  def reset(self):
    for b in self.diagram_blocks: b._reset()
    assert len(self.sink_blocks) > 0,\
           "Oops, block_diagram() didn't define any sinks!"
    assert self.run_queue.qsize() > 0,\
           "Oops, resetting block_diagram produced no runnable blocks"

  # position instrumentation windows in columns, starting at the upper left
  def get_position(self,wsize):
      dsize = wx.DisplaySize()   # size of display

      if self.next_position[1] + wsize[1] > dsize[1] and self.next_position[1] != self.TOP:
          # move to next column
          self.next_position = (self.next_position[0] + self.maxw + self.GAP,self.TOP)
          self.maxw = 0
      
      result = self.next_position
      self.next_position = (self.next_position[0],
                            self.next_position[1] + wsize[1] + self.GAP)
      self.maxw = max(self.maxw,wsize[0])
      return result

  def add_block(self,block):
      self.diagram_blocks.append(block)

  def add_sink(self,block):
      self.sink_blocks.append(block)

  def add_control(self,label,control):
      self.control_frame.AddControl(label,control)

  def run(self):
      self.MainLoop()

  # called by blocks when they reach end-of-stream and will no longer
  # be scheduled to run
  def reached_eos(self,block):
      # keep track of sinks that haven't reached eos!
      # Note: we can't simple call remove() since that uses == to
      # to test for equality and that operator is overridden for
      # dsp_blocks!
      for i in xrange(len(self.sink_blocks)):
          if block is self.sink_blocks[i]:
              del self.sink_blocks[i]
              break

  def block_diagram(self):
      assert False,'Oops, no block_diagram method defined'

__all__ = [
    'testbench'
    ]
