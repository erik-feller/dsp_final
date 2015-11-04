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

import blocks
import wx

class radio_buttons(wx.Panel):
    def __init__(self,value_list,callback,value=None,
                 title='radio_buttons',position=None,size=wx.DefaultSize,
                 orientation=wx.HORIZONTAL):
        wx.Panel.__init__(self,blocks.current_testbench.control_frame,style=wx.RAISED_BORDER)
        self.value_list = value_list
        self.callback = callback

        sizer = wx.BoxSizer(orientation)
        self.SetSizer(sizer)
        style = wx.RB_GROUP
        if value is None: value = value_list[0][0]
        for i in xrange(len(value_list)):
            label,v = value_list[i]
            button = wx.RadioButton(self,-1,label,style=style)
            button.SetValue(False)
            button.myindex = i
            style = 0  # subsequent buttons belong to the same group as the first
            if value == label or value == v:
                button.SetValue(True)
                self.callback(v)
            sizer.Add(button,0,wx.RIGHT,10)
            self.Bind(wx.EVT_RADIOBUTTON,self.OnRadio,button)

        blocks.current_testbench.add_control(title,self)

    def OnRadio(self,event):
        button = event.GetEventObject()
        self.callback(self.value_list[button.myindex][1])

class slider(wx.Panel):
    def __init__(self,range,callback,
                 title='slider',
                 quantization=None,scale=1,offset=0,
                 value=None,
                 orientation=wx.SL_HORIZONTAL,   # or wx.SL_VERTICAL
                 autoticks=wx.SL_AUTOTICKS,      # or 0
                 labels=wx.SL_LABELS,            # of 0
                 ):
        wx.Panel.__init__(self,blocks.current_testbench.control_frame,style=wx.RAISED_BORDER)
        self.callback = callback
        self.quantization = quantization
        self.scale = scale
        self.offset = offset

        vsizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(vsizer)
        self.slider = wx.Slider(self,minValue=range[0],maxValue=range[1],
                           style=orientation|autoticks|labels)
        vsizer.Add(self.slider,flag=wx.EXPAND|wx.ALIGN_CENTER_VERTICAL,proportion=1)
        self.Bind(wx.EVT_SLIDER,self.OnSlider,self.slider)

        if value:
            value = min(max(value,range[0]),range[1])
        else:
            value = range[0]
        self.SetValue(value)

        blocks.current_testbench.add_control(title,self)

    def SetValue(self,value):
        if self.quantization:
            value = round(float(value)/self.quantization)*self.quantization
        self.slider.SetValue(value)
        self.callback(value*self.scale + self.offset)

    def OnSlider(self,event):
        self.SetValue(self.slider.GetValue())

"""
class spinner(wx.Panel):
    def __init__(self,range,callback,
                 quantization=None,
                 value=None
                 ):
        wx.Panel.__init__(self,blocks.current_testbench.control_frame,style=wx.RAISED_BORDER)
        if value:
            value = min(max(value,range[0]),range[1])
        else:
            value = range[0]
        self.callback = callback
        self.quantization = quantization

        vsizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(vsizer)
        self.spinner = wx.SpinCtrl(self,initial=value,min=range[0],max=range[1],
                                   style=wx.SP_ARROW_KEYS)
        vsizer.Add(self.spinner,flag=wx.EXPAND|wx.ALIGN_CENTER_VERTICAL,proportion=1)
        self.Bind(wx.EVT_SPINCTRL,self.OnSpinner,self.spinner)

    def OnSpinner(self,event):
        value = self.spinner.GetValue()
        if self.quantization:
            value = round(float(value)/self.quantization)*self.quantization
            self.spinner.SetValue(value)
        self.callback(value)
"""
