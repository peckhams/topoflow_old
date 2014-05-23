
#  April 24, 2009
#  S.D. Peckham

import wx

#----------------------------------------------------------------
class Plot_Window(wx.Window):
    
    def __init__(self, parent, title, labels):
        wx.Window.__init__(self, parent)
        self.title  = title
        self.labels = labels
        self.data   = [0.0] * len(labels)

        self.Init_Buffer()

        self.Bind(wx.EVT_SIZE,  self.On_Resize)
        self.Bind(wx.EVT_PAINT, self.On_Paint)

    #   __init__()
    #------------------------------------------------------------
    def On_Resize(self, event)
    
        self.InitBuffer()  # (need a new buffer)

    #   On_Resize()
    #------------------------------------------------------------
    def On_Paint(self, event)
    
        dc = wx.BufferedPaintDC(self, self.buffer)

    #   On_Paint()
    #------------------------------------------------------------
    def Init_Buffer(self)
    
        w, h = self.GetClientSize()
        self.buffer = wx.EmptyBitmap(w, h)
        dc = wx.BufferedDC(wx.ClientDC(self), self.buffer)
        self.Draw_Graph(dc)

    #   Init_Buffer()
    #------------------------------------------------------------
    def Get_Data(self)
    
        return self.data

    #   Get_Data()
    #------------------------------------------------------------





    
