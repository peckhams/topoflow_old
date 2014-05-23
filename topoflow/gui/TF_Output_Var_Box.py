#!/usr/bin/env python

#  April 30, 2009
#  S.D. Peckham

import wx

#-------------------------------------------------------------
#   class TF_Output_Var_Box
#         __init__()
#         Timestep_Box()
#         On_Check_Box()

#-------------------------------------------------------------
class TF_Output_Var_Box(wx.Panel):

    #-----------------------------------------------------
    # Notes:  This class is for creating "output variable
    #         panels" within TopoFlow output dialogs.
    #         Initial settings, labels, etc. are read
    #         from the XML file provided
    #         Default is wx.ALIGN_LEFT for labels & text
    #-----------------------------------------------------
    def __init__(self, parent=None, id=-1, \
                 main_frame=None, data=None):
 
        wx.Panel.__init__(self, parent)
        self.SetBackgroundColour('Light Blue')
        
        #------------------------------------------------
        # Saving parent allows collected values to be
        # stored in parent frame before this one closes.
        #------------------------------------------------
        self.main_frame = main_frame
        self.parent     = parent
        self.data       = data       ###########
        #--------------------------
        self.vgap       = 10
        self.hgap       = 6
        self.text_width = 240
        
        #-------------------------------------------
        #  Create sizer box for all input var info
        #  (This provides a frame and "box label".)
        #-------------------------------------------
        box_label = "Output variables to save:"
        vbox  = wx.StaticBox(self, -1, box_label)
        sizer = wx.StaticBoxSizer(vbox, wx.VERTICAL)
        
        #---------------------------------------------
        #  Create another sizer box for rows of info
        #---------------------------------------------
        #  Use "vgap=0" for most compact dialog
        #---------------------------------------------
        header = ["Variable:", "Output filename:", "Units:"]
        nh = len(header)
        fg_sizer = wx.FlexGridSizer(cols=nh, hgap=self.hgap, vgap=0)

        # Should this be self.text_boxes ??
        ####################################        
        parent.text_boxes = []
        ####################################
        
        #-------------------------
        # Add the column headers
        #-------------------------
        for row in range(nh):
            L1 = wx.StaticText(self, -1, header[row])
            fg_sizer.Add(L1, 0, wx.ALL, self.hgap)

        if (data == None):
            msg = 'No data passed to TF_Output_Var_Box!'
            wx.MesageBox(msg, caption='SORRY,')
            return  #################################
        
        #-----------------------------------------------
        # Create a row in the dialog for each variable
        # using labels, types, values and units.
        #-----------------------------------------------
        for row in range(len(data.var_names)):      
            row_ID = (5000 + row)   #####
            ### row_ID = 'tf_row_' + str(row)  # (ID can't be string.)
            cbox  = wx.CheckBox(self, row_ID, data.var_names[row])
            ## check some boxes by default ???
            self.Bind(wx.EVT_CHECKBOX, self.On_Check_Box, cbox)      
            #----------------------------------------------------------            
            text  = wx.TextCtrl(self, -1, data.var_values[row],
                                size=(self.text_width,-1))
            ####################################
            parent.text_boxes.append( text )
            ####################################
            #----------------------------------------------------------
            units_str = "[" + data.var_units[row] + "]"
            # units_str = data.var_units[row]
            ustr  = wx.StaticText(self, -1, units_str)
            #----------------------------------------------------------
            fg_sizer.Add(cbox,  1, wx.EXPAND|wx.ALL, self.hgap)
            fg_sizer.Add(text,  1, wx.EXPAND|wx.ALL, self.hgap)
            fg_sizer.Add(ustr,  1, wx.EXPAND|wx.ALL, self.hgap)
        
        #---------------------------------
        # Add fg_sizer to the main sizer
        #---------------------------------
        sizer.Add(fg_sizer, 1, wx.EXPAND|wx.ALL, 5)

        #----------------------------------------
        # Add a sampling time box to main sizer
        #----------------------------------------
        pad_row1 = wx.StaticText(self, -1, " ")
        time_box = self.Timestep_Box()
        #--------------------------------------------
        # NB! Don't change proportion to 1 here !!
        #--------------------------------------------
        proportion = 0
        sizer.Add(pad_row1, proportion, wx.EXPAND|wx.ALL, 5)
        sizer.Add(time_box, proportion, wx.EXPAND|wx.ALL, 5)

        #--------------------
        # Set sizer and fit
        #--------------------
        self.SetSizer(sizer)
        
    #   __init__()
    #----------------------------------------------------------------
    def Timestep_Box(self):

        #---------------------------------------------
        #  Create sizer box for the process timestep
        #---------------------------------------------
        timestep = self.data.timestep
        unit_str  = "[" + timestep.units + "]"
        L1    = wx.StaticText(self, -1, timestep.label + ":")
        text  = wx.TextCtrl(self,   -1, timestep.value)
        L2    = wx.StaticText(self, -1, unit_str)
        #-------------------------------------------------------        
        box   = wx.BoxSizer(wx.HORIZONTAL | wx.ALIGN_LEFT)
        proportion = 0  # (do not use 1)
        box.Add((self.hgap, self.hgap), proportion)
        box.Add(L1)
        box.Add((self.hgap, self.hgap), proportion)
        box.Add(text)
        box.Add((self.hgap, self.hgap), proportion)
        box.Add(L2)

        return box
    
    #   Timestep_Box()
    #----------------------------------------------------------------
    def On_Check_Box(self, event):

        #########################################################
        # We need the state of the "top-most" parent here !!
        #########################################################
        if (self.main_frame == None):
            print 'Sorry:  The TopoFlow main dialog is missing,'
            print 'so checkbox choices cannot be saved.'
            print ' '
            return
        
        #-----------------------------------------
        #  Event handler for the Type droplists.
        #-----------------------------------------
        print 'IsChecked  =', event.IsChecked()
        # print 'event.Id      =', event.Id
        # print 'event.GetId() =', event.GetId()

        #---------------------------------------
        # Need to know the row of the droplist
        #---------------------------------------
        row = (event.Id - 5000)
        # Idea:  row_ID = 'tf_row_' + str(row)
        # row = int(event.Id[7:])  # (ID can't be string ?)
        print 'row =', row
        print ' '

        #----------------------------------
        # Change choice in internal state
        #----------------------------------
        self.data.out_choices[row] = event.IsChecked()

        #--------------------------------
        # Save choice in parent's state
        #--------------------------------
        # self.main_frame....
        
    #   On_Check_Box()
    #----------------------------------------------------------------

        
