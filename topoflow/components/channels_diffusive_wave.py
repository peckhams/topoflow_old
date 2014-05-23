
## Copyright (c) 2001-2013, Scott D. Peckham
##
## October 2012   (CSDMS Standard Names and BMI)
## January, May, July 2009
## May 2010 (changes to unit_test() and read_cfg_file()

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "diffusive wave" channel flow component
#          and related functions.  It inherits from the channels
#          "base class" in "channels_base.py".
#-----------------------------------------------------------------------
#
#  class channels_component
#
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (defined in channels_base.py)
#      get_output_var_names()    # (defined in channels_base.py)
#      get_var_name()            # (defined in channels_base.py)
#      get_var_units()           # (defined in channels_base.py)
#      ------------------------
#      update_velocity()
#
#-----------------------------------------------------------------------

import numpy as np

from topoflow.components import channels_base

#-----------------------------------------------------------------------
class channels_component(channels_base.channels_component):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'Channels_Diffusive_Wave',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #------------------------------------------------------
        'comp_name':          'ChannelsDiffWave',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Channels_Diffusive_Wave.cfg.in',
        'cfg_extension':      '_channels_diffusive_wave.cfg',
        'cmt_var_prefix':     '/ChannelsDiffWave/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Channels_Diffusive_Wave.xml',
        'dialog_title':       'Channels: Diffusive Wave Parameters',
        'time_units':         'seconds' }

    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        #-----------------------------------------------------------
        # This is done in channels_base.set_computed_input_vars()
        #-----------------------------------------------------------
        # self.KINEMATIC_WAVE = False 
        # self.DIFFUSIVE_WAVE = True 
        # self.DYNAMIC_WAVE   = False

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

    #   get_attribute() 
    #-------------------------------------------------------------------
    def update_velocity(self):

        #---------------------------------------------------------
        # Notes: Compute u from d and S_bed.  (7/13/05 version)

        #        cv        = channel vars structure
        #        cv.nval   = Manning's n values (grid)
        #        cv.z0_val = z0 roughness values (grid)
        #        cv.width  = channel bottom widths (grid)
        #        cv.angle  = channel bank angles (grid)

        #        Could use slopes in cv also, but S_bed has
        #        been modified from those values to impose a
        #        minimum slope that is nonzero.

        #        Rh = hydraulic radius (trapezoid here)

        #        S = S_bed  for KINEMATIC_WAVE option.
        #        S = S_free for DIFFUSIVE_WAVE_option.
        #---------------------------------------------------------

        #----------------------------
        # Update free surface slope
        #----------------------------
        self.update_free_surface_slope()   # (called S_free)

        #-------------------------------------------
        # Compute velocity using S_free vs. S_bed
        #-------------------------------------------
        # NB! This involves computing sqrt(S_free)
        # so disallow "backflow" this way ?
        #-------------------------------------------
        self.S_free = np.maximum(self.S_free, 0.0)        ############
        
        #------------------------
        # Use Manning's formula
        #------------------------
        if (self.MANNING):    
            self.u = self.manning_formula(self.S_free)
        
        #--------------------------------------
        # Use the Logarithmic Law of the Wall
        #--------------------------------------
        if (self.LAW_OF_WALL):    
            self.u = self.law_of_the_wall(self.S_free)

        #----------------------------------------
        # Allow negative velocity (backflow) ??
        # But to which child pixel ??
        #----------------------------------------
        # wn  = np.where( S_free < 0 )
        # nwn = np.size( wn[0] )
        # self.S_free = np.abs(self.S_free)
        # if (nwn > 0): self.u[wn] = -1.0 * self.u[wn]
        
    #    update_velocity()                       
    #-------------------------------------------------------------------


         
