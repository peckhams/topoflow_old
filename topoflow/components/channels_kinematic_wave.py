
## Copyright (c) 2001-2010, Scott D. Peckham
##
## October 2012   (CSDMS Standard Names and BMI)
## January, May, July 2009
## May 2010 (changes to unit_test() and read_cfg_file()

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "kinematic wave" channel flow component
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

import numpy

from topoflow.components import channels_base

#-----------------------------------------------------------------------
class channels_component(channels_base.channels_component):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'Channels_Kinematic_Wave',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #------------------------------------------------------
        'comp_name':          'ChannelsKinWave',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Channels_Kinematic_Wave.cfg.in',
        'cfg_extension':      '_channels_kinematic_wave.cfg',
        'cmt_var_prefix':     '/ChannelsKinWave/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Channels_Kinematic_Wave.xml',
        'dialog_title':       'Channels: Kinematic Wave Parameters',
        'time_units':         'seconds' }
    
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        #-----------------------------------------------------------
        # This is done in channels_base.set_computed_input_vars()
        #-----------------------------------------------------------
        # self.KINEMATIC_WAVE = True
        # self.DIFFUSIVE_WAVE = False
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

        #        nval   = Manning's n values (grid)
        #        z0_val = z0 roughness values (grid)
        #        width  = channel bottom widths (grid)
        #        angle  = channel bank angles (grid)

        #        Could use slopes in cp also, but S_bed has
        #        been modified from those values to impose a
        #        minimum slope that is nonzero.

        #        Rh = hydraulic radius (trapezoid here)

        #        S = S_bed  for KINEMATIC_WAVE option.
        #        S = S_free for DIFFUSIVE_WAVE_option.
        #---------------------------------------------------------
            
        #------------------------
        # Use Manning's formula
        #------------------------
        if (self.MANNING):    
            self.u = self.manning_formula(self.S_bed)
        
        #--------------------------------------
        # Use the Logarithmic Law of the Wall
        #--------------------------------------
        if (self.LAW_OF_WALL):    
            self.u = self.law_of_the_wall(self.S_bed)

        # print '(umin, umax) =', self.u.min(), self.u.max()
        
    #   update_velocity()                       
    #-------------------------------------------------------------------


         
