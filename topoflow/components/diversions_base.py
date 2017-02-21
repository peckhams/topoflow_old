
# (2/3/13) read_source_data() and read_sink_data() get "dt" from
#          source_file or sink_file vs. channels comp, but what
#          about canals ? Can't use "channel_dt" as before.
# (2/3/13) "vol@channel" is now obtained as a reference.

## NB! CFG file does not have "n_steps" so it is set
##     to 1 in BMI_base.py.  Maybe it should be
##     set to max of "nd_max" values that appear in
##     read_source_data(), read_sink_data(), etc. (11/14/11)
##============================================================
#
#  Copyright (c) 2001-2014, Scott D. Peckham
#
#  Feb. 2017. Changes to internal variable names.
#             Cleanup & testing with Test_Plane_Canal data.
#  Sep 2014.  New standard names and BMI updates and testing.
#  Nov 2013.  Converted TopoFlow to a Python package.
#  Feb 2013.  Adapted to use EMELI framework.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  May 2010.  Changes to initialize() and read_cfg_file().
#  Jul 2009.  Updates.
#  May 2009.  Updates.
#  Jan 2009.  Converted from IDL to Python with I2PY.
#
#---------------------------------------------------------------------
# Notes:  Maybe replace "dur_sums" approach with the same approach
#         now used in precip.py ??

#         Make sure volume in caller gets updated correctly.
#---------------------------------------------------------------------
#
#  class diversions_component:  (inherits from BMI_base.py)
# 
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#---------------------------------
#      read_input_files()
#      read_source_data()
#      read_sink_data()
#      read_canal_data()
#----------------------------
#      update_sources()
#      update_sinks()
#      update_canals()
#
#--------------------------------------------------------------------

import numpy as np
import glob
import os

from topoflow.utils import BMI_base
from topoflow.utils import cfg_files as cfg

#---------------------------------------------------------------------
class diversions_component( BMI_base.BMI_component ):
 
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print ' '
            print 'Diversions component: Initializing...'
        
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file 
        #-----------------------------------------------
        # self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()  # (need this, 5/19/10)
        self.initialize_basin_vars()  # (5/14/10)
        
        ############################################################   
        # With new framework approach, we can't request the time
        # step for a specific component as "dt" (due to conflicts).
        ############################################################   
        # The Diversions component "dt" should probably match that
        # of the Channels component.  (Must be named "dt" vs.
        # "canal_dt".)
        ############################################################
        self.initialize_time_vars()
        
        #-----------------------------------------------------
        # These are used by the Channels component (9/22/14)
        # to tell if diversions are available and "on".
        #-----------------------------------------------------
        if not(self.use_canals):
            self.n_canals = self.initialize_scalar( 0, dtype='int32')
        if not(self.use_sinks):
            self.n_sinks = self.initialize_scalar( 0, dtype='int32')
        if not(self.use_sources):
            self.n_sources = self.initialize_scalar( 0, dtype='int32')

        #----------------------------------
        # Return if component is disabled
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print 'Diversions component: Disabled in CFG file.'
            self.n_canals          = self.initialize_scalar( 0, dtype='int32')
            self.n_sinks           = self.initialize_scalar( 0, dtype='int32')
            self.n_sources         = self.initialize_scalar( 0, dtype='int32')
            #-----------------------------------------------------------------------
            self.sinks_x           = self.initialize_scalar( 0, dtype='float64')
            self.sinks_y           = self.initialize_scalar( 0, dtype='float64')
            self.sinks_Q           = self.initialize_scalar( 0, dtype='float64')
            #-----------------------------------------------------------------------
            self.sources_x         = self.initialize_scalar( 0, dtype='float64')
            self.sources_y         = self.initialize_scalar( 0, dtype='float64')
            self.sources_Q         = self.initialize_scalar( 0, dtype='float64')
            #-----------------------------------------------------------------------
            self.canals_in_x       = self.initialize_scalar( 0, dtype='float64')
            self.canals_in_y       = self.initialize_scalar( 0, dtype='float64')
            self.canals_in_Q_fraction = self.initialize_scalar( 0, dtype='float64')
            self.canals_out_Q      = self.initialize_scalar( 0, dtype='float64')
            self.canals_out_x      = self.initialize_scalar( 0, dtype='float64')
            self.canals_out_y      = self.initialize_scalar( 0, dtype='float64')
            #-----------------------------------------------------------------------
            self.DONE   = True
            self.status = 'initialized'  # (OpenMI 2.0 convention)
            return

        #----------------------------------------
        # Initialize all time-related variables
        #----------------------------------------
        self.initialize_time_vars()
        
        #-----------------------------------------------
        # Read from files as needed to initialize vars 
        #-----------------------------------------------
        # source_files and sink_files have their own "dt"
        # which will override the default above. (2/3/13)
        #--------------------------------------------------
        self.read_input_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #--------------------------------------------------------------------------
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        if (self.comp_status == 'Disabled'):
            return
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #-----------------------------------------------------
        # Update info from all sources, sinks and diversions
        #-----------------------------------------------------
        # print '### Calling update_sources()...'
        self.update_sources()
        # print '### Calling update_sinks()...'
        self.update_sinks()
        # print '### Calling update_canals()...'
        self.update_canals()

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time()
        self.status = 'updated'  # (OpenMI 2.0 convention)
        
    #   update()
    #--------------------------------------------------------------------------
    def finalize(self):
        
        self.status = 'finalizing'  # (OpenMI)
##        if (self.comp_status == 'Enabled'):
##            self.close_input_files()
##            self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Diversions component')
        ## print 'Diversions component: Finished.'
        
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #---------------------------------------------------------------    
        # Note: The initialize() method calls initialize_config_vars()
        #       (in BMI_base.py), which calls this method at the end.
        #       But read_input_files() has not been called yet.
        #       See initialize_computed_vars().
        #--------------------------------------------------------------
        self.use_sources = (self.use_sources == 'Yes')
        self.use_sinks   = (self.use_sinks   == 'Yes')
        self.use_canals  = (self.use_canals  == 'Yes')

    #   set_computed_input_vars()
    #--------------------------------------------------------------------------
    def read_input_files(self):
    
        self.source_file = self.in_directory + self.source_file
        self.sink_file   = self.in_directory + self.sink_file
        self.canal_file  = self.in_directory + self.canal_file
 
        self.read_source_data()
        self.read_sink_data()
        self.read_canal_data()

    #   read_input_files()
    #--------------------------------------------------------------------------                
    def read_source_data(self):
        
        pass  # (See diversions_fraction_method.py)
   
    #   read_source_data()
    #--------------------------------------------------------------------------
    def read_sink_data(self):

        pass  # (See diversions_fraction_method.py)
   
    #   read_sink_data()

    #--------------------------------------------------------------------------
    def read_canal_data(self):

        pass  # (See diversions_fraction_method.py)

    #   read_canal_data()
    #--------------------------------------------------------------------------
#     def update_sources(self):
# 
#         if not(self.use_sources): return     
#         # n = size(self.source_IDs)
#         n = self.source_IDs.size
#         
#         for k in xrange(n):
#             dvol = ( self.dt * self.sources_Q[self.time_index, k] )
#             self.vol.flat[self.source_IDs[k]] += dvol
#         
#     #   update_sources()
#     #--------------------------------------------------------------------------
#     def update_sinks(self):
# 
#         if not(self.use_sinks): return     
#         # n = size(self.sink_IDs)
#         n = self.sink_IDs.size
# 
#         for k in xrange(n):
#             dvol = (self.dt * self.sinks_Q[self.time_index, k])
#             self.vol.flat[self.sink_IDs[k]] -= dvol
#             self.vol = np.maximum(self.vol, 0.0)
#                   
#     #   update_sinks()
#     #--------------------------------------------------------------------------
#     def update_canals(self):
# 
#         #------------------------------------------------------------------
#         # Note: Q_canals is same at upstream and downstream ends, but the
#         #       downstream end lags the upstream end by the travel time
#         #       from in_ID to out_ID.  As a result, the duration and Q
#         #       vector for the downstream end are computed from those of
#         #       the upstream end, and the travel time, td, as:
#         #           Q_out   = [0,  Q_in]
#         #           dur_out = [td, dur_in]
#         #           dur_sum_out = [0, dur_sum_in] + td
#         #------------------------------------------------------------------
#         if not(self.use_canals): return
#         n = np.size( self.canal_in_IDs )
# 
#         #--------------------------------
#         # Process upstream IDs as sinks
#         #--------------------------------
#         for k in xrange(n):   
#             dvol = (self.dt * self.Q_canals_in[self.time_index, k])
#             self.vol.flat[self.canal_in_IDs[k]] -= dvol
#             self.vol = np.maximum( self.vol, 0.0 )
#         
#         #------------------------------------
#         # Process downstream IDs as sources
#         # Must account for time lag in Q.
#         #------------------------------------
#         for k in xrange(n):
#             #--------------------------------------------
#             # Compute canals_out_Q & dur_sum_canals_out
#             #--------------------------------------------
#             dur_sum_canals_out = np.array([0.0, self.dur_sum_canals_in[:,k]]) + self.canal_t_vals[k]
#             canals_out_Q       = np.array([0.0, self.Q_canals_in[:,k]])
#   
#             dvol = (self.dt * canals_out_Q[self.time_index])
#             self.vol[self.canal_out_IDs[k]] += dvol
# 
#     #   update_canals()
    #--------------------------------------------------------------------------


    

