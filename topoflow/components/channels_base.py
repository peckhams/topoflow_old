
## See "d_bankfull" in update_flow_depth()  ######## (2/21/13)

## See "(5/13/10)" for a temporary fix.

## Copyright (c) 2001-2013, Scott D. Peckham
##
## Jan 2013. Shared scalar doubles are now 0D numpy arrays.
##           This makes them mutable and allows components with
##           a reference to them to see them change.
##           So far:  Q_outlet, Q_peak, Q_min...
##
## Jan 2013. Revised handling of input/output names.
##
## Oct 2012. CSDMS Standard Names and BMI)
##
## May 2010. Changes to initialize() and read_cfg_file()
##
## May 2012. Commented out diversions.update() for now.  #######
##
## May 2012. Shared scalar doubles are now 1-element 1D numpy arrays.
##           This makes them mutable and allows components with
##           a reference to them to see them change.
##           So far:  Q_outlet, Q_peak, Q_min...
##
## Mar 2010. Changed codes to code, widths to width,
##           angles to angle, nvals to nval, z0vals to z0val,
##           slopes to slope (for GUI tools and consistency
##           across all process components)
##
## May, July, August 2009
##
## January 2009   (converted from IDL)

#-----------------------------------------------------------------------
#  NB!     In the CFG file, change MANNING and LAW_OF_WALL flags to
#          a single string entry like "friction method".   #########
#-----------------------------------------------------------------------
#  Notes:  Set self.u in manning and law_of_wall functions ??
#          Update friction factor in manning() and law_of_wall() ?
#          Double check how Rh is used in law_of_the_wall().

#          d8_flow has "flow_grids", but this one has "codes".
#          Make sure values are not stored twice.
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "base class" for channelized flow
#          components as well as functions used by most or
#          all channel flow methods.  The methods of this class
#          (especially "update_velocity") should be over-ridden as
#          necessary for different methods of modeling channelized
#          flow.  See channels_kinematic_wave.py,
#          channels_diffusive_wave.py and channels_dynamic_wave.py.
#-----------------------------------------------------------------------
#
#  class channels_component
#
#      ## get_attribute()        # (defined in each channel component)
#      get_input_var_names()     # (5/15/12)
#      get_output_var_names()    # (5/15/12)
#      get_var_name()            # (5/15/12)
#      get_var_units()           # (5/15/12)
#-----------------------------
#      set_constants()
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()   # (5/11/10)
#-----------------------------
#      CCA Port Related
#-----------------------------
#      get_cca_port_info()         # (OBSOLETE)
#      embed_child_components()    # (OBSOLETE)
#      add_child_ports()           # (OBSOLETE)
#      initialize_ports()          # (OBSOLETE)
#----------------------------------
#      initialize_d8_vars()          ########
#      initialize_computed_vars()
#      initialize_outlet_values()
#      initialize_peak_values()
#      initialize_min_and_max_values()  # (2/3/13)
#-------------------------------------
#      update_R()
#      update_R_integral()
#      update_discharge()
#      update_flow_volume()
#      update_flow_depth()
#      update_free_surface_slope()
#      update_trapezoid_Rh()
#      update_velocity()            # (override as needed)
#      update_velocity_on_edges()
#      update_outlet_values()
#      update_peak_values()         # (at the main outlet)
#      update_Q_out_integral()      # (moved here from basins.py)
#      update_mins_and_maxes()      # (don't add into update())
#      check_flow_depth()
#      check_flow_velocity()
#----------------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#----------------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#----------------------------------
#      manning_formula()
#      law_of_the_wall()
#      print_status_report()
#      remove_bad_slopes() 

#  Functions:               # (stand-alone versions of these)
#      Trapezoid_Rh()
#      Manning_Formula()
#      Law_of_the_Wall()
    
#-----------------------------------------------------------------------

import numpy as np
import os, os.path

from topoflow.utils import BMI_base

# from topoflow.utils import d8_base

from topoflow.utils import file_utils  ###
from topoflow.utils import model_input
from topoflow.utils import model_output
from topoflow.utils import ncgs_files  ###
from topoflow.utils import ncts_files  ###
from topoflow.utils import rtg_files   ###
from topoflow.utils import text_ts_files   ###
from topoflow.utils import tf_d8_base as d8_base
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
class channels_component( BMI_base.BMI_component ):

    #-----------------------------------------------------------
    # Notes: This is first called without the prefix argument.
    #        The grid filenames below are set using the data
    #        prefix by the Read_Run_Info_Panel routine.

    # NB!    (3/16/07) dt is now set to zero initially, so
    #        that Read_Run_Info_Panel can determine if user
    #        has set the value (e.g. by loading a saved value).
    #        If unset, then it is set by Read_Run_Info_Panel.
    #        If GUI is not used, it is set in the input file.
    #-----------------------------------------------------------
    # Should we use "land_snow" or just "snow" (current).
    #------------------------------------------------------------
    # For GW, we currently use:
    #   land_water__baseflow_emergence_rate
    # but should we use:
    #   land_subsurface_to_surface_water__baseflow_seepage_rate
    #------------------------------------------------------------    
    _input_var_names = [
        'atmosphere_water__liquid_equivalent_precipitation_rate', # (P)
        'glacier__melt_rate',                           # (MR)
        ## 'land_surface__elevation',
        ## 'land_surface__slope',
        ## 'land_subsurface_to_surface_water__baseflow_seepage_rate', # GW   
        'land_water__baseflow_emergence_rate',              # (GW)
        'land_water__evaporation_rate',                     # (ET)
        'land_water__infiltration_rate',                    # (IN)
        'snow__melt_rate' ]                                 # (SM)
        
    #----------------------------------
    # Maybe add these out_vars later.
    #----------------------------------
    #  ['time_sec', 'time_min' ]
    
    _output_var_names = [
        'channel_bed__manning_coefficient',                    # nval
        'channel_bed__max_over_domain_of_manning_coefficient', # nval_max
        'channel_bed__max_over_domain_of_roughness_length',    # z0val_max
        'channel_bed__min_over_domain_of_manning_coefficient', # nval_min
        'channel_bed__min_over_domain_of_roughness_length',    # z0val_min
        'channel_bed__roughness_length',                       # z0val
        'channel_bed_surface__slope',                          # S_bed
        'channel_centerline__straight_sinuosity',              # sinu
        'channel_cross_section__hydraulic_radius',             # Rh
        'channel_cross_section_trapezoid__bank_angle',         # angle    ####
        'channel_cross_section_trapezoid__bottom_width',       # width    ####
##        'channel_cross_section_water__depth',                  # d
##        'channel_cross_section_water__initial_depth',          # d0       ####
##        'channel_cross_section_water__speed',                  # u
        'channel_model__time_step',                            # dt
        'channel_water__depth',                                # d
        'channel_water__initial_depth',                        # d0       ####
        'channel_water__speed',                                # u
        'channel_water__discharge',                            # Q
        'channel_water__friction_factor',                      # f
        # 'channel_water__froude_number',
        'channel_water__volume',                               # vol
        'channel_water_model__time_step',                      # dt
        'channel_water_surface__slope',                        # S_free
        'model_grid_cell__area',                               # da
        'model__time_step',                                    # dt
        'land_water__runoff_rate',                                 # R
        'watershed_outlet_water__depth',                           # d_outlet
        'watershed_outlet_water__discharge',                       # Q_outlet
        'watershed_outlet_water__friction_factor',                 # f_outlet
        'watershed_outlet_water__time_integral_of_discharge',      # vol_Q
        'watershed_outlet_water__max_over_time_of_depth',          # d_peak
        'watershed_outlet_water__max_over_time_of_discharge',      # Q_peak
        'watershed_outlet_water__max_over_time_of_speed',          # u_peak
        'watershed_outlet_water__speed',                           # u_outlet
        'watershed_outlet_water__time_of_max_of_depth',            # Td_peak
        'watershed_outlet_water__time_of_max_of_discharge',        # T_peak
        'watershed_outlet_water__time_of_max_of_speed',            # Tu_peak
        'watershed_water__area_time_integral_of_runoff_rate',      # vol_R
        #-----------------------------------------------------
        # These might only be available at the end of run ??
        #-----------------------------------------------------
        'watershed_water__max_over_domain_of_depth',        # d_max
        'watershed_water__max_over_domain_of_discharge',    # Q_max
        'watershed_water__max_over_domain_of_speed',        # u_max
        'watershed_water__min_over_domain_of_depth',        # d_min
        'watershed_water__min_over_domain_of_discharge',    # Q_min
        'watershed_water__min_over_domain_of_speed' ]       # u_min   
        
    _var_name_map = {
        'atmosphere_water__liquid_equivalent_precipitation_rate': 'P',
        'glacier__melt_rate':                                  'MR',
        ## 'land_surface__elevation':                          'DEM',
        ## 'land_surface__slope':                              'S_bed',
        'land_water__baseflow_emergence_rate':                 'GW',
        'land_water__evaporation_rate':                        'ET',
        'land_water__infiltration_rate':                       'IN',
        'land_water__runoff_rate':                             'R',
        'snow__melt_rate':                                     'SM',
        #--------------------------------------------------------------------
        'channel_bed__manning_coefficient':                    'nval',
        'channel_bed__max_over_domain_of_manning_coefficient': 'nval_max',
        'channel_bed__max_over_domain_of_roughness_length':    'z0val_max',
        'channel_bed__min_over_domain_of_manning_coefficient': 'nval_min',
        'channel_bed__min_over_domain_of_roughness_length':    'z0val_min',
        'channel_bed__roughness_length':                       'z0val',
        'channel_bed_surface__slope':                          'S_bed',
        'channel_centerline__straight_sinuosity':              'sinu',
        'channel_cross_section__hydraulic_radius':             'Rh',
        'channel_cross_section_trapezoid__bank_angle':         'angle',   ####
        'channel_cross_section_trapezoid__bottom_width':       'width',   ####
##        'channel_cross_section_water__depth':                  'd',
##        'channel_cross_section_water__initial_depth':          'd0',    ####
##        'channel_cross_section_water__speed':                  'u',
        'channel_model__time_step':                            'dt',
        'channel_water__depth':                                'd',
        'channel_water__initial_depth':                        'd0',
        'channel_water__speed':                                'u',
        'channel_water__discharge':                            'Q',
        'channel_water__friction_factor':                      'f',
        # 'channel_water__froude_number':                      'Fr',
        'channel_water__volume':                               'vol',
        'channel_water_model__time_step':                      'dt',
        'channel_water_surface__slope':                        'S_free',
        'model_grid_cell__area':                               'da',  ### model?
        'model__time_step':                                    'dt',  ### model?
        'watershed_outlet_water__depth':                       'd_outlet',
        'watershed_outlet_water__discharge':                   'Q_outlet',
        'watershed_outlet_water__friction_factor':             'f_outlet',
        'watershed_outlet_water__max_over_time_of_depth':      'd_peak',
        'watershed_outlet_water__max_over_time_of_discharge':  'Q_peak',
        'watershed_outlet_water__max_over_time_of_speed':      'u_peak',
        'watershed_outlet_water__speed':                       'u_outlet',
        'watershed_outlet_water__time_integral_of_discharge':  'vol_Q',
        'watershed_outlet_water__time_of_max_of_depth':        'Td_peak',
        'watershed_outlet_water__time_of_max_of_discharge':    'T_peak',
        'watershed_outlet_water__time_of_max_of_speed':        'Tu_peak',
        'watershed_water__area_time_integral_of_runoff_rate':  'vol_R',
        'watershed_water__max_over_domain_of_depth':        'd_max',
        'watershed_water__max_over_domain_of_discharge':    'Q_max',
        'watershed_water__max_over_domain_of_speed':        'u_max',
        'watershed_water__min_over_domain_of_depth':        'd_min',
        'watershed_water__min_over_domain_of_discharge':    'Q_min',
        'watershed_water__min_over_domain_of_speed':        'u_min' }

    #------------------------------------------------
    # Create an "inverse var name map"
    # inv_map = dict(zip(map.values(), map.keys()))
    #------------------------------------------------
##    _inv_var_name_map = dict( zip(_var_name_map.values(),
##                                  _var_name_map.keys() ) )

    _var_units_map = {
        'atmosphere_water__liquid_equivalent_precipitation_rate': 'm s-1',
        'glacier__melt_rate':                                  'm s-1',
        ## 'land_surface__elevation':                          'm',
        ## 'land_surface__slope':                              '1',
        'land_water__baseflow_emergence_rate':                 'm s-1',
        'land_water__evaporation_rate':                        'm s-1',
        'land_water__infiltration_rate':                       'm s-1',
        'land_water__runoff_rate':                             'm s-1',
        'snow__melt_rate':                                     'm s-1',
        #-------------------------------------------------------------------
        'channel_bed__manning_coefficient':                    'm-1/3 s',
        'channel_bed__max_over_domain_of_manning_coefficient': 'm-1/3 s',
        'channel_bed__max_over_domain_of_roughness_length':    'm',
        'channel_bed__min_over_domain_of_manning_coefficient': 'm-1/3 s',
        'channel_bed__min_over_domain_of_roughness_length':    'm',
        'channel_bed__roughness_length':                       'm',
        'channel_bed_surface__slope':                          '1',
        'channel_centerline__straight_sinuosity':              '1',
        'channel_cross_section__hydraulic_radius':             'm',
        'channel_cross_section_trapezoid__bank_angle':         'rad', # CHECKED
        'channel_cross_section_trapezoid__bottom_width':       'm',
##        'channel_cross_section_water__depth':                  'm',
##        'channel_cross_section_water__initial_depth':          'm',
##        'channel_cross_section_water__speed':                  'm s-1',
        'channel_model__time_step':                            's',
        'channel_water__depth':                                'm',
        'channel_water__initial_depth':                        'm',
        'channel_water__speed':                                'm s-1',
        'channel_water__discharge':                            'm3 s-1',
        'channel_water__friction_factor':                      '1',
        # 'channel_water__froude_number':                      '1',
        'channel_water__volume':                               'm3',
        'channel_water_model__time_step':                      's',
        'channel_water_surface__slope':                        '1',
        'model_grid_cell__area':                               'm2', ### model?
        'model__time_step':                                    's',  ### model?
        'watershed_outlet_water__depth':                       'm',
        'watershed_outlet_water__discharge':                   'm3',
        'watershed_outlet_water__friction_factor':             '1',
        'watershed_outlet_water__time_integral_of_discharge':  'm3',
        'watershed_outlet_water__max_over_time_of_depth':       'm',
        'watershed_outlet_water__max_over_time_of_discharge':   'm3 s-1',
        'watershed_outlet_water__max_over_time_of_speed':       'm s-1',
        'watershed_outlet_water__speed':                        'm s-1',
        'watershed_outlet_water__time_of_max_of_depth':         'min',
        'watershed_outlet_water__time_of_max_of_discharge':     'min',
        'watershed_outlet_water__time_of_max_of_speed':         'min',
        'watershed_water__area_time_integral_of_runoff_rate':'m3',
        'watershed_water__max_over_domain_of_depth':            'm',
        'watershed_water__max_over_domain_of_discharge':        'm3 s-1',
        'watershed_water__max_over_domain_of_speed':            'm s-1',
        'watershed_water__min_over_domain_of_depth':            'm',
        'watershed_water__min_over_domain_of_discharge':        'm3 s-1',
        'watershed_water__min_over_domain_of_speed':            'm s-1' }
        
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )
        
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------   
        return self._input_var_names
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
 
        return self._output_var_names
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):
            
        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]
   
    #   get_var_units()
    #-------------------------------------------------------------------
##    def get_var_type(self, long_var_name):
##
##        #---------------------------------------
##        # So far, all vars have type "double",
##        # but use the one in BMI_base instead.
##        #---------------------------------------
##        return 'double'
##    
##    #   get_var_type()
    #-------------------------------------------------------------------
    def set_constants(self):

        #------------------------
        # Define some constants
        #------------------------
        self.g          = np.float64(9.81)    # (gravitation const.)
        self.aval       = np.float64(0.476)   # (integration const.)
        self.kappa      = np.float64(0.408)   # (von Karman's const.)
        self.law_const  = np.sqrt(self.g) / self.kappa
        self.one_third  = np.float64(1.0) / 3.0        
        self.two_thirds = np.float64(2.0) / 3.0
        self.deg_to_rad = np.pi / 180.0
        
    #   set_constants()
    #-------------------------------------------------------------------
    def initialize(self, cfg_prefix=None, mode="nondriver",
                   SILENT=False): 

        if not(SILENT):
            print ' '
            print 'Channels component: Initializing...'
        
        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_prefix = cfg_prefix

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()           # (12/7/09)
        # print 'CHANNELS calling initialize_config_vars()...'
        self.initialize_config_vars()   
        # print 'CHANNELS calling read_grid_info()...'
        self.read_grid_info()
        #print 'CHANNELS calling initialize_basin_vars()...'
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        # print 'CHANNELS calling initialize_time_vars()...'
        self.initialize_time_vars()
        
        #----------------------------------
        # Has component been turned off ?
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print 'Channels component: Disabled.'
            self.SAVE_Q_GRIDS  = False   # (It is True by default.)
            self.SAVE_Q_PIXELS = False   # (It is True by default.)
            self.DONE = True
            self.status = 'initialized'  # (OpenMI 2.0 convention) 
            return

##        print '################################################'
##        print 'min(d0), max(d0) =', self.d0.min(), self.d0.max()
##        print '################################################'
        
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        # Can't move read_input_files() to start of
        # update(), since initial values needed here.
        #---------------------------------------------
        # print 'CHANNELS calling open_input_files()...'
        self.open_input_files()
        print 'CHANNELS calling read_input_files()...'
        self.read_input_files()

        #-----------------------
        # Initialize variables
        #-----------------------
        print 'CHANNELS calling initialize_d8_vars()...'
        self.initialize_d8_vars()  # (depend on D8 flow grid)
        print 'CHANNELS calling initialize_computed_vars()...'
        self.initialize_computed_vars()

        #--------------------------------------------------
        # (5/12/10) I think this is obsolete now.
        #--------------------------------------------------
        # Make sure self.Q_ts_file is not NULL (12/22/05)
        
        # This is only output file that is set by default
        # and is still NULL if user hasn't opened the
        # output var dialog for the channel process.
        #--------------------------------------------------
##        if (self.SAVE_Q_PIXELS and (self.Q_ts_file == '')):    
##            self.Q_ts_file = (self.case_prefix + '_0D-Q.txt')       

        ###############################################################
        ###############################################################
        #  If this component is running in stand-alone mode,
        #  then it must initialize the process modules that
        #  it is going to use.  If it is being called by another
        #  component, should caller pass ports in via "set_value" ??
        ###############################################################
        ###############################################################

        ############################################
        # Not needed by new framework, 5/17/12.
        ############################################
        ## print 'CHANNELS calling initialize_required_components()...'
##        self.initialize_required_components( mode )
        
        self.open_output_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention) 
        
    #   initialize()
    #-------------------------------------------------------------------
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        #---------------------------------------------
        # Note that u and d from previous time step
        # must be used on RHS of the equations here.
        #---------------------------------------------
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #-------------------------------------------------------
        # There may be times where we want to call this method
        # even if component is not the driver.  But note that
        # the TopoFlow driver also makes this same call.
        #-------------------------------------------------------
        if (self.mode == 'driver'):
            self.print_time_and_value(self.Q_outlet, 'Q_out', '[m^3/s]')
                                      ### interval=0.5)  # [seconds]

        # For testing (5/19/12)
        # self.print_time_and_value(self.Q_outlet, 'Q_out', '[m^3/s]  CHANNEL')
            
        ## DEBUG = True
        DEBUG = False
 
        #-------------------------
        # Update computed values
        #-------------------------
        if (DEBUG): print '#### Calling update_R()...'
        self.update_R()
        if (DEBUG): print '#### Calling update_R_integral()...'
        self.update_R_integral()
        if (DEBUG): print '#### Calling update_discharge()...'
        self.update_discharge()
        if (DEBUG): print '#### Calling update_flow_volume()...'
        self.update_flow_volume()
        if (DEBUG): print '#### Calling update_flow_depth()...'
        self.update_flow_depth()
        if not(self.DYNAMIC_WAVE):
            if (DEBUG): print '#### Calling update_trapezoid_Rh()...'
            self.update_trapezoid_Rh()
            # print 'Rhmin, Rhmax =', self.Rh.min(), self.Rh.max()a
        if (DEBUG): print '#### Calling update_velocity()...'
        self.update_velocity()
        self.update_velocity_on_edges()     # (set to zero)
##        print 'Rmin, Rmax =', self.R.min(), self.R.max()
##        print 'Qmin,  Qmax =',  self.Q.min(), self.Q.max()
##        print 'umin,  umax =',  self.u.min(), self.u.max()
##        print 'dmin,  dmax =',  self.d.min(), self.d.max()
##        print 'nmin,  nmax =',  self.nval.min(), self.nval.max()
##        print 'Rhmin, Rhmax =', self.Rh.min(), self.Rh.max()
##        print 'Smin,  Smax =',  self.S_bed.min(), self.S_bed.max()
        if (DEBUG): print '#### Calling update_outlet_values()...'
        self.update_outlet_values()
        if (DEBUG): print '#### Calling update peak values()...'
        self.update_peak_values()
        if (DEBUG): print '#### Calling update_Q_out_integral()...'
        self.update_Q_out_integral()

        #---------------------------------------------
        # This takes extra time and is now done
        # only at the end, in finalize().  (8/19/13)
        #---------------------------------------------
        # But then "topoflow_driver" doesn't get
        # correctly updated values for some reason.
        #---------------------------------------------
        ## self.update_mins_and_maxes()
        
        #------------------------
        # Check computed values
        #------------------------
        D_OK = self.check_flow_depth()
        U_OK = self.check_flow_velocity()
        OK   = (D_OK and U_OK)

        #-------------------------------------------
        # Read from files as needed to update vars 
        #-----------------------------------------------------
        # NB! This is currently not needed for the "channel
        # process" because values don't change over time and
        # read_input_files() is called by initialize().
        #-----------------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        if (DEBUG): print '#### Calling write_output_files()...'
        self.write_output_files()
        ## self.write_output_files( time_seconds )

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        if (DEBUG): print '#### Calling update_time()'
        self.update_time( dt )
        
        if (OK):
            self.status = 'updated'  # (OpenMI 2.0 convention)
        else:
            self.status = 'failed'
            self.DONE   = True
            
    #   update()   
    #-------------------------------------------------------------------
    def finalize(self):

        #---------------------------------------------------
        # We can compute mins and maxes in the final grids
        # here, but the framework will not then pass them
        # to any component (e.g. topoflow_driver) that may
        # need them.
        #---------------------------------------------------
        REPORT = True
        self.update_mins_and_maxes( REPORT=REPORT )  ## (2/6/13)
        self.print_final_report(comp_name='Channels component')
        
        self.status = 'finalizing'  # (OpenMI)
        self.close_input_files()    # TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'   # (OpenMI)

        #---------------------------
        # Release all of the ports
        #----------------------------------------
        # Make this call in "finalize()" method
        # of the component's CCA Imple file
        #----------------------------------------
        # self.release_cca_ports( d_services )
        
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        cfg_extension = self.get_attribute( 'cfg_extension' ).lower()
        # cfg_extension = self.get_cfg_extension().lower()
        self.KINEMATIC_WAVE = ("kinematic" in cfg_extension)
        self.DIFFUSIVE_WAVE = ("diffusive" in cfg_extension)
        self.DYNAMIC_WAVE   = ("dynamic"   in cfg_extension)

        ##########################################################
        # (5/17/12) If MANNING, we need to set z0vals to -1 so
        # they are always defined for use with new framework.
        ##########################################################
        if (self.MANNING):
            if (self.nval != None):
                self.nval = np.float64( self.nval )  #### 10/9/10, NEED
                self.nval_min = self.nval.min()
                self.nval_max = self.nval.max()
            #-----------------------------------
            self.z0val     = np.float64(-1)
            self.z0val_min = np.float64(-1)
            self.z0val_max = np.float64(-1)
            
        if (self.LAW_OF_WALL):
            if (self.z0val != None):
                self.z0val = np.float64( self.z0val )  #### (10/9/10)
                self.z0val_min = self.z0val.min()
                self.z0val_max = self.z0val.max()
            #-----------------------------------
            self.nval      = np.float64(-1)
            self.nval_min  = np.float64(-1)
            self.nval_max  = np.float64(-1)
            
        #-------------------------------------------
        # These currently can't be set to anything
        # else in the GUI, but need to be defined.
        #-------------------------------------------
        self.code_type  = 'Grid'
        self.slope_type = 'Grid'

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = np.maximum(self.save_grid_dt,   self.dt)
        self.save_pixels_dt = np.maximum(self.save_pixels_dt, self.dt)
        
        #---------------------------------------------------
        # This is now done in CSDMS_base.read_config_gui()
        # for any var_name that starts with "SAVE_".
        #---------------------------------------------------
        # self.SAVE_Q_GRID = (self.SAVE_Q_GRID == 'Yes')
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
##    def get_cca_port_info(self):
##
##        self.cca_port_names = ['meteorology', 'snow', 'evap',
##                               'infil', 'satzone', 'diversions',
##                               'ice']
##        self.cca_port_short = ['mp', 'sp', 'ep', 'ip', 'gp',
##                               'dp', 'iip']
##        self.cca_port_type  = "IRFPort"
##        self.cca_project    = "edu.csdms.models"
##
##    #   get_cca_port_info()
##    #-------------------------------------------------------------------
##    def embed_child_components(self):
##
##        #------------------------------------------------
##        # Instantiate and embed "process components"
##        # in the place of the CCA ports.
##        #------------------------------------------------
##        # But how do we choose a given process method
##        # such as "kinematic wave" ??
##        #------------------------------------------------
##        import met_base
##        import snow_degree_day
##        import ET_priestley_taylor
##        import infil_green_ampt
##        import GW_darcy_layers
##        import diversions_fraction_method
##        import ice_base
##
##        self.mp = met_base.met_component()
##        self.sp = snow_degree_day.snow_component()
##        self.ep = ET_priestley_taylor.ET_component()
##        self.ip = infil_green_ampt.infil_component()
##        self.gp = GW_darcy_layers.GW_component()
##        self.dp = diversions_fraction_method.diversions_component()
##        self.iip= ice_base.ice_component()   ######
##        
##    #   embed_child_components()
##    #-------------------------------------------------------------------
##    def add_child_ports(self):
##
##        self.add_child_port('sp', 'mp')
##        #----------------------------------
##        self.add_child_port('ep', 'cp', SELF=True)
##        self.add_child_port('ep', 'sp')
##        self.add_child_port('ep', 'ip')
##        self.add_child_port('ep', 'gp')
##        self.add_child_port('ep', 'mp')
##        #----------------------------------
##        self.add_child_port('ip', 'mp')
##        self.add_child_port('ip', 'sp')
##        self.add_child_port('ip', 'ep')
##        self.add_child_port('ip', 'gp')
##        #----------------------------------
##        self.add_child_port('gp', 'cp', SELF=True)
##        self.add_child_port('gp', 'ip')
##        #----------------------------------
##        self.add_child_port('dp', 'cp', SELF=True)
##        #----------------------------------
##        self.add_child_port('mp', 'sp')
##        ## self.add_child_port('mp', 'pp')   ##### in future ?
##        #----------------------------------
##        self.add_child_port('iip', 'mp')
##        self.add_child_port('iip', 'sp')
##        
##    #   add_child_ports()
##    #-------------------------------------------------------------------
##    def initialize_ports(self):
##
##        #-------------------------------------------------
##        # Initialize the process objects/components
##        # This is also where output files are opened.
##        #-------------------------------------------------
##        DEBUG = True      
##        if (self.mp.get_status() != 'initialized'):      # met vars + precip 
##            self.mp.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized METEOROLOGY.\n'
##            
####        if (self.pp.get_status() != 'initialized'):        # precip vars
####            self.pp.initialize( cfg_prefix=self.cfg_prefix ) 
##
##        if (self.sp.get_status() != 'initialized'):        # snow vars         
##            self.sp.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized SNOW.\n'
##            
##        ## GW must get intialized before ET (9/25/09)
##        if (self.gp.get_status() != 'initialized'):        # GW vars
##            self.gp.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized SATZONE.\n'
##            
##        if (self.ep.get_status() != 'initialized'):        # ET vars        
##            self.ep.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized EVAP.\n'
##
##        if (self.dp.get_status() != 'initialized'):        # diversions
##            self.dp.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized DIVERSIONS.\n'
##            
##        if (self.ip.get_status() != 'initialized'):        # infil vars
##            self.ip.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized INFIL.\n'
##            
##        if (self.iip.get_status() != 'initialized'):        # ice vars
##            self.iip.initialize( cfg_prefix=self.cfg_prefix )
##            if (DEBUG): print 'CHANNELS component initialized ICE.\n'
##            
##
##    #   initialize_ports()
    #-------------------------------------------------------------------
    def initialize_d8_vars(self):

        #---------------------------------------------
        # Compute and store a variety of (static) D8
        # flow grid variables.  Embed structure into
        # the "channel_base" component.
        #---------------------------------------------
        self.d8 = d8_base.d8_component()
        ###############################################
        # (5/13/10)  Do next line here for now, until
        # the d8 cfg_file includes static prefix.
        # Same is done in GW_base.py.
        ###############################################
        # tf_d8_base.read_grid_info() also needs
        # in_directory to be set. (10/27/11)
        ###############################################
        self.d8.site_prefix  = self.site_prefix
        self.d8.in_directory = self.in_directory
 
        self.d8.initialize( cfg_prefix=self.cfg_prefix, 
                            SILENT=self.SILENT,
                            REPORT=self.REPORT )
        
        ## self.code = self.d8.code   # Don't need this.
        
        #-------------------------------------------      
        # We'll need this once we shift from using
        # "tf_d8_base.py" to the new "d8_base.py"
        #-------------------------------------------
        # self.d8.update(self.time, SILENT=False, REPORT=True)

    #   initialize_d8_vars()
    #-------------------------------------------------------------
    def initialize_computed_vars(self):
              
        #-----------------------------------------------
        # Convert bank angles from degrees to radians. 
        #-----------------------------------------------
        self.angle = self.angle * self.deg_to_rad  # [radians]
        
        #------------------------------------------------
        # 8/29/05.  Multiply ds by (unitless) sinuosity
        # Orig. ds is used by subsurface flow
        #------------------------------------------------
        # NB!  We should also divide slopes in S_bed by
        # the sinuosity, as now done here.
        #----------------------------------------------------
        # NB!  This saves a modified version of ds that
        #      is only used within the "channels" component.
        #      The original "ds" is stored within the
        #      topoflow model component and is used for
        #      subsurface flow, etc.
        #----------------------------------------------------
        ### self.d8.ds_chan = (self.sinu * ds)
        ### self.ds = (self.sinu * self.d8.ds)
        self.d8.ds = (self.sinu * self.d8.ds)  ### USE LESS MEMORY

        ###################################################
        ###################################################
        ### S_bed = (S_bed / self.sinu)     #*************
        self.slope = (self.slope / self.sinu)
        self.S_bed  = self.slope
        ###################################################
        ###################################################
        
        #---------------------------
        # Initialize spatial grids
        #-----------------------------------------------
        # NB!  It is not a good idea to initialize the
        # water depth grid to a nonzero scalar value.
        #-----------------------------------------------
        print 'Initializing u, f, d grids...'
        self.u = np.zeros([self.ny, self.nx], dtype='Float64')
        self.f = np.zeros([self.ny, self.nx], dtype='Float64')
        self.d = np.zeros([self.ny, self.nx], dtype='Float64') + self.d0

        #########################################################
        # Add this on (2/3/13) so make the TF driver happy
        # during its initialize when it gets reference to R.
        # But in "update_R()", be careful not to break the ref.
        # "Q" may be subject to the same issue.
        #########################################################
        self.Q = np.zeros([self.ny, self.nx], dtype='Float64')
        self.R = np.zeros([self.ny, self.nx], dtype='Float64')

        #---------------------------------------
        # These are used to check mass balance
        #---------------------------------------
        self.vol_R = self.initialize_scalar( 0, dtype='float64')
        self.vol_Q = self.initialize_scalar( 0, dtype='float64')
        
        #-------------------------------------------
        # Make sure all slopes are valid & nonzero
        # since otherwise flow will accumulate
        #-------------------------------------------
        if (self.KINEMATIC_WAVE):    
            self.remove_bad_slopes()      #(3/8/07. Only Kin Wave case)
        
        #----------------------------------------
        # Initial volume of water in each pixel
        #----------------------------------------
        ## self.vol = zeros([self.ny, self.nx], dtype='Float64')
        A_cross  = self.d * (self.width + (self.d * np.tan(self.angle)))
        self.vol = A_cross * self.d8.ds   #[m^3]

        self.initialize_outlet_values()
        self.initialize_peak_values()
        self.initialize_min_and_max_values()  ## (2/3/13)

        #########################################
        # Maybe save all refs in a dictionary
        # called "self_values" here ? (2/19/13)
        # Use a "reverse" var_name mapping?
        # inv_map = dict(zip(map.values(), map.keys()))
        #########################################
        
##        w  = np.where( self.width <= 0 )
##        nw = np.size( w[0] )   # (This is correct for 1D or 2D.)
##        if (nw > 0):
##            print 'WARNING:'
##            print 'Number of locations where width==0 =', nw
##            if (nw < 10):
##                print 'locations =', w
##            print ' '

    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def initialize_outlet_values(self):

        #---------------------------------------------------
        # Note:  These are retrieved and used by TopoFlow
        #        for the stopping condition.  TopoFlow
        #        receives a reference to these, but in
        #        order to see the values change they need
        #        to be stored as mutable, 1D numpy arrays.
        #---------------------------------------------------
        # Note:  Q_last is internal to TopoFlow.
        #---------------------------------------------------        
        # self.Q_outlet = self.Q[ self.outlet_ID ]
        self.Q_outlet = self.initialize_scalar(0, dtype='float64')
        self.u_outlet = self.initialize_scalar(0, dtype='float64')
        self.d_outlet = self.initialize_scalar(0, dtype='float64')
        self.f_outlet = self.initialize_scalar(0, dtype='float64')
          
    #   initialize_outlet_values()  
    #-------------------------------------------------------------------
    def initialize_peak_values(self):

        #-------------------------
        # Initialize peak values
        #-------------------------
        self.Q_peak  = self.initialize_scalar(0, dtype='float64')
        self.T_peak  = self.initialize_scalar(0, dtype='float64')
        self.u_peak  = self.initialize_scalar(0, dtype='float64')
        self.Tu_peak = self.initialize_scalar(0, dtype='float64') 
        self.d_peak  = self.initialize_scalar(0, dtype='float64')
        self.Td_peak = self.initialize_scalar(0, dtype='float64')

    #   initialize_peak_values()
    #-------------------------------------------------------------------
    def initialize_min_and_max_values(self):

        #-------------------------------
        # Initialize min & max values
        # (2/3/13), for new framework.
        #-------------------------------
        v = 1e6
        self.Q_min = self.initialize_scalar(v,  dtype='float64')
        self.Q_max = self.initialize_scalar(-v, dtype='float64')
        self.u_min = self.initialize_scalar(v,  dtype='float64')
        self.u_max = self.initialize_scalar(-v, dtype='float64')
        self.d_min = self.initialize_scalar(v,  dtype='float64')
        self.d_max = self.initialize_scalar(-v, dtype='float64')

    #   initialize_min_and_max_values() 
    #-------------------------------------------------------------------
    # def update_excess_rainrate(self):
    def update_R(self):

        #----------------------------------------
        # Compute the "excess rainrate", R.
        # Each term must have same units: [m/s]
        # Sum = net gain/loss rate over pixel.
        #----------------------------------------------------
        # R can be positive or negative.  If negative, then
        # water is removed from the surface at rate R until
        # surface water is consumed.
        #--------------------------------------------------------------
        # P  = precip_rate   [m/s]  (converted by read_input_data()).
        # SM = snowmelt rate [m/s]
        # GW = seep rate     [m/s]  (water_table intersects surface)
        # ET = evap rate     [m/s]
        # IN = infil rate    [m/s]
        # MR = icemelt rate  [m/s]

        #------------------------------------------------------------
        # Use refs to other comp vars from new framework. (5/18/12)
        #------------------------------------------------------------         
        P  = self.P
        SM = self.SM
        GW = self.GW
        ET = self.ET
        IN = self.IN
        MR = self.MR
        #--------------------------------------------------------------        
##        P  = self.get_port_data('P',  self.mp,  'METEOROLOGY')
##        SM = self.get_port_data('SM', self.sp,  'SNOW')
##        GW = self.get_port_data('GW', self.gp,  'SATZONE')
##        ET = self.get_port_data('ET', self.ep,  'EVAP')
##        IN = self.get_port_data('IN', self.ip,  'INFIL')
##        MR = self.get_port_data('MR', self.iip, 'ICE')

        
##        if (self.DEBUG):
##            print 'At time:', self.time_min, ', P =', P, '[m/s]'

##        ########################
##        # ONLY FOR TESTING
##        ########################
##        try:
##            H = self.get_port_data('H', self.iip) ## (ice depth   [m/s])
##        except:
##            H = np.float64(0)
##            print 'ERROR: Could not get "H" from ICE port.'
            
        #-------------------------------------------------------
        # (10/1/09)  This type conversion should not be needed
        # because it is done by get_port_data(), which is a
        # "private/local" vs. interface function/method.
        #-------------------------------------------------------
##        P  = np.float64(P)
##        SM = np.float64(SM)
##        GW = np.float64(GW)
##        ET = np.float64(ET)
##        IN = np.float64(IN)
##        MR = np.float64(MR)
##        H  = np.float64(H)

        #--------------
        # For testing
        #--------------        
##        print '(Pmin,  Pmax)  =', P.min(),  P.max()
##        print '(SMmin, SMmax) =', SM.min(), SM.max()
##        print '(GWmin, GWmax) =', GW.min(), GW.max()
##        print '(ETmin, ETmax) =', ET.min(), ET.max()
##        print '(INmin, INmax) =', IN.min(), IN.max()
##        print '(MRmin, MRmax) =', MR.min(), MR.max()
##        # print '(Hmin,  Hmax)  =', H.min(), H.max()
##        print ' '
        
        self.R = (P + SM + GW + MR) - (ET + IN)
            
    #   update_R()
    #-------------------------------------------------------------------
    def update_R_integral(self):

        #-----------------------------------------------
        # Update mass total for R, sum over all pixels
        #-----------------------------------------------   
        volume = np.double(self.R * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_R += (volume * self.rti.n_pixels)
        else:
            self.vol_R += np.sum(volume)

    #   update_R_integral()           
    #-------------------------------------------------------------------  
    def update_discharge(self):

        #---------------------------------------------------------
        # The discharge grid, Q, gives the flux of water _out_
        # of each grid cell.  This entire amount then flows
        # into one of the 8 neighbor grid cells, as indicated
        # by the D8 flow code. The update_flow_volume() function
        # is called right after this one in update() and uses
        # the Q grid.
        #---------------------------------------------------------
        # 7/15/05.  The cross-sectional area of a trapezoid is
        # given by:    Ac = d * (w + (d * tan(theta))),
        # where w is the bottom width.  If we were to
        # use: Ac = w * d, then we'd have Ac=0 when w=0.
        # We also need angle units to be radians.
        #---------------------------------------------------------

        #-----------------------------
        # Compute the discharge grid
        #-----------------------------
        L2 = self.d * np.tan(self.angle)
        Ac = self.d * (self.width + L2)
        ### self.Q = np.float64(self.u * Ac)
        self.Q[:] = self.u * Ac   ## (2/19/13, in place)

        #--------------
        # For testing
        #--------------  
##        print '(umin,   umax)  =', self.u.min(), self.u.max()
##        print '(d0min, d0max)  =', self.d0.min(), self.d0.max()
##        print '(dmin,   dmax)  =', self.d.min(), self.d.max()
##        print '(amin,   amax)  =', self.angle.min(), self.angle.max()
##        print '(wmin,   wmax)  =', self.width.min(), self.width.max()
##        print '(Qmin,   Qmax)  =', self.Q.min(),  self.Q.max()
##        print '(L2min,  L2max) =', L2.min(), L2.max()
##        print '(Qmin,   Qmax)  =', self.Q.min(),  self.Q.max()
        
        #--------------------------------------
        # Update Q and vol due to diversions
        #--------------------------------------
        # print 'Calling diversions.update() in update_discharge()...'

        ############################################
        ############################################  
        # What to do about new framework ??
        ############################################
        ############################################        
        ######## self.dp.update( self.time_sec )
        ############################################
        ############################################


        # print 'Finished with call to diversions.update()...'

        ############################################
        ############################################        
        ## print 'Updating DIVERSIONS...'
        ## print ' '
        
        #--------------
        # For testing
        #--------------
        # print 'dmin, dmax =', self.d.min(), self.d.max()
        # print 'umin, umax =', self.u.min(), self.u.max()
        # print 'Qmin, Qmax =', self.Q.min(), self.Q.max()
        # print ' '        
        # print 'u(outlet) =', self.u[self.outlet_ID]
        # print 'Q(outlet) =', self.Q[self.outlet_ID]  ########
          
        #----------------------------------------------------
        # Wherever depth is less than z0, assume that water
        # is not flowing and set u and Q to zero.
        # However, we also need (d gt 0) to avoid a divide
        # by zero problem, even when numerators are zero.
        #----------------------------------------------------
        # FLOWING = (d > (z0/aval))
        #*** FLOWING[self.d8.noflow_IDs] = False    ;******
        # u = (u * FLOWING)
        # Q = (Q * FLOWING)
        # d = np.maximum(d, 0.0)    ;(allow depths lt z0, if gt 0.)

    #   update_discharge()
    #-------------------------------------------------------------------
    def update_flow_volume(self):

        #-----------------------------------------------------------
        # Notes: This function must be called after
        #        update_discharge() and update_diversions().
        #-----------------------------------------------------------        
        # Notes: Q   = surface discharge  [m^3/s]
        #        R   = excess precip. rate  [m/s]
        #        da  = pixel area  [m^2]
        #        dt  = channel flow timestep  [s]
        #        vol = total volume of water in pixel [m^3]
        #        v2  = temp version of vol
        #        w1  = IDs of pixels that...
        #        p1  = IDs of parent pixels that...
        #-----------------------------------------------------------
        dt = self.dt  # [seconds]

        #----------------------------------------------------
        # Add contribution (or loss ?) from excess rainrate
        #----------------------------------------------------
        self.vol += (self.R * self.da) * dt
    
        #-----------------------------------------
        # Add contributions from neighbor pixels
        #-------------------------------------------------------------
        # Each grid cell passes flow to *one* downstream neighbor.
        # Note that multiple grid cells can flow toward a given grid
        # cell, so a grid cell ID may occur in d8.p1 and d8.p2, etc.
        #-------------------------------------------------------------
        # (2/16/10)  RETEST THIS.  Before, a copy called "v2" was
        # used but this doesn't seem to be necessary.
        #-------------------------------------------------------------        
        if (self.d8.p1_OK):    
            self.vol[ self.d8.p1 ] += (dt * self.Q[self.d8.w1])
        if (self.d8.p2_OK):    
            self.vol[ self.d8.p2 ] += (dt * self.Q[self.d8.w2])
        if (self.d8.p3_OK):    
            self.vol[ self.d8.p3 ] += (dt * self.Q[self.d8.w3])
        if (self.d8.p4_OK):    
            self.vol[ self.d8.p4 ] += (dt * self.Q[self.d8.w4])
        if (self.d8.p5_OK):    
            self.vol[ self.d8.p5 ] += (dt * self.Q[self.d8.w5])
        if (self.d8.p6_OK):    
            self.vol[ self.d8.p6 ] += (dt * self.Q[self.d8.w6])
        if (self.d8.p7_OK):    
            self.vol[ self.d8.p7 ] += (dt * self.Q[self.d8.w7])
        if (self.d8.p8_OK):    
            self.vol[ self.d8.p8 ] += (dt * self.Q[self.d8.w8])

        #----------------------------------------------------
        # Subtract the amount that flows out to D8 neighbor
        #----------------------------------------------------
        self.vol -= (self.Q * dt)
   
        #--------------------------------------------------------
        # While R can be positive or negative, the surface flow
        # volume must always be nonnegative. This also ensures
        # that the flow depth is nonnegative.  (7/13/06)
        #--------------------------------------------------------
        ## self.vol = np.maximum(self.vol, 0.0)
        ## self.vol[:] = np.maximum(self.vol, 0.0)  # (2/19/13)
        np.maximum( self.vol, 0.0, self.vol )  # (in place)
        
    #   update_flow_volume
    #-------------------------------------------------------------------
    def update_flow_depth(self):

        #-----------------------------------------------------------
        # Notes: 7/18/05.  Modified to use the equation for volume
        #        of a trapezoidal channel:  vol = Ac * ds, where
        #        Ac=d*[w + d*tan(t)], and to solve the resulting
        #        quadratic (discarding neg. root) for new depth, d.

        #        8/29/05.  Now original ds is used for subsurface
        #        flow and there is a ds_chan which can include a
        #        sinuosity greater than 1.  This may be especially
        #        important for larger pixel sizes.

        #        Removed (ds > 1) here which was only meant to
        #        avoid a "divide by zero" error at pixels where
        #        (ds eq 0).  This isn't necessary since the
        #        Flow_Lengths function in utils_TF.pro never
        #        returns a value of zero.
        #----------------------------------------------------------
        d = self.d   # (local alias)

        #----------------------------------------------------------
        # Commented this out on (2/18/10) because it doesn't
        #           seem to be used anywhere now.  Checked all
        #           of the Channels components.
        #----------------------------------------------------------        
        # self.d_last = self.d.copy()
        
        #------------------------------------------------------
        # (2/18/10) New code to deal with case where the flow
        #           depth exceeds a bankfull depth.
        #           For now, d_bankfull is hard-coded.
        #
        #           CHANGE Manning's n here, too?
        #------------------------------------------------------
        width = self.width  ###
        angle = self.angle
        SCALAR_ANGLES = (np.size(angle) == 1)
        #-------------------------------------
        d_bankfull = 4.0  # [meters]
        ################################
        w_overbank = np.where( d > d_bankfull )
        n_overbank = np.size( w_overbank[0] )
        if (n_overbank != 0):
            width[ w_overbank ] = self.d8.dw[ w_overbank ]
            if not(SCALAR_ANGLES): angle[w_overbank] = 0.0
    
        #------------------------------------------------------
        # (2/18/10) New code to deal with case where the top
        #           width exceeds the grid cell width, dw.
        #------------------------------------------------------            
        top_width  = width + (2.0 * d * np.sin(self.angle))
        w_bad      = np.where(top_width > self.d8.dw)
        n_bad      = np.size(w_bad[0])
        if (n_bad != 0):
            width[w_bad] = self.d8.dw[w_bad]
            if not(SCALAR_ANGLES): angle[w_bad] = 0.0

        #----------------------------------
        # Is "angle" a scalar or a grid ?
        #----------------------------------
        if (SCALAR_ANGLES):
            if (angle == 0.0):    
                d = self.vol / (width * self.d8.ds)
            else:
                term1 = 2.0 * np.tan(angle)
                arg   = 2.0 * term1 * self.vol / self.d8.ds
                arg  += width**(2.0)
                d     = (np.sqrt(arg) - width) / term1
        else:
            #-----------------------------------------------------
            # Pixels where angle is 0 must be handled separately
            #-----------------------------------------------------
            wz   = np.where( angle == 0 )
            nwz  = np.size( wz[0] )
            wzc  = np.where( angle != 0 )
            nwzc = np.size( wzc[0] )
            
            if (nwz != 0):
                A_top = width[wz] * self.d8.ds[wz]
                ## A_top = self.width[wz] * self.d8.ds_chan[wz]            
                d[wz] = self.vol[wz] / A_top
            
            if (nwzc != 0):    
                term1  = 2.0 * np.tan(angle[wzc])
                arg    = 2.0 * term1 * self.vol[wzc] / self.d8.ds[wzc]
                arg   += width[wzc]**(2.0)
                d[wzc] = (np.sqrt(arg) - width[wzc]) / term1

        #------------------------------------------
        # Set depth values on edges to zero since
        # they become spikes (no outflow) 7/15/06
        #------------------------------------------    
        d[ self.d8.noflow_IDs ] = 0.0

        #------------------------------------------------
        # 4/19/06.  Force flow depth to be positive ?
        #------------------------------------------------
        # This seems to be needed with the non-Richards
        # infiltration routines when starting with zero
        # depth everywhere, since all water infiltrates
        # for some period of time.  It also seems to be
        # needed more for short rainfall records to
        # avoid a negative flow depth error.
        #------------------------------------------------
        # 7/13/06.  Still needed for Richards method
        #------------------------------------------------
        ## self.d = np.maximum(d, 0.0)
        np.maximum(d, 0.0, self.d)  # (2/19/13, in place)
        
    #   update_flow_depth
    #-------------------------------------------------------------------
    def update_free_surface_slope(self):

        #-----------------------------------------------------------
        # Notes:  It is assumed that the flow directions don't
        #         change even though the free surface is changing.
        #-----------------------------------------------------------
        delta_d     = (self.d - self.d[self.d8.parent_IDs])
        self.S_free[:] = self.S_bed + (delta_d / self.d8.ds)
        
        #--------------------------------------------
        # Don't do this; negative slopes are needed
        # to decelerate flow in dynamic wave case
        # and for backwater effects.
        #--------------------------------------------
        # Set negative slopes to zero
        #------------------------------
        ###  self.S_free = np.maximum(self.S_free, 0)

    #   update_free_surface_slope()
    #-------------------------------------------------------------------
    def update_trapezoid_Rh(self):

        #-------------------------------------------------------------
        # Notes: Compute the hydraulic radius of a trapezoid that:
        #          (1) has a bed width of wb >= 0 (0 for triangular)
        #          (2) has a bank angle of theta (0 for rectangular)
        #          (3) is filled with water to a depth of d.
        #        The units of wb and d are meters.  The units of
        #        theta are assumed to be degrees and are converted.
        #-------------------------------------------------------------
        # NB!    wb should never be zero, so PW can never be 0,
        #        which would produce a NaN (divide by zero).
        #-------------------------------------------------------------
        #        See Notes for TF_Tan function in utils_TF.pro
        #            AW = d * (wb + (d * TF_Tan(theta_rad)) )
        #-------------------------------------------------------------
        d  = self.d        # (local synonyms)
        wb = self.width   # (trapezoid bottom width)
        
        theta_rad = (self.angle * self.deg_to_rad)

##        print 'dmin, dmax =', d.min(),  d.max()
##        print 'wmin, wmax =', wb.min(), wb.max()
##        print 'amin, amax =', theta_rad.min(), theta_rad.max()
        
        #---------------------------------------------------------
        # Compute hydraulic radius grid for trapezoidal channels
        #---------------------------------------------------------        
        AW = d * (wb + (d * np.tan(theta_rad)) )      
        PW = wb + (np.float64(2.0) * d / np.cos(theta_rad) )

        #---------------------------------------------------
        # At noflow_IDs (e.g. edges) PW may be zero
        # so do this to avoid "divide by zero". (10/29/11)
        #---------------------------------------------------
        PW[ self.d8.noflow_IDs ] = np.float64(1)
        Rh = (AW / PW)
        # w = np.where(PW == 0)
        # print 'In update_trapezoid_Rh():'
        # print '   PW = 0 at', w[0].size, 'cells'

        #------------------------------------
        # Force edge pixels to have Rh = 0.
        # This will make u = 0 there also.
        #------------------------------------
        Rh[ self.d8.noflow_IDs ] = np.float64(0)        
##        w  = np.where(wb <= 0)
##        nw = np.size(w[0])
##        if (nw > 0): Rh[w] = np.float64(0)
        
        self.Rh = Rh

    #   update_trapezoid_Rh()
    #-------------------------------------------------------------------
    def update_velocity(self):

        #---------------------------------------------------------
        # Note: Do nothing now unless this method is overridden
        #       by a particular method of computing velocity.
        #---------------------------------------------------------
        print "Warning: update_velocity() method is inactive."
        
        # print 'KINEMATIC WAVE =', self.KINEMATIC_WAVE
        # print 'DIFFUSIVE WAVE =', self.DIFFUSIVE_WAVE
        # print 'DYNAMIC WAVE   =', self.DYNAMIC_WAVE

    #   update_velocity()
    #-------------------------------------------------------------------
    def update_velocity_on_edges(self):

        #---------------------------------
        # Force edge pixels to have u=0.
        #----------------------------------------
        # Large slope around 1 flows into small
        # slope & leads to a negative velocity.
        #----------------------------------------
        self.u[ self.d8.noflow_IDs ] = np.float64(0)
        
    #   update_velocity_on_edges()
    #-------------------------------------------------------------
    def update_outlet_values(self):
        
        #-------------------------------------------------
        # Save computed values at outlet, which are used
        # by the TopoFlow driver.
        #-----------------------------------------------------
        # Note that Q_outlet, etc. are defined as 0D numpy
        # arrays to make them "mutable scalars" (i.e.
        # this allows changes to be seen by other components
        # who have a reference.  To preserver the reference,
        # however, we must use fill() to assign a new value.
        #-----------------------------------------------------
        Q_outlet = self.Q[ self.outlet_ID ]
        u_outlet = self.u[ self.outlet_ID ]
        d_outlet = self.d[ self.outlet_ID ]
        f_outlet = self.f[ self.outlet_ID ]
    
        self.Q_outlet.fill( Q_outlet )
        self.u_outlet.fill( u_outlet )
        self.d_outlet.fill( d_outlet )
        self.f_outlet.fill( f_outlet )
        
##        self.Q_outlet.fill( self.Q[ self.outlet_ID ] )
##        self.u_outlet.fill( self.u[ self.outlet_ID ] )
##        self.d_outlet.fill( self.d[ self.outlet_ID ] )
##        self.f_outlet.fill( self.f[ self.outlet_ID ] )
        
##        self.Q_outlet = self.Q[ self.outlet_ID ]
##        self.u_outlet = self.u[ self.outlet_ID ]
##        self.d_outlet = self.d[ self.outlet_ID ]
##        self.f_outlet = self.f[ self.outlet_ID ]
        
##        self.Q_outlet = self.Q.flat[self.outlet_ID]
##        self.u_outlet = self.u.flat[self.outlet_ID]
##        self.d_outlet = self.d.flat[self.outlet_ID]
##        self.f_outlet = self.f.flat[self.outlet_ID]
        
    #   update_outlet_values()
    #-------------------------------------------------------------
    def update_peak_values(self):

        if (self.Q_outlet > self.Q_peak):    
            self.Q_peak.fill( self.Q_outlet )
            self.T_peak.fill( self.time_min )      # (time to peak)
        #---------------------------------------
        if (self.u_outlet > self.u_peak):
            self.u_peak.fill( self.u_outlet )
            self.Tu_peak.fill( self.time_min )
        #---------------------------------------
        if (self.d_outlet > self.d_peak):    
            self.d_peak.fill(  self.d_outlet )
            self.Td_peak.fill( self.time_min )
            
##        if (self.Q_outlet > self.Q_peak):    
##            self.Q_peak  = self.Q_outlet
##            self.T_peak  = self.time_min      # (time to peak)
##        #-----------------------------------
##        if (self.u_outlet > self.u_peak):
##            self.u_peak  = self.u_outlet
##            self.Tu_peak = self.time_min
##        #-----------------------------------
##        if (self.d_outlet > self.d_peak):    
##            self.d_peak  = self.d_outlet
##            self.Td_peak = self.time_min

    #   update_peak_values()
    #-------------------------------------------------------------
    def update_Q_out_integral(self):

        #--------------------------------------------------------
        # Note: Renamed "volume_out" to "vol_Q" for consistency
        # with vol_P, vol_SM, vol_IN, vol_ET, etc. (5/18/12)
        #--------------------------------------------------------
        ## self.vol_Q += (self.Q_outlet * self.dt)
        self.vol_Q += (self.Q_outlet * self.dt)  ## Experiment: 5/19/12.
        ## self.vol_Q += (self.Q[self.outlet_ID] * self.dt)
        
    #   update_Q_out_integral()
    #-------------------------------------------------------------
    def update_mins_and_maxes(self, REPORT=False):

        #--------------------------------------
        # Get mins and max over entire domain
        #--------------------------------------
##        Q_min = self.Q.min()
##        Q_max = self.Q.max()
##        #---------------------
##        u_min = self.u.min()
##        u_max = self.u.max()        
##        #---------------------
##        d_min = self.d.min()
##        d_max = self.d.max()
        
        #--------------------------------------------
        # Exclude edges where mins are always zero.
        #--------------------------------------------
        nx = self.nx
        ny = self.ny
        Q_min = self.Q[1:(ny - 2)+1,1:(nx - 2)+1].min()
        Q_max = self.Q[1:(ny - 2)+1,1:(nx - 2)+1].max()
        #-------------------------------------------------
        u_min = self.u[1:(ny - 2)+1,1:(nx - 2)+1].min()
        u_max = self.u[1:(ny - 2)+1,1:(nx - 2)+1].max()        
        #-------------------------------------------------
        d_min = self.d[1:(ny - 2)+1,1:(nx - 2)+1].min()
        d_max = self.d[1:(ny - 2)+1,1:(nx - 2)+1].max()

        #-------------------------------------------------
        # (2/6/13) This preserves "mutable scalars" that
        # can be accessed as refs by other components.
        #-------------------------------------------------
        if (Q_min < self.Q_min):
            self.Q_min.fill( Q_min )
        if (Q_max > self.Q_max):
            self.Q_max.fill( Q_max )
        #------------------------------
        if (u_min < self.u_min):
            self.u_min.fill( u_min )
        if (u_max > self.u_max):
            self.u_max.fill( u_max )
        #------------------------------
        if (d_min < self.d_min):
            self.d_min.fill( d_min )
        if (d_max > self.d_max):
            self.d_max.fill( d_max )
        
        #-------------------------------------------------
        # (2/6/13) This preserves "mutable scalars" that
        # can be accessed as refs by other components.
        #-------------------------------------------------        
##        self.Q_min.fill( np.minimum( self.Q_min, Q_min ) )
##        self.Q_max.fill( np.maximum( self.Q_max, Q_max ) )
##        #---------------------------------------------------
##        self.u_min.fill( np.minimum( self.u_min, u_min ) )
##        self.u_max.fill( np.maximum( self.u_max, u_max ) )
##        #---------------------------------------------------
##        self.d_min.fill( np.minimum( self.d_min, d_min ) )
##        self.d_max.fill( np.maximum( self.d_max, d_max ) )

        #-------------------------------------------------
        # (2/6/13) This preserves "mutable scalars" that
        # can be accessed as refs by other components.
        #-------------------------------------------------        
##        self.Q_min.fill( min( self.Q_min, Q_min ) )
##        self.Q_max.fill( max( self.Q_max, Q_max ) )
##        #---------------------------------------------------
##        self.u_min.fill( min( self.u_min, u_min ) )
##        self.u_max.fill( max( self.u_max, u_max ) )
##        #---------------------------------------------------
##        self.d_min.fill( min( self.d_min, d_min ) )
##        self.d_max.fill( max( self.d_max, d_max ) )
        
        #----------------------------------------------
        # (2/6/13) This produces "immutable scalars".
        #----------------------------------------------
##        self.Q_min = self.Q.min()
##        self.Q_max = self.Q.max()
##        self.u_min = self.u.min()
##        self.u_max = self.u.max()
##        self.d_min = self.d.min()
##        self.d_max = self.d.max()

        if (REPORT):
            print 'In channels_base.update_mins_and_maxes():'
            print '(dmin, dmax) =', self.d_min, self.d_max
            print '(umin, umax) =', self.u_min, self.u_max
            print '(Qmin, Qmax) =', self.Q_min, self.Q_max
            print ' '
            
    #   update_mins_and_maxes()
    #-------------------------------------------------------------------
    def check_flow_depth(self):

        OK = True
        d  = self.d
        dt = self.dt
        nx = self.nx   #################
        
        #---------------------------------
        # All all flow depths positive ?
        #---------------------------------
        wbad = np.where( np.logical_or( d < 0.0, np.logical_not(np.isfinite(d)) ))
        nbad = np.size( wbad[0] )       
        if (nbad == 0):    
            return OK

        OK = False
        dmin = d[wbad].min()
        star_line = '*******************************************'
        
        msg = [ star_line, \
               'ERROR: Simulation aborted.', ' ', \
               'Negative depth found: ' + str(dmin), \
               'Time step may be too large.', \
               'Time step:      ' + str(dt) + ' [s]', ' ']
        for k in xrange(len(msg)):
            print msg[k]
        
        #-------------------------------------------
        # If not too many, print actual velocities
        #-------------------------------------------
        if (nbad < 30):          
            brow = wbad[0][0]
            bcol = wbad[1][0]
##            badi = wbad[0]
##            bcol = (badi % nx)
##            brow = (badi / nx)
            crstr = str(bcol) + ', ' + str(brow)

            msg = ['(Column, Row):  ' + crstr, \
                   'Flow depth:     ' + str(d[brow, bcol])]
            for k in xrange(len(msg)):
                print msg[k]

        print star_line 
        print ' '
        return OK

    #   check_flow_depth
    #-------------------------------------------------------------------
    def check_flow_velocity(self):

        OK = True
        u  = self.u
        dt = self.dt
        nx = self.nx
        
        #--------------------------------
        # Are all velocities positive ?
        #--------------------------------
        wbad = np.where( np.logical_or( u < 0.0, np.logical_not(np.isfinite(u)) ))
        nbad = np.size( wbad[0] )
        if (nbad == 0):    
            return OK

        OK = False
        umin = u[wbad].min()
        star_line = '*******************************************'
        msg = [ star_line, \
               'ERROR: Simulation aborted.', ' ', \
               'Negative or NaN velocity found: ' + str(umin), \
               'Time step may be too large.', \
               'Time step:      ' + str(dt) + ' [s]', ' ']
        for k in xrange(len(msg)):
            print msg[k]

        #-------------------------------------------
        # If not too many, print actual velocities
        #-------------------------------------------
        if (nbad < 30):
            brow = wbad[0][0]
            bcol = wbad[1][0]
##            badi = wbad[0]
##            bcol = (badi % nx)
##            brow = (badi / nx)
            crstr = str(bcol) + ', ' + str(brow)

            msg = ['(Column, Row):  ' + crstr, \
                   'Velocity:       ' + str(u[brow, bcol])]
            for k in xrange(len(msg)):
                print msg[k]

        print star_line
        print ' '
        return OK

            
##        umin = u[wbad].min()
##        badi = wbad[0]
##        bcol = (badi % nx)
##        brow = (badi / nx)
##        crstr = str(bcol) + ', ' + str(brow)
##        msg = np.array([' ', \
##                     '*******************************************', \
##                     'ERROR: Simulation aborted.', ' ', \
##                     'Negative velocity found: ' + str(umin), \
##                     'Time step may be too large.', ' ', \
##                     '(Column, Row):  ' + crstr, \
##                     'Velocity:       ' + str(u[badi]), \
##                     'Time step:      ' + str(dt) + ' [s]', \
##                     '*******************************************', ' '])
##        for k in xrange( np.size(msg) ):
##            print msg[k]

##        return OK                          


    #   check_flow_velocity
    #-------------------------------------------------------------------  
    def open_input_files(self):

        # This doesn't work, because file_unit doesn't get full path. (10/28/11)
        # start_dir = os.getcwd()
        # os.chdir( self.in_directory )

        # print '### start_dir =', start_dir
        # print '### in_directory =', self.in_directory

        in_files = ['slope_file', 'nval_file', 'z0val_file',
                    'width_file', 'angle_file', 'sinu_file', 'd0_file']
        self.prepend_directory( in_files, INPUT=True )

        # self.slope_file = self.in_directory + self.slope_file
        # self.nval_file  = self.in_directory + self.nval_file
        # self.z0val_file = self.in_directory + self.z0val_file
        # self.width_file = self.in_directory + self.width_file
        # self.angle_file = self.in_directory + self.angle_file
        # self.sinu_file  = self.in_directory + self.sinu_file
        # self.d0_file    = self.in_directory + self.d0_file

        #self.code_unit = model_input.open_file(self.code_type,  self.code_file)
        self.slope_unit = model_input.open_file(self.slope_type, self.slope_file)
        if (self.MANNING):
            self.nval_unit  = model_input.open_file(self.nval_type,  self.nval_file)
        if (self.LAW_OF_WALL):
            self.z0val_unit = model_input.open_file(self.z0val_type, self.z0val_file)
        self.width_unit = model_input.open_file(self.width_type, self.width_file)
        self.angle_unit = model_input.open_file(self.angle_type, self.angle_file)
        self.sinu_unit  = model_input.open_file(self.sinu_type,  self.sinu_file)
        self.d0_unit    = model_input.open_file(self.d0_type,    self.d0_file)

        # os.chdir( start_dir )

    #   open_input_files()        
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #---------------------------------------------------
        # The flow codes are always a grid, size of DEM.
        #---------------------------------------------------
        # NB! model_input.py also has a read_grid() function.
        #---------------------------------------------------        
        rti = self.rti
##        print 'Reading D8 flow grid (in CHANNELS)...'
##        self.code = rtg_files.read_grid(self.code_file, rti,
##                                        RTG_type='BYTE')
##        print ' '
        
        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        slope = model_input.read_next(self.slope_unit, self.slope_type, rti)
        if (slope != None): self.slope = slope
        
        # If EOF was reached, hopefully numpy's "fromfile"
        # returns None, so that the stored value will be
        # the last value that was read.

        if (self.MANNING):
            nval = model_input.read_next(self.nval_unit, self.nval_type, rti)
            if (nval != None):
                self.nval     = nval
                self.nval_min = nval.min()
                self.nval_max = nval.max()
                
        if (self.LAW_OF_WALL):
            z0val = model_input.read_next(self.z0val_unit, self.z0val_type, rti)
            if (z0val != None):
                self.z0val     = z0val
                self.z0val_min = z0val.min()
                self.z0val_max = z0val.max()
        
        width = model_input.read_next(self.width_unit, self.width_type, rti)
        if (width != None): self.width = width
        
        angle = model_input.read_next(self.angle_unit, self.angle_type, rti)
        if (angle != None): self.angle = angle

        sinu = model_input.read_next(self.sinu_unit, self.sinu_type, rti)
        if (sinu != None): self.sinu = sinu
        
        d0 = model_input.read_next(self.d0_unit, self.d0_type, rti)
        if (d0 != None): self.d0 = d0

        ## code = model_input.read_grid(self.code_unit, \
        ##                            self.code_type, rti, dtype='UInt8')
        ## if (code != None): self.code = code

    #   read_input_files()     
    #-------------------------------------------------------------------  
    def close_input_files(self):

        # if not(self.slope_unit.closed):
        # if (self.slope_unit != None):

        #-------------------------------------------------
        # NB!  self.code_unit was never defined as read.
        #-------------------------------------------------
        # if (self.code_type != 'scalar'): self.code_unit.close()

        if (self.slope_type != 'Scalar'): self.slope_unit.close()
        if (self.MANNING):
            if (self.nval_type != 'Scalar'): self.nval_unit.close()
        if (self.LAW_OF_WALL):
           if (self.z0val_type != 'Scalar'): self.z0val_unit.close()
        if (self.width_type != 'Scalar'): self.width_unit.close()
        if (self.angle_type != 'Scalar'): self.angle_unit.close()
        if (self.sinu_type  != 'Scalar'): self.sinu_unit.close()
        if (self.d0_type    != 'Scalar'): self.d0_unit.close()
    
##        if (self.slope_file != ''): self.slope_unit.close()
##        if (self.MANNING):
##            if (self.nval_file  != ''): self.nval_unit.close()
##        if (self.LAW_OF_WALL):
##           if (self.z0val_file != ''): self.z0val_unit.close()
##        if (self.width_file != ''): self.width_unit.close()
##        if (self.angle_file != ''): self.angle_unit.close()
##        if (self.sinu_file  != ''): self.sinu_unit.close()
##        if (self.d0_file    != ''): self.d0_unit.close()

    #   close_input_files()       
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.Q_gs_file = (self.out_directory + self.Q_gs_file)
        self.u_gs_file = (self.out_directory + self.u_gs_file)
        self.d_gs_file = (self.out_directory + self.d_gs_file) 
        self.f_gs_file = (self.out_directory + self.f_gs_file) 
        #--------------------------------------------------------
        self.Q_ts_file = (self.out_directory + self.Q_ts_file)
        self.u_ts_file = (self.out_directory + self.u_ts_file) 
        self.d_ts_file = (self.out_directory + self.d_ts_file) 
        self.f_ts_file = (self.out_directory + self.f_ts_file) 

    #   update_outfile_names()     
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_nio()
        self.update_outfile_names()

##        print 'self.SAVE_Q_GRIDS =', self.SAVE_Q_GRIDS
##        print 'self.SAVE_U_GRIDS =', self.SAVE_U_GRIDS
##        print 'self.SAVE_D_GRIDS =', self.SAVE_D_GRIDS
##        print 'self.SAVE_F_GRIDS =', self.SAVE_F_GRIDS
##        #---------------------------------------------------
##        print 'self.SAVE_Q_PIXELS =', self.SAVE_Q_PIXELS
##        print 'self.SAVE_U_PIXELS =', self.SAVE_U_PIXELS
##        print 'self.SAVE_D_PIXELS =', self.SAVE_D_PIXELS
##        print 'self.SAVE_F_PIXELS =', self.SAVE_F_PIXELS

        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_Q_GRIDS):   
            model_output.open_new_gs_file( self, self.Q_gs_file, self.rti,
                                           var_name='Q',
                                           long_name='volumetric_discharge',
                                           units_name='m^3/s')
            
        if (self.SAVE_U_GRIDS):    
            model_output.open_new_gs_file( self, self.u_gs_file, self.rti,
                                           var_name='u',
                                           long_name='mean_channel_flow_velocity',
                                           units_name='m/s')
        
        if (self.SAVE_D_GRIDS):    
            model_output.open_new_gs_file( self, self.d_gs_file, self.rti,
                                           var_name='d',
                                           long_name='max_channel_flow_depth',
                                           units_name='m')

        if (self.SAVE_F_GRIDS):    
            model_output.open_new_gs_file( self, self.f_gs_file, self.rti,
                                           var_name='f',
                                           long_name='friction_factor',
                                           units_name='none')
            
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_Q_PIXELS):  
            model_output.open_new_ts_file( self, self.Q_ts_file, IDs,
                                           var_name='Q',
                                           long_name='volumetric_discharge',
                                           units_name='m^3/s')
                                          
        if (self.SAVE_U_PIXELS):
            model_output.open_new_ts_file( self, self.u_ts_file, IDs,
                                           var_name='u',
                                           long_name='mean_channel_flow_velocity',
                                           units_name='m/s')
                                          
        if (self.SAVE_D_PIXELS):    
            model_output.open_new_ts_file( self, self.d_ts_file, IDs,
                                           var_name='d',
                                           long_name='max_channel_flow_depth',
                                           units_name='m')
            
        if (self.SAVE_F_PIXELS):    
            model_output.open_new_ts_file( self, self.f_ts_file, IDs,
                                           var_name='f',
                                           long_name='friction_factor',
                                           units_name='none')
        
    #   open_output_files()
    #-------------------------------------------------------------------  
    def write_output_files(self, time_seconds=None):

        #---------------------------------------------------------
        # Notes:  This function was written to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read by read_cfg_file().
        #
        #         read_cfg_file() makes sure that all of
        #         the "save_dts" are larger than or equal to the
        #         process dt.
        #---------------------------------------------------------
        
        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
##        if ((self.time_index % self.grid_save_step) == 0):
##             self.save_grids()
##        if ((self.time_index % self.pixel_save_step) == 0):
##             self.save_pixel_values()
        
    #   write_output_files()
    #-------------------------------------------------------------------  
    def close_output_files(self):

        if (self.SAVE_Q_GRIDS):  model_output.close_gs_file( self, 'Q')   
        if (self.SAVE_U_GRIDS):  model_output.close_gs_file( self, 'u')  
        if (self.SAVE_D_GRIDS):  model_output.close_gs_file( self, 'd')   
        if (self.SAVE_F_GRIDS):  model_output.close_gs_file( self, 'f')
        #---------------------------------------------------------------
        if (self.SAVE_Q_PIXELS): model_output.close_ts_file( self, 'Q')   
        if (self.SAVE_U_PIXELS): model_output.close_ts_file( self, 'u')    
        if (self.SAVE_D_PIXELS): model_output.close_ts_file( self, 'd')    
        if (self.SAVE_F_PIXELS): model_output.close_ts_file( self, 'f')
        
    #   close_output_files()              
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #-----------------------------------
        # Save grid stack to a netCDF file
        #---------------------------------------------
        # Note that add_grid() methods will convert
        # var from scalar to grid now, if necessary.
        #---------------------------------------------        
        if (self.SAVE_Q_GRIDS):
            model_output.add_grid( self, self.Q, 'Q', self.time_min )
            
        if (self.SAVE_U_GRIDS):
            model_output.add_grid( self, self.u, 'u', self.time_min )
            
        if (self.SAVE_D_GRIDS):
            model_output.add_grid( self, self.d, 'd', self.time_min )

        if (self.SAVE_F_GRIDS):
            model_output.add_grid( self, self.f, 'f', self.time_min )     

    #   save_grids()
    #-------------------------------------------------------------------  
    def save_pixel_values(self):   ##### save_time_series_data(self)  #######
        
        IDs  = self.outlet_IDs
        time = self.time_min       #####

        #-------------
        # New method
        #-------------
        if (self.SAVE_Q_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Q, 'Q', IDs )
                    
        if (self.SAVE_U_PIXELS):
            model_output.add_values_at_IDs( self, time, self.u, 'u', IDs )
            
        if (self.SAVE_D_PIXELS):
            model_output.add_values_at_IDs( self, time, self.d, 'd', IDs )
            
        if (self.SAVE_F_PIXELS):
            model_output.add_values_at_IDs( self, time, self.f, 'f', IDs )
        
    #   save_pixel_values()
    #-------------------------------------------------------------------
    def manning_formula(self, S):

        #---------------------------------------------------------
        # Notes: R = (A/P) = hydraulic radius [m]
        #        N = Manning's roughness coefficient
        #            (usually in the range 0.012 to 0.035)
        #        S = bed slope (assumed equal to friction slope)

        #        R,S, and N may be 2D arrays.

        #        If length units are all *feet*, then an extra
        #        factor of 1.49 must be applied.  If units are
        #        meters, no such factor is needed.

        #        Note that Q = Ac * u, where Ac is cross-section
        #        area.  For a trapezoid, Ac does not equal w*d.
        #---------------------------------------------------------
        ##  if (N == None): N = np.float64(0.03)
        
        u = (self.Rh ** self.two_thirds) * np.sqrt(S) / self.nval
        
        #------------------------------
        # Add a hydraulic jump option
        # for when u gets too big ??
        #------------------------------
        
        #--------------------------------------
        # Option to return friction factor, f
        # (Correct this for Rh vs. d.)
        # **** See update_velocity_dynamic()
        #--------------------------------------
        #*** self.f = (self.g * self.d * self.S) / (u^2d)
        
        return u
    
    #   manning_formula()
    #-------------------------------------------------------------------
    def law_of_the_wall(self, S):

        #---------------------------------------------------------
        # Notes: u  = flow velocity  [m/s]
        #        d  = flow depth [m]
        #        z0 = roughness height
        #        S  = bed slope (assumed equal to friction slope)

        #        g     = 9.81 = gravitation constant [m/s^2]
        #        kappa = 0.41 = von Karman's constant
        #        aval  = 0.48 = integration constant

        #        sqrt(g)/kappa = 7.6393d
        #        smoothness = (aval / z0) * d
        #        f = (kappa / alog(smoothness))^2d
        #        tau_bed = rho_w * f * u^2 = rho_w * g * d * S

        #        d, S, and z0 can be arrays.

        #        To make default z0 correspond to default
        #        Manning's n, can use this approximation:
        #        z0 = a * (2.34 * sqrt(9.81) * n / kappa)^6d
        #        For n=0.03, this gives: z0 = 0.011417
        #        However, for n=0.3, it gives: z0 = 11417.413
        #        which is 11.4 km!  So the approximation only
        #        holds within some range of values.
        #--------------------------------------------------------
##        if (self.z0val == None):    
##            self.z0val = np.float64(0.011417)   # (about 1 cm)

        smoothness = (self.aval / self.z0val) * self.d
          
        #-----------------------------
        # Make sure (smoothness > 1)
        #-----------------------------
        smoothness = np.maximum(smoothness, np.float64(1.1))

        u = self.law_const * np.sqrt(self.Rh * S) * np.log(smoothness)
                   
        ## d = self.d
        ## u = np.float64(7.676) * np.sqrt(d * S) * np.log(smoothness)
        
        #------------------------------
        # Add a hydraulic jump option
        # for when u gets too big ??
        #------------------------------
        
        #--------------------------------------
        # Option to return friction factor, f
        # **** See update_velocity_dynamic()
        #--------------------------------------
        ## self.f = (self.kappa / np.log(smoothness)) ** 2
        
        return u
    
    #   law_of_the_wall()
    #-------------------------------------------------------------------
    def print_status_report(self): 

        #----------------------------------------------------
        # Wherever depth is less than z0, assume that water
        # is not flowing and set u and Q to zero.
        # However, we also need (d gt 0) to avoid a divide
        # by zero problem, even when numerators are zero.
        #----------------------------------------------------
        # FLOWING = (d > (z0/aval))
        #*** FLOWING[noflow_IDs] = False    ;******
        
        wflow    = np.where( FLOWING != 0 )
        n_flow   = np.size( wflow[0] )
        n_pixels = self.rti.n_pixels
        percent  = np.float64(100.0) * (np.float64(n_flow) / n_pixels)
        fstr = ('%5.1f' % percent) + '%'
        # fstr = idl_func.string(percent, format='(F5.1)').strip() + '%'
        print ' Percentage of pixels with flow = ' + fstr
        print ' '

        self.update_mins_and_maxes(REPORT=True)
 
        wmax  = np.where(self.Q == self.Q_max)
        nwmax = np.size(wmax[0])
        print ' Max(Q) occurs at: ' + str( wmax[0] )
        #print,' Max attained at ', nwmax, ' pixels.'
        print ' '
        print '-------------------------------------------------'

    #   print_status_report()         
    #-------------------------------------------------------------------
    def remove_bad_slopes(self, FLOAT=False):

        #------------------------------------------------------------
        # Notes: The main purpose of this routine is to find
        #        pixels that have nonpositive slopes and replace
        #        then with the smallest value that occurs anywhere
        #        in the input slope grid.  For example, pixels on
        #        the edges of the DEM will have a slope of zero.

        #        With the Kinematic Wave option, flow cannot leave
        #        a pixel that has a slope of zero and the depth
        #        increases in an unrealistic manner to create a
        #        spike in the depth grid.

        #        It would be better, of course, if there were
        #        no zero-slope pixels in the DEM.  We could use
        #        an "Imposed gradient DEM" to get slopes or some
        #        method of "profile smoothing".

        #        It is possible for the flow code to be nonzero
        #        at a pixel that has NaN for its slope. For these
        #        pixels, we also set the slope to our min value.

        #        7/18/05. Broke this out into separate procedure.
        #------------------------------------------------------------

        #-----------------------------------
        # Are there any "bad" pixels ?
        # If not, return with no messages.
        #-----------------------------------  
        wb = np.where(np.logical_or((self.slope <= 0.0), \
                              np.logical_not(np.isfinite(self.slope))))
        nbad = np.size(wb[0])
        print 'size(slope) =', np.size(self.slope)
        print 'size(wb) =', nbad
        
        wg = np.where(np.invert(np.logical_or((self.slope <= 0.0), \
                                     np.logical_not(np.isfinite(self.slope)))))
        ngood = np.size(wg[0])
        if (nbad == 0) or (ngood == 0):
            return
        
        #---------------------------------------------
        # Find smallest positive value in slope grid
        # and replace the "bad" values with smin.
        #---------------------------------------------
        print '-------------------------------------------------'
        print 'WARNING: Zero or negative slopes found.'
        print '         Replacing them with smallest slope.'
        print '         Use "Profile smoothing tool" instead.'
        S_min = self.slope[wg].min()
        S_max = self.slope[wg].max()
        print '         min(S) = ' + str(S_min)
        print '         max(S) = ' + str(S_max)
        print '-------------------------------------------------'
        print ' '
        self.slope[wb] = S_min
        
        #--------------------------------
        # Convert data type to double ?
        #--------------------------------
        if (FLOAT):    
            self.slope = np.float32(self.slope)
        else:    
            self.slope = np.float64(self.slope)
        
    #   remove_bad_slopes
    #-------------------------------------------------------------------

#-------------------------------------------------------------------
def Trapezoid_Rh(d, wb, theta):

    #-------------------------------------------------------------
    # Notes: Compute the hydraulic radius of a trapezoid that:
    #          (1) has a bed width of wb >= 0 (0 for triangular)
    #          (2) has a bank angle of theta (0 for rectangular)
    #          (3) is filled with water to a depth of d.
    #        The units of wb and d are meters.  The units of
    #        theta are assumed to be degrees and are converted.
    #-------------------------------------------------------------
    # NB!    wb should never be zero, so PW can never be 0,
    #        which would produce a NaN (divide by zero).
    #-------------------------------------------------------------
    #        See Notes for TF_Tan function in utils_TF.pro
    #            AW = d * (wb + (d * TF_Tan(theta_rad)) )
    #-------------------------------------------------------------    
    theta_rad = (theta * np.pi / 180.0)
    
    AW = d * (wb + (d * np.tan(theta_rad)) )      
    PW = wb + (np.float64(2) * d / np.cos(theta_rad) )
    Rh = (AW / PW)

    w  = np.where(wb <= 0)
    nw = np.size(w[0])
    
    return Rh

#   Trapezoid_Rh()
#-------------------------------------------------------------------
def Manning_Formula(Rh, S, nval):

    #---------------------------------------------------------
    # Notes: R = (A/P) = hydraulic radius [m]
    #        N = Manning's roughness coefficient
    #            (usually in the range 0.012 to 0.035)
    #        S = bed slope (assumed equal to friction slope)

    #        R,S, and N may be 2D arrays.

    #        If length units are all *feet*, then an extra
    #        factor of 1.49 must be applied.  If units are
    #        meters, no such factor is needed.

    #        Note that Q = Ac * u, where Ac is cross-section
    #        area.  For a trapezoid, Ac does not equal w*d.
    #---------------------------------------------------------
    ##  if (N == None): N = np.float64(0.03)

    two_thirds = np.float64(2) / 3.0
    
    u = (Rh ** two_thirds) * np.sqrt(S) / nval
    
    #------------------------------
    # Add a hydraulic jump option
    # for when u gets too big ??
    #------------------------------
    
    return u

#   Manning_Formula()
#-------------------------------------------------------------------
def Law_of_the_Wall(d, Rh, S, z0val):

    #---------------------------------------------------------
    # Notes: u  = flow velocity  [m/s]
    #        d  = flow depth [m]
    #        z0 = roughness height
    #        S  = bed slope (assumed equal to friction slope)

    #        g     = 9.81 = gravitation constant [m/s^2]
    #        kappa = 0.41 = von Karman's constant
    #        aval  = 0.48 = integration constant

    #        sqrt(g)/kappa = 7.6393d
    #        smoothness = (aval / z0) * d
    #        f = (kappa / alog(smoothness))^2d
    #        tau_bed = rho_w * f * u^2 = rho_w * g * d * S

    #        d, S, and z0 can be arrays.

    #        To make default z0 correspond to default
    #        Manning's n, can use this approximation:
    #        z0 = a * (2.34 * sqrt(9.81) * n / kappa)^6d
    #        For n=0.03, this gives: z0 = 0.011417
    #        However, for n=0.3, it gives: z0 = 11417.413
    #        which is 11.4 km!  So the approximation only
    #        holds within some range of values.
    #--------------------------------------------------------
##        if (self.z0val == None):    
##            self.z0val = np.float64(0.011417)   # (about 1 cm)

    #------------------------
    # Define some constants
    #------------------------
    g          = np.float64(9.81)    # (gravitation const.)
    aval       = np.float64(0.476)   # (integration const.)
    kappa      = np.float64(0.408)   # (von Karman's const.)
    law_const  = np.sqrt(g) / kappa
        
    smoothness = (aval / z0val) * d
      
    #-----------------------------
    # Make sure (smoothness > 1)
    #-----------------------------
    smoothness = np.maximum(smoothness, np.float64(1.1))

    u = law_const * np.sqrt(Rh * S) * np.log(smoothness)
    
    #------------------------------
    # Add a hydraulic jump option
    # for when u gets too big ??
    #------------------------------
    
    return u

#   Law_of_the_Wall()
#-------------------------------------------------------------------                 
