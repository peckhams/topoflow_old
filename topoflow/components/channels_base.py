
## See "d_bankfull" in update_flow_depth()  ######## (2/21/13)

## See "(5/13/10)" for a temporary fix.
#------------------------------------------------------------------------
#  Copyright (c) 2001-2014, Scott D. Peckham
#
#  Sep 2014.  Wrote new update_diversions().
#             New standard names and BMI updates and testing.
#  Nov 2013.  Converted TopoFlow to a Python package.
#  Feb 2013.  Adapted to use EMELI framework.
#  Jan 2013.  Shared scalar doubles are now 0D numpy arrays.
#             This makes them mutable and allows components with
#             a reference to them to see them change.
#             So far:  Q_outlet, Q_peak, Q_min...
#  Jan 2013.  Revised handling of input/output names.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  May 2012.  Commented out diversions.update() for now.  #######
#  May 2012.  Shared scalar doubles are now 1-element 1D numpy arrays.
#             This makes them mutable and allows components with
#             a reference to them to see them change.
#             So far:  Q_outlet, Q_peak, Q_min...
#  May 2010.  Changes to initialize() and read_cfg_file()
#  Mar 2010.  Changed codes to code, widths to width,
#             angles to angle, nvals to nval, z0vals to z0val,
#             slopes to slope (for GUI tools and consistency
#             across all process components)
#  Aug 2009.  Updates.
#  Jul 2009.  Updates.
#  May 2009.  Updates.
#  Jan 2009.  Converted from IDL.

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
#  NOTES:  update_free_surface_slope() is called by the
#          update_velocity() methods of channels_diffusive_wave.py
#          and channels_dynamic_wave.py.
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
#----------------------------------
#      initialize_d8_vars()          ########
#      initialize_computed_vars()
#      initialize_diversion_vars()      # (9/22/14)
#      initialize_outlet_values()
#      initialize_peak_values()
#      initialize_min_and_max_values()  # (2/3/13)
#-------------------------------------
#      update_R()
#      update_R_integral()
#      update_discharge()
#      update_diversions()          # (9/22/14)
#      update_flow_volume()
#      update_flow_depth()
#      update_free_surface_slope()
#      update_shear_stress()        # (9/9/14, depth-slope product)
#      update_shear_speed()         # (9/9/14)
#      update_trapezoid_Rh()
#      update_friction_factor()     # (9/9/14)
#----------------------------------
#      update_velocity()            # (override as needed)
#      update_velocity_on_edges()
#      update_froude_number()       # (9/9/14)
#----------------------------------
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
#      bundle_output_files()        # (9/21/14. Not used yet)
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
    # Note: rainfall_volume_flux *must* be liquid-only precip.
    #-----------------------------------------------------------        
    _input_var_names = [
        'atmosphere_water__rainfall_volume_flux',          # (P_rain)
        'glacier_ice__melt_volume_flux',                   # (MR)
        ## 'land_surface__elevation',
        ## 'land_surface__slope',
        'land_surface_water__baseflow_volume_flux',        # (GW)
        'land_surface_water__evaporation_volume_flux',     # (ET)
        'soil_surface_water__infiltration_volume_flux',    # (IN)
        'snowpack__melt_volume_flux',                      # (SM)
        'water-liquid__mass-per-volume_density' ]           # (rho_H2O)
        #------------------------------------------------------------------
#         'canals__count',                                   # n_canals
#         'canals_entrance__x_coordinate',                   # canals_in_x
#         'canals_entrance__y_coordinate',                   # canals_in_y
#         'canals_entrance_water__volume_fraction',          # Q_canals_fraction
#         'canals_exit__x_coordinate',                       # canals_out_x
#         'canals_exit__y_coordinate',                       # canals_out_y
#         'canals_exit_water__volume_flow_rate',             # Q_canals_out
#         'sinks__count',                                    # n_sinks
#         'sinks__x_coordinate',                             # sinks_x
#         'sinks__y_coordinate',                             # sinks_y
#         'sinks_water__volume_flow_rate',                   # Q_sinks
#         'sources__count',                                  # n_sources
#         'sources__x_coordinate',                           # sources_x
#         'sources__y_coordinate',                           # sources_y
#         'sources_water__volume_flow_rate' ]                # Q_sources
        
    #----------------------------------
    # Maybe add these out_vars later.
    #----------------------------------
    #  ['time_sec', 'time_min' ]
    
    _output_var_names = [
        'basin_outlet_water_flow__half_of_fanning_friction_factor',        # f_outlet
        'basin_outlet_water_x-section__mean_depth',                        # d_outlet
        'basin_outlet_water_x-section__peak_time_of_depth',                # Td_peak
        'basin_outlet_water_x-section__peak_time_of_volume_flow_rate',     # T_peak
        'basin_outlet_water_x-section__peak_time_of_volume_flux',          # Tu_peak
        'basin_outlet_water_x-section__time_integral_of_volume_flow_rate', # vol_Q
        'basin_outlet_water_x-section__time_max_of_mean_depth',            # d_peak
        'basin_outlet_water_x-section__time_max_of_volume_flow_rate',      # Q_peak
        'basin_outlet_water_x-section__time_max_of_volume_flux',           # u_peak
        'basin_outlet_water_x-section__volume_flow_rate',                  # Q_outlet
        'basin_outlet_water_x-section__volume_flux',                       # u_outlet
         #--------------------------------------------------
        'canals_entrance_water__volume_flow_rate',                         # Q_canals_in 
         #-------------------------------------------------- 
        'channel_bottom_surface__slope',                           # S_bed 
        'channel_bottom_water_flow__domain_max_of_log_law_roughness_length',  # z0val_max
        'channel_bottom_water_flow__domain_min_of_log_law_roughness_length',  # z0val_min
        'channel_bottom_water_flow__log_law_roughness_length',                # z0val
        'channel_bottom_water_flow__magnitude_of_shear_stress',    # tau
        'channel_bottom_water_flow__shear_speed',                  # u_star
        'channel_centerline__sinuosity',                           # sinu
        'channel_water__volume',                                   # vol
        'channel_water_flow__froude_number',                       # froude
        'channel_water_flow__half_of_fanning_friction_factor',     # f
        'channel_water_flow__domain_max_of_manning_n_parameter',   # nval_max
        'channel_water_flow__domain_min_of_manning_n_parameter',   # nval_min
        'channel_water_flow__manning_n_parameter',                 # nval
        'channel_water_surface__slope',                            # S_free
        #---------------------------------------------------
        # These might only be available at the end of run.
        #---------------------------------------------------
        'channel_water_x-section__domain_max_of_mean_depth',       # d_max
        'channel_water_x-section__domain_min_of_mean_depth',       # d_min
        'channel_water_x-section__domain_max_of_volume_flow_rate', # Q_max
        'channel_water_x-section__domain_min_of_volume_flow_rate', # Q_min
        'channel_water_x-section__domain_max_of_volume_flux',      # u_max
        'channel_water_x-section__domain_min_of_volume_flux',      # u_min
        #---------------------------------------------------------------------    
        'channel_water_x-section__hydraulic_radius',               # Rh
        'channel_water_x-section__initial_mean_depth',             # d0
        'channel_water_x-section__mean_depth',                     # d
        'channel_water_x-section__volume_flow_rate',               # Q  
        'channel_water_x-section__volume_flux',                    # u
        'channel_water_x-section__wetted_area',                    # A_wet
        'channel_water_x-section__wetted_perimeter',               # P_wet
        ## 'channel_water_x-section_top__width',                   # (not used)
        'channel_x-section_trapezoid_bottom__width',               # width
        'channel_x-section_trapezoid_side__flare_angle',           # angle
        'land_surface_water__runoff_volume_flux',                  # R  
        'land_surface_water__domain_time_integral_of_runoff_volume_flux', # vol_R     
        'model__time_step',                                        # dt
        'model_grid_cell__area' ]                                  # da
        
    _var_name_map = {
        'atmosphere_water__rainfall_volume_flux':              'P_rain',
        'glacier_ice__melt_volume_flux':                       'MR',
        ## 'land_surface__elevation':                          'DEM',
        ## 'land_surface__slope':                              'S_bed',
        'land_surface_water__baseflow_volume_flux':            'GW',
        'land_surface_water__evaporation_volume_flux':         'ET',
        'soil_surface_water__infiltration_volume_flux':        'IN',
        'snowpack__melt_volume_flux':                          'SM',
        'water-liquid__mass-per-volume_density':               'rho_H2O',
        #------------------------------------------------------------------------
        'basin_outlet_water_flow__half_of_fanning_friction_factor':'f_outlet',
        'basin_outlet_water_x-section__mean_depth':                'd_outlet',
        'basin_outlet_water_x-section__peak_time_of_depth':            'Td_peak',
        'basin_outlet_water_x-section__peak_time_of_volume_flow_rate': 'T_peak',
        'basin_outlet_water_x-section__peak_time_of_volume_flux':      'Tu_peak',
        'basin_outlet_water_x-section__volume_flow_rate':            'Q_outlet',
        'basin_outlet_water_x-section__volume_flux':                 'u_outlet',
        'basin_outlet_water_x-section__time_integral_of_volume_flow_rate': 'vol_Q',
        'basin_outlet_water_x-section__time_max_of_mean_depth':      'd_peak',
        'basin_outlet_water_x-section__time_max_of_volume_flow_rate':'Q_peak',
        'basin_outlet_water_x-section__time_max_of_volume_flux':     'u_peak',
        #--------------------------------------------------------------------------
        'canals_entrance_water__volume_flow_rate':                 'Q_canals_in', 
        #--------------------------------------------------------------------------    
        'channel_bottom_surface__slope':                           'S_bed',
        'channel_bottom_water_flow__domain_max_of_log_law_roughness_length': 'z0val_max',
        'channel_bottom_water_flow__domain_min_of_log_law_roughness_length': 'z0val_min',
        'channel_bottom_water_flow__log_law_roughness_length':     'z0val',
        'channel_bottom_water_flow__magnitude_of_shear_stress':    'tau',
        'channel_bottom_water_flow__shear_speed':                  'u_star',
        'channel_centerline__sinuosity':                           'sinu',
        'channel_water__volume':                                   'vol',
        'channel_water_flow__domain_max_of_manning_n_parameter':   'nval_max',
        'channel_water_flow__domain_min_of_manning_n_parameter':   'nval_min',
        'channel_water_flow__froude_number':                       'froude',
        'channel_water_flow__half_of_fanning_friction_factor':     'f',
        'channel_water_flow__manning_n_parameter':                 'nval',
        'channel_water_surface__slope':                            'S_free',
        #-----------------------------------------------------------------------
        'channel_water_x-section__domain_max_of_mean_depth':       'd_max',
        'channel_water_x-section__domain_min_of_mean_depth':       'd_min',
        'channel_water_x-section__domain_max_of_volume_flow_rate': 'Q_max',
        'channel_water_x-section__domain_min_of_volume_flow_rate': 'Q_min',
        'channel_water_x-section__domain_max_of_volume_flux':      'u_max',
        'channel_water_x-section__domain_min_of_volume_flux':      'u_min',
        #-----------------------------------------------------------------------      
        'channel_water_x-section__hydraulic_radius':               'Rh',
        'channel_water_x-section__initial_mean_depth':             'd0',
        'channel_water_x-section__mean_depth':                     'd',
        'channel_water_x-section__volume_flow_rate':               'Q',                
        'channel_water_x-section__volume_flux':                    'u',
        'channel_water_x-section__wetted_area':                    'A_wet',
        'channel_water_x-section__wetted_perimeter':               'P_wet',
        ## 'channel_water_x-section_top__width':                   # (not used)
        'channel_x-section_trapezoid_bottom__width':               'width',   ####
        'channel_x-section_trapezoid_side__flare_angle':           'angle',   ####
        'land_surface_water__domain_time_integral_of_runoff_volume_flux': 'vol_R',
        'land_surface_water__runoff_volume_flux':                  'R',
        'model__time_step':                                        'dt',
        'model_grid_cell__area':                                   'da',
        #------------------------------------------------------------------
        'canals__count':                          'n_canals',
        'canals_entrance__x_coordinate':          'canals_in_x',
        'canals_entrance__y_coordinate':          'canals_in_y',
        'canals_entrance_water__volume_fraction': 'Q_canals_fraction',
        'canals_exit__x_coordinate':              'canals_out_x',
        'canals_exit__y_coordinate':              'canals_out_y',
        'canals_exit_water__volume_flow_rate':    'Q_canals_out',
        'sinks__count':                           'n_sinks',
        'sinks__x_coordinate':                    'sinks_x',
        'sinks__y_coordinate':                    'sinks_y',
        'sinks_water__volume_flow_rate':          'Q_sinks',
        'sources__count':                         'n_sources',
        'sources__x_coordinate':                  'sources_x',
        'sources__y_coordinate':                  'sources_y',
        'sources_water__volume_flow_rate':        'Q_sources' }


    #------------------------------------------------
    # Create an "inverse var name map"
    # inv_map = dict(zip(map.values(), map.keys()))
    #------------------------------------------------
##    _long_name_map = dict( zip(_var_name_map.values(),
##                               _var_name_map.keys() ) )

    _var_units_map = {
        'atmosphere_water__rainfall_volume_flux':              'm s-1',
        'glacier_ice__melt_volume_flux':                       'm s-1',
        ## 'land_surface__elevation':                          'm',
        ## 'land_surface__slope':                              '1',
        'land_surface_water__baseflow_volume_flux':            'm s-1',
        'land_surface_water__evaporation_volume_flux':         'm s-1',
        'soil_surface_water__infiltration_volume_flux':        'm s-1',
        'snowpack__melt_volume_flux':                          'm s-1',
        'water-liquid__mass-per-volume_density':               'kg m-3',
        #--------------------------------------------------------------------------- 
        'basin_outlet_water_flow__half_of_fanning_friction_factor':        '1',       
        'basin_outlet_water_x-section__mean_depth':                        'm',
        'basin_outlet_water_x-section__peak_time_of_depth':                'min',
        'basin_outlet_water_x-section__peak_time_of_volume_flow_rate':     'min',
        'basin_outlet_water_x-section__peak_time_of_volume_flux':          'min',        
        'basin_outlet_water_x-section__time_integral_of_volume_flow_rate': 'm3',
        'basin_outlet_water_x-section__time_max_of_mean_depth':            'm',
        'basin_outlet_water_x-section__time_max_of_volume_flow_rate':      'm3 s-1',
        'basin_outlet_water_x-section__time_max_of_volume_flux':           'm s-1',
        'basin_outlet_water_x-section__volume_flow_rate':                  'm3',
        'basin_outlet_water_x-section__volume_flux':                       'm s-1',
        #---------------------------------------------------------------------------
        'canals_entrance_water__volume_flow_rate':                 'm3 s-1', 
        #---------------------------------------------------------------------------  
        'channel_bottom_surface__slope':                           '1',
        'channel_bottom_water_flow__domain_max_of_log_law_roughness_length':  'm',
        'channel_bottom_water_flow__domain_min_of_log_law_roughness_length':  'm',
        'channel_bottom_water_flow__log_law_roughness_length':     'm',
        'channel_bottom_water_flow__magnitude_of_shear_stress':    'kg m-1 s-2',
        'channel_bottom_water_flow__shear_speed':                  'm s-1',
        'channel_centerline__sinuosity':                           '1',    
        'channel_water__volume':                                   'm3', 
        'channel_water_flow__froude_number':                       '1',
        'channel_water_flow__half_of_fanning_friction_factor':     '1',               
        'channel_water_flow__manning_n_parameter':                 'm-1/3 s',
        'channel_water_flow__domain_max_of_manning_n_parameter':   'm-1/3 s',
        'channel_water_flow__domain_min_of_manning_n_parameter':   'm-1/3 s',
        'channel_water_surface__slope':                            '1',
        #--------------------------------------------------------------------
        'channel_water_x-section__domain_max_of_mean_depth':       'm',
        'channel_water_x-section__domain_min_of_mean_depth':       'm',
        'channel_water_x-section__domain_max_of_volume_flow_rate': 'm3 s-1',
        'channel_water_x-section__domain_min_of_volume_flow_rate': 'm3 s-1',
        'channel_water_x-section__domain_max_of_volume_flux':      'm s-1',
        'channel_water_x-section__domain_min_of_volume_flux':      'm s-1',
        #--------------------------------------------------------------------
        'channel_water_x-section__hydraulic_radius':               'm',
        'channel_water_x-section__initial_mean_depth':             'm',
        'channel_water_x-section__mean_depth':                     'm',
        'channel_water_x-section__volume_flow_rate':               'm3 s-1',
        'channel_water_x-section__volume_flux':                    'm s-1',
        'channel_water_x-section__wetted_area':                    'm2',
        'channel_water_x-section__wetted_perimeter':               'm',
        'channel_x-section_trapezoid_bottom__width':               'm',
        'channel_x-section_trapezoid_side__flare_angle':           'rad', # CHECKED 
        'land_surface_water__domain_time_integral_of_runoff_volume_flux': 'm3',  
        'land_surface_water__runoff_volume_flux':                  'm s-1',        
        'model__time_step':                                        's',
        'model_grid_cell__area':                                   'm2',
        #------------------------------------------------------------------
        'canals__count':                          '1',
        'canals_entrance__x_coordinate':          'm',
        'canals_entrance__y_coordinate':          'm',
        'canals_entrance_water__volume_fraction': '1',
        'canals_exit__x_coordinate':              'm',
        'canals_exit__y_coordinate':              'm',
        'canals_exit_water__volume_flow_rate':    'm3 s-1',
        'sinks__count':                           '1',
        'sinks__x_coordinate':                    'm',
        'sinks__y_coordinate':                    'm',
        'sinks_water__volume_flow_rate':          'm3 s-1',
        'sources__count':                         '1',
        'sources__x_coordinate':                  'm',
        'sources__y_coordinate':                  'm',
        'sources_water__volume_flow_rate':        'm3 s-1' }

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
##        return 'float64'
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
    def initialize(self, cfg_file=None, mode="nondriver", SILENT=False): 

        if not(SILENT):
            print ' '
            print 'Channels component: Initializing...'
        
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
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

        # Attempt to define default values for running channels
        # components in standalone mode. (@mdpiper, 8/19/15)
        try:
            self.P_rain
        except AttributeError:
            self.P_rain = np.float64(0.0)
        try:
            self.SM
        except AttributeError:
            self.SM = np.float64(0.0)
        try:
            self.ET
        except AttributeError:
            self.ET = np.float64(0.0)
        try:
            self.GW
        except AttributeError:
            self.GW = np.float64(0.0)
        try:
            self.IN
        except AttributeError:
            self.IN = np.float64(0.0)
        try:
            self.MR
        except AttributeError:
            self.MR = np.float64(0.0)
        try:
            self.rho_H2O
        except AttributeError:
            self.rho_H2O = np.float64(1000.0)  # met_base.py:525
        try:
            self.S_free
        except AttributeError:
            self.S_free = np.zeros(self.slope.shape)  # DEM shape from rti file

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
        if (DEBUG): print '#### Calling update_diversions()...'
        self.update_diversions()
        if (DEBUG): print '#### Calling update_flow_volume()...'
        self.update_flow_volume()
        if (DEBUG): print '#### Calling update_flow_depth()...'
        self.update_flow_depth()
        #-----------------------------------------------------------------
        if not(self.DYNAMIC_WAVE):
            if (DEBUG): print '#### Calling update_trapezoid_Rh()...'
            self.update_trapezoid_Rh()
            # print 'Rhmin, Rhmax =', self.Rh.min(), self.Rh.max()a
        #-----------------------------------------------------------------
        # (9/9/14) Moved this here from update_velocity() methods.
        #-----------------------------------------------------------------        
        if not(self.KINEMATIC_WAVE):
            if (DEBUG): print '#### Calling update_free_surface_slope()...' 
            self.update_free_surface_slope()
        if (DEBUG): print '#### Calling update_shear_stress()...'
        self.update_shear_stress()
        if (DEBUG): print '#### Calling update_shear_speed()...'
        self.update_shear_speed()  
        #-----------------------------------------------------------------
        # Must update friction factor before velocity for DYNAMIC_WAVE.
        #-----------------------------------------------------------------        
        if (DEBUG): print '#### Calling update_friction_factor()...'
        self.update_friction_factor()      
        #-----------------------------------------------------------------          
        if (DEBUG): print '#### Calling update_velocity()...'
        self.update_velocity()
        self.update_velocity_on_edges()     # (set to zero)
        if (DEBUG): print '#### Calling update_froude_number()...'
        self.update_froude_number()
        #-----------------------------------------------------------------
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

        #---------------------------------------------------------------    
        # Note: The initialize() method calls initialize_config_vars()
        #       (in BMI_base.py), which calls this method at the end.
        #--------------------------------------------------------------
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
        
        #--------------------------------------------------         
        # D8 component builds its cfg filename from these  
        #--------------------------------------------------      
        self.d8.site_prefix  = self.site_prefix
        self.d8.in_directory = self.in_directory
        self.d8.initialize( cfg_file=None,
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

        #---------------------------------------------------
        # Initialize new grids. Is this needed?  (9/13/14)
        #---------------------------------------------------
        self.tau    = np.zeros([self.ny, self.nx], dtype='Float64')
        self.u_star = np.zeros([self.ny, self.nx], dtype='Float64')
        self.froude = np.zeros([self.ny, self.nx], dtype='Float64')
                        
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
        #-----------------------------------------------------------
        # Note: angles were read as degrees & converted to radians
        #-----------------------------------------------------------
        L2         = self.d * np.tan(self.angle)
        self.A_wet = self.d * (self.width + L2)
        self.P_wet = self.width + (np.float64(2) * self.d / np.cos(self.angle) )
        self.vol   = self.A_wet * self.d8.ds   # [m3]

        #-------------------------------------------------------        
        # Note: depth is often zero at the start of a run, and
        # both width and then P_wet are also zero in places.
        # Therefore initialize Rh as shown.
        #-------------------------------------------------------
        self.Rh = np.zeros([self.ny, self.nx], dtype='Float64')
        ## self.Rh = self.A_wet / self.P_wet   # [m]
        ## print 'P_wet.min() =', self.P_wet.min()
        ## print 'width.min() =', self.width.min()
       
        ## self.initialize_diversion_vars()    # (9/22/14)
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
    #-------------------------------------------------------------
    def initialize_diversion_vars(self):

        #-----------------------------------------
        # Compute source IDs from xy coordinates
        #-----------------------------------------
        source_rows     = np.int32( self.sources_y / self.ny )
        source_cols     = np.int32( self.sources_x / self.nx )
        self.source_IDs = (source_rows, source_cols)
        ## self.source_IDs = (source_rows * self.nx) + source_cols
   
        #---------------------------------------
        # Compute sink IDs from xy coordinates
        #---------------------------------------
        sink_rows     = np.int32( self.sinks_y / self.ny )
        sink_cols     = np.int32( self.sinks_x / self.nx )
        self.sink_IDs = (sink_rows, sink_cols)
        ## self.sink_IDs = (sink_rows * self.nx) + sink_cols
        
        #-------------------------------------------------
        # Compute canal entrance IDs from xy coordinates
        #-------------------------------------------------
        canal_in_rows     = np.int32( self.canals_in_y / self.ny )
        canal_in_cols     = np.int32( self.canals_in_x / self.nx )
        self.canal_in_IDs = (canal_in_rows, canal_in_cols)
        ## self.canal_in_IDs = (canal_in_rows * self.nx) + canal_in_cols
        
        #---------------------------------------------
        # Compute canal exit IDs from xy coordinates
        #---------------------------------------------
        canal_out_rows     = np.int32( self.canals_out_y / self.ny )
        canal_out_cols     = np.int32( self.canals_out_x / self.nx )
        self.canal_out_IDs = (canal_out_rows, canal_out_cols)
        ## self.canal_out_IDs = (canal_out_rows * self.nx) + canal_out_cols

        #--------------------------------------------------
        # This will be computed from Q_canal_fraction and
        # self.Q and then passed back to Diversions
        #--------------------------------------------------
        self.Q_canals_in = np.array( self.n_sources, dtype='float64' )

    #   initialize_diversion_vars()
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
        P  = self.P_rain  # (This is now liquid-only precip. 9/14/14)
        SM = self.SM
        GW = self.GW
        ET = self.ET
        IN = self.IN
        MR = self.MR
        
##        if (self.DEBUG):
##            print 'At time:', self.time_min, ', P =', P, '[m/s]'

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
        #------------------------------------------------------ 
        # A_wet is initialized in initialize_computed_vars().
        # A_wet is updated in update_trapezoid_Rh().
        #------------------------------------------------------     
        ### self.Q = np.float64(self.u * A_wet)
        self.Q[:] = self.u * self.A_wet   ## (2/19/13, in place)

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
    def update_diversions(self):

        #--------------------------------------------------------------    
        # Note: The Channel component requests the following input
        #       vars from the Diversions component by including
        #       them in its "get_input_vars()":
        #       (1) Q_sources, Q_sources_x, Q_sources_y
        #       (2) Q_sinks,   Q_sinks_x, Q_sinks_y
        #       (3) Q_canals_out, Q_canals_out_x, Q_canals_out_y
        #       (4) Q_canals_fraction, Q_canals_in_x, Q_canals_in_y.
        
        #       source_IDs are computed from (x,y) coordinates during
        #       initialize().
        #
        #       Diversions component needs to get Q_canals_in from the
        #       Channel component.
        #--------------------------------------------------------------
        # Note: This *must* be called after update_discharge() and
        #       before update_flow_volume().
        #--------------------------------------------------------------
        # Note: The Q grid stores the volume flow rate *leaving* each
        #       grid cell in the domain.  For sources, an extra amount
        #       is leaving the cell which can flow into its D8 parent
        #       cell.  For sinks, a lesser amount is leaving the cell
        #       toward the D8 parent.
        #--------------------------------------------------------------
        # Note: It is not enough to just update Q and then call the
        #       update_flow_volume() method.  This is because it
        #       won't update the volume in the channels in the grid
        #       cells that the extra discharge is leaving from.
        #--------------------------------------------------------------
        # If a grid cell contains a "source", then an additional Q
        # will flow *into* that grid cell and increase flow volume.
        #-------------------------------------------------------------- 

        #-------------------------------------------------------------         
        # This is not fully tested but runs.  However, the Diversion
        # vars are still computed even when Diversions component is
        # disabled. So it slows things down somewhat.
        #-------------------------------------------------------------              
        return
        ########################
        ########################
        
        #----------------------------------------            
        # Update Q and vol due to point sources
        #----------------------------------------
        ## if (hasattr(self, 'source_IDs')): 
        if (self.n_sources > 0): 
			self.Q[ self.source_IDs ]   += self.Q_sources
			self.vol[ self.source_IDs ] += (self.Q_sources * self.dt)

        #--------------------------------------            
        # Update Q and vol due to point sinks
        #--------------------------------------
        ## if (hasattr(self, 'sink_IDs')):
        if (self.n_sinks > 0): 
			self.Q[ self.sink_IDs ]   -= self.Q_sinks
			self.vol[ self.sink_IDs ] -= (self.Q_sinks * self.dt)
 
        #---------------------------------------            
        # Update Q and vol due to point canals
        #---------------------------------------    
        ## if (hasattr(self, 'canal_in_IDs')):
        if (self.n_canals > 0):   
			#-----------------------------------------------------------------
			# Q grid was just modified.  Apply the canal diversion fractions
			# to compute the volume flow rate into upstream ends of canals.
			#-----------------------------------------------------------------
			Q_canals_in = self.Q_canals_fraction * self.Q[ self.canal_in_IDs ]
			self.Q_canals_in = Q_canals_in

			#----------------------------------------------------        
			# Update Q and vol due to losses at canal entrances
			#----------------------------------------------------
			self.Q[ self.canal_in_IDs ]   -= Q_canals_in
			self.vol[ self.canal_in_IDs ] -= (Q_canals_in * self.dt)        

			#-------------------------------------------------       
			# Update Q and vol due to gains at canal exits.
			# Diversions component accounts for travel time.
			#-------------------------------------------------        
			self.Q[ self.canal_out_IDs ]   += self.Q_canals_out
			self.vol[ self.canal_out_IDs ] += (self.Q_canals_out * self.dt)    
        
    #   update_diversions()
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
        # Contributions over entire grid cell from rainfall,
        # snowmelt, icemelt and baseflow (minus losses from
        # evaporation and infiltration) are assumed to flow
        # into the channel within the grid cell.
        # Note that R is allowed to be negative.
        #----------------------------------------------------        
        self.vol += (self.R * self.da) * dt   # (in place)
    
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
        self.vol -= (self.Q * dt)  # (in place)
   
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
        #        Modified to avoid double where calls, which
        #        reduced cProfile run time for this method from
        #        1.391 to 0.644.  (9/23/14)
        #----------------------------------------------------------
        # Commented this out on (2/18/10) because it doesn't
        #           seem to be used anywhere now.  Checked all
        #           of the Channels components.
        #----------------------------------------------------------        
        # self.d_last = self.d.copy()

        #-----------------------------------        
        # Make some local aliases and vars
        #-----------------------------------------------------------
        # Note: angles were read as degrees & converted to radians
        #-----------------------------------------------------------
        d = self.d
        width = self.width  ###
        angle = self.angle
        SCALAR_ANGLES = (np.size(angle) == 1)
        
        #------------------------------------------------------
        # (2/18/10) New code to deal with case where the flow
        #           depth exceeds a bankfull depth.
        #           For now, d_bankfull is hard-coded.
        #
        #           CHANGE Manning's n here, too?
        #------------------------------------------------------
        d_bankfull = 4.0  # [meters]
        ################################
        wb = (self.d > d_bankfull)  # (array of True or False)
        self.width[ wb ]  = self.d8.dw[ wb ]
        if not(SCALAR_ANGLES):
            self.angle[ wb ] = 0.0
     
#         w_overbank = np.where( d > d_bankfull )
#         n_overbank = np.size( w_overbank[0] )
#         if (n_overbank != 0):
#             width[ w_overbank ] = self.d8.dw[ w_overbank ]
#             if not(SCALAR_ANGLES): angle[w_overbank] = 0.0

        #------------------------------------------------------
        # (2/18/10) New code to deal with case where the top
        #           width exceeds the grid cell width, dw.
        #------------------------------------------------------            
        top_width = width + (2.0 * d * np.sin(self.angle))
        wb = (top_width > self.d8.dw)  # (array of True or False)
        self.width[ wb ] = self.d8.dw[ wb ]
        if not(SCALAR_ANGLES):
            self.angle[ wb ] = 0.0
                    
#         wb = np.where(top_width > self.d8.dw)
#         nb = np.size(w_bad[0])
#         if (nb != 0):
#             width[ wb ] = self.d8.dw[ wb ]
#             if not(SCALAR_ANGLES): angle[ wb ] = 0.0

        #----------------------------------
        # Is "angle" a scalar or a grid ?
        #----------------------------------
        if (SCALAR_ANGLES):
            if (angle == 0.0):    
                d = self.vol / (width * self.d8.ds)
            else:
                denom = 2.0 * np.tan(angle)
                arg   = 2.0 * denom * self.vol / self.d8.ds
                arg  += width**(2.0)
                d     = (np.sqrt(arg) - width) / denom
        else:
            #-----------------------------------------------------
            # Pixels where angle is 0 must be handled separately
            #-----------------------------------------------------
            w1 = ( angle == 0 )  # (arrays of True or False)
            w2 = np.invert( w1 )
            #-----------------------------------
            A_top = width[w1] * self.d8.ds[w1]          
            d[w1] = self.vol[w1] / A_top
            #-----------------------------------               
            denom  = 2.0 * np.tan(angle[w2])
            arg    = 2.0 * denom * self.vol[w2] / self.d8.ds[w2]
            arg   += width[w2]**(2.0)
            d[w2] = (np.sqrt(arg) - width[w2]) / denom

            #-----------------------------------------------------
            # Pixels where angle is 0 must be handled separately
            #-----------------------------------------------------
#             wz   = np.where( angle == 0 )
#             nwz  = np.size( wz[0] )
#             wzc  = np.where( angle != 0 )
#             nwzc = np.size( wzc[0] )
#             
#             if (nwz != 0):
#                 A_top = width[wz] * self.d8.ds[wz]
#                 ## A_top = self.width[wz] * self.d8.ds_chan[wz]            
#                 d[wz] = self.vol[wz] / A_top
#             
#             if (nwzc != 0):    
#                 term1  = 2.0 * np.tan(angle[wzc])
#                 arg    = 2.0 * term1 * self.vol[wzc] / self.d8.ds[wzc]
#                 arg   += width[wzc]**(2.0)
#                 d[wzc] = (np.sqrt(arg) - width[wzc]) / term1

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

        #-------------------------------------------------        
        # Find where d <= 0 and save for later (9/23/14)
        #-------------------------------------------------
        self.d_is_pos = (self.d > 0)
        self.d_is_neg = np.invert( self.d_is_pos )
        
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
    def update_shear_stress(self):

        #--------------------------------------------------------
        # Notes: 9/9/14.  Added so shear stress could be shared.
        #        This uses the depth-slope product.
        #--------------------------------------------------------
        if (self.KINEMATIC_WAVE):
        	slope = self.S_bed
        else:
            slope = self.S_free
        self.tau[:] = self.rho_H2O * self.g * self.d * slope
               
    #   update_shear_stress()
    #-------------------------------------------------------------------
    def update_shear_speed(self):

        #--------------------------------------------------------
        # Notes: 9/9/14.  Added so shear speed could be shared.
        #--------------------------------------------------------
        self.u_star[:] = np.sqrt( self.tau / self.rho_H2O )
               
    #   update_shear_speed()
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
        # NB!    wb should never be zero, so P_wet can never be 0,
        #        which would produce a NaN (divide by zero).
        #-------------------------------------------------------------
        #        See Notes for TF_Tan function in utils_TF.pro
        #            AW = d * (wb + (d * TF_Tan(theta_rad)) )
        #-------------------------------------------------------------
        # 9/9/14.  Bug fix.  Angles were already in radians but
        #          were converted to radians again.
        #--------------------------------------------------------------

        #---------------------------------------------------------
        # Compute hydraulic radius grid for trapezoidal channels
        #-----------------------------------------------------------
        # Note: angles were read as degrees & converted to radians
        #-----------------------------------------------------------
        d     = self.d        # (local synonyms)
        wb    = self.width    # (trapezoid bottom width)
        L2    = d * np.tan( self.angle )          
        A_wet = d * (wb + L2)      
        P_wet = wb + (np.float64(2) * d / np.cos(self.angle) )

        #---------------------------------------------------
        # At noflow_IDs (e.g. edges) P_wet may be zero
        # so do this to avoid "divide by zero". (10/29/11)
        #---------------------------------------------------
        P_wet[ self.d8.noflow_IDs ] = np.float64(1)
        Rh = (A_wet / P_wet)
        #--------------------------------
        # w = np.where(P_wet == 0)
        # print 'In update_trapezoid_Rh():'
        # print '   P_wet= 0 at', w[0].size, 'cells'

        #------------------------------------
        # Force edge pixels to have Rh = 0.
        # This will make u = 0 there also.
        #------------------------------------
        Rh[ self.d8.noflow_IDs ] = np.float64(0)        
##        w  = np.where(wb <= 0)
##        nw = np.size(w[0])
##        if (nw > 0): Rh[w] = np.float64(0)
        
        self.Rh[:]    = Rh
        self.A_wet[:] = A_wet   ## (Now shared: 9/9/14)
        self.P_wet[:] = P_wet   ## (Now shared: 9/9/14)

        #---------------
        # For testing
        #--------------
##        print 'dmin, dmax =', d.min(),  d.max()
##        print 'wmin, wmax =', wb.min(), wb.max()
##        print 'amin, amax =', self.angle.min(), self.angle.max()

    #   update_trapezoid_Rh()
    #-------------------------------------------------------------------
    def update_friction_factor(self):    

        #----------------------------------------    
        # Note:  Added on 9/9/14 to streamline.
        #----------------------------------------------------------
        # Note:  f  = half of the Fanning friction factor
        #        d  = flow depth [m]
        #        z0 = roughness length
        #        S  = bed slope (assumed equal to friction slope)
        #        g  = 9.81 = gravitation constant [m/s^2]
        #---------------------------------------------------------       
        #        For law of the wall:
        #        kappa = 0.41 = von Karman's constant
        #        aval  = 0.48 = integration constant

        #        law_const  = sqrt(g)/kappa = 7.6393d
        #        smoothness = (aval / z0) * d
        #        f = (kappa / alog(smoothness))^2d
        #        tau_bed = rho_w * f * u^2 = rho_w * g * d * S

        #        d, S, and z0 can be arrays.

        #        To make default z0 correspond to default
        #        Manning's n, can use this approximation:
        #        z0 = a * (2.34 * sqrt(9.81) * n / kappa)^6d
        #        For n=0.03, this gives: z0 = 0.011417
        #########################################################
        #        However, for n=0.3, it gives: z0 = 11417.413
        #        which is 11.4 km!  So the approximation only
        #        holds within some range of values.
        #--------------------------------------------------------

        ###############################################################
        # cProfile:  This method took: 0.369 secs for topoflow_test()
        ###############################################################            
        #--------------------------------------
        # Find where (d <= 0).  g=good, b=bad
        #-------------------------------------- 
        wg = self.d_is_pos
        wb = self.d_is_neg
#         wg = ( self.d > 0 )
#         wb = np.invert( wg )
        
        #-----------------------------
        # Compute f for Manning case
        #-----------------------------------------
		# This makes f=0 and du=0 where (d <= 0)
		#-----------------------------------------
        if (self.MANNING):
            n2 = self.nval ** np.float64(2)  
            self.f[ wg ] = self.g * (n2[wg] / (self.d[wg] ** self.one_third))
            self.f[ wb ] = np.float64(0)
 
        #---------------------------------
        # Compute f for Law of Wall case
        #---------------------------------
        if (self.LAW_OF_WALL):
            #------------------------------------------------
            # Make sure (smoothness > 1) before taking log.
            # Should issue a warning if this is used.
            #------------------------------------------------
            smoothness = (self.aval / self.z0val) * self.d
            np.maximum(smoothness, np.float64(1.1), smoothness)  # (in place)
            self.f[wg] = (self.kappa / np.log(smoothness[wg])) ** np.float64(2)
            self.f[wb] = np.float64(0)

        ##############################################################
        # cProfile:  This method took: 0.93 secs for topoflow_test()
        ##############################################################        
#         #--------------------------------------
#         # Find where (d <= 0).  g=good, b=bad
#         #-------------------------------------- 
#         wg = np.where( self.d > 0 )
#         ng = np.size( wg[0])
#         wb = np.where( self.d <= 0 )
#         nb = np.size( wb[0] )
# 
#         #-----------------------------
#         # Compute f for Manning case
#         #-----------------------------------------
# 		  # This makes f=0 and du=0 where (d <= 0)
# 		  #-----------------------------------------
#         if (self.MANNING):
#             n2 = self.nval ** np.float64(2)  
#             if (ng != 0):
#                 self.f[wg] = self.g * (n2[wg] / (self.d[wg] ** self.one_third))
#             if (nb != 0):
#                 self.f[wb] = np.float64(0)
# 
#         #---------------------------------
#         # Compute f for Law of Wall case
#         #---------------------------------
#         if (self.LAW_OF_WALL):
#             #------------------------------------------------
#             # Make sure (smoothness > 1) before taking log.
#             # Should issue a warning if this is used.
#             #------------------------------------------------
#             smoothness = (self.aval / self.z0val) * self.d
#             np.maximum(smoothness, np.float64(1.1), smoothness)  # (in place)
#             ## smoothness = np.maximum(smoothness, np.float64(1.1))
#             if (ng != 0):
#                 self.f[wg] = (self.kappa / np.log(smoothness[wg])) ** np.float64(2)
#             if (nb != 0):
#                 self.f[wb] = np.float64(0)                       

        #---------------------------------------------
        # We could share the Fanning friction factor
        #---------------------------------------------
        ### self.fanning = (np.float64(2) * self.f)

    #   update_friction_factor()       
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
    #-------------------------------------------------------------------
    def update_froude_number(self):

        #----------------------------------------------------------
        # Notes: 9/9/14.  Added so Froude number could be shared.
        # This use of wg & wb reduced cProfile time from:
        # 0.644 sec to: 0.121.  (9/23/14)
        #----------------------------------------------------------
        # g = good, b = bad
        #-------------------- 
        wg = self.d_is_pos
        wb = self.d_is_neg

        self.froude[ wg ] = self.u[wg] / np.sqrt( self.g * self.d[wg] )       
        self.froude[ wb ] = np.float64(0)
               
    #   update_froude_number()
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
        if (angle != None):
            #-----------------------------------------------
            # Convert bank angles from degrees to radians. 
            #-----------------------------------------------
            self.angle = angle * self.deg_to_rad  # [radians]
            ### self.angle = angle  # (before 9/9/14)

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
    def bundle_output_files(self):    

        ###################################################
        # NOT READY YET. Need "get_long_name()" and a new
        # version of "get_var_units".  (9/21/14)
        ###################################################
                
        #-------------------------------------------------------------       
        # Bundle the output file info into an array for convenience.
        # Then we just need one open_output_files(), in BMI_base.py,
        # and one close_output_files().  Less to maintain. (9/21/14)
        #-------------------------------------------------------------        
        # gs = grid stack, ts = time series, ps = profile series.
        #-------------------------------------------------------------
        self.out_files = [
        {var_name:'Q', 
        save_gs:self.SAVE_Q_GRIDS,  gs_file:self.Q_gs_file,
        save_ts:self.SAVE_Q_PIXELS, ts_file:self.Q_ts_file, 
        long_name:get_long_name('Q'), units_name:get_var_units('Q')}, 
        #-----------------------------------------------------------------
        {var_name:'u',
        save_gs:self.SAVE_U_GRIDS,  gs_file:self.u_gs_file,
        save_ts:self.SAVE_U_PIXELS, ts_file:self.u_ts_file,
        long_name:get_long_name('u'), units_name:get_var_units('u')},
        #-----------------------------------------------------------------
        {var_name:'d',
        save_gs:self.SAVE_D_GRIDS,  gs_file:self.d_gs_file,
        save_ts:self.SAVE_D_PIXELS, ts_file:self.d_ts_file,
        long_name:get_long_name('d'), units_name:get_var_units('d')}, 
        #-----------------------------------------------------------------
        {var_name:'f',
        save_gs:self.SAVE_F_GRIDS,  gs_file:self.f_gs_file,
        save_ts:self.SAVE_F_PIXELS, ts_file:self.f_ts_file,
        long_name:get_long_name('f'), units_name:get_var_units('f')} ]
                        
    #   bundle_output_files
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_netcdf()
        self.update_outfile_names()
        ## self.bundle_output_files()
        

##        print 'self.SAVE_Q_GRIDS =', self.SAVE_Q_GRIDS
##        print 'self.SAVE_U_GRIDS =', self.SAVE_U_GRIDS
##        print 'self.SAVE_D_GRIDS =', self.SAVE_D_GRIDS
##        print 'self.SAVE_F_GRIDS =', self.SAVE_F_GRIDS
##        #---------------------------------------------------
##        print 'self.SAVE_Q_PIXELS =', self.SAVE_Q_PIXELS
##        print 'self.SAVE_U_PIXELS =', self.SAVE_U_PIXELS
##        print 'self.SAVE_D_PIXELS =', self.SAVE_D_PIXELS
##        print 'self.SAVE_F_PIXELS =', self.SAVE_F_PIXELS

#         IDs = self.outlet_IDs
#         for k in xrange( len(self.out_files) ):
# 			#--------------------------------------
# 			# Open new files to write grid stacks
# 			#--------------------------------------
#             if (self.out_files[k].save_gs):
#                 model_output.open_new_gs_file( self, self.out_files[k], self.rti )
# 			#--------------------------------------
# 			# Open new files to write time series
# 			#--------------------------------------
#             if (self.out_files[k].save_ts):
#                 model_output.open_new_ts_file( self, self.out_files[k], IDs )
                                                          
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
    def manning_formula(self):

        #---------------------------------------------------------
        # Notes: R = (A/P) = hydraulic radius [m]
        #        N = Manning's roughness coefficient
        #            (usually in the range 0.012 to 0.035)
        #        S = bed slope or free slope

        #        R,S, and N may be 2D arrays.

        #        If length units are all *feet*, then an extra
        #        factor of 1.49 must be applied.  If units are
        #        meters, no such factor is needed.

        #        Note that Q = Ac * u, where Ac is cross-section
        #        area.  For a trapezoid, Ac does not equal w*d.
        #---------------------------------------------------------
        if (self.KINEMATIC_WAVE):
            S = self.S_bed
        else:
        	S = self.S_free

        u = (self.Rh ** self.two_thirds) * np.sqrt(S) / self.nval
        
        #--------------------------------------------------------
        # Add a hydraulic jump option for when u gets too big ?
        #--------------------------------------------------------
          
        return u
    
    #   manning_formula()
    #-------------------------------------------------------------------
    def law_of_the_wall(self):

        #---------------------------------------------------------
        # Notes: u  = flow velocity  [m/s]
        #        d  = flow depth [m]
        #        z0 = roughness length
        #        S  = bed slope or free slope

        #        g     = 9.81 = gravitation constant [m/s^2]
        #        kappa = 0.41 = von Karman's constant
        #        aval  = 0.48 = integration constant

        #        law_const  = sqrt(g)/kappa = 7.6393d
        #        smoothness = (aval / z0) * d
        #        f = (kappa / alog(smoothness))^2d
        #        tau_bed = rho_w * f * u^2 = rho_w * g * d * S

        #        d, S, and z0 can be arrays.

        #        To make default z0 correspond to default
        #        Manning's n, can use this approximation:
        #        z0 = a * (2.34 * sqrt(9.81) * n / kappa)^6d
        #        For n=0.03, this gives: z0 = 0.011417
        #########################################################
        #        However, for n=0.3, it gives: z0 = 11417.413
        #        which is 11.4 km!  So the approximation only
        #        holds within some range of values.
        #--------------------------------------------------------
        if (self.KINEMATIC_WAVE):
            S = self.S_bed
        else:
        	S = self.S_free

        smoothness = (self.aval / self.z0val) * self.d
          
        #------------------------------------------------
        # Make sure (smoothness > 1) before taking log.
        # Should issue a warning if this is used.
        #------------------------------------------------
        smoothness = np.maximum(smoothness, np.float64(1.1))

        u = self.law_const * np.sqrt(self.Rh * S) * np.log(smoothness)
        
        #--------------------------------------------------------
        # Add a hydraulic jump option for when u gets too big ?
        #--------------------------------------------------------
  
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
