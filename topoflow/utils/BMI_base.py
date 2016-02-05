#
#  We should use "initialize_scalar()" for all scalar assignments.
#  See the Notes for that method.  Search for "np.float64(0".
#      
#  Copyright (c) 2009-2014, Scott D. Peckham
#
#  Sep 2014. New initialize_basin_vars(), using outlets.py.
#            Removed obsolete functions.
#
#  Jan 2013. Added "initialize_scalar()" method.
#
#  Feb 2012. Complete BMI, starting from CSDMS_base.py
#            This now takes the place of CSDMS_base.py.
#
#  Nov 2011. Cleanup and conversion to BMI function names.
#
#  May 2010. initialize_config_vars(), cleanup, etc.
#
#  Aug 2009. Created, as CSDMS_base.py.
#
#-----------------------------------------------------------------------
#
#  Notes:  This file defines a "base class" with a BMI (Basic Model
#          Interface) for CSDMS "process" components.  These methods
#          allow a component to "fit into" a CMI harness (IMPL file)
#          that allows it to be used in the CSDMS/CCA framework.
#
#          The BMI interface is designed to be completely framework
#          independent, so this class contains no methods related to
#          CCA framework concepts such as ports.
#
#          Some "private" utility methods are defined at the end.
#
#          Need to figure out UDUNITS and use in get_var_units().
#
#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class BMI_component
#
#      __init__()
#
#      -------------------------------
#      BMI methods to get model info
#      -------------------------------
#      get_status()
#      get_attribute()
#      set_attribute()  (Experimental: not yet BMI)

#      ------------------------------
#      BMI methods to get grid info
#      ------------------------------
#      get_grid_shape()
#      get_grid_spacing()
#      get_grid_lower_left_corner()
#      get_grid_attribute()               ## NEW ADDITION TO BMI ??
#      read_grid_info()                   ## (Not part of BMI)
#
#      ----------------------------------
#      BMI methods to get variable info
#      ----------------------------------
#      get_input_var_names()
#      get_output_var_names()
#      -------------------------
#      get_var_name()                     # (override)
#      get_var_units()                    # (override)
#      get_var_rank()
#      get_var_type()
#      get_var_state()  or "mode()"       ############### "static" or "dynamic"  ###########
#      -------------------------
#      get_values()                       # (9/22/14)
#      set_values()                       # (9/22/14)
#      get_values_at_indices()
#      set_values_at_indices()
#
#      ------------------------------
#      BMI methods to get time info
#      ------------------------------
#      get_time_step()
#      get_time_units()
#      get_time()                        ## NEW ADDITION TO BMI ??
#      get_start_time()
#      get_current_time()
#      get_end_time()
#
#      --------------------------------------
#      BMI methods for fine-grained control
#      --------------------------------------
#      initialize           # (template)
#      update()             # (template)
#      finalize()
#      ------------------
#      run_model()          # (not required for BMI)
#      check_finished()     # (not part of BMI)
#
#      -----------------------------
#      More Time-related (not BMI)
#      -----------------------------
#      initialize_time_vars()
#      update_time()
#      print_time_and_value()
#      get_run_time_string()
#      print_run_time()
#
#      -------------------------------
#      Convenience methods (not BMI)
#      -------------------------------
#      print_final_report()          # (6/30/10)
#      print_traceback()             # (10/10/10)
#      -------------------------
#      read_config_file()            # (5/17/10, 5/9/11)
#      initialize_config_vars()      # (5/6/10)
#      set_computed_input_vars       # (5/6/10) over-ridden by each comp.
#      initialize_basin_vars()       # (9/19/14) New version that uses outlets.py.
#      initialize_basin_vars0()
#      -------------------------
#      prepend_directory()           # (may not work yet)
#      check_directories()
#      -------------------------
#      initialize_scalar()           # (2/5/13, for ref passing)
#      is_scalar()
#      is_vector()
#      is_grid()
#
#-----------------------------------------------------------------------

import numpy as np
import os
import sys
import time
import traceback        # (10/10/10)

#--------------------------------------------
# (5/14/10. Can't be here because basins.py
# has import BMI_base at top (circular).
# See initialize_basin_vars() below.
#--------------------------------------------
# import basins
   
## import cfg_files as cfg   # (not used)

import outlets          ## (9/19/14)
import pixels
import rti_files

#---------------------------------------------
# Experiment.  9/19/14
#---------------------------------------------
# from topoflow.utils import basins
# from topoflow.utils import cfg_files as cfg
# from topoflow.utils import pixels
# from topoflow.utils import rti_files

#-----------------------------------------------------------------------
def unit_test():

    c = BMI_component()
    print 'Instantiated BMI component.'

#   unit_test()
#-----------------------------------------------------------------------
class BMI_component:

    def __init__(self):

        ######################################################
        # These should be obsolete now.
        # self.CCA              = tf_utils.TF_Use_CCA()
        # self.USE_GUI_SETTINGS = False
        ######################################################
        
        ## self.DEBUG    = True
        self.DEBUG       = False
        self.SKIP_ERRORS = False
        self.SILENT      = True   # (new default: 11/16/11) 
        self.REPORT      = False
        self.DONE        = False
        self.status      = 'created'   # (OpenMI 2.0 conventions)

        self.in_directory     = None
        self.out_directory    = None
        self.site_prefix      = None
        self.case_prefix      = None
        self.cfg_prefix       = None       ###### (9/18/14)
        self.comp_status      = 'Enabled'
        
        # NB! This probably won't work here, because a driver
        #     may be instantiated later that then changes the
        #     current working directory.
##        self.cfg_directory  = os.getcwd()
##        self.cfg_prefix     = 'Case5'
        
    #   __init__()
    #-------------------------------------------------------------------
    def get_status(self):

        #-----------------------------------------------------
        # Notes: Return component status as a string.  The
        #        possible return values are from OpenMI 2.0:
        #
        #           created, initializing, initialized,
        #           updating, updated, finalizing, finalized,
        #           failed (could add "stopped").
        #-----------------------------------------------------
        return self.status

    #   get_status()
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        #----------------------------------------------------------
        # Note: This method must be overridden by each component.
        #----------------------------------------------------------
        att_map = {
            'model_name':         'Model_Name',
            'version':            '1.0',
            'author_name':        'Scott D. Peckham',
            'grid_type':          'uniform',
            'time_step_type':     'fixed',
            'step_method':        'explicit',
            #------------------------------------------------------
            'comp_name':          'Model_Name',
            'model_family':       'TopoFlow',
            'cfg_template_file':  'Model_Name.cfg.in',
            'cfg_extension':      '_Model_Name.cfg',
            'cmt_var_prefix':     '/ModelName/Input/Var/',
            'gui_xml_file':       '/home/csdms/cca/model_name/1.0/src/share/cmt/gui/Model_Name.xml',
            'dialog_title':       'ProcessName: Component Parameters',
            'time_units':         'seconds' }

        try:
            return att_map[ att_name.lower() ]
        except:
            print '######################################################'
            print ' ERROR: Could not find BMI attribute: ' + att_name
            print '######################################################'
            print ' '

    #   get_attribute()
    #-------------------------------------------------------------------
##    def set_attribute(self, att_name, att_string):
##
##        #-------------------------------------------------------
##        # Note: This is experimental, for setting in_directory
##        #       and out_directory, etc. (8/20/13)
##        #       Equivalent to: self.att_name = att_string.
##        #-------------------------------------------------------
##        setattr(self, att_name, att_string)
##
##    #   set_attribute()
    #-------------------------------------------------------------------

    def get_component_name(self):
        """Name of the component.

        Returns
        -------
        str
            The name of the component.
        """
        return self.get_attribute('model_name')

    def get_grid_rank(self, grid_id):
        """Get number of dimensions of the computational grid.

        Parameters
        ----------
        grid_id : int
          A grid identifier.

        Returns
        -------
        int
          Rank of the grid.
        """
        return len(self.get_grid_shape(grid_id))

    def get_grid_size(self, grid_id):
        """Get number of elements in the computational grid.

        Parameters
        ----------
        grid_id : int
          A grid identifier.

        Returns
        -------
        int
          Size of the grid.
        """
        return np.prod(self.get_grid_shape(grid_id))

    def get_grid_type(self, grid_id):
        """Get the grid type as a string.

        Parameters
        ----------
        grid_id : int
          A grid identifier.

        Returns
        -------
        str
          Type of grid as a string.
        """
        return self.get_attribute('grid_type')

    #-------------------------------------------------------------------
    def get_grid_shape(self, grid_id):

        #-------------------------------------------------------
        # Note: This assumes same grid info for all var_names.
        #-------------------------------------------------------
        # Note: We may want to reverse the order of the
        #       returned values in this method and the next 2.
        #       Ordering hasn't been established for BMI yet?
        #-------------------------------------------------------
        if not(hasattr( self, 'grid_info' )):
            self.read_grid_info()

        info  = self.grid_info
        shape = np.array([info.ny, info.nx])
        ## shape = np.array( [info.ncols, info.nrows, 0] )
        
        return shape
    
    #   get_grid_shape()
    #-------------------------------------------------------------------
    def get_grid_spacing(self, grid_id):

        #-------------------------------------------------------
        # Note: This assumes same grid info for all var_names.
        #--------------------------------------------------------
        # Note: xres and yres could have angular units like
        #       arcseconds, if (info.pixel_geom == 0).  In
        #       that case, xres (or dx) will vary with latitude
        #       and this method isn't really appropriate.
        #--------------------------------------------------------
        if not(hasattr( self, 'grid_info' )):
            self.read_grid_info()

        info = self.grid_info
        spacing = np.array([info.yres, info.xres])
        
        return spacing

    #   get_grid_spacing()
    #-------------------------------------------------------------------
    def get_grid_origin(self, grid_id):

        #-------------------------------------------------------
        # Note: This assumes same grid info for all var_names.
        #-------------------------------------------------------        
        if not(hasattr( self, 'grid_info' )):
            self.read_grid_info()

        info = self.grid_info
        corner = np.array([info.y_south_edge, info.x_west_edge])
        
        return corner

    #   get_grid_lower_left_corner()
    #-------------------------------------------------------------------
    def get_grid_attribute(self, long_var_name, att_name):

        amap = {'x_units':'meters', 'y_units':'meters',
                'z_units':'meters', 'ellipsoid':'None',
                'datum':'None', 'projection':'None'}
                ## 'utm_zone':'None'}

        try:
            return amap[ att_name.lower() ]
        except:
            print '######################################################'
            print ' ERROR: Could not find grid attribute: ' + att_name
            print '######################################################'
            print ' '            

    #   get_grid_attribute()
    #-------------------------------------------------------------------
    def read_grid_info(self):

        #-------------------------------------------------------
        # Note: This isn't a BMI method, but is currently used
        #       by all of the TopoFlow and Erode components.
        #       Notice that it uses "rti_files" and "pixels".
        #-------------------------------------------------------
        # Read grid info from an RTI file that is
        # in the current working directory.
        #------------------------------------------
        if (self.DEBUG):
            print 'Process component: Reading grid info...'
        
        self.grid_info_file = (self.in_directory +
                               self.site_prefix + '.rti')
        info = rti_files.read_info( self.grid_info_file )
        if (info == None):
            #-------------------------------------------------------
            # The RTI file is usually in "in_directory" and
            # uses "site_prefix", but when it needs to be
            # created, as for Erode (see erode_base.py), then
            # it uses "out_directory" and "case_prefix". (2/17/13)
            #-------------------------------------------------------
            print '### In BMI_base.read_grid_info():'
            print '### out_directory =', self.out_directory
            print ' '
            
            self.grid_info_file = (self.out_directory +
                                   self.case_prefix + '.rti')
            info = rti_files.read_info( self.grid_info_file )

        # print '##### In BMI_base.read_grid_info():'
        # print '##### in_directory   =', self.in_directory
        # print '##### grid_info_file =', self.grid_info_file
        
        #----------------------
        # Convenient synonyms
        #-----------------------------------------------------
        # Note that "info" has additional derived attributes
        # such as: n_pixels, bpe, grid_size and SWAP_ENDIAN.
        #-----------------------------------------------------
        self.rti = info
        self.nx  = info.ncols
        self.ny  = info.nrows

        #------------------------------------------------
        # (2/17/12) Start to phase out "rti" in favor
        #           of "grid_info".  See methods below.
        #------------------------------------------------
        self.grid_info    = self.rti
        self.grid_info.nx = self.rti.ncols
        self.grid_info.ny = self.rti.nrows
        
        #------------------------------------------------
        # Get grid cell areas, "da", which is either a
        # scalar (if same for all grid cells) or a grid
        # with default units of "m^2".
        #------------------------------------------------
        self.da = pixels.get_da( info )

    #   read_grid_info()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #-------------------------------------------------------------
        # Note: This method must be overridden by each model.
        #-------------------------------------------------------------
        # Note: There may be a way to retrieve this list auto-
        #       matically using code similar to self.has_variable().
        #-------------------------------------------------------------
        # Note:  NumPy dtype of returned array will be "|Sn",
        #        where n is the length of the longest string.
        #-------------------------------------------------------------
        items = ['None']
        return np.array( items )   # (string array vs. list)
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):

        #------------------------------------------------------
        # Note: This method must be overridden by each model.
        #------------------------------------------------------
        items = ['None']
        return np.array( items )   # (string array vs. list)
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):

        #------------------------------------------------------
        # Note: This method must be overridden by each model.
        #------------------------------------------------------
        # Define this map just once in "__init__()"  ??
        #-------------------------------------------------
        name_map = {
            'model__time_step':'dt' }

        return name_map[ long_var_name ]
    
##            #--------------------------
##            # Shared erode model vars
##            #--------------------------
##            'derivative_wrt_time_of_land_surface_elevation':'dz_dt',
##            'land_surface_d8_contributing_area':'A',
##            'land_surface_d8_slope':'S',
##            'land_surface_elevation':'DEM',
##            
##            #--------------------------
##            # Shared meteorology vars
##            #--------------------------
##            'air_density':'rho_air',
##            'air_relative_humidity':'RH',
##            'air_temperature':'T_air',
##            'air_thermal_capacity':'Cp_air',   # (CF uses thermal vs. heat, seems better)
##            'air_turbulent_boundary_layer_roughness_length':'z0',
##            'basin_cumulative_lwe_precipitated_water_volume':'vol_P',
##            'land_surface_air_pressure':'p0',   # (or just "surface_air_pressure" ??)
##            'land_surface_net_downward_energy_flux':'Q_sum',  # (not just SW+LW)
##            'land_surface_net_downward_longwave_flux':'Qn_LW',   # (clear-sky in metadata)
##            'land_surface_net_downward_shortwave_flux':'Qn_SW',
##            'land_surface_net_latent_heat_flux':'Qe',   ##########
##            'land_surface_temperature':'T_surf',
##            'lwe_max_precipitation_rate':'P_max',
##            'lwe_precipitation_rate':'P',
##            'water_density':'rho_H2O',
##            'wind_reference_height':'z',
##            'wind_reference_height_speed':'uz',  ########### ???
##            
##            #----------------------
##            # Shared channel vars
##            #----------------------
##            'basin_cumulative_runoff_water_volume':'vol_R',  ## (runoff vs. excess_rainrate ??)
##            'basin_outlet_water_discharge':'Q_outlet',       ####
##            'basin_outlet_water_mean_depth':'d_outlet',      ####           
##            'basin_outlet_water_mean_speed':'u_outlet',      ####
##            'channel_bed_max_roughness_length':'z0val_max',
##                ## 'max_of_channel_bed_roughness_length':'z0val_max',  #("Function" operator)
##            'channel_bed_max_manning_roughness_parameter':'nval_max',
##            'channel_bed_min_roughness_length':'z0val_min',
##            'channel_bed_min_manning_roughness_parameter':'nval_min',
##            'channel_x-section_hydraulic_radius':'Rh',
##            'channel_outgoing_sediment_discharge':'Qs',  #### check if "mass_flux", etc.
##            'channel_outgoing_water_discharge':'Q',
##            'channel_outgoing_peak_water_discharge':'Q_peak',
##            'channel_outgoing_peak_water_discharge_time':'T_peak',   ### (time as "quantity suffix")          
##                 ###'channel_incoming_water_discharge'   (if we need it)
##            'channel_reach_total_water_volume':'vol',
##                 ## channel_water_in_cell_volume
##            'channel_water_mean_depth':'d',   # (in_cell_mean_depth ??)
##                 ## 'mean_of_channel_water_depth':'d',  (mean => section_mean)
##            'channel_water_mean_speed':'u',
##                 ## 'mean_of_channel_water_speed':'d',
##            'channel_water_peak_mean_depth':'d_peak',
##            'channel_water_peak_mean_depth_time':'Td_peak',
##            'channel_water_peak_mean_speed':'u_peak',
##            'channel_water_peak_mean_speed_time':'Tu_peak',
##            
##            #---------------------------
##            # Shared infiltration vars
##            #---------------------------
##            'basin_cumulative_infiltrated_water_volume':'vol_IN',
##            'basin_cumulative_saturated_zone_infiltrated_water_volume':'vol_Rg',   ########
##            'land_surface_water_infiltration_rate':'IN',   ### (downward_flow_rate ? remove "land"?)
##            'subsurface_water_downward_flow_rate':'v',
##                ## 'ground_water_downward_flow_rate':'v',  ######
##            'water_table_recharge_rate':'Rg',  #####
##                ## 'saturated_zone_infiltration_rate':'Rg',
##            
##            #-----------------------------
##            # Shared saturated-zone vars
##            #-----------------------------
##            'basin_cumulative_subsurface_to_surface_seeped_water_volume':'vol_GW',
##            'ground_water_table_elevation':'h_table',    # (insert the word "surface" after table?)
##            'land_surface_elevation':'elev',  ######
##            'soil_layer_0_porosity':'qs[0]',
##            'soil_layer_0_thickness':'th[0,:,:]',
##            'soil_layer_0_saturated_thickness':'y[0,:,:]',            
##            'soil_layer_1_porosity':'qs[1]',
##            'soil_layer_1_thickness':'th[1,:,:]',
##            'soil_layer_1_saturated_thickness':'y[1,:,:]',            
##            'soil_layer_2_porosity':'qs[2]',
##            'soil_layer_2_thickness':'th[2,:,:]',
##            'soil_layer_2_saturated_thickness':'y[2,:,:]',
##            'soil_layer_porosity':'qs',
##            'soil_layer_thickness':'th',
##            'soil_layer_saturated_thickness':'y',  ######### (all soil layers;  how to specify an index?)
##            'subsurface_to_surface_water_seepage_rate':'GW',   ################
##                 ### land_surface_water_baseflow_seepage_rate
##                 ### subsurface_water_seepage_rate
##            'surface_soil_layer_porosity':'qs[0]',             # (same as soil_layer_0_porosity above)  ####
##            'surface_soil_layer_saturated_thickness':'y[0,:,:]',  # (same as soil_layer_0_wetted_thickness above)  ####
##            
##            #--------------------------
##            # Shared evaporation vars
##            #--------------------------
##            'basin_cumulative_evaporated_water_volume':'vol_ET',    ## (use "total" vs. "cumulative" ??)
##            'land_surface_water_evaporation_rate':'ET',
##            'land_surface_water_potential_evaporation_rate':'PET',  #### (not used)
##            
##            #-------------------
##            # Shared snow vars
##            #-------------------                   
##            'basin_cumulative_snow_meltwater_volume':'vol_SM',  # (meltwater better than melted_water)
##            'surface_snow_density':'rho_snow',
##            'surface_snow_depth':'h_snow',
##            'surface_snow_melt_rate':'SM',
##            'surface_swe_depth':'h_swe',
##            
##            #------------------
##            # Shared ice vars
##            #------------------
##            'basin_cumulative_ice_meltwater_volume':'vol_MR',
##            'surface_ice_density':'rho_ice', ########
##            'surface_ice_depth':'H',         ########  (change to "h_ice" ???)
##            'surface_ice_melt_rate':'MR',    ########  (land_ice_basal_melt_rate ??)          
##            #--------------------------------
##            'time_step_size':'dt'}
##
##        return name_map[ long_var_name ]
        
    #   get_var_name()
    #-------------------------------------------------------------------

    def get_var_grid(self, long_var_name):
        """Get grid identifier for a variable name.

        Parameters
        ----------
        long_var_name : str
            An input or output variable name as a CSDMS standard name.

        Returns
        -------
        int
            The grid identifier.
        """
        return 0

    def get_var_itemsize(self, long_var_name):
        """The memory use of each array element in bytes.

        Parameters
        ----------
        long_var_name : str
            An input or output variable name as a CSDMS standard name.

        Returns
        -------
        int
            Item size in bytes.
        """
        var_name = self.get_var_name(long_var_name)
        return getattr(self, var_name).itemsize

    def get_var_units(self, long_var_name):
        #-------------------------------------------------
        # Define this map just once in "__init__()"  ??
        #-------------------------------------------------
        units_map = { 'Q':     'm3/s',
                      'Qs':    'm3/s', ########
                      'A':     'km2',
                      'S':     'none',  ## (or 'm/m')
                      'z':     'm',
                      'dz_dt': 'm/yr',
                      'DEM':   'm',
                      #--------------------------
                      'T_air':     'degrees_C',
                      'T_surf':    'degrees_C',
                      'rho_H2O':   'kg/m3',
                      'rho_snow':  'kg/m3',
                      'Q_sum':     'W/m2',
                      'h_swe':     'm',
                      'h_snow':    'm' }

        var_name = self.get_var_name( long_var_name )
        
        return units_map[ var_name ]
    
    #   get_var_units()
    #-------------------------------------------------------------------
    def get_var_rank(self, long_var_name):

        var_name = self.get_var_name( long_var_name )  # (2/20/12)
        
        exec("rank = np.rank(self." + var_name + ")")

        ### print '######## rank(' + var_name + ') =', rank
        
        return rank
    
        ## return np.int32( rank )  ###### (need this ??)

    #   get_var_rank()
    #-------------------------------------------------------------------
    def get_var_type(self, long_var_name):

        #--------------------------------------------------------
        # Notes: I think we should use the same type names that
        #        NumPy "ndarrays" store as "dtype" because they
        #        are unambiguous.  We assume here that all vars
        #        are NumPy types with dtype defined. e.g.
        #        uint8, int16, int32, float32, float64, etc.
        #--------------------------------------------------------
        var_name = self.get_var_name( long_var_name )  # (2/20/12)

        try:
            exec( "dtype = self." + var_name + ".dtype" )
        except:
            dtype = 'unknown'
        return str(dtype)       # (need str() here)

        #----------------------------------------------------------------
##        # This should also work.
##        exec( "HAS_DTYPE = hasattr( self." + var_name + ", 'dtype')" )
##        if (HAS_DTYPE):
##            exec( "dtype = self." + var_name + ".dtype" )
##        else:
##            dtype = 'unknown'   ###############
##        return dtype
    
    #   get_var_type() 
    #-------------------------------------------------------------------     

    def get_values(self, long_var_name):

        #------------------------------------------------------- 
        # Note: The value returned by getattr() will have rank
        #       and data type that goes with long_var_name.
        #------------------------------------------------------- 
        var_name = self.get_var_name( long_var_name )

        try:
            return getattr(self, var_name)   ## (2/19/13)

            #--------------------------------------------------            
            # Return 0 as default if attribute doesn't exist.
            #--------------------------------------------------
            # return getattr(self, var_name, 0)
            
            #-----------------------------------
            # Using exec works, but is slower.
            #-----------------------------------
            ## exec("result = self." + var_name)
            ## return result

            #------------------------------------
            # This doesn't work as a one-liner.
            #------------------------------------
            ## exec("return self." + var_name)

            #----------------------------
            # This breaks the reference.
            #----------------------------
            ## return np.float64(result)

        except:
            print 'ERROR in BMI_base.get_values()'
            print '    for var_name =', var_name
            print '    Returning 0.'
 
            #-----------------------------------------           
            # We could also call get_var_rank() and
            # return an array of the right rank, but
            # a 0D array seems better.
            #-----------------------------------------
            ## value = np.float64(0)   # (scalar)
            dtype = self.get_var_type( long_var_name )
            value = np.array(0, dtype=dtype)  # (0D array)
    
            #---------------------------------------------        
            # Save new variable with this name into self
            # to avoid repeating this message later.
            #---------------------------------------------
            setattr(self, var_name, value)
            return value
        
    #   get_values()
    #-------------------------------------------------------------------
    def set_values(self, long_var_name, value):

        #---------------------------------------------------------------
        # Notes: The "var_name" string cannot contain a ".". (5/17/12)
        #---------------------------------------------------------------
        # (2/7/13) We are now using 0D numpy arrays as a way to
        # produce "mutable scalars" that allow a component with a
        # reference to the scalar to see changes to its value.
        # But we can't apply np.float64() to the value as we did
        # before or it destroys the reference.
        # See BMI_base.initialize_scalar() for more information.
        #--------------------------------------------------------------- 
        var_name = self.get_var_name( long_var_name )
        setattr( self, var_name, value )  ## (2/19/13)
         
    #   set_values()
    #-------------------------------------------------------------------
    def get_values_at_indices(self, long_var_name, IDs):

        #---------------------------------------------------------
        # Note: This function was causing a "segmentation fault
        #       in gui-backend.sh" error message when trying to
        #       run TopoFlow through the CMT (in CCA framework).
        #       Solution was to use np.array, as shown.
        #       (2/18/10)
        #---------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #
        #        Use of getattr() has not been tested. (2/19/13)
        #---------------------------------------------------------
        var_name = self.get_var_name( long_var_name )
        
        try:
            var_IDs_name = var_name + '.flat[IDs]'
            result = getattr(self, var_IDs_name)  ## (2/19/13)
            return np.array(result, dtype='float64')
        except:
            print 'ERROR in BMI_base.get_values_at_indices().'
            print '    Returning zeros.'
            dtype = self.get_var_type( long_var_name )
            return np.zeros(len(IDs), dtype=dtype)
        
    #   get_values_at_indices()
    #-------------------------------------------------------------------
    def set_values_at_indices(self, long_var_name, IDs, values):

        #--------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #--------------------------------------------------------
        var_name = self.get_var_name( long_var_name )
        var_IDs_name = var_name + '.flat[IDs]'
        setattr( self, var_IDs_name, values )  ## (2/19/13)
        
    #   set_values_at_indices()          
    #-------------------------------------------------------------------

    def get_value(self, var_name):
        """Get a copy of the values of the given variable.

        This is a getter for the model, used to access the model's
        current state. It returns a *copy* of a model variable, with
        the return type, size and rank dependent on the variable.

        Parameters
        ----------
        var_name : str
          An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        array_like
          The value of a model variable.
        """
        return self.get_values(var_name)

    def set_value(self, var_name, src):
        """Specify a new value for a model variable.

        This is the setter for the model, used to change the model's
        current state. It accepts, through *src*, a new value for a
        model variable, with the type, size and rank of *src*
        dependent on the variable.

        Parameters
        ----------
        var_name : str
          An input or output variable name, a CSDMS Standard Name.
        src : array_like
          The new value for the specified variable.
        """
        if src.shape == (1,):
            src = np.array(src[0])
        self.set_values(var_name, src)

    def get_value_at_indices(self, var_name, indices):
        """Get values at particular indices.

        Parameters
        ----------
        var_name : str
          An input or output variable name, a CSDMS Standard Name.
        indices : array_like
          The indices into the variable array.

        Returns
        -------
        array_like
            Value of the model variable at the given location.
        """
        return self.get_values_at_indices(var_name, indices)

    def set_value_at_indices(self, var_name, indices, src):
        """Specify a new value for a model variable at particular indices.

        Parameters
        ----------
        var_name : str
          An input or output variable name, a CSDMS Standard Name.
        indices : array_like
          The indices into the variable array.
        src : array_like
          The new value for the specified variable.
        """
        self.set_values_at_indices(var_name, indices, src)

    #-------------------------------------------------------------------
    # BMI methods to get time-related info
    #-------------------------------------------------------------------
    def get_time_step(self):

        return np.float64( self.dt )
    
    #   get_time_step()
    #-------------------------------------------------------------------        
    def get_time_units(self):

        #-----------------------------------------------------------
        # Notes: The get_attribute() method should always provide
        #        a model's time units.  However, some components
        #        have an initialize_time_vars() method that saves
        #        a units string in self.time_units.  These strings
        #        are lower-case, e.g. 'seconds', 'minutes',
        #        or 'years'.  They should conform to UDUNITS.
        #-----------------------------------------------------------
        unit_str = self.get_attribute( 'time_units' )
        return unit_str
        ## return self.time_units  
    
    #   get_time_units()
    #-------------------------------------------------------------------
    def get_time(self, when='current'):

        #----------------------------------------------------------
        # Note:  The "when" argument can be: current, start, end.
        #----------------------------------------------------------
        # Notes: Most TF and Erode components do not specify
        #        a particular start time.  An exception is the
        #        TF Meteorology component, which we could check
        #        for here and process differently.
        #----------------------------------------------------------
        if (when == 'current'):
            return self.time
        elif (when == 'start'):
            return np.float64( 0 )
        else:
            pass  # (continue below)

        #--------------------------------------------------
        # Does model have a "stop_time" attribute ?
        # Even if it does, it may not be honored exactly.
        #--------------------------------------------------
        if (hasattr( self, 'stop_time' )):
            return np.float64( self.stop_time )

        #---------------------------------
        # Can we compute the stop_time ?
        #---------------------------------
        dt_type    = self.get_attribute( 'time_step_type' )
        COMPUTABLE = (dt_type == 'fixed') and (hasattr(self, 'n_steps'))
        
        if (hasattr( self, 'stop_code' )):
            if (self.stop_code != 0):
                COMPUTABLE = False  # (overrides COMPUTABLE above)

        if (COMPUTABLE):
            return np.float64( self.n_steps * self.dt )
        else:  
            print '##############################################'
            print ' ERROR: Unable to compute model stop_time.'
            print '##############################################'
            print ' '
            return np.float64( -1 )
        
    #   get_time()
    #-------------------------------------------------------------------
    def get_start_time(self):

        #--------------------------------------------------------
        # Notes: Most TF and Erode components do not specify
        #        a particular start time.  An exception is the
        #        TF Meteorology component, which we could check
        #        for here and process differently.
        #--------------------------------------------------------
        return np.float64( 0  )
    
    #   get_start_time()
    #-------------------------------------------------------------------
    def get_current_time(self):

        return self.time
    
    #   get_current_time()
    #-------------------------------------------------------------------
    def get_end_time(self):

        #--------------------------------------------------
        # Does model have a "stop_time" attribute ?
        # Even if it does, it may not be honored exactly.
        #--------------------------------------------------
        if (hasattr( self, 'stop_time' )):
            return np.float64( self.stop_time )

        #---------------------------------
        # Can we compute the stop_time ?
        #---------------------------------
        dt_type    = self.get_attribute( 'time_step_type' )
        COMPUTABLE = (dt_type == 'fixed') and (hasattr(self, 'n_steps'))
        
        if (hasattr( self, 'stop_code' )):
            if (self.stop_code != 0):
                COMPUTABLE = False  # (overrides COMPUTABLE above)

        if (COMPUTABLE):
            return np.float64( self.n_steps * self.dt )
        else:  
            print '##############################################'
            print ' ERROR: Unable to compute model stop_time.'
            print '##############################################'
            print ' '
            return np.float64( -1 )           
            
    #   get_end_time()
    #-------------------------------------------------------------------
    # BMI methods for fine-grained control of model
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver"):

        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        # Will use cfg_file from above.
        #-----------------------------------------------
        ## self.set_constants()
        self.initialize_config_vars()  # calls check_directories().
        self.read_grid_info()
        self.initialize_basin_vars()
        self.initialize_time_vars()

        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

##        self.open_output_files()

        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
##    def update(self, time_seconds=None):
##
##        self.status = 'updating'  # (OpenMI 2.0 convention)
##
##        #-----------------------------------
##        # Update computed values here with
##        # a series of "update_var()" calls
##        #----------------------------------
##        # self.update_var1()
##        # self.update_var2()
##
##        #------------------------
##        # Update internal clock
##        #------------------------
##        self.update_time()
##
##        #-------------------------------
##        # Check for NaNs, etc. in var1
##        #-------------------------------    
##        # self.check_var1()
##
##        #------------------------------------------
##        # Read next infil vars from input files ?
##        #------------------------------------------
##        self.read_input_files()
##
##        #----------------------------------------------
##        # Write user-specified data to output files ?
##        #----------------------------------------------
##        self.write_output_files(time_seconds)
##        self.status = 'updated'  # (OpenMI 2.0 convention)
##        
##    #   update()
    #-------------------------------------------------------------------

    def update(self):
        """Advance model state by one time step."""
        pass

    def update_frac(self, time_frac):
        """Update model by a fraction of a time step.

        Parameters
        ----------
        time_frac : float
            Fraction of a time step.
        """
        time_step = self.get_time_step()
        self.dt = time_frac * time_step
        if self.dt > 0.0:
            self.update()
        self.dt = time_step

    def update_until(self, then):
        """Advance model state until the given time.

        Parameters
        ----------
        then : float
            Time to run model until.
        """
        n_steps = (then - self.get_current_time()) / self.get_time_step()

        for _ in xrange(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))

    def finalize(self):

        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        self.close_input_files()    ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)
        
    #   finalize()
    #-------------------------------------------------------------------
    def run_model(self, cfg_directory=None, cfg_prefix=None,
                  n_steps=5):

        #--------------------------------------------------
        # NOTE: This method is not called from the go()
        #       method in the IMPL file, but it is often
        #       called by a unit_test().
        #--------------------------------------------------

        #--------------------------------------------------
        # All components, including this one (the driver)
        # will look in the CWD for their CFG file.
        #--------------------------------------------------
        if (cfg_directory != None):
            os.chdir( cfg_directory )
        self.cfg_prefix = cfg_prefix

        #-----------------------------------------
        # Initialize the model run (driver mode)
        #---------------------------------------------
        # This will set in_directory, out_directory,
        # site_prefix and case_prefix
        #---------------------------------------------
        self.initialize( cfg_prefix=cfg_prefix, mode='driver' )
            
        #----------------------------------------------------------- 
        # Note:  If the number of timesteps is specified in a
        #        component's CFG file and is then saved by
        #        "read_cfg_file()" as "n_steps" then we should
        #        honor that setting.  Otherwise we use the n_steps
        #        argument.
        #-----------------------------------------------------------
        if not(hasattr( self, 'stop_code' )):
           self.stop_code = 0
        #--------------------------------------
        if (hasattr(self, 'n_steps')):
            ## print 'NUMBER OF STEPS =', self.n_steps  ####
            n_steps = self.n_steps
        # else:
            # (Use the n_steps argument.)

        #-------------------------------------------        
        # Note: __init__() sets self.DONE to False
        #-------------------------------------------
        while not(self.DONE):
            if not(self.SKIP_ERRORS):
                #-------------------------------------------
                # Exceptions will not be caught and
                # Python error messages will be displayed.
                #-------------------------------------------
                if (self.DEBUG):
                    print 'time_index =', self.time_index
                self.update()
                # self.update( -1 )  # (BMI later: use own dt)
            else:   
                try:
                    self.update()
                    # self.update( -1 )  # (BMI later: use own dt)
                except:
                    print 'ERROR in run_model() method at:'
                    print '   time_index =', self.time_index
                    self.status = 'failed'
                    self.DONE = True

            #------------------------------------------------
            # If the model has set self.DONE = True, then
            # stop, even if we haven't reached n_steps yet.
            #------------------------------------------------
            self.check_finished()   # (see below: 2/13/12)

        #-------------------------
        # Finalize the model run
        #-------------------------
        self.finalize()

    #   run_model()
    #-------------------------------------------------------------------
    def check_finished(self):

        #---------------------------------------------------------
        # Note: If self.DONE has already been set to True by
        #       another function or component, this function
        #       preserves that setting (see below).
        #---------------------------------------------------------
        #       TINY_DZ can occur either because dt required for
        #       stability is really small or because we have
        #       converged to a steady-state landscape.
        #---------------------------------------------------------
        #       Moved here from erode_base.py on 2/13/12.
        #---------------------------------------------------------

        #---------------------------------------------------
        # If component isn't the driver, return. (2/20/12)
        #---------------------------------------------------
        ### if (self.mode == 'nondriver'):
        if (self.mode != 'driver'):
            return
        
        if (self.stop_code == 0):
            #---------------------
            # Stop after n_steps
            #---------------------
            TIMES_UP = (self.time_index >= self.n_steps)
        elif (self.stop_code == 1):
            #-----------------------
            # Stop after stop_time 
            #-----------------------
            TIMES_UP = (self.time >= self.stop_time)
        elif (self.stop_code == 2):
            #-----------------------------------------
            # Stop if "steady-state", but only works
            # as written here for global timesteps. 
            #-----------------------------------------
            TIMES_UP = (self.dz_max < self.dz_tolerance)

        self.DONE = (self.DONE or TIMES_UP)

    #   check_finished()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def initialize_time_vars(self, units='seconds'):

        #------------------
        # Start the clock
        #------------------
        self.start_time = time.time()
        
        #--------------------------------
        # Initialize the time variables
        #--------------------------------
        self.time_units = units.lower()
        self.time_index = np.int32(0)
        self.time       = self.initialize_scalar(0, dtype='float64') 
        self.DONE       = False

        #-------------------------------------------
        # (2/5/12) Set default stopping method:
        # (0=n_steps, 1=stop_time, 2=steady-state)
        #-------------------------------------------
        # Don't do this in set_constants().
        #-------------------------------------------
        if not(hasattr(self, 'stop_code')):
            self.stop_code = 0

        #---------------------------------------
        # Make sure n_steps is defined for new
        # use in update_time().  Was not set
        # for diversions_base. (11/14/11)
        #---------------------------------------
        if not(hasattr(self, 'n_steps')):
            self.n_steps = 1
 
        #--------------------------
        # Time conversion factors
        #--------------------------
        self.sec_per_year = np.float64(365) * 24 * 3600
        self.min_per_year = np.float64(365) * 24 * 60
        
        #-------------------------------------------
        # For backward compatibility with TopoFlow
        #------------------------------------------------
        # Use 1D array (mutable) vs. scalar (immutable)
        # to permit passing as a reference.  (2/4/13)
        #------------------------------------------------
        self.time_sec = self.initialize_scalar(0, dtype='float64')
        self.time_min = self.initialize_scalar(0, dtype='float64')
            
        #--------------------------------------------
        # For print_time_and_value() function below
        #--------------------------------------------
        # Substract 100 seconds so we'll always
        # print values at time zero. (6/29/10)
        #--------------------------------------------
        self.last_print_time = time.time() - 100.0
        
##        self.last_check_time  = time.time()  # (for user interrupt)
##        self.last_plot_time   = np.float32(0.0)   ### CHECK ###
        
    #   initialize_time_vars()
    #-------------------------------------------------------------------
    def update_time(self, dt=-1):

        #--------------------------------------------------
        # Note: The BMI update() method has a dt argument
        #       that is passed to this routine.
        #--------------------------------------------------
        
        #---------------------
        # Increment the time
        #---------------------
        self.time_index += 1
        if (dt == -1):
            self.time += self.dt  # (use same units as dt)
        else:
            self.time += dt

        #------------------------------------------------
        # Compute and store the current time in seconds
        # (time_sec) and in minutes (time_min).
        #------------------------------------------------
        if (self.time_units == 'seconds'):
            self.time_sec = self.time         # [seconds]
        elif (self.time_units == 'minutes'):
            self.time_sec = self.time * np.float64(60)
            ## self.time_min = self.time
        elif (self.time_units == 'hours'):
            self.time_sec = self.time * np.float64(3600)
        elif (self.time_units == 'days'):
            self.time_sec = self.time * np.float64(3600) * 24
        elif (self.time_units == 'years'):
            #-----------------------------------
            # Used by GC2D and Erode (12/4/09)
            #-----------------------------------
            self.time_sec = self.time * self.sec_per_year  ####

        #------------------------------------------
        # Compute & store current time in minutes
        #------------------------------------------
        self.time_min = self.time_sec / np.float64(60)   # [minutes]

        #---------------------------------------------------
        # Are we DONE yet?  BMI_base.py contains a basic
        # version of check_finished(), but components that
        # inherit from BMI_base.py may override it.
        #---------------------------------------------------
        # Moved here for new framework approach. (5/18/12)
        #---------------------------------------------------
        self.check_finished()
                
    #   update_time()
    #-------------------------------------------------------------------
    def print_time_and_value(self, var, var_name='Q_out',
                             units_name='[m^3/s]',
                             interval=2.0,
                             PRINT_INDEX=False):

        #-----------------------------------------
        # (8/2/10) Print message about interval.
        #-----------------------------------------
        if (self.time_index == 0):
            print 'Will print values every', interval, 'seconds.'
            
        #---------------------------------------------------
        # Note: Print the model time, in minutes, and the
        #       current value of "var", at the specified
        #       real-time "interval" (in seconds).
        #---------------------------------------------------
        # Note: Plotting hydrograph at same interval is
        #       generally too infrequent.
        #---------------------------------------------------
        elapsed_time = (time.time() - self.last_print_time)
        if (elapsed_time > interval):
            if (self.time_units == 'seconds'):
                cur_time = self.time_min
                time_units_str = ' [min]'
            else:
                cur_time = self.time
                time_units_str = ' [' + self.time_units + ']' 
            time_str = 'Time = ' + ("%10.2f" % cur_time)
            time_str = time_str + time_units_str
            #-------------------------------------------------
            var_str  = var_name + ' = ' + ("%10.5f" % var)
            var_str  = var_str  + ' ' + units_name          
            #-------------------------------------------------      
            print (time_str + ',  ' + var_str)
            #-----------------------------------------------------
            if (PRINT_INDEX):
                index = (self.time_index + 1)  # (starts at 0)
                print 'n =', index, 'of', self.n_steps
            #-----------------------------------------------------                
            self.last_print_time = time.time()

    #   print_time_and_value()
    #-------------------------------------------------------------------
    def get_run_time_string(self, proc_name='component',
                            sec_digits=4, seconds=None,
                            SILENT=None):
                            ### SUB_PROCESS=False) 

        #------------------------------------------------------
        # If "seconds" argument is only provided for testing.
        # You can provide this value to make sure that the
        # minuts, hours, days, etc. are computed correctly.
        #------------------------------------------------------
        if (seconds == None):    
            finish  = time.time()
            seconds = (finish - self.start_time)

        #-------------------------------------------------
        # Could later add this to help gauge performance
        # for all models.  It is currently coded into
        # erode_d8_local.finalize(). But not all models
        # have self.time in "years". (10/4/11)
        #-------------------------------------------------
##        finish         = time.time()
##        run_time_secs  = (finish - self.start_time)
##        run_time_hours = (run_time_secs / 3600.0)
##        sim_time_years = self.time
##        years_per_hour = (sim_time_years / run_time_hours)
##        print ' '
##        print 'Years per hour =', years_per_hour
        
        #----------------------------------
        # Compute minutes, hours and days
        #----------------------------------
        dec_part  = (seconds % np.float32(1.0))     #(Save decimal part)
        days      = np.int32(seconds) / np.int32(86400)
        secs_left = np.int32(seconds) % np.int32(86400)
        hours     = (secs_left / np.int32(3600))
        secs_left = (secs_left % np.int32(3600))
        minutes   = (secs_left / np.int32(60))
        seconds   = (secs_left % np.int32(60))
        #-----------------------------------------
        #hours     = long(seconds)  /  3600L
        #secs_left = long(seconds) mod 3600L
        #minutes   = (secs_left  /  60L)
        #seconds   = (secs_left mod 60L)
        
        #----------------------------
        # Construct the time string
        #----------------------------
        time_string = ''
        #--------------------------------------------------------
        if (days > 0):    
            if (days > 1):    
                e0 = ' days, '
            else:    
                e0 = ' day, '
            time_string += str(days) + e0
        #--------------------------------------------------------
        if (hours > 0):    
            if (hours > 1):    
                e1 = ' hours, '
            else:    
                e1 = ' hour, '
            time_string += str(hours) + e1
        #--------------------------------------------------------
        if (minutes > 0):    
            if (minutes > 1):    
                e2 = ' minutes, '
            else:    
                e2 = ' minute, '
            time_string += str(minutes) + e2
        
        #-----------------------------------------
        # Default is 4 digits after the decimal.
        #-----------------------------------------
        dec_pastr = ('.' + str(dec_part)[2:2+sec_digits])
        time_string += str(seconds) + dec_pastr + ' seconds.'

        return time_string
            
##            if (SUB_PROCESS):    
##                PART1 = '>> '
##            else:    
##                PART1 = ''
##            print (PART1 + 'Run time for ' + procname + ' = ')
##            print (PART1 + time_string)
##            print ' '
     
    #   get_run_time_string()
    #-------------------------------------------------------------------
    def print_run_time(self, proc_name='component',
                       sec_digits=4, seconds=None,
                       SILENT=None):
                       ### SUB_PROCESS=False) 


        print ('Run time for ' + proc_name + ' = ')
        print  self.get_run_time_string()
        print ' '
     
    #   print_run_time()
    #-------------------------------------------------------------------    
    def print_final_report(self, comp_name='BMI component',
                           mode='nondriver'):

        if (mode == 'nondriver'):
            print comp_name + ': Finished.'
            return

        if not(hasattr( self, 'in_directory' )):
            return   # (2/20/12)
        
        #-------------------
        # Print the report
        #-------------------
        hline = ''.ljust(60, '-')
        print hline
        print comp_name
        print time.asctime()
        print ' '
        print 'Input directory:      ' + self.in_directory
        print 'Output directory:     ' + self.out_directory
        print 'Site prefix:          ' + self.site_prefix
        print 'Case prefix:          ' + self.case_prefix
        print ' '

        #-----------------------------------
        # Construct sinulation time string
        #-----------------------------------
        sim_units    = ' [' + self.time_units + ']'
        sim_time_str = str(self.time) + sim_units
        
        #----------------------------
        # Construct run time string
        #----------------------------
        run_time_str = self.get_run_time_string() 

        print 'Simulated time:      ' + sim_time_str
        print 'Program run time:    ' + run_time_str
        print ' '
        print 'Number of timesteps: ' + str(self.time_index)
        print 'Process timestep:    ' + str(self.dt) + sim_units
        print 'Number of columns:   ' + str(self.nx)
        print 'Number of rows:      ' + str(self.ny)
        print ' '
        print 'Finished. (' + self.case_prefix + ')'
        print ' '

        ## finish_str = ': Finished. (' + self.case_prefix + ')'
        ## print finish_str
        ## print comp_name + finish_str
        
    #   print_final_report()
    #-------------------------------------------------------------------    
    def print_traceback(self, caller_name='TopoFlow'):

        print '################################################'
        print ' ERROR encountered in ' + caller_name
        print '       Please check your input parameters.'
        print '################################################'
        
        traceback.print_exc()
        
    #   print_traceback()
    #-------------------------------------------------------------------
    def read_config_file(self):

        #----------------------------------------------------------
        # Notes: This version reads configuration settings from
        #        a new type of CFG file that only has var_name,
        #        value and type; more like key-value. (5/9/11)
        #        The GUI is no longer generated using a CFG file.
        #----------------------------------------------------------
        #        All of the TopoFlow and Erode components now use
        #        the same type of CFG file and all of them use
        #        this method to read settings from the file.
        #----------------------------------------------------------
        if not(self.SILENT):
            print 'Reading config file into component state.'

        # print '######## In read_config_file(): cfg_file =', self.cfg_file
        # print '######## In read_config_file(): cfg_prefix =', self.cfg_prefix

        #############################################################
        # No longer needed; cfg_file is passed to BMI.initialize()
        # and saved as self.cfg_file. (9/17/14)
        #############################################################
        #---------------------------
        # Get name of the cfg_file
        #---------------------------
#         cfg_extension = self.get_attribute( 'cfg_extension' )  # (10/26/11)
#         cfg_directory = (os.getcwd() + os.sep)
#         file_name     = (self.cfg_prefix + cfg_extension)
#         self.cfg_file = (cfg_directory + file_name) 

        #------------------------
        # Does CFG file exist ?
        #------------------------
        if (self.DEBUG):
            print 'cfg_file =', self.cfg_file
        if not(os.path.exists(self.cfg_file)):
            print 'WARNING: cfg_file not found:'
            print '         ' + self.cfg_file
            return

        #-----------------------------
        # Open CFG file to read data
        #-----------------------------
        cfg_unit = open( self.cfg_file, 'r' )
        last_var_name = ''

        #-----------------------------------------
        # Save user input into component's state
        #--------------------------------------------------
        # Recall that a "blank line", with just a (hidden)
        # newline character will not be null and will
        # have len(line) = 1.
        #--------------------------------------------------
        while (True):
            line  = cfg_unit.readline()
            if (line == ''):
                break                  # (reached end of file)
            
            COMMENT = (line[0] == '#')
            #--------------------------------------------
            # Using "|" as a delimiter means we can use
            # " ", "," or "=" in the filenames.
            #--------------------------------------------
            # Added ".strip()" on 10/25/11.
            #--------------------------------------------
            words   = line.split('|')  # (split on equals)
            if (len(words) == 4) and not(COMMENT):
                var_name = words[0].strip()
                value    = words[1].strip()
                var_type = words[2].strip()
                ## help_str = words[3]
                READ_SCALAR   = False
                READ_FILENAME = False

                # For debugging
                # print 'var_name, value, var_type =', var_name, value, var_type

                #----------------------------------------------
                # Does var_name end with an array subscript ?
                #----------------------------------------------
                p1 = var_name.rfind('[')
                p2 = var_name.rfind(']')
                if (p1 > 0) and (p2 > p1):
                    var_base  = var_name[:p1]
                    subscript = var_name[p1:p2+1]
                    var_name_file_str = var_base + '_file' + subscript
                else:
                    var_base = var_name
                    var_name_file_str = var_name + '_file'

                #--------------------------------------------
                # Update var_type based on droplist setting
                #--------------------------------------------
                if (last_var_name.startswith(var_base + '_type')):
                    exec( "type_choice = self." + last_var_name )
                    if (type_choice.lower() == 'scalar'):
                        #--------------------------------------------------
                        # It seems that things will work as long as the
                        # "type" and "value" fields in the GUI CFG file
                        # are consistent.  Don't change var_type here.
                        #
                        # Otherwise get this message:
                        # "Mismatch with value found in typemap
                        #  (requested type String, actual type Double)."
                        #--------------------------------------------------
                        exec( "self." + var_name_file_str + " = ''")
                        READ_SCALAR = True
                        ## var_type = 'float64'
                    else:
                        exec( "self." + var_name + " = 0.0")
                        READ_FILENAME = True
                        ## var_type = 'string'

                #-----------------------------------           
                # Read a value of type "var_type"
                #-----------------------------------
                # Convert scalars to numpy scalars
                #-----------------------------------
                if (var_type in ['float64', 'float']):
                    value = np.float64( value )

                    #------------------------
                    # For testing (5/18/12)
                    #------------------------
##                    print 'var_name =', var_name
##                    print 'var_type =', var_type
##                    print 'value    =', value
##                    print '---------------------------------'
                    
                    exec( "self." + var_name + " = value" )
                elif (var_type in ['long', 'int']):
                    value = np.int32( value )
                    exec( "self." + var_name + " = value" )
                elif (var_type == 'string'):
                    #-----------------------------------------
                    # Get the value string for this var_name
                    #----------------------------------------------------
                    # If string has a placeholder filename prefix, then
                    # expand it here.  Need to use original "var_name"
                    # without appending "_file" until assignment.
                    #----------------------------------------------------
                    # case_str = '<case_prefix>'
                    # site_str = '<site_prefix>'
                    case_str = '[case_prefix]'
                    site_str = '[site_prefix]'
                    #---------------------------------
                    s = value
                    if (s[:13] == case_str):
                        value_str = (self.case_prefix + s[13:])
                    elif (s[:13] == site_str):
                        value_str = (self.site_prefix  + s[13:])
                    else:
                        value_str = s

                    #-----------------------------------------------
                    # If var_name starts with "SAVE_" and value is
                    # Yes or No, then convert to Python boolean.
                    #-----------------------------------------------
                    if (var_name[:5] == 'SAVE_'):
                        VALUE_SET = True
                        if (s.lower() == 'yes'):
                            exec( "self." + var_name + " = True" )
                        elif (s.lower() == 'no'):
                            exec( "self." + var_name + " = False" )
                        else:
                            VALUE_SET = False
                    else:
                        VALUE_SET = False
                    #----------------------------------------------------------
                    if not(VALUE_SET):
                        if (READ_FILENAME):
                            exec( "self." + var_name_file_str + " = value_str" )
                        elif (READ_SCALAR):
                            exec( "self." + var_name + " = np.float64(value_str)")
                        else:
                            exec( "self." + var_name + " = value_str" )
                else:
                    print 'ERROR in BMI_base.read_config_file().'
                    print '   Unsupported data type = ' + var_type + '.'
                    print ' '

                last_var_name = var_name

    #   read_config_file()
    #-------------------------------------------------------------------
    def initialize_config_vars(self):
   
        # print '## At start of initialize_config_vars(): cfg_prefix =', self.cfg_prefix
        if (self.cfg_prefix == None):
                self.cfg_prefix = self.case_prefix  # (10/25/11)
             
        #---------------------------------------
        # Read input variables from a CFG file
        #---------------------------------------
        # CFG filename was built and saved in
        # the initialize() method before now.
        #---------------------------------------
        # print '#### CALLING read_config_file()...'
        self.read_config_file()
        # print '#### AFTER read_config_file():'
        # print '#### in_directory  =', self.in_directory
        # print '#### out_directory =', self.out_directory
        # print ' '

        #-------------------------------------
        # Check the directories and prefixes
        #-------------------------------------
        # Defaults are set in __init__() and
        # some may have just been read in.
        #-------------------------------------
        # print '#### CALLING check_directories()...'
        self.check_directories()  # (Moved up: 10/25/11)
        
        #--------------------------------------------
        # Let driver set CWD to its in_directory ??
        #---------------------------------------------
        # if (self.mode == "driver"):
        #    if (hasattr(self, 'in_directory')):
        #        if (self.in_directory != None):
        #            os.chdir( self.in_directory )

        #--------------------------------------------------
        # (10/27/11) Maybe set to out_directory instead?
        # First try not setting it at all.  All filenames
        # should be fully qualified paths.
        #--------------------------------------------------
        # if (self.mode == "driver"):
        #     if (hasattr(self, 'out_directory')):
        #         if (self.out_directory != None):
        #             os.chdir( self.out_directory )

        #-------------------------------------------------------
        # These comments apply to idea of using "sync_files".
        #-------------------------------------------------------
        # If component is in "driver" mode, then set its
        # in_directory as current working directory (CWD).
        #
        # There are several reasons for this:
        #
        # (1) Each component's "read_cfg_file()" method
        #     needs to know where to find its CFG file.
        #     We _could_ store the paths to all of a project's
        #     CFG files in the sync_file (created by driver).
        #     However, assuming that all of the CFG files are
        #     in the CWD is simpler and keeps them together.
        #
        # (2) If a user is running an example w/o the GUI,
        #     then it is simplest if the CFG files needed for
        #     the example are stored together with the input
        #     data for the example in a central place on
        #     beach.  It is likely that the user won't have
        #     write permission there, which is why we need a
        #     separate output directory for model output.
        #     We also want to avoid the need for users to
        #     create copies of example CFG files in their own
        #     directories. (But they will have read permission
        #     and can make copies if they want.)
        #
        # (3) Component CFG files represent model input and
        #     may contain input parameters and/or the names
        #     of input files, such as initial-value grids.
        #     If a user is not running an example, then they
        #     will need to create an appropriate set of input
        #     files using whatever tools we give them.
        #
        # (4) Except when running examples (or using someone
        #     else's input files) the directories for input
        #     and output will typically be the same.
        #
        #-------------------------------------------------------
        # The driver also creates a "sync file" to share
        # the site_prefix, case_prefix, etc. with
        # the nondrivers.  This needs to be a name that is
        # not already in use and may be stored in "/tmp".
        #
        # self.mode is set to "driver" or "nondriver" in the
        # initialize() method.
        #-------------------------------------------------------
##        if (self.mode == "driver"):
##            os.chdir( self.in_directory )
##            ## self.write_sync_file()
##        else:
##            #------------------------------------------------
##            # Use information from sync file, such as
##            # site_prefix and case_prefix, to override
##            # those of a component in "nondriver" mode.
##            #------------------------------------------------
##            # The driver passes one argument, the name of
##            # the sync_file to the initialize() method of
##            # each nondriver component.
##            #------------------------------------------------            
##            ## self.read_sync_file()
##            pass
        
        #--------------------------------------------
        # Set any input variables that are computed
        #--------------------------------------------
        # print '#### CALLING set_computed_input_vars()...'
        self.set_computed_input_vars()
        
        #-----------------------------------------------------
        # Not all components need this, so don't do it here.
        #-----------------------------------------------------
        # self.read_grid_info()   # (stores rti in self, adds "da")
        
    #   initialize_config_vars()
    #-------------------------------------------------------------------
    def set_computed_input_vars( self):

        pass
    
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def initialize_basin_vars( self ):

        #------------------------------------------------------------
        # This saves outlet_file, outlet_IDs, outlet_ID, n_outlets,
        # basin_area and basin_relief into self.  (9/19/14)
        #------------------------------------------------------------
        outlets.read_outlet_file( self )

    #   initialize_basin_vars()    
    #-------------------------------------------------------------------
    def initialize_basin_vars0( self ):

        #------------------------------------------------------------
        # Notes: Most of the TopoFlow and Erode components have the
        #        ability to save values at monitored pixels.  This
        #        method supports this by embedding an instance of
        #        "basins_component" within each component.
        #------------------------------------------------------------
        # Notes: The routine BMI_base.read_grid_info() looks first
        #        in "in_directory" and then in "out_directory" for
        #        the RTI file, which supports TopoFlow and Erode.
        #------------------------------------------------------------
        ## from topoflow.utils import basins
        ## import basins
        
        self.bp = basins.basins_component()

##        print 'self.cfg_prefix  =', self.cfg_prefix
##        print 'self.site_prefix =', self.site_prefix    ##########
##        print 'self.case_prefix =', self.case_prefix

        #-----------------------------------------------------
        # Copy all of this from the "host" component, since
        # "basin components" don't have CFG files.\
        #-----------------------------------------------------
        # Added "out_directory" and "case_prefix" for Erode
        # on (11/5/13) to fix a bug.
        #-----------------------------------------------------
        self.bp.site_prefix   = self.site_prefix
        self.bp.case_prefix   = self.case_prefix
        self.bp.in_directory  = self.in_directory
        self.bp.out_directory = self.out_directory

        ## This isn't actually used by basins.initialize.
        cfg_file = (self.in_directory + self.site_prefix + '.rti')
        
        self.bp.initialize( cfg_file=cfg_file,
                            SILENT=not(self.DEBUG) )

        #-------------------------------
        # Store the outlet IDs in self
        #-------------------------------
        outlet_IDs = self.bp.outlet_IDs
        outlet_ID  = outlet_IDs[0]
        self.outlet_IDs = (outlet_IDs / self.nx, outlet_IDs % self.nx)
        self.outlet_ID  = (outlet_ID  / self.nx, outlet_ID  % self.nx)
        
##        self.outlet_IDs = outlet_IDs   # (long-int calendar indices)
##        self.outlet_ID  = outlet_ID

        #--------------------------------------------------
        # Before 5/14/10, get outlet_IDs from the
        # basins port of a Basins component.  Also,
        # every initialize() called "store_outlet_IDs()".
        #--------------------------------------------------
        # outlet_IDs = self.bp.get_vector_long('outlet_IDs')
        
    #   initialize_basin_vars0()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def prepend_directory(self, file_list, INPUT=True):

        #-----------------------------------------------------------
        # (11/14/11) Call from a component's open_input_file() 
        # method something like this:
        #   self.prepend_directory( ['slope_file', 'width_file'] ) 
        #-----------------------------------------------------------
        if (INPUT):
            dir_part = " = self.in_directory + "
        else:
            dir_part = " = self.out_directory + "

        for file_str in file_list:
            self_part = "self." + file_str
            exec( 'filename = ' + self_part )
            if (filename != ''):
                exec( self_part + dir_part + self_part ) 

    #   prepend_directory
    #-------------------------------------------------------------------
    def check_directories(self):

        #----------------------------------------
        # Note:  Defaults are set in __init_().
        #----------------------------------------
        if (self.in_directory == None):
            self.in_directory = os.getcwd() + os.sep
            ## self.in_directory = self.cfg_directory
        #-----------------------------------------------
        if (self.out_directory == None):
            self.out_directory = os.getcwd() + os.sep
        #-----------------------------------------------
        if (self.site_prefix == None):
            self.site_prefix = self.cfg_prefix
        #-----------------------------------------------
        if (self.case_prefix == None):
            self.case_prefix = self.cfg_prefix

##        print 'self.in_directory  =', self.in_directory
##        print 'self.out_directory =', self.out_directory

        #------------------------------------------------------------
        # In CFG files, the input directory is often set to ".",
        # which indicates the directory that contains the CFG file.
        # However, once "os.chdir()" is called, "." will expand to
        # something else.  Also, we want to avoid calling
        # "os.chdir()", because it creates problems for finding
        # package paths.  So to address these issues, we expand
        # the "." to the full CFG_file directory. (9/21/14)
        #------------------------------------------------------------
        if (self.cfg_file != None):
            cfg_directory = os.path.dirname(os.path.realpath(self.cfg_file))
            ## print 'cfg_directory =', cfg_directory
            self.cfg_directory = cfg_directory
            if (self.in_directory[0] == '.'):
                self.in_directory = self.cfg_directory
                
        #------------------------------------------------------
        # Expand path abbreviations: "." and "..", but NOT
        # "~" (11/5/13)
        #------------------------------------------------------
        # Note: os.path.expanduser() does not expand the
        # relative paths ".", ".." or "./". (2/13/12)
        #-----------------------------------------------------
        # Note: This removes trailing path separator !!
        #------------------------------------------------------
##        self.in_directory  = os.path.realpath( self.in_directory  )
##        self.out_directory = os.path.realpath( self.out_directory )

        #-----------------------------------------------
        # Expand path abbreviations like "~" (5/19/10)
        #---------------------------------------------------
        # Note that os.path.expanduser() does not
        # expand the relative paths "." or "./". (2/13/12)
        #------------------------------------------------------
        # Note: This does NOT remove trailing path separator.
        #------------------------------------------------------
        self.in_directory  = os.path.expanduser( self.in_directory  )
        self.out_directory = os.path.expanduser( self.out_directory )
        
        #--------------------------------------------------
        # Add trailing separator to directory, if missing
        # Note: Must come AFTER os.path.realpath calls.
        #--------------------------------------------------
        # (self.in_directory != ''):
        if (self.in_directory[-1] != os.sep):
            self.in_directory += os.sep
        #----------------------------------------
        # (self.out_directory != ''):
        if (self.out_directory[-1] != os.sep):
            self.out_directory += os.sep
   
    #   check_directories()
    #-------------------------------------------------------------------
    def initialize_scalar(self, value=0.0, dtype='float64'):

        #--------------------------------------------------------
        # If we create scalars as 0D numpy arrays, then:
        # - return var is MUTABLE and other components with
        #   a reference to it WILL see it change.
        # - Values will print without square brackets.
        # - Trying to subscript it will generate the error:
        #   "IndexError: 0-d arrays can't be indexed"
        # - We CAN compare this type of scalar to others with
        #   "<", ">", "==", etc.
        # - We can get the data type with ".dtype".
        # - The data type will be "numpy.ndarray".
        # - numpy.rank(x) will return 0.
        #--------------------------------------------------------
        # In order to preserve the reference (mutable scalar),
        # assignments must be done carefully.
        #
        # - Direct assignment or replacement:
        #   x=5 will BREAK the reference, but
        #   x.fill(5) is OK.
        #   x[:]=5 will generate this error:
        #      "ValueError: cannot slice a 0-d array"
        # 
        # - Incrementing the value:
        #   x = (x+1) will BREAK the reference, but
        #   x += 1 is OK.
        #--------------------------------------------------------
        # Recall that most numpy array operators have an
        # optional third argument that allows "in-place" calcs.
        #--------------------------------------------------------        
        return np.array(value, dtype=dtype)
    
        #--------------------------------------------------------
        # If we create scalars like this, then:
        # - return var is MUTABLE and other components with
        #   a reference to it WILL see it change.
        # - Values will print WITH square brackets unless we
        #   ask for element 0 "[0]".
        # - We CAN compare this type of scalar to others with
        #    "<", ">", "==", etc.
        # - We can get the data type with ".dtype".
        # - The data type will be "numpy.ndarray".
        # - numpy.rank(x) will return 1.
        # - They are really 1D arrays with one element.
        #--------------------------------------------------------
        ## return np.array([value], dtype=dtype)

        #---------------------------------------------------------
        # If we create scalars like this, then:
        # -  return var is IMMUTABLE and other components with
        #    a reference to it WILL NOT see it change.
        # - We CAN compare this type of scalar to others with
        #    "<", ">", "==", etc.
        # - We can get the data type with ".dtype".
        # - The data type will be "numpy.float64".
        # - numpy.rank(x) will return 0.
        #---------------------------------------------------------
        # return np.float64(0)
        
    #   initialize_scalar()
    #-------------------------------------------------------------------
    # These are for convenience;  not part of BMI.
    #-------------------------------------------------------------------
    def is_scalar(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #-------------------------------------------------      
        exec("n = np.rank(self." + var_name + ")")       
        return (n == 0)
    
    #   is_scalar()
    #-------------------------------------------------------------------
    def is_vector(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------     
        exec("n = np.rank(self." + var_name + ")")       
        return (n == 1)
    
    #   is_vector()
    #-------------------------------------------------------------------
    def is_grid(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------ 

        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------         
##        exec("type_str = str(type(self." + var_name + "))")
##        p1 = type_str.find("ndarray")
##        p2 = type_str.find("float")
##        if (p1 == -1) and (p2 == -1):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------        
##        if ("ndarray" not in type_str) and \
##           ("float" not in type_str):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------------        
        exec("n = np.rank(self." + var_name + ")")
        return (n == 2)

    #   is_grid()
    #-------------------------------------------------------------------

    




    
    
