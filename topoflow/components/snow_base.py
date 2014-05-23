
## Copyright (c) 2001-2013, Scott D. Peckham
## January 2009  (converted from IDL)
## May, July, August 2009
## May 2010 (changes to initialize() and read_cfg_file()

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "base class" for snowmelt
#          components as well as functions used by most or
#          all snowmelt methods.  The methods of this class
#          should be over-ridden as necessary for different
#          methods of modeling snowmelt.

#          update_snow_vars() in precip.py sets values here.
#-----------------------------------------------------------------------
#
#  class snow_component    (inherits from CSDMS_base.py)
#
#      set_constants()
#      -----------------------
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      --------------------------
#      check_input_types()
#      initialize_computed_vars()
#      ----------------------------
#      update_meltrate()
#      update_SM_integral()
#      update_depth()
#      update_swe()
#      -----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      -----------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()

#  Functions:
#      Initial_Cold_Content()
#      Max_Meltrate()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.utils import BMI_base
from topoflow.utils import model_input
from topoflow.utils import model_output

#-----------------------------------------------------------------------
class snow_component( BMI_base.BMI_component ):


    #------------------------------------------------------------
    # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are currently
    #        hardwired (not adjustable with GUI).

    #        Cp_snow is from NCAR CSM Flux Coupler web page

    #        h_snow = snow depth, h0_snow = initial depth
    
    #        hs = snow depth
    #        sw = snow water equivalent (depth)
    #        mr = snow melt rate
    #        cc = cold content
    #-------------------------------------------------------------------
    def set_constants(self):

        #-----------------------------------
        # Constants not changeable by user
        #-----------------------------------
        self.Cp_snow  = np.float64( 2090.0 )

        #--------------------------------------
        # Not a constant; read from CFG file.
        #--------------------------------------
        ## self.rho_snow = np.float64(300)
        ## self.rho_H2O  = np.float64(1000)  # (See initialize() method.)
        
    #   set_constants()         
    #-------------------------------------------------------------------
    def initialize(self, cfg_prefix=None, mode="nondriver",
                   SILENT=False):

        #---------------------------------------------------------
        # Notes:  Need to make sure than h_swe matches h_snow ?
        #         User may have entered incompatible values.
        #---------------------------------------------------------
        # (3/14/07) If the Energy Balance method is used for ET,
        # then we must initialize and track snow depth even if
        # there is no snowmelt method because the snow depth
        # affects the ET rate.  Otherwise, return to caller.
        #---------------------------------------------------------
        if not(SILENT):
            print ' '
            print 'Snow component: Initializing...'
            
        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_prefix = cfg_prefix
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------       
        self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars()

        ############################################
        # Commented out for new method, 3/15/12.
        ############################################
        #-----------------------------------------------------------
        # Need this, even if (method==0) for mp.update_snow_vars()
        #-----------------------------------------------------------
##        self.initialize_required_components(mode)  # (NEED BEFORE NEXT LINE)
##        self.rho_H2O = self.mp.get_scalar_double('rho_H2O')
   
        if (self.comp_status == 'Disabled'):
            #########################################
            #  DOUBLE CHECK THIS; SEE NOTES ABOVE
            #########################################
               ####### and (ep.method != 2):  ??????
            if not(SILENT):
                print 'Snow component: Disabled.'
            self.h_snow = self.initialize_scalar(0, dtype='float64')
            self.h_swe  = self.initialize_scalar(0, dtype='float64')
            self.SM     = self.initialize_scalar(0, dtype='float64')
            self.vol_SM = self.initialize_scalar(0, dtype='float64') # [m3]
            self.DONE   = True
            self.status = 'initialized'
            return
 
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #---------------------------
        # Initialize computed vars
        #---------------------------
        self.check_input_types()  # (maybe not used yet)
        self.initialize_computed_vars()  # (h_snow, h_swe, etc.)

        self.open_output_files()
        self.status = 'initialized'        
        
    #   initialize()
    #-------------------------------------------------------------------
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        #----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        #----------------------------------------------------------
        
        #-------------------------------------------------
        # Note: self.SM already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI)
        
        #-------------------------
        # Update computed values 
        #-------------------------
        self.update_meltrate()    # (meltrate = SM)
        self.update_SM_integral()
        self.enforce_max_meltrate()
        self.update_depth()
        self.update_swe()
        
        #-----------------------------------------
        # Read next snow vars from input files ?
        #-------------------------------------------
        # Note that read_input_files() is called
        # by initialize() and these values must be
        # used for "update" calls before reading
        # new ones.
        #-------------------------------------------
        if (self.time_index > 0):
            self.read_input_files()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        self.write_output_files()
        ## self.write_output_files( time_seconds )

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI)
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI)   
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Snow component')

        #---------------------------
        # Release all of the ports
        #----------------------------------------
        # Make this call in "finalize()" method
        # of the component's CCA Imple file
        #----------------------------------------
        # self.release_cca_ports( port_names, d_services )
        
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = np.maximum(self.save_grid_dt,   self.dt)
        self.save_pixels_dt = np.maximum(self.save_pixels_dt, self.dt)
        
    #   set_computed_input_vars()        
    #-------------------------------------------------------------------
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing snow meltrate.
        #----------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #----------------------------------------------------        
        are_scalars = np.array([
##                          self.pp.is_scalar('rate'),
##                          self.pp.is_scalar('duration'),
                          #----------------------------------
##                          self.mp.is_scalar('P'),
##                          self.mp.is_scalar('rho_H2O'),
##                          self.mp.is_scalar('rho_air'),
##                          self.mp.is_scalar('Cp_air'),
                          #----------------------------------
                          self.is_scalar('P'),
                          self.is_scalar('rho_H2O'),
                          self.is_scalar('rho_air'),
                          self.is_scalar('Cp_air'),
                          #----------------------------------
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Cp_snow'),
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe') ])

        self.ALL_SCALARS = np.all(are_scalars)
  
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #------------------------------------------
        # If T_air or precip are grids, then make
        # sure that h_snow and h_swe are grids
        #------------------------------------------
        T_IS_GRID = self.is_grid('T_air')
        P_IS_GRID = self.is_grid('P')
        ## P_IS_GRID = self.is_grid('rate')
##        T_IS_GRID = self.mp.is_grid('T_air')
##        P_IS_GRID = self.mp.is_grid('P')
##        ## P_IS_GRID = self.pp.is_grid('rate')

        H0_SNOW_IS_SCALAR = self.is_scalar('h0_snow')
        H0_SWE_IS_SCALAR  = self.is_scalar('h0_swe') 

        #------------------------------------------------------
        # If h0_snow or h0_swe are scalars, the use of copy()
        # here requires they were converted to numpy scalars.
        # Using copy() may not be necessary for scalars.
        #------------------------------------------------------
        h_snow = self.h0_snow.copy()    # [meters]
        h_swe  = self.h0_swe.copy()     # [meters]
        
        if (T_IS_GRID or P_IS_GRID):
            if (H0_SNOW_IS_SCALAR):
                self.h_snow = h_snow + zeros([self.ny, self.nx], dtype='Float64')
            else:
                self.h_snow = h_snow  # (is already a grid)
            #------------------------------------------------
            if (H0_SWE_IS_SCALAR):
                self.h_swe = h_swe + zeros([self.ny, self.nx], dtype='Float64')
            else:
                self.h_swe = h_swe    # (is already a grid)              
        else:
            self.h_snow = h_snow      # (both are scalars and that's OK)
            self.h_swe  = h_swe

        self.SM     = self.initialize_scalar( 0, dtype='float64')
        self.vol_SM = self.initialize_scalar( 0, dtype='float64') # (m3)
        
        #----------------------------------------------------
        # Initialize the cold content of snowpack (2/21/07)
        #----------------------------------------------------
        # Should we rename T_surf to T_snow for clarity ??
        #----------------------------------------------------
        T_surf   = self.T_surf  # (2/3/13, new framework)
        #----------------------------------------------------
##        T_surf   = self.get_port_data('T_surf', self.mp)
        self.Ecc = Initial_Cold_Content(self.h0_snow, T_surf, \
                                        self.rho_snow, self.Cp_snow)
        
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_meltrate(self):
    
        #--------------------------------------------------
        # Note: We don't need to update any variables if
        #       the snowmelt method is None.  But we need
        #       to make sure that self.SM = 0.0.
        #       This "method" will be over-ridden by a
        #       particular snowmelt method.
        #--------------------------------------------------
        print "WARNING: 'update_meltrate' method is inactive."
       
    #   update_meltrate()
    #-------------------------------------------------------------------
    def update_SM_integral(self):

        #------------------------------------------------
        # Update mass total for SM, sum over all pixels
        #------------------------------------------------   
        volume = np.float64(self.SM * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_SM += (volume * self.rti.n_pixels)
        else:
            self.vol_SM += np.sum(volume)  #### np.sum vs. sum ???
            
    #   update_SM_integral()
    #-------------------------------------------------------------------
    def enforce_max_meltrate(self):

        #-----------------------------------
        # Make sure that meltrate does not
        # exceed max possible meltrate
        #-----------------------------------     
        SM_max  = Max_Meltrate(self.h_snow, self.rho_H2O,
                               self.rho_snow, self.dt)
        self.SM = np.minimum(self.SM, SM_max)

        #------------------------------------------------------
        # Make sure meltrate is positive, while we're at it ?
        # Is already done by "Energy-Balance" component.
        #------------------------------------------------------
        self.SM = np.maximum(self.SM, np.float64(0))
   
    #   enforce_max_meltrate()
    #-------------------------------------------------------------------
##    def update_depth_from_snowfall(self):
##
##        #---------------------------------------------------
##        # Note:  See update_snow_vars() in precip_base.py.
##        #---------------------------------------------------
##        
##    #   update_depth_from_snowfall()     
    #-------------------------------------------------------------------
    def update_depth(self):

        #--------------------------------------------------------
        # Note: When precipitation falls as snow, then we need
        #       to update the snow depth at the precipitation
        #       timestep.  However, the snow depth, h_snow, is
        #       maintained and updated as part of the "snow
        #       process" component, using the "snow process"
        #       timestep.  What is the best solution?
        #--------------------------------------------------------

        #------------------------------------------
        # Increase snow depth due to falling snow
        #---------------------------------------------
        # This is currently done by pp.update() at
        # the precip timestep vs. the snow timestep.
        #---------------------------------------------
        
        #-------------------------------------
        # Decrease snow depth due to melting
        #-------------------------------------   
        ratio = (self.rho_H2O / self.rho_snow)
        dh    = self.SM * ratio * self.dt
        self.h_snow = np.maximum((self.h_snow - dh), np.float64(0))
        
    #   update_depth() 
    #-------------------------------------------------------------------
    def update_swe(self):

        #---------------------------------------------
        # If P or T_air is a grid, then we must have
        # h_swe and h_snow be grids.  This is set
        # up at start of Route_Flow.
        #---------------------------------------------
##        P = self.P          # (2/3/13, new framework)
##        T_air = self.T_air  # (2/3/13, new framework)
        #---------------------------------------------       
##        P     = self.get_port_data('P',     self.pp)
##        T_air = self.get_port_data('T_air', self.mp)
##        
##        #------------------------------------------------
##        # Increase snow water equivalent due to snowfall
##        #------------------------------------------------
##        # (3/14/07) New method that works regardless
##        # of whether P and T are scalars or grids.
##        #---------------------------------------------
##        dh1_swe     = (P * (T_air <= 0)) * self.main_dt
##        self.h_swe += dh1_swe
        
        #------------------------------------------------
        # Decrease snow water equivalent due to melting
        #------------------------------------------------
        dh2_swe    = self.SM * self.dt
        self.h_swe = np.maximum((self.h_swe - dh2_swe), np.float64(0))
        
    #   update_swe()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        self.c0_file       = self.in_directory + self.c0_file
        self.T0_file       = self.in_directory + self.T0_file
        self.rho_snow_file = self.in_directory + self.rho_snow_file
        self.h0_snow_file  = self.in_directory + self.h0_snow_file
        self.h0_swe_file   = self.in_directory + self.h0_swe_file

        self.c0_unit       = model_input.open_file(self.c0_type,       self.c0_file)
        self.T0_unit       = model_input.open_file(self.T0_type,       self.T0_file)
        self.rho_snow_unit = model_input.open_file(self.rho_snow_type, self.rho_snow_file)
        self.h0_snow_unit  = model_input.open_file(self.h0_snow_type,  self.h0_snow_file)
        self.h0_swe_unit   = model_input.open_file(self.h0_swe_type,   self.h0_swe_file)

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        c0 = model_input.read_next(self.c0_unit, self.c0_type, rti)
        if (c0 != None): self.c0 = c0

        T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        if (T0 != None): self.T0 = T0

        rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        if (rho_snow != None): self.rho_snow = rho_snow

        h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        if (h0_snow != None): self.h0_snow = h0_snow
        
        h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        if (h0_swe != None): self.h0_swe = h0_swe
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.c0_type       != 'Scalar'): self.c0_unit.close()        
        if (self.T0_type       != 'Scalar'): self.T0_unit.close()
        if (self.rho_snow_type != 'Scalar'): self.rho_snow_unit.close()
        if (self.h0_snow_type  != 'Scalar'): self.h0_snow_unit.close()
        if (self.h0_swe_type   != 'Scalar'): self.h0_swe_unit.close()
        
##        if (self.c0_file       != ''): self.c0_unit.close()        
##        if (self.T0_file       != ''): self.T0_unit.close()
##        if (self.rho_snow_file != ''): self.rho_snow_unit.close()
##        if (self.h0_snow_file  != ''): self.h0_snow_unit.close()
##        if (self.h0_swe_file   != ''): self.h0_swe_unit.close()

    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.mr_gs_file = (self.out_directory + self.mr_gs_file)
        self.hs_gs_file = (self.out_directory + self.hs_gs_file)
        self.sw_gs_file = (self.out_directory + self.sw_gs_file)
        self.cc_gs_file = (self.out_directory + self.cc_gs_file)
        #---------------------------------------------------------
        self.mr_ts_file = (self.out_directory + self.mr_ts_file)
        self.hs_ts_file = (self.out_directory + self.hs_ts_file)
        self.sw_ts_file = (self.out_directory + self.sw_ts_file)
        self.cc_ts_file = (self.out_directory + self.cc_ts_file)

        
##        self.mr_gs_file = (self.case_prefix + '_2D-SMrate.rts')
##        self.hs_gs_file = (self.case_prefix + '_2D-hsnow.rts')
##        self.sw_gs_file = (self.case_prefix + '_2D-hswe.rts')
##        self.cc_gs_file = (self.case_prefix + '_2D-Ecc.rts')
##        #-----------------------------------------------------------
##        self.mr_ts_file = (self.case_prefix + '_0D-SMrate.txt')
##        self.hs_ts_file = (self.case_prefix + '_0D-hsnow.txt')
##        self.sw_ts_file = (self.case_prefix + '_0D-hswe.txt')
##        self.cc_ts_file = (self.case_prefix + '_0D-Ecc.txt')

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_nio()
        self.update_outfile_names()
        
        #----------------------------------
        # Open files to write grid stacks
        #----------------------------------
        if (self.SAVE_MR_GRIDS):
            model_output.open_new_gs_file( self, self.mr_gs_file, self.rti,
                                           ## var_name='MR',
                                           var_name='mr',
                                           long_name='snow_meltrate',
                                           units_name='m/s')
            
        if (self.SAVE_HS_GRIDS):
            model_output.open_new_gs_file( self, self.hs_gs_file, self.rti,
                                           ## var_name='h_snow',
                                           var_name='hs',
                                           long_name='snow_depth',
                                           units_name='m')
            
        if (self.SAVE_SW_GRIDS):
            model_output.open_new_gs_file( self, self.sw_gs_file, self.rti,
                                           ## var_name='SWE',
                                           var_name='sw',
                                           long_name='snow_water_equivalent',
                                           units_name='m')
            
        if (self.SAVE_CC_GRIDS):
            model_output.open_new_gs_file( self, self.cc_gs_file, self.rti,
                                           ## var_name='SCC',
                                           var_name='cc',
                                           long_name='snow_cold_content',
                                           units_name='J/m^2')

        #---------------------------------------
        # Open text files to write time series
        #---------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_MR_PIXELS):
            model_output.open_new_ts_file( self, self.mr_ts_file, IDs,
                                           ## var_name='MR',
                                           var_name='mr',
                                           long_name='snow_meltrate',
                                           units_name='m/s')

        if (self.SAVE_HS_PIXELS):
            model_output.open_new_ts_file( self, self.hs_ts_file, IDs,
                                           ## var_name='h_snow',
                                           var_name='hs',
                                           long_name='snow_depth',
                                           units_name='m')

        if (self.SAVE_SW_PIXELS):
            model_output.open_new_ts_file( self, self.sw_ts_file, IDs,
                                           ## var_name='SWE',
                                           var_name='sw',
                                           long_name='snow_water_equivalent',
                                           units_name='m')
            
        if (self.SAVE_CC_PIXELS):
            model_output.open_new_ts_file( self, self.cc_ts_file, IDs,
                                           ## var_name='SCC',
                                           var_name='cc',
                                           long_name='snow_cold_content',
                                           units_name='J/m^2')
            
    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

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
    
        if (self.SAVE_MR_GRIDS): model_output.close_gs_file( self, 'mr')   
        if (self.SAVE_HS_GRIDS): model_output.close_gs_file( self, 'hs')   
        if (self.SAVE_SW_GRIDS): model_output.close_gs_file( self, 'sw')   
        if (self.SAVE_CC_GRIDS): model_output.close_gs_file( self, 'cc')
        #-----------------------------------------------------------------        
        if (self.SAVE_MR_PIXELS): model_output.close_ts_file( self, 'mr')  
        if (self.SAVE_HS_PIXELS): model_output.close_ts_file( self, 'hs')   
        if (self.SAVE_SW_PIXELS): model_output.close_ts_file( self, 'sw')   
        if (self.SAVE_CC_PIXELS): model_output.close_ts_file( self, 'cc')
        
    #-------------------------------------------------------------------  
    def save_grids(self):
     
        if (self.SAVE_MR_GRIDS):
            model_output.add_grid( self, self.SM, 'mr', self.time_min )
            
        if (self.SAVE_HS_GRIDS):
            model_output.add_grid( self, self.h_snow, 'hs', self.time_min )
            
        if (self.SAVE_SW_GRIDS):
            model_output.add_grid( self, self.h_swe, 'sw', self.time_min )

        if (self.SAVE_CC_GRIDS):
            model_output.add_grid( self, self.Ecc, 'cc', self.time_min )

    #   save_grids()     
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min   ###
        
        if (self.SAVE_MR_PIXELS):
            model_output.add_values_at_IDs( self, time, self.SM, 'mr', IDs )
            
        if (self.SAVE_HS_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_snow, 'hs', IDs )
            
        if (self.SAVE_SW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_swe, 'sw', IDs )
            
        if (self.SAVE_CC_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Ecc, 'cc', IDs )

    #   save_pixel_values()
#---------------------------------------------------------------------
#---------------------------------------------------------------------
def Initial_Cold_Content(h_snow, T_snow, rho_snow, Cp_snow):

    #----------------------------------------------------------------
    #NOTES:  This function is used to initialize the cold content
    #        of a snow pack (in Initialize_Snow_Vars in route.pro).
    #        The cold content has units of [J/m^2] (_NOT_ [W/m^2]).
    #        It is an energy (per unit area) threshold (or deficit)
    #        that must be overcome before melting of snow can occur.
    #        Cold content changes over time as the snowpack warms or
    #        cools, but must always be non-negative.  See the Notes
    #        for the Energy_Balance_Meltrate function.

    #        Caller sets T_snow argument to T_surf.

    #        K_snow is between 0.063 and 0.71  [W/m/deg_C]
    #        All of the Q's have units of W/m^2 = J/(m^2 s).

    #        12/22/05.  Removed T0 from argument list and
    #        set to zero here.
    #---------------------------------------------------------------
  
    #-------------------------------------------
    # Compute initial cold content of snowpack
    #-------------------------------------------
    T0   = np.float64(0)
    Ecc0 = (rho_snow * Cp_snow) * h_snow * (T0 - T_snow)
    
    return Ecc0

#   Initial_Cold_Content
#---------------------------------------------------------------------
def Max_Meltrate(h_snow, rho_H2O, rho_snow, snow_dt):

    #--------------------------------------------------------
    #NOTES:  This function returns the max possible meltrate
    #        for snow, assuming that all snow (given by snow
    #        depth) is melted in the time interval, snow_dt.
    #        Snow meltrates should never exceed this value.

    #        h_snow, rho_H2O, and rho_snow are pointers.
    #        snow_dt is real-valued.
    #--------------------------------------------------------  
    M_max = h_snow * (rho_H2O / rho_snow) / snow_dt    #[m/s]
    
    return M_max

#   Max_Meltrate
#---------------------------------------------------------------------
