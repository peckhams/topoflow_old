
## Copyright (c) 2009-2013, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## May, July, August 2009
## May 2010 (changes to unit_test() and read_cfg_file()

#################################################
#  NB!  See note at "THIS MAY BE COSTLY"
#################################################

#-----------------------------------------------------------------------
#  NOTES:  This file defines a Priestley-Taylor ET component
#          and related functions.  It inherits from the ET
#          "base class" in "evap_base.py".
#-----------------------------------------------------------------------
#
#  class evap_component
#
#      get_attribute()             # (10/26/11)
#      get_input_var_names()       # (5/15/12)
#      get_output_var_names()      # (5/15/12)
#      get_var_name()              # (5/16/12), Bolton
#      get_var_units()             # (5/16/12), Bolton
#      ------------------------
#      check_input_types()
#      update_ET_rate()
#      ------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#
#  Functions:
#      Priestley_Taylor_ET_Rate
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import evap_base

from topoflow.utils import model_input

#-----------------------------------------------------------------------
class evap_component( evap_base.evap_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Evaporation_Priestley_Taylor',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #------------------------------------------------------
        'comp_name':          'EvapPriestleyTaylor',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Evap_Priestley_Taylor.cfg.in',
        'cfg_extension':      '_evap_priestley_taylor.cfg',
        'cmt_var_prefix':     '/EvapPriestleyTaylor/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Evap_Priestley_Taylor.xml',
        'dialog_title':       'Evaporation: Priestley-Taylor Parameters',
        'time_units':         'seconds' }

    #----------------------------------------------------------------
    # Note that the "meteorology" component uses the following to
    # compute Q_sum and Qe, but they aren't needed directly here:
    #     uz, z, z0_air, rho_air, Cp_air, Qn_SW, Qn_LW 
    #----------------------------------------------------------------
    # TopoFlow currently assumes that "soil_layer_0__porosity" is
    # the same as "soil_layer_0__saturated_water_content".
    #----------------------------------------------------------------
    
    #############################################################
    # We need to add "porosity" as a reference to theta_sat in
    # the approprate place.
    #
    # Check algorithms that add or remove water from the soil;
    # make sure to account for the *current* water content.
    #############################################################
    # Note that evap_base.py is inherited and contains additional
    # methods that need vars for top soil layer.  Note that
    # "model_top_layer" = "model_layer_0".
    #############################################################
    _input_var_names = [
        'air__temperature',                             # (meteorology)
        'channel_water__depth',                         # (channels)
        'land_surface__net_longwave_irradiation_flux',  # (meteorology)
        'land_surface__net_shortwave_irradiation_flux', # (meteorology)
        'land_surface__temperature',                    # (meteorology)
        'soil_model_top_layer__porosity',               # (satzone)
        'soil_model_top_layer__wetted_thickness',       # (satzone)
        'soil_water_table_surface__elevation' ]         # (satzone)
                 
         #-------------------------------------------------
         # These are currently obtained from the GUI/file
         # and are not obtained from other components.
         #-------------------------------------------------
##        'land_surface__elevation',                        # (GUI, DEM)
##        'soil__priestley_taylor_evaporation_coefficient'  # (GUI, alpha)
##        'soil__reference_depth_temperature',     # (GUI, T_soil_x)
##        'soil_surface__temperature',             # (GUI, T_surf)
##        'soil__temperature_reference_depth',     # (GUI, soil_x)
##        'soil__thermal_conductivity' :           # (GUI, K_soil)
        
        #----------------------------------------------------
        # These could be added in the future; not used yet.
        #----------------------------------------------------
##        'soil_model_top_layer__saturated_water_content', # (satzone comp)
##        'land_surface_water_potential_evaporation_rate':'PET' }

    _output_var_names = [
        'land_water__evaporation_rate',  # (ET)
        'land_water__area_time_integral_of_evaporation_rate',  # (vol_ET)
        'model__time_step', # (dt)
        'soil_surface__conduction_energy_flux' ]  # (Qc)
        #-----------------------------------------------------
        # These are read from GUI/file, but can be returned.
        #-----------------------------------------------------       
        #'land_surface__elevation',
        #'soil__priestley_taylor_evaporation_coefficient',
        #'soil__reference_depth_temperature',
        ## 'soil_surface__temperature',
        #'soil__temperature_reference_depth',
        #'soil__thermal_conductivity' ]
        
    #----------------------------------------------------------------
    # Should we use "ponded_water__depth" or "surface_water__depth"
    # instead of "channel_water__depth" in this case ?
    #----------------------------------------------------------------
    # Should we use "soil_surface__temperature" or
    # "land_surface__temperature" here ?  (Both, for now.)
    #----------------------------------------------------------------   
    _var_name_map = {   
        'air__temperature' :                            'T_air',
        'channel_water__depth' :                        'depth',
        'land_surface__net_longwave_irradiation_flux':  'Qn_LW',
        'land_surface__net_shortwave_irradiation_flux': 'Qn_SW',
        'land_surface__temperature':                    'T_surf',
        'soil_model_top_layer__porosity':               'p0',
        'soil_model_top_layer__wetted_thickness' :      'y0',
        'soil_water_table_surface__elevation' :         'h_table',
        'snow__depth' :                                 'h_snow',
        #----------------------------------------------------------------
        'land_water__evaporation_rate' :                      'ET',
        'land_water__area_time_integral_of_evaporation_rate': 'vol_ET',
        'model__time_step':                                   'dt',
        'soil_surface__conduction_energy_flux' :              'Qc',   # (computed)
        #-----------------------------------------------------
        # These are read from GUI/file, but can be returned.
        #-----------------------------------------------------       
        'land_surface__elevation' :                       'DEM',
        'soil__priestley_taylor_evaporation_coefficient': 'alpha',
        'soil__reference_depth_temperature' :             'T_soil_x',
        # 'soil_surface__temperature' :                   'T_surf',    # (from met)
        'soil__temperature_reference_depth':              'soil_x',
        'soil__thermal_conductivity' :                    'K_soil' }   # (thermal !)
        
    #------------------------------------------------
    # What is the correct unit string for "deg_C" ?
    #------------------------------------------------
    _var_units_map = {
        'air__temperature' :                            'deg_C',
        'channel_water__depth' :                        'm',        
        'land_surface__net_longwave_irradiation_flux':  'W m-2',
        'land_surface__net_shortwave_irradiation_flux': 'W m-2',
        'land_surface__temperature':                    'deg_C',
        'soil_model_top_layer__porosity':               '1',      ### CHECK
        'soil_model_top_layer__wetted_thickness' :      'm',
        'soil_water_table_surface__elevation' :         'm',
        'snow__depth' :                                 'm',
        #------------------------------------------------------------ 
        'land_water__evaporation_rate' :         'm s-1',
        'land_water__area_time_integral_of_evaporation_rate': 'm3',
        'model__time_step' :                     's',
        'soil_surface__conduction_energy_flux' : 'W m-2',
        #-----------------------------------------------------
        # These are read from GUI/file, but can be returned.
        #-----------------------------------------------------
        'land_surface__elevation' :                       'm',
        'soil__priestley_taylor_evaporation_coefficient': '1',     ### CHECK
        'soil__reference_depth_temperature' :             'deg_C',
        # 'soil_surface__temperature' :                   'deg_C',
        'soil__temperature_reference_depth':              'm',
        'soil__thermal_conductivity' :                    'W m-1 deg_C-1]' } 
        
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )
            
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

    #   get_attribute()
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
    def check_input_types(self):

        #--------------------------------------------------------
        # As of 7/9/10, Qn_SW and Qn_LW are computed internally
        # from other vars, including slope and aspect grids.
        # So they'll always be grids and so will self.ET
        # unless PRECIP_ONLY = True.
        #--------------------------------------------------------
        # Need to get a boolean value, PRECIP_ONLY from the
        # "meteorology" component.  This doesn't work:
        # PRECIP_ONLY = self.get_port_data('PRECIP_ONLY', 'mp')
        # But the method just below does work.  Should we add
        # a new IRF port method called "get_boolean()" ??
        # Note that "self.mp" here is an embedded CCA port.
        #--------------------------------------------------------
        # (2/5/13) Modified for use with new framework.
        #--------------------------------------------------------        
##        PRECIP_ONLY = self.mp.get_scalar_long('PRECIP_ONLY')
##        Qn_SW_IS_SCALAR = (PRECIP_ONLY == 1)
##        Qn_LW_IS_SCALAR = (PRECIP_ONLY == 1)
##        if (self.DEBUG):
##            print 'In ET component: PRECIP_ONLY =', PRECIP_ONLY
##        ## Qn_SW_IS_SCALAR = self.mp.is_scalar('Qn_SW')
##        ## Qn_LW_IS_SCALAR = self.mp.is_scalar('Qn_LW')
        
        are_scalars = np.array([
                         self.is_scalar('T_soil_x'),
                         self.is_scalar('soil_x'),
                         self.is_scalar('K_soil'),
                         self.is_scalar('alpha'),
                         #-------------------------------
                         self.is_scalar('Qn_SW'),  # @met
                         self.is_scalar('Qn_LW'),  # @met
                         self.is_scalar('T_air'),  # @met
                         self.is_scalar('T_surf'), # @met
                         #---------------------------------
                         self.is_scalar('depth'),  # d@chan
                         #---------------------------------
                         self.is_scalar('h_table') ]) # satzone
                         #-------------------------------
##                         Qn_SW_IS_SCALAR,
##                         Qn_LW_IS_SCALAR,
##                         self.mp.is_scalar('T_air'),
##                         self.mp.is_scalar('T_surf'),
##                         #-------------------------------
##                         self.cp.is_scalar('d'),
##                         #-------------------------------
##                         self.gp.is_scalar('h_table') ])

        self.ALL_SCALARS = np.all(are_scalars)

        ## self.ALL_SCALARS = False
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def update_ET_rate(self):

        #--------------------------------------------------------------
        # Notes: Qet   = energy used for ET of water from surface
        #        Qn_SW = net shortwave irradiation flux (solar)
        #        Qn_LW = net longwave irradiation flux (air, surface)
        #        Qh    = sensible heat flux from turbulent convection
        #                between snow surface and air
        #        Qc    = energy transferred from surface to subsurface

        #        All of the Q's have units of [W/m^2].

        #        T_air    = air temperature [deg_C]
        #        T_surf   = soil temp at the surface [deg_C]
        #        T_soil_x = soil temp at depth of x meters [deg_C]

        #        Ks   = thermal conductivity of soil [W/m/deg_C]
        #        Ks = 0.45   ;[W/m/deg_C] (thawed soil; moisture content
        #                     near field capacity)
        #        Ks = 1.0    ;[W/m/deg_C] (frozen soil)

        #        alpha = evaporation parameter
        #        alpha = 0.95   ;(average found by Rouse)
        #        alpha = 1.26   ;(Jackson et al. (1996), at saturation)

        #        Modification of alpha:  alpha = (a1 * R) + a2
        #        R  = 1.0d   ;(equals 1 for saturation; R in [0,1])
        #        a1 = 1.0d   ;(accounts for moisture content of soil)
        #        a2 = 0.2d   ;(accounts for vegetation effect)
        #--------------------------------------------------------------
        Qn_SW  = self.Qn_SW   # (2/3/13, new framework)
        Qn_LW  = self.Qn_LW   # (2/3/13, new framework)
        T_air  = self.T_air   # (2/3/13, new framework)
        T_surf = self.T_surf  # (2/3/13, new framework)
        #------------------------------------------------
##        Qn_SW  = self.get_port_data('Qn_SW',  self.mp)
##        Qn_LW  = self.get_port_data('Qn_LW',  self.mp)
##        T_air  = self.get_port_data('T_air',  self.mp)
##        T_surf = self.get_port_data('T_surf', self.mp)
        #------------------------------------------------
        Q_net  = Qn_SW + Qn_LW
        
        #---------------------------------------------
        # Compute the conductive energy between the
        # surface and subsurface using Fourier's law
        #---------------------------------------------
        # soil_x is converted from [cm] to [m] when
        # it is read from the GUI and then stored
        #---------------------------------------------
        # In Qet formula, the constant 0.011 has
        # units of 1/[deg_C] to cancel T_air units.
        #---------------------------------------------  
        Qc  = self.K_soil * (self.T_soil_x - T_surf) / self.soil_x
        Qet = self.alpha * (np.float64(0.406) + (np.float64(0.011) * T_air)) * (Q_net - Qc)
        self.Qc = Qc  ## (2/3/13)
        
        #-----------------------------------
        # Convert ET energy to a loss rate
        #------------------------------------------
        # Lf = latent heat of fusion [J/kg]
        # Lv = latent heat of vaporization [J/kg]
        # ET = (Qet / (rho_w * Lv))
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lv    = -2500000d   ;[J/kg]
        # So (rho_w * Lv) = -2.5e+9  [J/m^3]
        #-------------------------------------
        ET = (Qet / np.float64(2.5E+9))  #[m/s]  (A loss, but returned as positive.)
        
        self.ET = np.maximum(ET, np.float64(0))

        ##########################################
        #  THIS MAY BE COSTLY.  BETTER WAY OR
        #  ALLOW ET TO BE A SCALAR ??
        ##########################################
        if (np.size(self.ET) == 1):
            self.ET += np.zeros((self.ny, self.nx), dtype='Float64')
            
    #   update_ET_rate()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #----------------------------------------------------
        # Note: Priestley-Taylor method needs alpha but the
        #       energy balance method doesn't. (2/5/13)
        #----------------------------------------------------
        self.alpha_file    = self.in_directory + self.alpha_file
        self.K_soil_file   = self.in_directory + self.K_soil_file
        self.soil_x_file   = self.in_directory + self.soil_x_file
        self.T_soil_x_file = self.in_directory + self.T_soil_x_file

        self.alpha_unit    = model_input.open_file(self.alpha_type,    self.alpha_file)
        self.K_soil_unit   = model_input.open_file(self.K_soil_type,   self.K_soil_file)
        self.soil_x_unit   = model_input.open_file(self.soil_x_type,   self.soil_x_file)
        self.T_soil_x_unit = model_input.open_file(self.T_soil_x_type, self.T_soil_x_file)
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti
        
        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        alpha = model_input.read_next(self.alpha_unit, self.alpha_type, rti)
        if (alpha != None): self.alpha = alpha

        K_soil = model_input.read_next(self.K_soil_unit, self.K_soil_type, rti)
        if (K_soil != None): self.K_soil = K_soil

        soil_x = model_input.read_next(self.soil_x_unit, self.soil_x_type, rti)
        if (soil_x != None): self.soil_x = soil_x

        T_soil_x = model_input.read_next(self.T_soil_x_unit, self.T_soil_x_type, rti)
        if (T_soil_x != None): self.T_soil_x = T_soil_x
        
    #   read_input_files()        
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.alpha_type    != 'Scalar'): self.alpha_unit.close()        
        if (self.K_soil_type   != 'Scalar'): self.K_soil_unit.close()
        if (self.soil_x_type   != 'Scalar'): self.soil_x_unit.close()
        if (self.T_soil_x_type != 'Scalar'): self.T_soil_x_unit.close()
        
##        if (self.alpha_file    != ''): self.alpha_unit.close()        
##        if (self.K_soil_file   != ''): self.K_soil_unit.close()
##        if (self.soil_x_file   != ''): self.soil_x_unit.close()
##        if (self.T_soil_x_file != ''): self.T_soil_x_unit.close()
        
    #   close_input_files()
    #-------------------------------------------------------------------


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------   
def Priestley_Taylor_ET_Rate(alpha, Ks, T_soil_x, soil_x, \
                             Qn_SW, Qn_LW, T_air, T_surf):

    #--------------------------------------------------------------
    # Notes: Qet   = energy used for ET of water from surface
    #        Qn_SW = net shortwave radiation flux (solar)
    #        Qn_LW = net longwave radiation flux (air, surface)
    #        Qh    = sensible heat flux from turbulent convection
    #                between snow surface and air
    #        Qc    = energy transferred from surface to subsurface

    #        All of the Q's have units of [W/m^2].

    #        T_air    = air temperature [deg_C]
    #        T_surf   = soil temp at the surface [deg_C]
    #        T_soil_x = soil temp at depth of x meters [deg_C]

    #        Ks   = thermal conductivity of soil [W/m/deg_C]
    #        Ks = 0.45   ;[W/m/deg_C] (thawed soil; moisture content
    #                     near field capacity)
    #        Ks = 1.0    ;[W/m/deg_C] (frozen soil)

    #        alpha = evaporation parameter
    #        alpha = 0.95   ;(average found by Rouse)
    #        alpha = 1.26   ;(Jackson et al. (1996), at saturation)

    #        Modification of alpha:  alpha = (a1 * R) + a2
    #        R  = 1.0d   ;(equals 1 for saturation; R in [0,1])
    #        a1 = 1.0d   ;(accounts for moisture content of soil)
    #        a2 = 0.2d   ;(accounts for vegetation effect)
    #--------------------------------------------------------------
    
    #---------------------------------------------
    # Compute the conductive energy between the
    # surface and subsurface using Fourier's law
    #---------------------------------------------
    # soil_x is converted from [cm] to [m] when
    # it is read from the GUI and then stored
    #---------------------------------------------
    # In Qet formula, the constant 0.011 has
    # units of 1/[deg_C] to cancel T_air units.
    #---------------------------------------------  
    Qc   = Ks * (T_soil_x - T_surf) / (soil_x)
    Qnet = Qn_SW + Qn_LW
    Qet  = alpha * (np.float64(0.406) + (np.float32(0.011) * T_air)) * (Qnet - Qc)
    
    #-----------------------------------
    # Convert ET energy to a loss rate
    #------------------------------------------
    # Lf = latent heat of fusion [J/kg]
    # Lv = latent heat of vaporization [J/kg]
    # ET = (Qet / (rho_w * Lv))
    #------------------------------------------
    # rho_w = 1000d       ;[kg/m^3]
    # Lv    = -2500000d   ;[J/kg]
    # So (rho_w * Lv) = -2.5e+9  [J/m^3]
    #-------------------------------------
    ET = (Qet / np.float32(2.5E+9))  #[m/s]  (A loss, but returned as positive.)
    
    return np.maximum(ET, np.float64(0))

#   Priestley_Taylor_ET_Rate
#-----------------------------------------------------------------------

