
## Copyright (c) 2001-2013, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## July, August 2009
## May 2010 (changes to unit_test() and read_cfg_file()

#-----------------------------------------------------------------------
#  NOTES:  This file defines a Degree-Day snowmelt component
#          and related functions.  It inherits from the snowmelt
#          "base class" in "snow_base.py".
#-----------------------------------------------------------------------
#
#  class snow_degree_day
#
#      get_attribute()          # (10/26/11)
#      get_input_var_names()    # (5/14/12)
#      get_output_var_names()   # (5/14/12)
#      get_var_name()           # (5/16/12, Bolton)
#      get_var_units()          # (5/16/12, Bolton)
#      -----------------------
#      check_input_types()
#      update_meltrate()
#
#-----------------------------------------------------------------------

import numpy as np

from topoflow.components import snow_base

#-----------------------------------------------------------------------
class snow_component( snow_base.snow_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Snowmelt_Degree_Day',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'SnowDegreeDay',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Snow_Degree_Day.cfg.in',
        'cfg_extension':      '_snow_degree_day.cfg',
        'cmt_var_prefix':     '/SnowDegreeDay/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Snow_Degree_Day.xml',
        'dialog_title':       'Snowmelt: Degree-Day Parameters',    # (Snowmelt ?)
        'time_units':         'seconds' }

    #----------------------------------------------------------
    # The only input variable that the Degree-Day method
    # needs from the Met component is "T_air".  The others
    # are used for consistent use of scalars vs. grids;
    # see "check_input_types()" below.
    #----------------------------------------------------------
    # Some input vars, like c0 and T0 are read from the CFG
    # file and not from other components.  They are therefore
    # included with the output_vars.
    #----------------------------------------------------------
    _input_var_names = [
        'air__density',
        'air__relative_humidity',
        'air__temperature',
        'air__thermal_capacity',
        'atmosphere_turbulent_boundary_layer__roughness_length',
        'atmosphere_water__liquid_equivalent_precipitation_rate',
        # 'land_surface_air__pressure',
        'land_surface__net_irradiation_flux',
        'land_surface__net_longwave_irradiation_flux',
        'land_surface__net_shortwave_irradiation_flux',
        'land_surface__temperature',        
        'water__density',
        'wind__speed_reference_height',
        'wind__reference_height_speed' ]

    #------------------------------------------------------------
    # Note: The main output of the Degree-Day method is SM.
    #       Functions in snow_base.py compute additional output
    #       vars, such as:
    #       update_SM_integral(), update_depth(), update_swe().
    #       They need things like: rho_H2O and rho_snow.
    #------------------------------------------------------------
    # Note: Cp_snow is a constant set in the "set_constants()"
    #       function in snow_base.py.
    #------------------------------------------------------------
    #       vol_SM was "basin_cumulative_snow_meltwater_volume"
    #------------------------------------------------------------
    _output_var_names = [
        'model__time_step',                        # dt        
        'snow__area_time_integral_of_melt_rate',   # vol_SM
        'snow__degree_day_coefficient',            # c0   (read from CFG)
        'snow__degree_day_threshold_temperature',  # T0   (read from CFG)
        'snow__density',                           # rho_snow
        'snow__depth',                             # h_snow
        'snow__initial_depth',                     # h0_snow
        'snow__initial_liquid_equivalent_depth',   # h0_swe
        'snow__liquid_equivalent_depth',           # h_swe
        'snow__melt_rate',                         # SM   (MR is used for ice)
        'snow__thermal_capacity' ]                 # Cp_snow
    
    _var_name_map = {
        'air__density': 'rho_air',
        'air__relative_humidity': 'RH',
        'air__temperature': 'T_air',
        'air__thermal_capacity': 'Cp_air',
        'atmosphere_turbulent_boundary_layer__roughness_length': 'z0_air',
        'atmosphere_water__liquid_equivalent_precipitation_rate': 'P',
        # 'land_surface_air__pressure': 'p0',   # (check if needed)
        'land_surface__net_irradiation_flux': 'Q_sum',
        'land_surface__net_longwave_irradiation_flux': 'Qn_LW',
        'land_surface__net_shortwave_irradiation_flux': 'Qn_SW',
        'land_surface__temperature': 'T_surf',        
        'water__density': 'rho_H2O',
        'wind__speed_reference_height': 'z',
        'wind__reference_height_speed': 'uz',
        #----------------------------------------------------------
        'model__time_step': 'dt',
        'snow__area_time_integral_of_melt_rate': 'vol_SM',
        'snow__degree_day_coefficient': 'c0',                  ########
        'snow__degree_day_threshold_temperature': 'T0',        ########
        'snow__density': 'rho_snow',
        'snow__depth': 'h_snow',
        'snow__initial_depth': 'h0_snow',
        'snow__initial_liquid_equivalent_depth': 'h0_swe',
        'snow__liquid_equivalent_depth': 'h_swe',
        'snow__melt_rate': 'SM',                # (MR is used for ice)
        'snow__thermal_capacity': 'Cp_snow' }
    
    #-----------------------------------------------------------------
    # Note: We need to be careful with whether units are C or K,
    #       for all "thermal" quantities (e.g. thermal_capacity).
    #-----------------------------------------------------------------       
    _var_units_map = {
        'air__density': 'kg m-3',
        'air__relative_humidity': '1',
        'air__temperature': 'C',                                  # (see Notes above)
        'air__thermal_capacity': 'J kg-1 C-1',                    # (see Notes above)
        'atmosphere_turbulent_boundary_layer__roughness_length': 'm',
        'atmosphere_water__liquid_equivalent_precipitation_rate': 'm s-1',
        # 'land_surface_air__pressure': 'mbar',
        'land_surface__net_irradiation_flux': 'W m-2',
        'land_surface__net_longwave_irradiation_flux': 'W m-2',
        'land_surface__net_shortwave_irradiation_flux': 'W m-2',
        'land_surface__temperature': 'C',        
        'water__density': 'kg m-3',
        'wind__speed_reference_height': 'm',
        'wind__reference_height_speed': 'm s-1',
        #----------------------------------------------------------
        'model__time_step': 's',
        'snow__area_time_integral_of_melt_rate': 'm3',
        'snow__degree_day_coefficient': 'mm day-1 C-1',     ##########
        'snow__degree_day_threshold_temperature': 'C',      ##########
        'snow__density': 'kg m-3',
        'snow__depth': 'm',
        'snow__initial_depth': 'm',
        'snow__initial_liquid_equivalent_depth': 'm',
        'snow__liquid_equivalent_depth': 'm',
        'snow__melt_rate': 'm s-1',
        'snow__thermal_capacity': 'J kg-1 C-1' }

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

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #--------------------------------------------------        
        are_scalars = np.array([
##                          self.pp.is_scalar('rate'),   # (not "rates")
##                          self.pp.is_scalar('duration'),
                          #-------------------------------
##                          self.mp.is_scalar('P'),
##                          self.mp.is_scalar('rho_H2O'),
##                          self.mp.is_scalar('rho_air'),
##                          self.mp.is_scalar('T_air'),
##                          self.mp.is_scalar('Cp_air'),
                          #-------------------------------
                          self.is_scalar('P'),
                          self.is_scalar('rho_H2O'),
                          self.is_scalar('rho_air'),
                          self.is_scalar('T_air'),
                          self.is_scalar('Cp_air'),
                          #-------------------------------
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Cp_snow'),
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe'),
                          self.is_scalar('c0'),
                          self.is_scalar('T0') ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def update_meltrate(self):

        #------------------------------------------------------------
        #Notes:  Arguments are assumed to be references to a scalar
        #        or grid so that:
        #        c0     = degree-day melt factor [mm/day/deg_C]
        #        T0     = threshold temperature [deg_C]
        #        M      = water equivalent of snowmelt [m/s]
        #        M_max  = max possible meltrate if all snow melts
        #        T_air  = air temperature [deg_C]

        #        Model must start when snow is isothermal.
        #        Cooling of the snowpack is not considered.

        #        This is a simple snowmelt model that can be used
        #        when there is insufficient data to use the energy
        #        balance method.  If c0 and T0 are established for a
        #        range of conditions, Kane et al. (1997) showed that
        #        this method gives comparable results to the energy
        #        balance method.

        #        86400000d = 1000 [mm/m] * 60 [sec/min] *
        #                    60 [min/hr] * 24 [hrs/day]

        #        rho_snow is not needed in formula here, but is
        #        needed to convert snowmelt to snow water
        #        equivalent (swe) in update_swe() in "snow_base.py".
        #-------------------------------------------------------------
        T_air = self.T_air  # (2/3/13, new framework)
        #-----------------------------------------------
##        T_air = self.get_port_data('T_air', self.mp)
        
        #------------------------------
        # Compute degree-day meltrate
        #------------------------------
        M = (self.c0 / np.float64(8.64E7)) * (T_air - self.T0)   #[m/s]

        # This is really an "enforce_min_meltrate()"
        self.SM = np.maximum(M, np.float64(0))
   
        #-------------------------------------------------------
        # Note: enforce_max_meltrate() method is always called
        #       by the base class to make sure that meltrate
        #       does not exceed the max possible.
        #-------------------------------------------------------
        #  self.enforce_max_meltrate()
        
    #   update_meltrate()
    #-------------------------------------------------------------------
    
