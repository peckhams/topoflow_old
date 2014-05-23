
## Copyright (c) 2001-2013, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## July, August 2009
## May 2010 (changes to unit_test() and read_cfg_file()

#-----------------------------------------------------------------------
#  NOTES:  This file defines an Energy-Balance snowmelt component
#          and related functions.  It inherits from the snowmelt
#          "base class" in "snow_base.py".
#-----------------------------------------------------------------------
#
#  class snow_energy_balance
#
#      get_attribute()          # (10/26/11)
#      get_input_var_names()    # (5/14/12)
#      get_output_var_names()   # (5/14/12)
#      get_var_name()           # (5/16/12, Bolton)
#      get_var_units()          # (5/16/12, Bolton)
#      ----------------------
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
        'model_name':         'TopoFlow_Snowmelt_Energy_Balance',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'SnowEnergyBalance',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Snow_Energy_Balance.cfg.in',
        'cfg_extension':      '_snow_energy_balance.cfg',
        'cmt_var_prefix':     '/SnowEnergyBalance/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Snow_Energy_Balance.xml',
        'dialog_title':       'Snowmelt: Energy Balance Parameters',  # (Snowmelt ?)
        'time_units':         'seconds' }

    #----------------------------------------------------------
    # The Energy-Balance method needs several vars from the
    # Met component, such as T_air, Q_sum.  The others
    # are used for consistent use of scalars vs. grids;
    # see "check_input_types()" below. (e.g. p0 and RH)
    #----------------------------------------------------------
    # Some input vars, like ****** are read from the CFG
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
        'land_surface_air__pressure',
        'land_surface__net_irradiation_flux',
        'land_surface__net_longwave_irradiation_flux',
        'land_surface__net_shortwave_irradiation_flux',
        'land_surface__temperature',        
        'water__density',
        'wind__speed_reference_height',
        'wind__reference_height_speed' ]
                          
    #------------------------------------------------------------
    # Note: The main output of the Energy-Balance method is SM.
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
        'land_surface_air__pressure': 'p0',
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
        'land_surface_air__pressure': 'mbar',
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
##                          self.pp.is_scalar('rate'),
##                          self.pp.is_scalar('duration'),
                          #--------------------------------
##                          self.mp.is_scalar('P'),
##                          self.mp.is_scalar('rho_H2O'),
##                          self.mp.is_scalar('rho_air'),
##                          self.mp.is_scalar('Cp_air'),
##                          self.mp.is_scalar('T_air'),
##                          self.mp.is_scalar('T_surf'),
##                          self.mp.is_scalar('RH'),
##                          self.mp.is_scalar('p0'),
##                          self.mp.is_scalar('uz'),
##                          self.mp.is_scalar('z'),
##                          self.mp.is_scalar('z0_air'),
##                          self.mp.is_scalar('Qn_SW'),
##                          self.mp.is_scalar('Qn_LW'),
                          #--------------------------------
                          self.is_scalar('P'),
                          self.is_scalar('rho_H2O'),
                          self.is_scalar('rho_air'),
                          self.is_scalar('Cp_air'),
                          self.is_scalar('T_air'),
                          self.is_scalar('T_surf'),
                          self.is_scalar('RH'),
                          self.is_scalar('p0'),
                          self.is_scalar('uz'),
                          self.is_scalar('z'),
                          self.is_scalar('z0_air'),
                          self.is_scalar('Qn_SW'),
                          self.is_scalar('Qn_LW'),
                          #--------------------------------
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Cp_snow'),
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe') ])


        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def update_meltrate(self):

        #------------------------------------------------------------
        # Notes: See notes in "met_base.py" for the method called
        #        "update_net_energy_flux()".

        #        This version uses "Q_sum" which is computed as a
        #        state variable for a meteorology component
        #        (e.g. met_base.py).

        #        Arguments are assumed to be references to a scalar
        #            or grid so that:

        #        M  = water equivalent of snowmelt [m/s]
        #        M_max  = max possible meltrate if all snow melts
        #        T_air  = air temperature [deg_C]

        #        Model must start when snow is isothermal. (CHECK)
        #        Cooling of the snowpack is not considered.

        #        86400000d = 1000 [mm/m] * 60 [sec/min] *
        #                    60 [min/sec] * 24 [hrs/day]

        #        rho_snow is not needed in formula here, but is
        #        needed to convert snowmelt to water equivalent?
        #-------------------------------------------------------------
     
        #----------------------------------
        # Compute energy-balance meltrate   
        #------------------------------------------------------
        # Ecc is initialized with the Initial_Cold_Content
        # function by Initialize_Snow_Vars function (2/21/07)
        #------------------------------------------------------
        # The following pseudocode only works for scalars but
        # is otherwise equivalent to that given below and
        # clarifies the logic:
        #------------------------------------------------------
        #  if (Q_sum gt 0) then begin
        #      if ((Q_sum * dt) gt Ecc) then begin
        #          ;-------------------------------------------
        #          ; Snow is melting.  Use some of Q_sum to
        #          ; overcome Ecc, and remainder to melt snow
        #          ;-------------------------------------------
        #          Qm  = Q_sum - (Ecc/dt)
        #          Ecc = 0
        #          M   = (Qm / (rho_w * Lf))
        #      endif else begin
        #          ;------------------------------
        #          ; Snow is warming; reduce Ecc
        #          ;------------------------------
        #          Ecc = (Ecc - (Q_sum * dt))
        #          M   = 0d
        #      endelse
        #  endif else begin
        #      ;--------------------------------
        #      ; Snow is cooling; increase Ecc
        #      ;--------------------------------
        #      Ecc = Ecc - (Q_sum * dt)
        #      M   = 0d
        #  endelse
        #----------------------------------------------------------
        # Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc     # [W/m^2]
        #----------------------------------------------------------
        Q_sum = self.Q_sum  # (2/3/13, new framework)
        #----------------------------------------------------------
##        Q_sum = self.get_port_data('Q_sum', self.mp)    # [W/m^2]
        Qcc   = (self.Ecc / self.dt)                    # [W/m^2]
        Qm    = np.maximum((Q_sum - Qcc), float64(0))                   # [W/m^2]
        Ecc   = np.maximum((self.Ecc - (Q_sum * self.dt)), np.float64(0))  # [J/m^2]
        
        #print 'Qm = ', Qm
        #print ' '
        
        #-------------------------------------
        # Convert melt energy to a melt rate
        #------------------------------------------
        # Lf = latent heat of fusion [J/kg]
        # Lv = latent heat of vaporization [J/kg]
        # M  = (Qm/ (rho_w * Lf))
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lf    = 334000d     ;[J/kg = W*s/kg]
        # So (rho_w * Lf) = 3.34e+8  [J/m^3]
        #------------------------------------------
        M       = (Qm / np.float64(3.34E+8))   #[m/s]
        self.SM = np.maximum(M, np.float64(0))
   
        #-------------------------------------------------------
        # Note: enforce_max_meltrate() method is always called
        #       by the base class to make sure that meltrate
        #       does not exceed the max possible.
        #-------------------------------------------------------
        #  self.enforce_max_meltrate()
            
    #   update_meltrate()
    #-------------------------------------------------------------------
    
#-------------------------------------------------------------------------  
###-------------------------------------------------------------------------
##def Energy_Balance_Meltrate(Qn_SW, Qn_LW, T_air, T_surf, RH, p0, \
##                            uz, z, z0_air, rho_air, Cp_air, Ecc, \
##                            h_snow, rho_snow, Cp_snow, dt, \
##                            e_air, e_surf):  #(returned)
##
##    #-----------------------------------------------------------------
##    #Notes:  3/13/07.  This function used to have vapor pressure
##    #        arguments e_air and e_surf.  However, e_air is always
##    #        computed as a function of T_air and RH and e_surf is
##    #        computed as a function of T_surf (saturated vap. press.)
##    #        So it makes sense to remove these two arguments and add
##    #        RH (relative humidity).  This change only affects the
##    #        Latent_Heat_Flux function call, which calls a new
##    #        function called Vapor_Pressure.
##    #-----------------------------------------------------------------
##    #        Qm    = energy used to melt snowpack (if > 0)
##    #        Qn_SW = net shortwave radiation flux (solar)
##    #        Qn_LW = net longwave radiation flux (air, surface)
##    #        Qh    = sensible heat flux from turbulent convection
##    #                between snow surface and air
##    #        Qe    = latent heat flux from evaporation, sublimation,
##    #                and condensation
##    #        Qa    = energy advected by moving water (i.e. rainfall)
##    #                (ARHYTHM assumes this to be negligible; Qa=0.)
##    #        Qc    = energy flux via conduction from snow to soil
##    #                (ARHYTHM assumes this to be negligible; Qc=0.)
##    #        Ecc   = cold content of snowpack = amount of energy
##    #                needed before snow can begin to melt [J/m^2]
##
##    #        All Q's here have units of [W/m^2].
##    #        Are they all treated as positive quantities ?
##
##    #        rho_air  = density of air [kg/m^3]
##    #        rho_snow = density of snow [kg/m^3]
##    #        Cp_air   = specific heat of air [J/kg/deg_C]
##    #        Cp_snow  = heat capacity of snow [J/kg/deg_C]
##    #                 = ???????? = specific heat of snow
##    #        Kh       = eddy diffusivity for heat [m^2/s]
##    #        Ke       = eddy diffusivity for water vapor [m^2/s]
##    #        Lv       = latent heat of vaporization [J/kg]
##    #        Lf       = latent heat of fusion [J/kg]
##    #        ------------------------------------------------------
##    #        Dn       = bulk exchange coeff for the conditions of
##    #                   neutral atmospheric stability [m/s]
##    #        Dh       = bulk exchange coeff for heat
##    #        De       = bulk exchange coeff for vapor
##    #        ------------------------------------------------------
##    #        T_air    = air temperature [deg_C]
##    #        T_surf   = surface temperature [deg_C]
##    #        T_snow   = average snow temperature [deg_C]
##    #        RH       = relative humidity [none] (in [0,1])
##    #        e_air    = air vapor pressure at height z [mbar]
##    #        e_surf   = surface vapor pressure [mbar]
##    #        ------------------------------------------------------
##    #        h_snow   = snow depth [m]
##    #        z        = height where wind speed is uz [m]
##    #        uz       = wind speed at height z [m/s]
##    #        P0       = atmospheric pressure [mbar]
##    #        T0       = snow temperature when isothermal [deg_C]
##    #                   (This is usually 0.)
##    #        z0_air   = surface roughness length scale [m]
##    #                   (includes vegetation not covered by snow)
##    #                   (Values from page 1033: 0.0013, 0.02 [m])
##    #        kappa    = von Karman's constant [unitless] = 0.41
##    #        dt       = snowmelt timestep [seconds]
##    #----------------------------------------------------------------
##   
##    # FORWARD_FUNCTION Richardson_Number
##    
##    #---------------------------------
##    #Some required physical constants
##    #are defined in the functions:
##    #e.g. Lv, Lf
##    #---------------------------------
##    
##    #------------------------------
##    #Compute the Richardson number
##    #------------------------------
##    Ri = Richardson_Number(z, uz, T_air, T_surf)
##    
##    #-------------------------------------------------
##    #Compute bulk exchange coeffs (neutral stability)
##    #-------------------------------------------------
##    Dn = Bulk_Exchange_Coeff(uz, z, h_snow, z0_air, T_air, T_surf)
##    Dh = Dn
##    De = Dn
##    
##    #---------------------------
##    #Compute sensible heat flux
##    #---------------------------
##    Qh = Sensible_Heat_Flux(rho_air, Cp_air, Dh, T_air, T_surf)
##    #Formula:  Qh = rho_air * Cp_air * Dh * (T_air - T_surf)
##    #print,'Dh = ', Dh
##    #print,'Qh = ', Qh
##    
##    #-------------------------
##    #Compute latent heat flux
##    #-------------------------
##    Qe = Latent_Heat_Flux(rho_air, De, T_air, T_surf, RH, p0,
##                          e_air, e_surf)  #(these 2 returned)
##    #Formula:  Qe = rho_air * Lv * De * (0.662/p0) * (e_air - e_surf)
##    
##    #print,'Qe = ', Qe
##    
##    #-----------------------------
##    #Compute conduction heat flux
##    #-----------------------------
##    Qc = Conduction_Heat_Flux()
##    #Formula:  Qc = 0d
##    
##    #-----------------------------
##    #Compute advective heat flux
##    #-----------------------------
##    Qa = Advection_Heat_Flux()
##    #Formula:  Qa = 0d
##    
##    #---------------------------------
##    #Qn_SW, Qn_SW & Ecc are pointers,
##    #others are local variables
##    #----------------------------------------------------
##    #Ecc is initialized with the Initial_Cold_Content
##    #function by Initialize_Snow_Vars function (2/21/07)
##    #----------------------------------------------------
##    #The following pseudocode only works for scalars but
##    #is otherwise equivalent to that given below and
##    #clarifies the logic:
##    #----------------------------------------------------
##    #  if (Q_sum gt 0) then begin
##    #      if ((Q_sum * dt) gt Ecc) then begin
##    #          ;-----------------------------------------
##    #          ;Snow is melting.  Use some of Q_sum to
##    #          ;overcome Ecc, and remainder to melt snow
##    #          ;-----------------------------------------
##    #          Qm  = Q_sum - (Ecc/dt)
##    #          Ecc = 0
##    #          M   = (Qm / (rho_w * Lf))
##    #      endif else begin
##    #          ;----------------------------
##    #          ;Snow is warming; reduce Ecc
##    #          ;----------------------------
##    #          Ecc = (Ecc - (Q_sum * dt))
##    #          M   = 0d
##    #      endelse
##    #  endif else begin
##    #      ;------------------------------
##    #      ;Snow is cooling; increase Ecc
##    #      ;------------------------------
##    #      Ecc = Ecc - (Q_sum * dt)
##    #      M   = 0d
##    #  endelse
##    #-------------------------------------------------------
##    Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc    #[W/m^2]
##    Qcc  = (Ecc / dt)                            #[W/m^2]
##    Qm   = maximum((Q_sum - Qcc), float64(0))                      #[W/m^2]
##    Ecc  = maximum((Ecc - (Q_sum * dt)), float64(0))              #[J/m^2]
##    #print,'Qm = ', Qm
##    #print,' '
##    
##    #-----------------------------------
##    #Convert melt energy to a melt rate
##    #----------------------------------------
##    #Lf = latent heat of fusion [J/kg]
##    #Lv = latent heat of vaporization [J/kg]
##    #M  = (Qm/ (rho_w * Lf))
##    #----------------------------------------
##    #rho_w = 1000d       ;[kg/m^3]
##    #Lf    = 334000d     ;[J/kg = W*s/kg]
##    #So (rho_w * Lf) = 3.34e+8  [J/m^3]
##    #-------------------------------------
##    M = (Qm / float32(3.34E+8))   #[m/s]
##    
##    return maximum(M, float32(0.0))
##
###   Energy_Balance_Meltrate
###-------------------------------------------------------------------------
