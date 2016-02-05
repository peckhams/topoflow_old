
## Does "land_surface_air__latent_heat_flux" make sense? (2/5/13)

# Copyright (c) 2001-2014, Scott D. Peckham
#
#  Sep 2014.  Fixed sign error in update_bulk_richardson_number().
#             Ability to compute separate P_snow and P_rain.
#  Aug 2014.  New CSDMS Standard Names and clean up.
#  Nov 2013.  Converted TopoFlow to a Python package.
#
#  Jan 2013. Revised handling of input/output names.
#  Oct 2012. CSDMS Standard Names (version 0.7.9) and BMI.
#  May 2012. P is now a 1D array with one element and mutable,
#            so any comp with ref to it can see it change.
#  Jun 2010. update_net_shortwave_radiation(), etc.  
#  May 2010. Changes to initialize() and read_cfg_file().
#  Aug 2009
#  Jan 2009. Converted from IDL.
#
#-----------------------------------------------------------------------
# NOTES: This file defines a "base class" for meteorology
#        components as well as any functions used by most or
#        all meteorology methods.  The methods of this class
#        should be over-ridden as necessary for different
#        methods of modeling meteorology.
#-----------------------------------------------------------------------
# Notes: Do we ever need to distinguish between a surface
#        temperature and snow temperature (in the snow) ?
#        Recall that a separate T_soil_x variable is used
#        to compute Qc.
#
#        Cp_snow is from NCAR CSM Flux Coupler web page
#
#        rho_H2O is currently not adjustable with GUI. (still true?)
#    
#-----------------------------------------------------------------------    
#
#  class met_component     (inherits from BMI_base.py)
#
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (5/15/12)
#      get_output_var_names()    # (5/15/12)
#      get_var_name()            # (5/15/12)
#      get_var_units()           # (5/15/12)
#      ---------------------
#      set_constants()
#      initialize()
#      update()
#      finalize()
#      ----------------------------
#      set_computed_input_vars()
#      initialize_computed_vars()
#      ----------------------------
#      update_P_integral()
#      update_P_max()
#      update_P_rain()    # (9/14/14, new method)
#      update_P_snow()    # (9/14/14, new method)
#      ------------------------------------
#      update_bulk_richardson_number()
#      update_bulk_aero_conductance()
#      update_sensible_heat_flux()
#      update_saturation_vapor_pressure()
#      update_vapor_pressure()
#      update_dew_point()                    # (7/6/10)
#      update_precipitable_water_content()   # (7/6/10)
#      ------------------------------------
#      update_latent_heat_flux()
#      update_conduction_heat_flux()
#      update_advection_heat_flux()
#      ------------------------------------
#      update_julian_day()                   # (7/1/10)
#      update_net_shortwave_radiation()      # (7/1/10)
#      update_em_air()                       # (7/1/10)
#      update_net_longwave_radiation()       # (7/1/10)
#      update_net_energy_flux()              # ("Q_sum")
#      ------------------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#  Functions:
#      compare_em_air_methods()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import solar_funcs as solar

from topoflow.utils import BMI_base
from topoflow.utils import model_input
from topoflow.utils import model_output
from topoflow.utils import rtg_files

#-----------------------------------------------------------------------
class met_component( BMI_base.BMI_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Meteorology',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'Meteorology',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Meteorology.cfg.in',
        'cfg_extension':      '_meteorology.cfg',
        'cmt_var_prefix':     '/Meteorology/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Meteorology.xml',
        'dialog_title':       'Meteorology: Method 1 Parameters',
        'time_units':         'seconds' }

    #---------------------------------------------------------
    # Note that SWE = "snow water equivalent", but it really
    # just means "liquid_equivalent".
    #---------------------------------------------------------   
    _input_var_names = [
        'snowpack__z_mean_of_mass-per-volume_density', # rho_snow
        'snowpack__depth',                             # h_snow
        'snowpack__liquid-equivalent_depth',           # h_swe
        'snowpack__melt_volume_flux' ]                 # SM   (MR used for ice?)

    #-----------------------------------------------------------
    # albedo, emissivity and transmittance are dimensionless.
    #-----------------------------------------------------------
    # "atmosphere_aerosol_dust__reduction_of_transmittance" vs.
    # This TF parameter comes from Dingman, App. E, p. 604.
    #-----------------------------------------------------------
    # There is an Optical_Air_Mass function in solar_funcs.py.
    # However, this quantity is not saved in comp state.
    #
    # "optical_path_length_ratio" vs. "optical_air_mass" OR
    # "airmass_factor" OR "relative_airmass" OR
    # "relative_optical_path_length"
    #-----------------------------------------------------------
    # Our term "liquid_equivalent_precipitation" is widely
    # used on the Internet, with 374,000 Google hits.
    #--------------------------------------------------------------
    # Note: "bulk exchange coefficient" has 2460 Google hits.
    #       It is closely related to a "transfer coefficient"
    #       for mass, momentum or heat.  There are no CF
    #       Standard Names with "bulk", "exchange" or "transfer".
    #
    # Zhang et al. (2000) use "bulk exchange coefficient" in a
    # nonstandard way, with units of velocity vs. unitless.
    #
    # Dn = bulk exchange coeff for the conditions of
    #      neutral atmospheric stability [m/s]
    # Dh = bulk exchange coeff for heat  [m/s]
    # De = bulk exchange coeff for vapor [m/s]
    #---------------------------------------------------------------
    # Now this component uses T_air to break the liquid-equivalent
    # precip rate into separate P_rain and P_snow components.
    # P_rain is used by channel_base.update_R()
    # P_snow is used by snow_base.update_depth()
    #---------------------------------------------------------------
    _output_var_names = [
        # 'atmosphere__optical_path_length_ratio',                           # M_opt [1]  (in solar_funcs.py)
        # 'atmosphere__von_karman_constant',                                 # kappa
        'atmosphere_aerosol_dust__reduction_of_transmittance',               # dust_atten  ##### (from GUI)
        'atmosphere_air-column_water-vapor__liquid-equivalent_depth',        # W_p ("precipitable depth")
        'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor',         # canopy_factor
        'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor',          # cloud_factor
        'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance',   # De [m s-1], latent
        'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance', # Dh [m s-1], sensible
        'atmosphere_bottom_air__emissivity',                                 # em_air
        'atmosphere_bottom_air__mass-per-volume_density',                    # rho_air
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity',       # Cp_air
        'atmosphere_bottom_air__neutral_bulk_aerodynamic_conductance',       # Dn [m s-1], neutral
        'atmosphere_bottom_air__pressure',                                   # p0
        'atmosphere_bottom_air__temperature',                                # T_air
        'atmosphere_bottom_air_flow__bulk_richardson_number',                # Ri [1]
        'atmosphere_bottom_air_flow__log_law_roughness_length',              # z0_air
        'atmosphere_bottom_air_flow__reference-height_speed',                # uz
        'atmosphere_bottom_air_flow__speed_reference_height',                # z
        'atmosphere_bottom_air_land_net-latent-heat__energy_flux',           # Qe [W m-2]
        'atmosphere_bottom_air_land_net-sensible-heat__energy_flux',         # Qh [W m-2]
        'atmosphere_bottom_air_water-vapor__dew_point_temperature',          # T_dew
        'atmosphere_bottom_air_water-vapor__partial_pressure',               # e_air # (insert "reference_height" ??)
        'atmosphere_bottom_air_water-vapor__relative_saturation',            # RH
        'atmosphere_bottom_air_water-vapor__saturated_partial_pressure',     # e_sat_air         
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux',  # vol_P
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux',     # P_max
        'atmosphere_water__precipitation_leq-volume_flux',                        # P [m s-1]
        'atmosphere_water__rainfall_volume_flux',            # P_rain [m s-1] (liquid)      
        'atmosphere_water__snowfall_leq-volume_flux',        # P_snow [m s-1]
        'earth__standard_gravity_constant',                  # g   [m s-2]
        'land_surface__albedo',                              # albedo
        'land_surface__aspect_angle',                        # alpha  (from GUI)
        'land_surface__emissivity',                          # em_surf
        'land_surface__latitude',                            # lat_deg [degrees]
        'land_surface__longitude',                           # lon_deg [degrees]
        'land_surface__slope_angle',                         # beta  (from GUI)
        'land_surface__temperature',                         # T_surf   ### OR JUST "land__temperature"?
        # 'land_surface_air__temperature',                          # T_air
        'land_surface_air_water-vapor__partial_pressure',           # e_surf # (insert "reference_height" ??)
        'land_surface_air_water-vapor__saturated_partial_pressure', # e_sat_surf
        'land_surface_net-longwave-radiation__energy_flux',  # Qn_LW [W m-2]
        'land_surface_net-shortwave-radiation__energy_flux', # Qn_SW [W m-2]
        'land_surface_net-total-energy__energy_flux',        # Q_sum [W w-2]
        'model__time_step',                                  # dt
        'physics__stefan_boltzmann_constant',                # sigma  [W m-2 K-4]
        'physics__von_karman_constant',                      # kappa  [1]
        'water__mass-specific_latent_fusion_heat',           # Lf     [J kg-1]
        'water__mass-specific_latent_vaporization_heat',     # Lv     [J kg-1]
        'water-liquid__mass-per-volume_density' ]            # rho_H2O
        
    #-----------------------------------------
    # These are used only in solar_funcs.py
    # Later, create a Radiation component.
    #---------------------------------------------
    # Should we allow "day" as a base quantity ?
    # "day_length" is confusing.  Think about "date" also.
    # Maybe something like:
    #
    #    "earth__mean_solar_rotation_period"
    #    "earth__sidereal_rotation_period"
    #    "earth__stellar_rotation_period"   (relative to "fixed stars")
    #         maybe:  "earth__complete_rotation_period" ??
    #
    #    OR:
    #    "earth_mean_solar_day__duration"
    #    "earth_sidereal_day__duration"
    #    "earth_stellar_day__duration"
    #
    #    OR perhaps:
    #    "earth_mean_solar_day__rotation_period"
    #    "earth_sidereal_day__rotation_period"
    #    "earth_stellar_day__rotation_period"
    #
    #    "stellar rotation period" gives 84,500 Google hits.
    #    "solar_rotation_period" gives 41,100 Google hits.
    #    "sidereal_roation_period" gives 86,000 Google hits.
    #    "stellar day" gives 136,000 Google hits (but many unrelated).
    #
    #    NB! "stellar_rotation_period" is ambiguous since it is also
    #         used for the rotation period of a star.
    #
    #    "earth_mean_solar_day__hour_count"  ("standard_day" ?)
    #    "earth_sidereal_day__hour_count"
    #    "earth_sidereal_day__duration"
    #    "earth__rotation_period"   = "sidereal_day"
    #
    #    "earth_stellar_day__period"  ??
    #    "earth_stellar_day__duration" ??
    #
    #------------------------------------------------------------------
    # For "earth__rotation_rate", it seems this should be based on
    # the sidereal day (23.93 hours) instead of the mean solar day.
    #------------------------------------------------------------------
    # There are at least a few online sources that use both terms:
    # "equivalent latitude" and "equivalent longitude".  See:
    # "The Construction and Application of a Martian Snowpack Model".
    #------------------------------------------------------------------
    # Adopt the little-used term:  "topographic_sunrise" ?
    # Or maybe "illuminated_topography", or "local_sunrise" ??
    #------------------------------------------------------------------
    # For angle relations between the earth and the sun, should we
    # just use the adjective "solar" in the quantity name or include
    # sun in the object name?  We could also use terms like:
    #     earth_to_sun__declination_angle
    #     earth_to_sun__right_ascension_angle
    #
    #------------------------------------------------------------------
    # The adjective "local" in "earth_local_apparent_noon__time"
    # may be helpful in other contexts such as:
    # 'earth__local_longitude' and 'land_surface__local_elevation'.
    #------------------------------------------------------------------    
    # 'earth__autumnal_equinox_date',
    # 'earth__autumnal_equinox_time',
    # 'earth_axis__ecliptic_tilt_angle',  # tilt_angle
    # 'earth__julian_day_number',         ########
    # 'earth__julian_day_angle',
    # 'earth__local_apparent_noon_time'
    # 'earth__mean_radius', 
    # 'earth__mean_solar_day_duration',   # (exactly 24 hours)
    # 'earth_orbit__eccentricity',
    # 'earth_orbit__period',     # (one year)
    # 'earth__perihelion_julian_day',  ######
    # 'earth__rotation_period',        ######
    # 'earth__rotation_rate', # Omega       ###### What about Angular Velocity ?
    # 'earth__sidereal_day_duration',     # (one rotation = 23.934470 hours)
    # 'earth__solar_declination_angle',
    # 'earth__solar_hour_angle',
    # 'earth__solar_irradiation_constant',  ## (or "insolation_constant" ??)
    # 'earth__solar_right_ascension_angle',
    # 'earth__solar_vertical_angle',   (complement of zenith angle)
    # 'earth__solar_zenith_angle',
    # 'earth__stellar_day_duration',  # (relative to the "fixed stars")
    # 'earth__summer_solstice_date',
    # 'earth__summer_solstice_time',
    # 'earth__topographic_sunrise_equivalent_latitude',
    # 'earth__topographic_sunrise_equivalent_longitude',  (flat_lon + offset) 
    # 'earth__topographic_sunrise_equivalent_longitude_offset',
    # 'earth__topographic_sunrise_time',
    # 'earth__topographic_sunset_time',
    # 'earth_true_solar_noon___time',  #####
    #     'earth_clock__true_solar_noon_time'
    # 'earth__vernal_equinox_date',
    # 'earth__vernal_equinox_time',
    # 'earth__winter_solstice_date',
    # 'earth__winter_solstice_time',
    #
    # What about a "slope_corrected" or "topographic" version of K_dir ?
    #
    # 'land_surface__backscattered_shortwave_irradiation_flux', # K_bs
    # 'land_surface__diffuse_shortwave_irradiation_flux',       # K_dif
    # 'land_surface__direct_shortwave_irradiation_flux',        # K_dir
    # 'land_surface__global_shortwave_irradiation_flux',        # K_glob = K_dif + K_dir
    #------------------------------------------------------------------


    #------------------------------------------------------------------   
    # Maybe we should rename "z" to "z_ref" and "uz" to "uz_ref" ?
    #------------------------------------------------------------------   
    _var_name_map = {
        'snowpack__z_mean_of_mass-per-volume_density': 'rho_snow',
        'snowpack__depth': 'h_snow',
        'snowpack__liquid-equivalent_depth': 'h_swe',
        'snowpack__melt_volume_flux': 'SM',              # (MR is used for ice)
        #-----------------------------------------------------------------
        #'atmosphere__optical_path_length_ratio': 'M_opt',    # (in solar_funcs.py)
        # 'atmosphere__von_karman_constant': 'kappa',
        'atmosphere_aerosol_dust__reduction_of_transmittance': 'dust_atten',
        'atmosphere_air-column_water-vapor__liquid-equivalent_depth': 'W_p',   #########
        'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor': 'canopy_factor',
        'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor':  'cloud_factor',
        'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance':   'De',
        'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance': 'Dh',
        'atmosphere_bottom_air__emissivity': 'em_air',               
        'atmosphere_bottom_air__mass-per-volume_density': 'rho_air',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'Cp_air',
        'atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance':  'Dn',
        'atmosphere_bottom_air__pressure':                           'p0',
        'atmosphere_bottom_air__temperature':                        'T_air',
        'atmosphere_bottom_air_flow__bulk_richardson_number':        'Ri',        
        'atmosphere_bottom_air_flow__log_law_roughness_length':      'z0_air', ## (not "z0")
        'atmosphere_bottom_air_flow__reference-height_speed':        'uz',
        'atmosphere_bottom_air_flow__speed_reference_height':        'z',
        'atmosphere_bottom_air_land_net-latent-heat__energy_flux':   'Qe',
        'atmosphere_bottom_air_land_net-sensible-heat__energy_flux': 'Qh',
        'atmosphere_bottom_air_water-vapor__dew_point_temperature':  'T_dew',  
        'atmosphere_bottom_air_water-vapor__partial_pressure':       'e_air',
        'atmosphere_bottom_air_water-vapor__relative_saturation':    'RH',
        'atmosphere_bottom_air_water-vapor__saturated_partial_pressure': 'e_sat_air',
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux': 'vol_P', 
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux': 'P_max',       
        'atmosphere_water__precipitation_leq-volume_flux': 'P',
        'atmosphere_water__rainfall_volume_flux':          'P_rain',    
        'atmosphere_water__snowfall_leq-volume_flux':      'P_snow',
        'earth__standard_gravity_constant':                'g',
        'land_surface__albedo':                            'albedo',
        'land_surface__aspect_angle':                      'alpha',
        'land_surface__emissivity':                        'em_surf',
        'land_surface__latitude':                          'lat_deg',
        'land_surface__longitude':                         'lon_deg',
        'land_surface__slope_angle':                       'beta',
        'land_surface__temperature':                       'T_surf',
         # 'land_surface_air__temperature': 'T_surf',
        'land_surface_air_water-vapor__partial_pressure':           'e_surf',
        'land_surface_air_water-vapor__saturated_partial_pressure': 'e_sat_surf',        
        'land_surface_net-longwave-radiation__energy_flux':         'Qn_LW',
        'land_surface_net-shortwave-radiation__energy_flux':        'Qn_SW',
        'land_surface_net-total-energy__energy_flux':               'Q_sum',
        'model__time_step':                                         'dt',
        'physics__stefan_boltzmann_constant':                       'sigma',
        'physics__von_karman_constant':                             'kappa',
        'water__mass-specific_latent_fusion_heat':                  'Lf',
        'water__mass-specific_latent_vaporization_heat':            'Lv',
        'water-liquid__mass-per-volume_density': 'rho_H2O' }

    #-----------------------------------------------------------------
    # Note: The "update()" function calls several functions with the
    #       MBAR keyword set to get units of "mbar" vs. "kPa".
    #-----------------------------------------------------------------
    # Note: We need to be careful with whether units are C or K,
    #       for all "thermal" quantities (e.g. thermal_capacity).
    #-----------------------------------------------------------------
    # Note: ARHYTHM had 3 "bulk exchange coefficients" that are all
    #       equal and therefore have the same units of [m s-1].
    #       Double-check that this is what is intended.   ##########
    #-----------------------------------------------------------------
    # Note: "atmosphere_column_water__liquid_equivalent_depth" has
    #       units of "cm", as in Dingman's book.  Make sure it gets
    #       used correctly in equations.
    #-----------------------------------------------------------------
    # Note: slope_angle and aspect_angle have units of RADIANS.
    #       aspect_angle is measured CW from north.
    #       RT files ending in "_mf-angle.rtg" and "fd-aspect.rtg"
    #       contain aspect values.  The former are in [0, 2 Pi]
    #       while the latter are in [-Pi, Pi] and both measure
    #       CCW from due east.  They are converted for use here.
    #-----------------------------------------------------------------    
    _var_units_map = {
        'snowpack__z_mean_of_mass-per-volume_density':     'kg m-3',
        'snowpack__depth':                   'm',
        'snowpack__liquid-equivalent_depth': 'm',
        'snowpack__melt_volume_flux':        'm s-1',
        #-------------------------------------------------------------
        # 'atmosphere__optical_path_length_ratio': '1',
        # 'atmosphere__von_karman_constant': '1',
        'atmosphere_aerosol_dust__reduction_of_transmittance':        '1',
        'atmosphere_air-column_water-vapor__liquid-equivalent_depth': 'cm',           # (see Notes above)
        'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor':  '1',
        'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor':   '1',
        'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance': 'm s-1',   # (see Notes above)
        'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance': 'm s-1', # (see Notes above)  
        'atmosphere_bottom_air__emissivity': '1',                            
        'atmosphere_bottom_air__mass-per-volume_density': 'kg m-3',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1',   # (see Notes above)
        'atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance': 'm s-1',   # (see Notes above)
        'atmosphere_bottom_air__pressure':                               'mbar',
        'atmosphere_bottom_air__temperature':                            'deg_C',      # (see Notes above)
        'atmosphere_bottom_air_flow__bulk_richardson_number':            '1',
        'atmosphere_bottom_air_flow__log_law_roughness_length':          'm',
        'atmosphere_bottom_air_flow__reference-height_speed':            'm s-1',
        'atmosphere_bottom_air_flow__speed_reference_height':            'm',
        'atmosphere_bottom_air_land_net-latent-heat__energy_flux':       'W m-2',
        'atmosphere_bottom_air_land_net-sensible-heat__energy_flux':     'W m-2',
        'atmosphere_bottom_air_water-vapor__dew_point_temperature':      'deg_C',
        'atmosphere_bottom_air_water-vapor__partial_pressure':           'mbar', # (see Notes above)
        'atmosphere_bottom_air_water-vapor__relative_saturation':        '1',
        'atmosphere_bottom_air_water-vapor__saturated_partial_pressure': 'mbar',     # (see Notes above)
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux': 'm3',
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux': 'm s-1',        
        'atmosphere_water__precipitation_leq-volume_flux': 'm s-1',
        'atmosphere_water__rainfall_volume_flux': 'm s-1',      # (see Notes above)  
        'atmosphere_water__snowfall_leq-volume_flux': 'm s-1',  # (see Notes above)
        'earth__standard_gravity_constant': 'm s-2',
        'land_surface__albedo': '1',
        'land_surface__aspect_angle': 'radians',                      # (see Notes above)
        'land_surface__emissivity': '1',
        'land_surface__latitude': 'degrees',
        'land_surface__longitude': 'degrees',
        'land_surface__slope_angle': 'radians',
        'land_surface__temperature': 'deg_C',
        # 'land_surface_air__temperature': 'deg_C',
        'land_surface_air_water-vapor__partial_pressure':           'mbar',
        'land_surface_air_water-vapor__saturated_partial_pressure': 'mbar',
        'land_surface_net-longwave-radiation__energy_flux': 'W m-2',
        'land_surface_net-shortwave-radiation__energy_flux': 'W m-2',
        'land_surface_net-total-energy__energy_flux': 'W m-2',
        'model__time_step': 's',
        'physics__stefan_boltzmann_constant': 'W m-2 K-4',
        'physics__von_karman_constant': '1',
        'water__mass-specific_latent_fusion_heat': 'J kg-1',
        'water__mass-specific_latent_vaporization_heat': 'J kg-1',
        'water-liquid__mass-per-volume_density': 'kg m-3' }

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
##        return 'float64'
##    
##    #   get_var_type()
    #-------------------------------------------------------------------
    def set_constants(self):

        #---------------------------------
        # Define some physical constants
        #---------------------------------
        self.g        = np.float64(9.81)    # [m s-2, gravity]
        self.kappa    = np.float64(0.408)   # [1]  (von Karman)
        self.rho_H2O  = np.float64(1000)    # [kg m-3]
        self.rho_air  = np.float64(1.2614)  # [kg m-3]
        self.Cp_air   = np.float64(1005.7)  # [J kg-1 K-1]
        self.Lv       = np.float64(2500000) # [J kg-1] Latent heat of vaporiz.
        self.Lf       = np.float64(334000)  # [J kg-1 = W s kg-1], Latent heat of fusion
        self.sigma    = np.float64(5.67E-8) # [W m-2 K-4]  (Stefan-Boltzman constant)
        self.C_to_K   = np.float64(273.15)  # (add to convert deg C to K)

        self.twopi         = np.float64(2) * np.pi
        self.one_seventh   = np.float64(1) / 7
        self.hours_per_day = np.float64(24)
        self.secs_per_day  = np.float64(3600) * self.hours_per_day

        #---------------------------
        # See update_latent_heat()
        #-----------------------------------------------------------        
        # According to Dingman (2002, p. 273), constant should
        # be 0.622 instead of 0.662 (Zhang et al., 2000, p. 1002).
        # Is this constant actually the dimensionless ratio of
        # the molecular weight of water to that of dry air ?
        #-----------------------------------------------------------
        ## self.latent_heat_constant = np.float64(0.622)
        self.latent_heat_constant = np.float64(0.662)
        
        #----------------------------------------
        # Constants related to precip (9/24/09)
        #----------------------------------------
        self.mmph_to_mps = (np.float64(1) / np.float64(3600000))
        self.mps_to_mmph = np.float64(3600000)
        self.forever     = np.float64(999999999)  # [minutes]
        
        #------------------------------------------------
        # Only needed for method 1, where all rates and
        # durations are read as 1D arrays from GUI.
        # Method 1 may be removed in a future version.
        #------------------------------------------------
##        self.method1_rates     = None
##        self.method1_durations = None
##        self.method1_n_rates   = 0
        
    #   set_constants()             
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print ' '
            print 'Meteorology component: Initializing...'
            
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
                
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()     
        ## print '    Calling read_grid_info()...'
        self.read_grid_info()
        ## print '    Calling initialize_basin_vars()...'
        self.initialize_basin_vars()  # (5/14/10)
        #----------------------------------------------------
        # NB! This read_input_files() uses self.time_index.
        #     Also needs to be before "Disabled" test.
        #----------------------------------------------------
        ## print '    Calling initialize_time_vars()...'
        self.initialize_time_vars()

        #-------------------------------------------------
        # (5/19/12) This makes P "mutable", which allows
        # its updated values to be seen by any component
        # that has a reference to it.
        #-------------------------------------------------
        # Write a "initialize_computed_vars()" method?
        #-------------------------------------------------
        self.P      = self.initialize_scalar(0, dtype='float64')
        self.P_rain = self.initialize_scalar(0, dtype='float64')
        self.P_snow = self.initialize_scalar(0, dtype='float64')
                    
        #------------------------------------------------------
        # NB! "Sample steps" must be defined before we return
        #     Check all other process modules.
        #------------------------------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print 'Meteorology component: Disabled.'
            ## self.P        = np.float64(0)
            self.e_air    = self.initialize_scalar(0, dtype='float64')
            self.e_surf   = self.initialize_scalar(0, dtype='float64')
            self.em_air   = self.initialize_scalar(0, dtype='float64')
            self.Qn_SW    = self.initialize_scalar(0, dtype='float64')
            self.Qn_LW    = self.initialize_scalar(0, dtype='float64')
            self.Q_sum    = self.initialize_scalar(0, dtype='float64')
            self.Qc       = self.initialize_scalar(0, dtype='float64')
            self.Qa       = self.initialize_scalar(0, dtype='float64')
            self.DONE     = True
            self.status   = 'initialized'
            return

        #-----------------------------------------------
        # Read from files as needed to initialize vars 
        #-----------------------------------------------
        self.open_input_files()
        self.read_input_files()  # (initializes P)

        # Some output variables aren't defined until update() is called.
        # Initialize them here, instead. (@mdpiper, 9/8/15)
        try:
            self.Ri
        except AttributeError:
            self.Ri = np.float64(0.0)
        try:
            self.h_snow
        except AttributeError:
            self.h_snow = np.float64(0.0)
        
        ## self.check_input_types()  # (not needed so far)
        
        #-----------------------
        # Initialize variables
        #-----------------------
        ## print '    Calling initialize_computed_vars()...'
        self.initialize_computed_vars() # (after read_input_files)
        
        if not(self.PRECIP_ONLY):
            self.open_output_files() 
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
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
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #-------------------------------------------
        # Update computed values related to precip
        #-------------------------------------------
        self.update_P_integral()
        self.update_P_max()
        self.update_P_rain()
        self.update_P_snow()

        #-------------------------
        # Update computed values
        #-------------------------
        if not(self.PRECIP_ONLY):
            self.update_bulk_richardson_number()
            self.update_bulk_aero_conductance()
            self.update_sensible_heat_flux()
            self.update_saturation_vapor_pressure(MBAR=True)
            self.update_saturation_vapor_pressure(MBAR=True, SURFACE=True)  ########
            self.update_vapor_pressure(MBAR=True)
            self.update_dew_point() ###
            self.update_precipitable_water_content()  ###
            self.update_vapor_pressure(MBAR=True, SURFACE=True)   ########
            self.update_latent_heat_flux()      # (uses e_air and e_surf)
            self.update_conduction_heat_flux()
            self.update_advection_heat_flux()
            self.update_julian_day()
            self.update_net_shortwave_radiation()
            self.update_em_air()
            self.update_net_longwave_radiation()
            self.update_net_energy_flux()  # (at the end)
            
        #----------------------------------------
        # Read next met vars from input files ?
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
        if not(self.PRECIP_ONLY):
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
        if (self.comp_status == 'Enabled'):
            self.close_input_files()   ##  TopoFlow input "data streams"
            if not(self.PRECIP_ONLY):
                self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Meteorology component')
        
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

        #-----------------------------------------------
        # Convert precip rate units from mm/h to m/s ?
        #-----------------------------------------------
        # NB! read_input_files() does this for files.
        #-----------------------------------------------
        if (self.P_type == 'Scalar'):
##            print '######## self.P_type  =', self.P_type
##            print '######## type(self.P) =', type(self.P)
##            print '######## self.P       =', self.P
##            print '######## Converting scalar P from MMPH to MPS.'
            #-----------------------------------------------------
            # (2/7/13) Must use "*=" here to preserve reference.
            #-----------------------------------------------------
            self.P *= self.mmph_to_mps
            ## self.P = self.P * self.mmph_to_mps
            print 'Scalar rainrate set to:', self.P, ' [mmph]'
            
        #---------------------------------
        # Process the PRECIP_ONLY toggle
        #---------------------------------
        if not(hasattr(self, 'PRECIP_ONLY')):
            self.PRECIP_ONLY = False
        elif (self.PRECIP_ONLY.lower() == 'yes'):
            self.PRECIP_ONLY = True
        else:
            self.PRECIP_ONLY = False

        #---------------------------------------
        # Print info message about PRECIP_ONLY
        #---------------------------------------
        if (self.PRECIP_ONLY):
            print '-----------------------------------------'
            print ' NOTE: Since PRECIP_ONLY = True, output'
            print '       variables will not be computed'
            print '       or saved to files.'
            print '-----------------------------------------'
            print' '
            
        #----------------------------------------------------
        # Toggle to use SATTERLUND or BRUTSAERT methods
        # for computing e_air and em_air. (Not in GUI yet.)
        #----------------------------------------------------
        if not(hasattr(self, 'SATTERLUND')):
            self.SATTERLUND = False

        #---------------------------------------------
        # Convert GMT_offset from string to int
        # because GUI can't use ints in droplist yet
        #---------------------------------------------
        self.GMT_offset = np.int16( self.GMT_offset )

        #------------------------------------------------
        # Convert start_month from string to integer
        # January should be 1.  See solar.Julian_Day().
        #------------------------------------------------
        month_list = ['January', 'February', 'March', 'April',
                      'May', 'June', 'July', 'August', 'September',
                      'October', 'November', 'December']
        self.start_month = month_list.index( self.start_month ) + 1
                 
        #-------------------------------
        # Initialize some more toggles
        #-------------------------------                   
        if not(hasattr(self, 'SAVE_QSW_GRIDS')):
            self.SAVE_QSW_GRIDS = False
        if not(hasattr(self, 'SAVE_QLW_GRIDS')):
            self.SAVE_QLW_GRIDS = False
        #-------------------------------------------
        if not(hasattr(self, 'SAVE_QSW_PIXELS')):
            self.SAVE_QSW_PIXELS = False
        if not(hasattr(self, 'SAVE_QLW_PIXELS')):
            self.SAVE_QLW_PIXELS = False
   
        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = np.maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt = np.maximum(self.save_pixels_dt,  self.dt)
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #------------------------------------------------------
        # Note: Some of these require "self.rti", which is
        #       only stored by read_grid_info() after the
        #       set_computed_input_vars() function is called.
        #       So these parts can't go there.
        #------------------------------------------------------
        
        #---------------------------------------
        # Add self.in_directory to:
        #   slope_grid_file & aspect_grid_file
        #---------------------------------------
        self.slope_grid_file  = (self.in_directory + self.slope_grid_file) 
        self.aspect_grid_file = (self.in_directory + self.aspect_grid_file)
        
        #-------------------------------------------------
        # Read slope grid & convert to slope angle, beta
        # NB!  RT slope grids have NaNs on edges.
        #-------------------------------------------------
        slopes = rtg_files.read_grid( self.slope_grid_file, self.rti,
                                      RTG_type='FLOAT' )
        beta   = np.arctan( slopes )
        beta   = (self.twopi + beta) % self.twopi
        #---------------------------------------------
        w_nan = np.where( np.logical_not(np.isfinite(beta)) )
        n_nan = np.size(w_nan[0])
        if (n_nan != 0):    
            beta[ w_nan ] = np.float64(0)
        #------------------------------------------------------------------
        w_bad = np.where( np.logical_or( (beta < 0), (beta > np.pi / 2) ) )
        n_bad = np.size(w_bad[0])
        if (n_bad != 0):    
            msg = array(['ERROR:  Some slope angles are out of range.', ' '])
            for line in msg:
                print line
            ## result = GUI_Message(msg, INFO=True, TITLE='ERROR MESSAGE')
            return
        self.beta = beta  ######
        
        #------------------------------------------------------
        # Read aspect grid.  Alpha must be CW from north.
        # NB!  RT aspect grids have NaNs on edges.
        #---------------------------------------------------------
        # RT files ending in "_mf-angle.rtg" and "fd-aspect.rtg"
        # contain aspect values.  The former are in [0, 2 Pi]
        # while the latter are in [-Pi, Pi] and both measure
        # CCW from due east.
        #---------------------------------------------------------
        aspects = rtg_files.read_grid( self.aspect_grid_file, self.rti,
                                       RTG_type='FLOAT' )
        alpha   = (np.pi / 2) - aspects
        alpha   = (self.twopi + alpha) % self.twopi
        #-----------------------------------------------
        w_nan = np.where( np.logical_not( np.isfinite(alpha) ) )
        n_nan = np.size( w_nan[0] )
        if (n_nan != 0):    
            alpha[ w_nan ] = np.float64(0)
        self.alpha = alpha  ######

        #---------------------------        
        # Create lon and lat grids
        #---------------------------
        if (self.rti.pixel_geom == 0):
            self.lon_deg = solar.Longitude_Grid( self.rti )
            self.lat_deg = solar.Latitude_Grid( self.rti )

##            print 'Lon grid ='
##            print self.lon_deg
##            print 'Lat grid ='
##            print self.lat_deg
            
            #-----------------------------
            # Write grids to RTG files ?
            #-----------------------------
##            lon_file = (self.out_directory + self.site_prefix + '_lons.bin')
##            rtg_files.write_grid( self.lon_deg, lon_file, self.rti )
##            lat_file = (self.out_directory + self.site_prefix + '_lats.bin')
##            rtg_files.write_grid( self.lat_deg, lat_file, self.rti )
        else:    
            print 'SORRY: Cannot yet create lon and lat grids for'
            print '       this DEM because it uses UTM coordinates.'
            print '       Will use lat/lon for Denver, Colorado.'
            print ' '
            #--------------------------------------------
            # For now, use scalar values for Denver, CO
            #--------------------------------------------
            self.lon_deg = np.float64( -104.9841667 )
            self.lat_deg = np.float64( 39.7391667 )
            ## return

        #-------------------------------------------------
        # Initialize max precip rate with the first rate
        #------------------------------------------------
        # Note: Need this here because rate may be
        #       zero at the end of update_precip_rate()
        #------------------------------------------------
        # vol_P is used for mass balance check.
        #------------------------------------------------
        P_max = self.P.max()       # (after read_input_files)
        ## self.P_max = self.P.max()
        self.P_max = self.initialize_scalar( P_max, dtype='float64')
        self.vol_P = self.initialize_scalar( 0, dtype='float64')


        #----------------------------------------------------------
        # For using new framework which embeds references from
        # meteorology to snow, etc., these need to be defined
        # in the initialize step.  However, they will most likely
        # change from scalar to grid during update, so we need to
        # check that the reference isn't broken when the dtype
        # changes. (5/17/12)
        #---------------------------------------------------------- 
        # These depend on grids alpha and beta, so will be grids.
        #----------------------------------------------------------                
        self.Qn_SW  = np.zeros([self.ny, self.nx], dtype='float64')
        self.Qn_LW  = np.zeros([self.ny, self.nx], dtype='float64')
        self.Qn_tot = np.zeros([self.ny, self.nx], dtype='float64')
        self.Q_sum  = np.zeros([self.ny, self.nx], dtype='float64')
        #----------------------------------------------------------  
#         self.Qn_SW  = self.initialize_scalar( 0, dtype='float64')
#         self.Qn_LW  = self.initialize_scalar( 0, dtype='float64')
#         self.Qn_tot = self.initialize_scalar( 0, dtype='float64')
#         self.Q_sum  = self.initialize_scalar( 0, dtype='float64')
        #---------------------------------------------------------- 
        # These may be scalars or grids.
        #---------------------------------         
        self.Qe     = self.initialize_scalar( 0, dtype='float64')
        self.e_air  = self.initialize_scalar( 0, dtype='float64')
        self.e_surf = self.initialize_scalar( 0, dtype='float64')
        self.em_air = self.initialize_scalar( 0, dtype='float64')
        self.Qc     = self.initialize_scalar( 0, dtype='float64')
        self.Qa     = self.initialize_scalar( 0, dtype='float64')
         
        #------------------------------------
        # Initialize the decimal Julian day
        #------------------------------------
        self.julian_day = solar.Julian_Day( self.start_month,
                                            self.start_day,
                                            self.start_hour )
        ## print '    julian_day =', self.julian_day
        
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_P_integral(self):

        #---------------------------------------------------
        # Notes: This can be used for mass balance checks,
        #        such as now done by update_mass_totals()
        #        in topoflow.py.  The "dt" here should be
        #        TopoFlow's "main dt" vs. the process dt.
        
        #        dV[i] = P[i] * da[i] * dt, dV = sum(dV[i])
        #---------------------------------------------------
        if (self.DEBUG):
            print 'Calling update_P_integral()...'
            
        #------------------------------------------------
        # Update mass total for P, sum over all pixels
        #------------------------------------------------   
        volume = np.double(self.P * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_P += (volume * self.rti.n_pixels)
        else:
            self.vol_P += np.sum(volume)
        
    #   update_P_integral()
    #-------------------------------------------------------------------
    def update_P_max(self):

        if (self.DEBUG):
            print 'Calling update_P_max()...'
            
        #-----------------------------------------
        # Save the maximum precip. rate in [m/s]
        #-------------------------------------------
        # Must use "fill()" to preserve reference.
        #-------------------------------------------
        self.P_max.fill( np.maximum(self.P_max, self.P.max()) )
        ## self.P_max = np.maximum(self.P_max, self.P.max())

        ### print '##### P =', self.P
        
    #   update_P_max()  
    #-------------------------------------------------------------------
    def update_P_rain(self):

        #-----------------------------------------------------------
        # Note:  This routine is written so that it doesn't matter
        #        whether P and T_air are grids or scalars.
        #        For scalars: 1.5 * True = 1.5, 1.5 * False = 0.
        #        Here are the possible combinations for checking.
        #-----------------------------------------------------------
        # P       T_air     P_rain
        #----------------------------
        # scalar  scalar    scalar    
        # scalar  grid      grid
        # grid    scalar    grid
        # grid    grid      grid
        #----------------------------
        if (self.DEBUG):
            print 'Calling update_P_rain()...'

        #-------------------------------------------------
        # P_rain is the precip that falls as liquid that
        # can contribute to runoff production.
        #-------------------------------------------------
        # P_rain is used by channel_base.update_R.
        #-------------------------------------------------
        P_rain = self.P * (self.T_air > 0)
        if (np.rank( self.P_rain ) == 0):
            self.P_rain.fill( P_rain )   #### (mutable scalar)
        else:
            self.P_rain[:] = P_rain
  
        if (self.DEBUG):
            if (self.P_rain.max() > 0):
                print '   >> Rain is falling...'

        #--------------
        # For testing
        #--------------
        ## print 'shape(P)      =', shape(self.P)
        ## print 'shape(T_air)  =', shape(self.T_air)
        ## print 'shape(P_rain) =', shape(self.P_rain)
        ## print 'T_air         =', self.T_air

        #########################################
        #### Old note, to remember for later.
        #--------------------------------------------------
        # (2/7/13) We must use "*=" to preserve reference
        # if P is a "mutable scalar".
        #--------------------------------------------------
                             
    #   update_P_rain()
    #-------------------------------------------------------------------
    def update_P_snow(self):

        #----------------------------------------------------
        # Notes:  Rain and snow may fall simultaneously at
        #         different grid cells in the model domain.
        #----------------------------------------------------
        if (self.DEBUG):
            print 'Calling update_P_snow()...'
            
        #-------------------------------------------------
        # P_snow is the precip that falls as snow or ice
        # that contributes to the snow depth.  This snow
        # may melt to contribute to runoff later on.
        #-------------------------------------------------
        # P_snow is used by snow_base.update_depth.
        #-------------------------------------------------
        P_snow = self.P * (self.T_air <= 0)
        if (np.rank( self.P_snow ) == 0):
            self.P_snow.fill( P_snow )   #### (mutable scalar)
        else:
            self.P_snow[:] = P_snow

        if (self.DEBUG):
            if (self.P_snow.max() > 0):
                print '   >> Snow is falling...'
                     
    #   update_P_snow()
    #-------------------------------------------------------------------
    def update_bulk_richardson_number(self):

        if (self.DEBUG):
            print 'Calling update_bulk_richardson_number()...'

        #---------------------------------------------------------------
        # (9/6/14)  Found a typo in the Zhang et al. (2000) paper,
        # in the definition of Ri.  Also see Price and Dunne (1976).
        # We should have (Ri > 0) and (T_surf > T_air) when STABLE.
        # This also removes problems/singularities in the corrections
        # for the stable and unstable cases in the next function.      
        #---------------------------------------------------------------
        # Notes: Other definitions are possible, such as the one given
        #        by Dingman (2002, p. 599).  However, this one is the
        #        one given by Zhang et al. (2000) and is meant for use
        #        with the stability criterion also given there.
        #---------------------------------------------------------------
        #### top     = self.g * self.z * (self.T_air - self.T_surf)  # BUG.
        top     = self.g * self.z * (self.T_surf - self.T_air)
        bot     = (self.uz)**2.0 * (self.T_air + np.float64(273.15))
        self.Ri = (top / bot)

    #   update_bulk_richardson_number()
    #-------------------------------------------------------------------        
    def update_bulk_aero_conductance(self):

        if (self.DEBUG):
            print 'Calling update_bulk_aero_conductance()...'

        #----------------------------------------------------------------
        # Notes: Dn       = bulk exchange coeff for the conditions of
        #                   neutral atmospheric stability [m/s]
        #        Dh       = bulk exchange coeff for heat  [m/s]
        #        De       = bulk exchange coeff for vapor [m/s]
        #        h_snow   = snow depth [m]
        #        z0_air   = surface roughness length scale [m]
        #                   (includes vegetation not covered by snow)
        #        z        = height that has wind speed uz [m]
        #        uz       = wind speed at height z [m/s]
        #        kappa    = 0.408 = von Karman's constant [unitless]
        #        RI       = Richardson's number (see function)
        #----------------------------------------------------------------
        h_snow = self.h_snow  # (ref from new framework)
        
        #---------------------------------------------------
        # Compute bulk exchange coeffs (neutral stability)
        # using the logarithm "law of the wall".
        #-----------------------------------------------------
        # Note that "arg" = the drag coefficient (unitless).
        #-----------------------------------------------------        
        arg = self.kappa / np.log((self.z - h_snow) / self.z0_air)
        Dn  = self.uz * (arg)**2.0
        
        #-----------------------------------------------
        # NB! Dn could be a scalar or a grid, so this
        #     must be written to handle both cases.
        #     Note that WHERE can be used on a scalar:
        
        #     IDL> a = 1
        #     IDL> print, size(a)
        #     IDL> w = where(a ge 1, nw)
        #     IDL> print, nw
        #     IDL> a[w] = 2
        #     IDL> print, a
        #     IDL> print, size(a)
        #-----------------------------------------------

        ###########################################################
        #  NB!  If T_air and T_surf are both scalars, then next
        #       few lines won't work because we can't index the
        #       resulting empty "w" (even if T_air == T_surf).
        ###########################################################
##        w  = np.where(self.T_air != self.T_surf)
##        nw = np.size(w[0])
##        ## nw = np.size(w,0)  # (doesn't work if 2 equal scalars)
        #----------------------------------------------------------
        T_AIR_SCALAR  = (np.rank( self.T_air )  == 0)
        T_SURF_SCALAR = (np.rank( self.T_surf ) == 0)
        if (T_AIR_SCALAR and T_SURF_SCALAR):
            if (self.T_air == self.T_surf):  nw=1
            else: nw=0      
        else:
            w  = np.where(self.T_air != self.T_surf)
            nw = np.size(w[0])
        
        if (nw == 0):
            #--------------------------------------------
            # All pixels are neutral. Set Dh = De = Dn.
            #--------------------------------------------
            self.Dn = Dn
            self.Dh = Dn
            self.De = Dn
            return
        
        #-------------------------------------
        # One or more pixels are not neutral
        # so make a correction using RI
        #---------------------------------------------
        # NB!  RI could be a grid when Dn is a
        # scalar, and this will change Dn to a grid.
        #---------------------------------------------
        # Ri = Richardson_Number(z, uz, T_air, T_surf)
        #--------------------------------------------
        # Before 12/21/07.  Has bug if RI is a grid
        #--------------------------------------------
        # w_stable = where(*T_air gt *T_surf, n_stable)
        # if (n_stable ne 0) then begin
        #     Dn[w_stable] = Dn[w_stable]/(1d + (10d * RI))
        # endif
        # w_unstable = where(*T_air lt *T_surf, n_unstable)
        # if (n_unstable ne 0) then begin
        #----------------------------------------------
        # Multiplication and substraction vs. opposites
        # for the stable case.  Zhang et al. (2000)
        # Hopefully not just a typo.
        #----------------------------------------------
        #    Dn[w_unstable] = Dn[w_unstable]*(1d - (10d * self.Ri))
        # endif
        
        #-----------------
        # After 12/21/07
        #------------------------------------------------------------
        # If T_air, T_surf or uz is a grid, then Ri will be a grid.
        # This version makes only one call to WHERE, so its faster.
        #------------------------------------------------------------
        # Multiplication and substraction vs. opposites for the
        # stable case (Zhang et al., 2000); hopefully not a typo.
        # It plots as a smooth curve through Ri=0.
        #------------------------------------------------------------
        # (9/7/14)  Modified so that Dn is saved, but Dh = De.
        #------------------------------------------------------------        
        Dh = Dn.copy()   ### (9/7/14.  Save Dn also.)
        nD = np.size( Dh )
        nR = np.size( self.Ri )
        if (nR > 1):    
            #--------------------------
            # Case where RI is a grid
            #--------------------------
            ws = np.where( self.Ri > 0 )
            ns = np.size( ws[0] )
            wu = np.where( np.invert(self.Ri > 0) )
            nu = np.size( wu[0] )
            if (nD == 1):    
                #******************************************
                # Convert Dn to a grid here or somewhere
                # Should stop with an error message
                #******************************************
                dum = np.int16(0)
            if (ns != 0):
                #----------------------------------------------------------
                # If (Ri > 0), or (T_surf > T_air), then STABLE. (9/6/14)
                #----------------------------------------------------------   
                Dh[ws] = Dh[ws] / (np.float64(1) + (np.float64(10) * self.Ri[ws]))
            if (nu != 0):    
                Dh[wu] = Dh[wu] * (np.float64(1) - (np.float64(10) * self.Ri[wu]))
        else:    
            #----------------------------
            # Case where Ri is a scalar
            #--------------------------------
            # Works if Dh is grid or scalar
            #--------------------------------
            if (self.Ri > 0):    
                Dh = Dh / (np.float64(1) + (np.float64(10) * self.Ri))
            else:    
                Dh = Dh * (np.float64(1) - (np.float64(10) * self.Ri))

        #----------------------------------------------------
        # NB! We currently assume that these are all equal.
        #----------------------------------------------------
        self.Dn = Dn
        self.Dh = Dh
        self.De = Dh   ## (assumed equal)
        
    #   update_bulk_aero_conductance()
    #-------------------------------------------------------------------
    def update_sensible_heat_flux(self):

        #--------------------------------------------------------
        # Notes: All the Q's have units of W/m^2 = J/(m^2 s).
        #        Dh is returned by Bulk_Exchange_Coeff function
        #        and is not a pointer.
        #--------------------------------------------------------
        if (self.DEBUG):
            print 'Callilng update_sensible_heat_flux()...'
            
        #---------------------
        # Physical constants
        #---------------------
        # rho_air = 1.225d   ;[kg m-3, at sea-level]
        # Cp_air  = 1005.7   ;[J kg-1 K-1]
        
        #-----------------------------
        # Compute sensible heat flux
        #-----------------------------
        delta_T = (self.T_air - self.T_surf)
        self.Qh = (self.rho_air * self.Cp_air) * self.Dh * delta_T

    #   update_sensible_heat_flux()
    #-------------------------------------------------------------------
    def update_saturation_vapor_pressure(self, MBAR=False,
                                         SURFACE=False):

        if (self.DEBUG):
            print 'Calling update_saturation_vapor_pressure()...'

        #----------------------------------------------------------------
        #Notes:  Saturation vapor pressure is a function of temperature.
        #        T is temperature in Celsius.  By default, the method
        #        of Brutsaert (1975) is used.  However, the SATTERLUND
        #        keyword is set then the method of Satterlund (1979) is
        #        used.  When plotted, they look almost identical.  See
        #        the Compare_em_air_Method routine in Qnet_file.pro.
        #        Dingman (2002) uses the Brutsaert method.
        #        Liston (1995, EnBal) uses the Satterlund method.

        #        By default, the result is returned with units of kPa.
        #        Set the MBAR keyword for units of millibars.
        #        100 kPa = 1 bar = 1000 mbars
        #                => 1 kPa = 10 mbars
        #----------------------------------------------------------------
        #NB!     Here, 237.3 is correct, and not a misprint of 273.2.
        #        See footnote on p. 586 in Dingman (Appendix D).
        #----------------------------------------------------------------
        if (SURFACE):
##            if (self.T_surf_type in ['Scalar', 'Grid']):
##                return
            T = self.T_surf
        else:
##            if (self.T_air_type in ['Scalar', 'Grid']):
##                return
            T = self.T_air
        
        if not(self.SATTERLUND):    
            #------------------------------
            # Use Brutsaert (1975) method
            #------------------------------
            term1 = (np.float64(17.3) * T) / (T + np.float64(237.3))
            e_sat = np.float64(0.611) * np.exp(term1)        # [kPa]
        else:    
            #-------------------------------
            # Use Satterlund (1979) method     ############ DOUBLE CHECK THIS (7/26/13)
            #-------------------------------
            term1 = np.float64(2353) / (T + np.float64(273.15))
            e_sat = np.float64(10) ** (np.float64(11.4) - term1)   # [Pa]
            e_sat = (e_sat / np.float64(1000))   # [kPa]

        #-----------------------------------
        # Convert units from kPa to mbars?
        #-----------------------------------
        if (MBAR):    
            e_sat = (e_sat * np.float64(10))   # [mbar]

        if (SURFACE):
            self.e_sat_surf = e_sat
        else:
            self.e_sat_air  = e_sat

    #   update_saturation_vapor_pressure()
    #------------------------------------------------------------------- 
    def update_vapor_pressure(self, MBAR=False,
                              SURFACE=False):

        if (self.DEBUG):
            print 'Calling update_vapor_pressure()...'

        #---------------------------------------------------
        # Notes: T is temperature in Celsius
        #        RH = relative humidity, in [0,1]
        #             by definition, it equals (e / e_sat)
        #        e has units of kPa.
        #---------------------------------------------------
        if (SURFACE):
##            if (self.T_surf_type in ['Scalar', 'Grid']) and \
##               (self.RH_type in ['Scalar', 'Grid']):
##                return
            e_sat = self.e_sat_surf
        else:
##            if (self.T_air_type in ['Scalar', 'Grid']) and \
##               (self.RH_type in ['Scalar', 'Grid']):
##                return
            e_sat = self.e_sat_air
            
        e = (self.RH * e_sat)
        
        #-----------------------------------
        # Convert units from kPa to mbars?
        #-----------------------------------
        if (MBAR):    
            e = (e * np.float64(10))   # [mbar]

        if (SURFACE):
            self.e_surf = e
        else:
            self.e_air  = e
 
    #   update_vapor_pressure()
    #-------------------------------------------------------------------
    def update_dew_point(self):

        if (self.DEBUG):
            print 'Calling update_dew_point()...'

        #-----------------------------------------------------------
        # Notes:  The dew point is a temperature in degrees C and
        #         is a function of the vapor pressure, e_air.
        #         Vapor pressure is a function of air temperature,
        #         T_air, and relative humidity, RH.
        #         The formula used here needs e_air in kPa units.
        #         See Dingman (2002, Appendix D, p. 587).
        #-----------------------------------------------------------
        e_air_kPa = self.e_air / np.float64(10)   # [kPa]
        log_vp    = np.log( e_air_kPa )

        top = log_vp + np.float64(0.4926)
        bot = np.float64(0.0708) - (np.float64(0.00421) * log_vp)

        self.T_dew = (top / bot)    # [degrees C]
    
    #   update_dew_point()
    #-------------------------------------------------------------------
    def update_precipitable_water_content(self):

        if (self.DEBUG):
            print 'Calling update_precipitable_water_content()...'

        #------------------------------------------------------------
        # Notes:  W_p is precipitable water content in centimeters,
        #         which depends on air temp and relative humidity.
        #------------------------------------------------------------
        arg      = np.float64( 0.0614 * self.T_dew )
        self.W_p = np.float64(1.12) * np.exp( arg )  # [cm]

    #   update_precipitable_water_content()
    #-------------------------------------------------------------------
    def update_latent_heat_flux(self):

        if (self.DEBUG):
            print 'Calling update_latent_heat_flux()...'

        #--------------------------------------------------------
        # Notes:  Pressure units cancel out because e_air and
        #         e_surf (in numer) have same units (mbar) as
        #         p0 (in denom).
        #--------------------------------------------------------        
        # According to Dingman (2002, p. 273), constant should
        # be 0.622 instead of 0.662 (Zhang et al., 2000).
        #--------------------------------------------------------
        const   = self.latent_heat_constant
        factor  = (self.rho_air * self.Lv * self.De)
        delta_e = (self.e_air - self.e_surf)
        self.Qe = factor * delta_e * (const / self.p0)

    #   update_latent_heat_flux()
    #-------------------------------------------------------------------
    def update_conduction_heat_flux(self):

        if (self.DEBUG):
            print 'Calling update_conduction_heat_flux()...'

        #-----------------------------------------------------------------
        # Notes: The conduction heat flux from snow to soil for computing
        #        snowmelt energy, Qm, is close to zero.

        #        However, the conduction heat flux from surface and sub-
        #        surface for computing Qet is given by Fourier's Law,
        #        namely Qc = Ks(Tx - Ts)/x.

        #        All the Q's have units of W/m^2 = J/(m^2 s).
        #-----------------------------------------------------------------
        pass  # (initialized at start)

    #   update_conduction_heat_flux()
    #-------------------------------------------------------------------
    def update_advection_heat_flux(self):

        if (self.DEBUG):
            print 'Calling update_advection_heat_flux()...'

        #------------------------------------------------------
        # Notes: All the Q's have units of W/m^2 = J/(m^2 s).
        #------------------------------------------------------
        pass  # (initialized at start)
        
    #   update_advection_heat_flux()
    #-------------------------------------------------------------------
    def update_julian_day(self):

        if (self.DEBUG):
            print 'Calling update_julian_day()...'

        #----------------------------------
        # Update the *decimal* Julian day
        #----------------------------------
        self.julian_day += (self.dt / self.secs_per_day) # [days]
  
        #------------------------------------------
        # Compute the offset from True Solar Noon
        # clock_hour is in 24-hour military time
        # but it can have a decimal part.
        #------------------------------------------
        dec_part   = self.julian_day - np.int16(self.julian_day)
        clock_hour = dec_part * self.hours_per_day
        ## print '    Computing solar_noon...'
        solar_noon = solar.True_Solar_Noon( self.julian_day,
                                            self.lon_deg,
                                            self.GMT_offset )
        ## print '    Computing TSN_offset...'
        self.TSN_offset = (clock_hour - solar_noon)    # [hours]
    
    #   update_julian_day()
    #-------------------------------------------------------------------
    def update_net_shortwave_radiation(self):

        #---------------------------------------------------------
        # Notes:  If time is before local sunrise or after local
        #         sunset then Qn_SW should be zero.
        #---------------------------------------------------------
        if (self.DEBUG):
            print 'Calling update_net_shortwave_radiation()...'

        #--------------------------------
        # Compute Qn_SW for this time
        #--------------------------------
        Qn_SW = solar.Clear_Sky_Radiation( self.lat_deg,
                                           self.julian_day,
                                           self.W_p,
                                           self.TSN_offset,
                                           self.alpha,
                                           self.beta,
                                           self.albedo,
                                           self.dust_atten )

        if (np.rank( self.Qn_SW ) == 0):
            self.Qn_SW.fill( Qn_SW )   #### (mutable scalar)
        else:
            self.Qn_SW[:] = Qn_SW  # [W m-2]
        
    #   update_net_shortwave_radiation()
    #-------------------------------------------------------------------
    def update_em_air(self):

        if (self.DEBUG):
            print 'Calling update_em_air()...'

        #---------------------------------------------------------
        # NB!  The Brutsaert and Satterlund formulas for air
        #      emissivity as a function of air temperature are in
        #      close agreement; see compare_em_air_methods().
        #      However, we must pay close attention to whether
        #      equations require units of kPa, Pa, or mbar.
        #
        #             100 kPa = 1 bar = 1000 mbars
        #                => 1 kPa = 10 mbars
        #---------------------------------------------------------
        # NB!  Temperatures are assumed to be given with units
        #      of degrees Celsius and are converted to Kelvin
        #      wherever necessary by adding C_to_K = 273.15.
        #
        #      RH = relative humidity [unitless]
        #---------------------------------------------------------
        # NB!  I'm not sure about how F is added at end because
        #      of how the equation is printed in Dingman (2002).
        #      But it reduces to other formulas as it should.
        #---------------------------------------------------------
        T_air_K = self.T_air + self.C_to_K
        
        if not(self.SATTERLUND):
            #-----------------------------------------------------
            # Brutsaert (1975) method for computing emissivity
            # of the air, em_air.  This formula uses e_air with
            # units of kPa. (From Dingman (2002, p. 196).)
            # See notes for update_vapor_pressure().
            #-----------------------------------------------------
            e_air_kPa = self.e_air / np.float64(10)  # [kPa]
            F       = self.canopy_factor
            C       = self.cloud_factor
            term1   = (1.0 - F) * 1.72 * (e_air_kPa / T_air_K) ** self.one_seventh
            term2   = (1.0 + (0.22 * C ** 2.0))
            self.em_air  = (term1 * term2) + F
        else:
            #--------------------------------------------------------
            # Satterlund (1979) method for computing the emissivity
            # of the air, em_air, that is intended to "correct
            # apparent deficiencies in this formulation at air
            # temperatures below 0 degrees C" (see G. Liston)
            # Liston cites Aase and Idso(1978), Satterlund (1979)
            #--------------------------------------------------------
            e_air_mbar = self.e_air
            eterm  = np.exp(-1 * (e_air_mbar)**(T_air_K / 2016) )
            self.em_air = 1.08 * (1.0 - eterm)
 
        #--------------------------------------------------------------  
        # Can't do this yet.  em_air is always initialized scalar now
        # but may change to grid on assignment. (9/23/14)
        #--------------------------------------------------------------
#         if (np.rank( self.em_air ) == 0):
#             self.em_air.fill( em_air )   #### (mutable scalar)
#         else:
#             self.em_air[:] = em_air
           
    #   update_em_air()    
    #-------------------------------------------------------------------
    def update_net_longwave_radiation(self):

        #----------------------------------------------------------------
        # Notes: Net longwave radiation is computed using the
        #        Stefan-Boltzman law.  All four data types
        #        should be allowed (scalar, time series, grid or
        #        grid stack).
        #
        #        Qn_LW = (LW_in - LW_out)
        #        LW_in   = em_air  * sigma * (T_air  + 273.15)^4
        #        LW_out  = em_surf * sigma * (T_surf + 273.15)^4
        #
        #        Temperatures in [deg_C] must be converted to
        #        [K].  Recall that absolute zero occurs at
        #        0 [deg_K] or -273.15 [deg_C].
        #
        #----------------------------------------------------------------
        # First, e_air is computed as:
        #   e_air = RH * 0.611 * exp[(17.3 * T_air) / (T_air + 237.3)]
        # Then, em_air is computed as:
        #   em_air = (1 - F) * 1.72 * [e_air / (T_air + 273.15)]^(1/7) *
        #             (1 + 0.22 * C^2) + F
        #----------------------------------------------------------------
        if (self.DEBUG):
            print 'Calling update_net_longwave_radiation()...'

        #--------------------------------
        # Compute Qn_LW for this time
        #--------------------------------
        T_air_K  = self.T_air  + self.C_to_K
        T_surf_K = self.T_surf + self.C_to_K
        LW_in    = self.em_air  * self.sigma * (T_air_K)** 4.0
        LW_out   = self.em_surf * self.sigma * (T_surf_K)** 4.0
        LW_out   = LW_out + ((1.0 - self.em_surf) * LW_in)
               
        self.Qn_LW = (LW_in - LW_out)   # [W m-2]

        #--------------------------------------------------------------  
        # Can't do this yet.  Qn_LW is always initialized grid now
        # but will often be created above as a scalar. (9/23/14)
        #--------------------------------------------------------------
#         if (np.rank( self.Qn_LW ) == 0):
#             self.Qn_LW.fill( Qn_LW )   #### (mutable scalar)
#         else:
#             self.Qn_LW[:] = Qn_LW  # [W m-2]
        
    #   update_net_longwave_radiation()
    #-------------------------------------------------------------------
    def update_net_total_radiation(self):

        #-----------------------------------------------
        # Notes: Added this on 9/11/14.  Not used yet.
        #------------------------------------------------------------
        #        Qn_SW = net shortwave radiation flux (solar)
        #        Qn_LW = net longwave radiation flux (air, surface)
        #------------------------------------------------------------       
        if (self.DEBUG):
            print 'Calling update_net_total_radiation()...'
            
        Qn_tot = self.Qn_SW + self.Qn_LW   # [W m-2]

        if (np.rank( self.Qn_tot ) == 0):
            self.Qn_tot.fill( Qn_tot )   #### (mutable scalar)
        else:
            self.Qn_tot[:] = Qn_tot  # [W m-2]
                       
    #   update_net_total_radiation()
    #-------------------------------------------------------------------
    def update_net_energy_flux(self):

        if (self.DEBUG):
            print 'Calling update_net_energy_flux()...'

        #------------------------------------------------------
        # Notes: Q_sum is used by "snow_energy_balance.py".
        #------------------------------------------------------
        #        Qm    = energy used to melt snowpack (if > 0)
        #        Qn_SW = net shortwave radiation flux (solar)
        #        Qn_LW = net longwave radiation flux (air, surface)
        #        Qh    = sensible heat flux from turbulent convection
        #                between snow surface and air
        #        Qe    = latent heat flux from evaporation, sublimation,
        #                and condensation
        #        Qa    = energy advected by moving water (i.e. rainfall)
        #                (ARHYTHM assumes this to be negligible; Qa=0.)
        #        Qc    = energy flux via conduction from snow to soil
        #                (ARHYTHM assumes this to be negligible; Qc=0.)
        #        Ecc   = cold content of snowpack = amount of energy
        #                needed before snow can begin to melt [J m-2]

        #        All Q's here have units of [W m-2].
        #        Are they all treated as positive quantities ?

        #        rho_air  = density of air [kg m-3]
        #        rho_snow = density of snow [kg m-3]
        #        Cp_air   = specific heat of air [J kg-1 K-1]
        #        Cp_snow  = heat capacity of snow [J kg-1 K-1]
        #                 = ???????? = specific heat of snow
        #        Kh       = eddy diffusivity for heat [m2 s-1]
        #        Ke       = eddy diffusivity for water vapor [m2 s-1]
        #        Lv       = latent heat of vaporization [J kg-1]
        #        Lf       = latent heat of fusion [J kg-1]
        #        ------------------------------------------------------
        #        Dn       = bulk exchange coeff for the conditions of
        #                   neutral atmospheric stability [m/s]
        #        Dh       = bulk exchange coeff for heat
        #        De       = bulk exchange coeff for vapor
        #        ------------------------------------------------------
        #        T_air    = air temperature [deg_C]
        #        T_surf   = surface temperature [deg_C]
        #        T_snow   = average snow temperature [deg_C]
        #        RH       = relative humidity [unitless] (in [0,1])
        #        e_air    = air vapor pressure at height z [mbar]
        #        e_surf   = surface vapor pressure [mbar]
        #        ------------------------------------------------------
        #        h_snow   = snow depth [m]
        #        z        = height where wind speed is uz [m]
        #        uz       = wind speed at height z [m/s]
        #        p0       = atmospheric pressure [mbar]
        #        T0       = snow temperature when isothermal [deg_C]
        #                   (This is usually 0.)
        #        z0_air   = surface roughness length scale [m]
        #                   (includes vegetation not covered by snow)
        #                   (Values from page 1033: 0.0013, 0.02 [m])
        #        kappa    = von Karman's constant [unitless] = 0.41
        #        dt       = snowmelt timestep [seconds]
        #----------------------------------------------------------------
        Q_sum = self.Qn_SW + self.Qn_LW + self.Qh + \
                self.Qe + self.Qa + self.Qc    # [W m-2]

        if (np.rank( self.Q_sum) == 0):
            self.Q_sum.fill( Q_sum )   #### (mutable scalar)
        else:
            self.Q_sum[:] = Q_sum  # [W m-2]
            
    #   update_net_energy_flux()   
    #-------------------------------------------------------------------  
    def open_input_files(self):

        if (self.DEBUG):
            print 'Calling open_input_files()...'
       
        self.P_file      = self.in_directory + self.P_file
        self.T_air_file  = self.in_directory + self.T_air_file
        self.T_surf_file = self.in_directory + self.T_surf_file
        self.RH_file     = self.in_directory + self.RH_file
        self.p0_file     = self.in_directory + self.p0_file
        self.uz_file     = self.in_directory + self.uz_file
        self.z_file      = self.in_directory + self.z_file
        self.z0_air_file = self.in_directory + self.z0_air_file

        self.albedo_file        = self.in_directory + self.albedo_file
        self.em_surf_file       = self.in_directory + self.em_surf_file
        self.dust_atten_file    = self.in_directory + self.dust_atten_file
        self.cloud_factor_file  = self.in_directory + self.cloud_factor_file
        self.canopy_factor_file = self.in_directory + self.canopy_factor_file

        self.P_unit      = model_input.open_file(self.P_type,      self.P_file)
        self.T_air_unit  = model_input.open_file(self.T_air_type,  self.T_air_file)
        self.T_surf_unit = model_input.open_file(self.T_surf_type, self.T_surf_file)
        self.RH_unit     = model_input.open_file(self.RH_type,     self.RH_file)
        self.p0_unit     = model_input.open_file(self.p0_type,     self.p0_file)
        self.uz_unit     = model_input.open_file(self.uz_type,     self.uz_file)
        self.z_unit      = model_input.open_file(self.z_type,      self.z_file)
        self.z0_air_unit = model_input.open_file(self.z0_air_type, self.z0_air_file)
               
        #-----------------------------------------------
        # These are needed to compute Qn_SW and Qn_LW.
        #-----------------------------------------------
        self.albedo_unit        = model_input.open_file(self.albedo_type,
                                                        self.albedo_file)
        self.em_surf_unit       = model_input.open_file(self.em_surf_type,
                                                        self.em_surf_file)
        self.dust_atten_unit    = model_input.open_file(self.dust_atten_type,
                                                        self.dust_atten_file)
        self.cloud_factor_unit  = model_input.open_file(self.cloud_factor_type,
                                                        self.cloud_factor_file)
        self.canopy_factor_unit = model_input.open_file(self.canopy_factor_type,
                                                        self.canopy_factor_file)
        #----------------------------------------------------------------------------
        # Note: GMT_offset plus slope and aspect grids will be read separately.
        #----------------------------------------------------------------------------

##        self.Qn_SW_unit  = model_input.open_file(self.Qn_SW_type,  self.Qn_SW_file)
##        self.Qn_LW_unit  = model_input.open_file(self.Qn_LW_type,  self.Qn_LW_file)
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        if (self.DEBUG):
            print 'Calling read_input_files()...'

        rti = self.rti

        #--------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #--------------------------------------------------------
        # NB! read_next() returns None if TYPE arg is "Scalar".
        #--------------------------------------------------------
        P = model_input.read_next(self.P_unit, self.P_type, rti,
                                  factor=self.mmph_to_mps)

        if (P != None):
            ## print 'MET: (time,P) =', self.time, P
            
            ## if (self.P_type.lower() != 'scalar'):
            if (np.rank( self.P ) == 0):
                self.P.fill( P )  #### (2/7/13, mutable scalar)
            else:
                self.P = P

            if (self.DEBUG or (self.time_index == 0)):
                print 'In read_input_files():'
                print '   min(P) =', P.min() * self.mps_to_mmph, ' [mmph]'
                print '   max(P) =', P.max() * self.mps_to_mmph, ' [mmph]'
                print ' '
        else:
            #-----------------------------------------------
            # Either self.P_type is "Scalar" or we've read
            # all of the data in the rain_rates file.
            #-----------------------------------------------
            if (self.P_type.lower() != 'scalar'):
                #------------------------------------
                # Precip is unique in this respect.
                #--------------------------------------------------
                # 2/7/13. Note that we don't change P from grid
                # to scalar since that could cause trouble for
                # other comps that use P, so we just zero it out.
                #--------------------------------------------------                #
                self.P.fill( 0 )
                if (self.DEBUG):
                    print 'Reached end of file:', self.P_file
                    print '  P set to 0 by read_input_files().'
            elif (self.time_sec >= self.dt):
                self.P.fill( 0 )
                if (self.DEBUG):
                    print 'Reached end of scalar rainfall duration.'
                    print '  P set to 0 by read_input_files().'
                    ## print 'time_sec =', self.time_sec
                    ## print 'met dt   =', self.dt
                    
##        print '######### In met_base.read_input_files() #######'
##        print 'self.P_type =', self.P_type
##        print 'self.P      =', self.P

        ###############################################################
        # If any of these are scalars (read from a time series file)
        # then we'll need to use "fill()" method to prevent breaking
        # the reference to the "mutable scalar". (2/7/13)
        ###############################################################
        T_air = model_input.read_next(self.T_air_unit, self.T_air_type, rti)
        if (T_air != None): self.T_air = T_air

        T_surf = model_input.read_next(self.T_surf_unit, self.T_surf_type, rti)
        if (T_surf != None): self.T_surf = T_surf

        RH = model_input.read_next(self.RH_unit, self.RH_type, rti)
        if (RH != None): self.RH = RH

        p0 = model_input.read_next(self.p0_unit, self.p0_type, rti)
        if (p0 != None): self.p0 = p0

        uz = model_input.read_next(self.uz_unit, self.uz_type, rti)
        if (uz != None): self.uz = uz

        z = model_input.read_next(self.z_unit, self.z_type, rti)
        if (z != None): self.z = z

        z0_air = model_input.read_next(self.z0_air_unit, self.z0_air_type, rti)
        if (z0_air != None): self.z0_air = z0_air

        #----------------------------------------------------------------------------
        # These are needed to compute Qn_SW and Qn_LW.
        #----------------------------------------------------------------------------
        # Note: We could later write a version of read_next() that takes "self"
        #       and "var_name" as args and that uses "exec()".
        #----------------------------------------------------------------------------
        albedo = model_input.read_next(self.albedo_unit, self.albedo_type, rti)
        if (albedo != None): self.albedo = albedo

        em_surf = model_input.read_next(self.em_surf_unit, self.em_surf_type, rti)
        if (em_surf != None): self.em_surf = em_surf

        dust_atten = model_input.read_next(self.dust_atten_unit, self.dust_atten_type, rti)
        if (dust_atten != None): self.dust_atten = dust_atten

        cloud_factor = model_input.read_next(self.cloud_factor_unit, self.cloud_factor_type, rti)
        if (cloud_factor != None): self.cloud_factor = cloud_factor

        canopy_factor = model_input.read_next(self.canopy_factor_unit, self.canopy_factor_type, rti)
        if (canopy_factor != None): self.canopy_factor = canopy_factor

        #-------------------------------------------------------------
        # Compute Qsw_prefactor from cloud_factor and canopy factor.
        #-------------------------------------------------------------
        ## self.Qsw_prefactor = 
        
        #-------------------------------------------------------------
        # These are currently treated as input data, but are usually
        # generated by functions in Qnet_file.py.  Later on, we'll
        # provide the option to compute them "on the fly" with new
        # functions called "update_net_shortwave_radiation()" and
        # "update_net_longwave_radiation()", called from update().
        #-------------------------------------------------------------        
##        Qn_SW = model_input.read_next(self.Qn_SW_unit, self.Qn_SW_type, rti)
##        if (Qn_SW != None): self.Qn_SW = Qn_SW
##
##        Qn_LW = model_input.read_next(self.Qn_LW_unit, self.Qn_LW_type, rti)
##        if (Qn_LW != None): self.Qn_LW = Qn_LW
         
    #   read_input_files()
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.DEBUG):
            print 'Calling close_input_files()...'

        if (self.P_type      != 'Scalar'): self.P_unit.close()
        if (self.T_air_type  != 'Scalar'): self.T_air_unit.close()
        if (self.T_surf_type != 'Scalar'): self.T_surf_unit.close()
        if (self.RH_type     != 'Scalar'): self.RH_unit.close()
        if (self.p0_type     != 'Scalar'): self.p0_unit.close()
        if (self.uz_type     != 'Scalar'): self.uz_unit.close()
        if (self.z_type      != 'Scalar'): self.z_unit.close()
        if (self.z0_air_type != 'Scalar'): self.z0_air_unit.close()

        #---------------------------------------------------
        # These are needed to compute Qn_SW and Qn_LW.
        #---------------------------------------------------
        if (self.albedo_type        != 'Scalar'): self.albedo_unit.close()
        if (self.em_surf_type       != 'Scalar'): self.em_surf_unit.close()
        if (self.dust_atten_type    != 'Scalar'): self.dust_atten_unit.close()
        if (self.cloud_factor_type  != 'Scalar'): self.cloud_factor_unit.close()
        if (self.canopy_factor_type != 'Scalar'): self.canopy_factor_unit.close()
        
##        if (self.Qn_SW_type  != 'Scalar'): self.Qn_SW_unit.close()        
##        if (self.Qn_LW_type  != 'Scalar'): self.Qn_LW_unit.close()
        
##        if (self.P_file      != ''): self.P_unit.close()
##        if (self.T_air_file  != ''): self.T_air_unit.close()
##        if (self.T_surf_file != ''): self.T_surf_unit.close()
##        if (self.RH_file     != ''): self.RH_unit.close()
##        if (self.p0_file     != ''): self.p0_unit.close()
##        if (self.uz_file     != ''): self.uz_unit.close()
##        if (self.z_file      != ''): self.z_unit.close()
##        if (self.z0_air_file != ''): self.z0_air_unit.close()
##        #--------------------------------------------------------
##        if (self.Qn_SW_file  != ''): self.Qn_SW_unit.close()        
##        if (self.Qn_LW_file  != ''): self.Qn_LW_unit.close()
        
    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        if (self.DEBUG):
            print 'Calling update_outfile_names()...'
        
        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.ea_gs_file  = (self.out_directory + self.ea_gs_file  )
        self.es_gs_file  = (self.out_directory + self.es_gs_file  )
        self.Qsw_gs_file = (self.out_directory + self.Qsw_gs_file )
        self.Qlw_gs_file = (self.out_directory + self.Qlw_gs_file )
        self.ema_gs_file = (self.out_directory + self.ema_gs_file )
        #------------------------------------------------------------
        self.ea_ts_file  = (self.out_directory + self.ea_ts_file  )
        self.es_ts_file  = (self.out_directory + self.es_ts_file  )
        self.Qsw_ts_file = (self.out_directory + self.Qsw_ts_file )
        self.Qlw_ts_file = (self.out_directory + self.Qlw_ts_file )
        self.ema_ts_file = (self.out_directory + self.ema_ts_file )
        
##        self.ea_gs_file = (self.case_prefix + '_2D-ea.rts')
##        self.es_gs_file = (self.case_prefix + '_2D-es.rts')
##        #-----------------------------------------------------
##        self.ea_ts_file = (self.case_prefix + '_0D-ea.txt')
##        self.es_ts_file = (self.case_prefix + '_0D-es.txt')

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        if (self.DEBUG):
            print 'Calling open_output_files()...'
        model_output.check_netcdf()
        self.update_outfile_names()
        
        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_EA_GRIDS):
            model_output.open_new_gs_file( self, self.ea_gs_file, self.rti,
                                           ## var_name='e_air',
                                           var_name='ea',
                                           long_name='vapor_pressure_in_air',
                                           units_name='mbar')
            
        if (self.SAVE_ES_GRIDS):
            model_output.open_new_gs_file( self, self.es_gs_file, self.rti,
                                           ## var_name='e_surf',
                                           var_name='es',
                                           long_name='vapor_pressure_at_surface',
                                           units_name='mbar')
        if (self.SAVE_QSW_GRIDS):
            model_output.open_new_gs_file( self, self.Qsw_gs_file, self.rti,
                                           var_name='Qsw',
                                           long_name='net_shortwave_radiation',
                                           units_name='W/m^2')
            
        if (self.SAVE_QLW_GRIDS):
            model_output.open_new_gs_file( self, self.Qlw_gs_file, self.rti,
                                           var_name='Qlw',
                                           long_name='net_longwave_radiation',
                                           units_name='W/m^2')
        if (self.SAVE_EMA_GRIDS):
            model_output.open_new_gs_file( self, self.ema_gs_file, self.rti,
                                           var_name='ema',
                                           long_name='air_emissivity',
                                           units_name='none')
            
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs 
        if (self.SAVE_EA_PIXELS):
            model_output.open_new_ts_file( self, self.ea_ts_file, IDs,
                                           ## var_name='e_air',
                                           var_name='ea',
                                           long_name='vapor_pressure_in_air',
                                           units_name='mbar')

        if (self.SAVE_ES_PIXELS):
            model_output.open_new_ts_file( self, self.es_ts_file, IDs,
                                           ## var_name='e_surf',
                                           var_name='es',
                                           long_name='vapor_pressure_at_surface',
                                           units_name='mbar')

        if (self.SAVE_QSW_PIXELS):
            model_output.open_new_ts_file( self, self.Qsw_ts_file, IDs,
                                           var_name='Qsw',
                                           long_name='net_shortwave_radiation',
                                           units_name='W/m^2')

        if (self.SAVE_QLW_PIXELS):
            model_output.open_new_ts_file( self, self.Qlw_ts_file, IDs,
                                           var_name='Qlw',
                                           long_name='net_longwave_radiation',
                                           units_name='W/m^2')
            
        if (self.SAVE_EMA_PIXELS):
            model_output.open_new_ts_file( self, self.ema_ts_file, IDs,
                                           var_name='ema',
                                           long_name='air_emissivity',
                                           units_name='none')
            
    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        if (self.DEBUG):
            print 'Calling write_output_files()...'

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

    #  write_output_files()
    #-------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_EA_GRIDS):   model_output.close_gs_file( self, 'ea')
        if (self.SAVE_ES_GRIDS):   model_output.close_gs_file( self, 'es')
        if (self.SAVE_QSW_GRIDS):  model_output.close_gs_file( self, 'Qsw')
        if (self.SAVE_QLW_GRIDS):  model_output.close_gs_file( self, 'Qlw')
        if (self.SAVE_EMA_GRIDS):  model_output.close_gs_file( self, 'ema')
        #-------------------------------------------------------------------
        if (self.SAVE_EA_PIXELS):  model_output.close_ts_file( self, 'ea') 
        if (self.SAVE_ES_PIXELS):  model_output.close_ts_file( self, 'es') 
        if (self.SAVE_QSW_PIXELS): model_output.close_ts_file( self, 'Qsw') 
        if (self.SAVE_QLW_PIXELS): model_output.close_ts_file( self, 'Qlw')
        if (self.SAVE_EMA_PIXELS): model_output.close_ts_file( self, 'ema')
        
    #   close_output_files()        
    #-------------------------------------------------------------------  
    def save_grids(self):
       
        if (self.SAVE_EA_GRIDS):
            model_output.add_grid( self, self.e_air,  'ea', self.time_min )
            
        if (self.SAVE_ES_GRIDS):
            model_output.add_grid( self, self.e_surf, 'es', self.time_min )

        if (self.SAVE_QSW_GRIDS):
            model_output.add_grid( self, self.Qn_SW, 'Qsw', self.time_min )
            
        if (self.SAVE_QLW_GRIDS):
            model_output.add_grid( self, self.Qn_LW, 'Qlw', self.time_min )

        if (self.SAVE_EMA_GRIDS):
            model_output.add_grid( self, self.em_air, 'ema', self.time_min )
            
    #   save_grids()            
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min   ######
        
        if (self.SAVE_EA_PIXELS):
            model_output.add_values_at_IDs( self, time, self.e_air,  'ea', IDs )
            
        if (self.SAVE_ES_PIXELS):
            model_output.add_values_at_IDs( self, time, self.e_surf, 'es', IDs )

        if (self.SAVE_QSW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Qn_SW, 'Qsw', IDs )
            
        if (self.SAVE_QLW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Qn_LW, 'Qlw', IDs )

        if (self.SAVE_EMA_PIXELS):
            model_output.add_values_at_IDs( self, time, self.em_air, 'ema', IDs )
            
    #   save_pixel_values()
    #-------------------------------------------------------------------
#---------------------------------------------------------------------------------
def compare_em_air_methods():

    #--------------------------------------------------------------
    # Notes:  There are two different methods that are commonly
    #         used to compute the vapor pressure of air, e_air,
    #         and then the emissivity of air, em_air, for use in
    #         longwave radiation calculations.  This routine
    #         compares them graphically.
    #
    # NB!     This hasn't been tested since conversion from IDL.
    #-------------------------------------------------------------
    import matplotlib.pyplot
    
    T_air = np.arange(80, dtype='Float32') - np.float64(40)   #[Celsius]  (-40 to 40)
    RH  = np.float64(1.0)
    C2K = np.float64(273.15)
    
    #--------------------------
    # Brutsaert (1975) method
    #--------------------------
    term1   = (np.float64(17.3) * T_air) / (T_air + np.float64(237.3))   ######### DOUBLE CHECK THIS (7/26/13)
    e_air1  = RH * np.float64(0.611) * np.exp( term1 )  # [kPa]
    em_air1 = np.float64(1.72) * (e_air1 / (T_air + C2K)) ** (np.float64(1) / 7)
    
    #---------------------------
    # Satterlund (1979) method
    #----------------------------
    # NB! e_air has units of Pa
    #----------------------------
    term2   = np.float64(2353) / (T_air + C2K)
    e_air2  = RH * np.float64(10) ** (np.float64(11.40) - term2)   # [Pa]
    eterm   = np.exp(-np.float64(1) * (e_air2 / np.float64(100)) ** ((T_air + C2K) / np.float64(2016)))
    em_air2 = np.float64(1.08) * (np.float64(1) - eterm)
    
    #----------------------------
    # Plot the two e_air curves
    #--------------------------------
    # These two agree quite closely
    #--------------------------------
    matplotlib.pyplot.figure(figsize=(8, 6), dpi=80)
    matplotlib.pyplot.show()
    matplotlib.pyplot.plot(T_air, e_air1)
    matplotlib.pyplot.show()
    ## oplot(T_air, (e_air2 / np.float64(1000)), psym=-3)   # [Pa -> kPa]
    
    #-----------------------------
    # Plot the two em_air curves
    #--------------------------------------------------
    # These two don't agree very well for some reason
    #--------------------------------------------------
    matplotlib.pyplot.figure(figsize=(8, 6), dpi=80)
    matplotlib.pyplot.show()
    matplotlib.pyplot.plot(T_air, em_air1)
    matplotlib.pyplot.show()
    ## oplot(T_air, em_air2, psym=-3)
    
#   compare_em_air_Methods
#---------------------------------------------------------------------------------
    
        
