
## Copyright (c) 2001-2013, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## April, May, July, August 2009
## May 2010 (changes to initialize() and read_cfg_file()

#********************************************************
# NB!  update_water_balance() method is not ready yet!
#********************************************************

#-----------------------------------------------------------------------
#  Notes:  This file defines a "base class" for evaporation
#          components as well as functions used by most or
#          all evaporation methods.  The methods of this class
#          should be over-ridden as necessary for different
#          methods of modeling evaporation.
#-----------------------------------------------------------------------

# (5/7/09)  "Nested WHERE calls" work differently
# in numpy than in IDL. For a 1D array called "a":
#
#     >>> a = np.arange(11)-5
#     >>> w = np.where(a < 0)
#     >>> w2 = np.where(a[w] > -3)
#     >>> print a[w[w2]]  # (this gives an error)
#     >>> print a[w2]     # (this works)
#     >>> print a[w][w2]  # (this works, too)

# For a 2D array called "a":
#
#     >>> a = np.arange(9) - 4
#     >>> a = a.reshape(3,3)
#     >>> w = np.where(a < 0)
#     >>> w2 = np.where(a[w] > -3)
#     >>> print a[w[w2]]    # TypeError: tuple indices must be integers
#     >>> print a[w2]       # IndexError: index (3) out of range (0<=index<=2) in dimension 0
#     >>> a[w][w2] = 99     # No error, but this doesn't work.
#     >>> print a[w][w2]    # (this works)
#     >>> print a.flat[w2]  # (this works, same result as last line)
#     >>> a.flat[w2] = 99   # (this works)
#     >>> a.flat[w2] = [-2,-1]  # (this works)
#     >>> np.put(a, w2, 99)  # (this works)

#-----------------------------------------------------------------------
#
#  class evap_component    (inherits from CSDMS_base)
#
#      (see non-base components for BMI functions)
#
#      ------------------------
#      set_constants()
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      ---------------------------
#      get_cca_ports()           # (in CSDMS_base.py)
#      get_cca_port_info()       # (5/11/10)
#      embed_child_components()
#      add_child_ports()
#      initialize_ports()
#      release_cca_ports()       # (in CSDMS_base.py)
#      -----------------------------
#      check_input_types()
#      check_if_types_match()
#      initialize_computed_vars()
#      -----------------------------
#      update_Qc()                   # (not used yet)
#      update_ET_rate()
#      update_ET_integral()
#      update_water_balance()
#      ------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()     #####
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.utils import BMI_base
from topoflow.utils import cfg_files as cfg
from topoflow.utils import model_input
from topoflow.utils import model_output

#-----------------------------------------------------------------------
class evap_component( BMI_base.BMI_component):

    #-------------------------------------------------------------------
    def set_constants(self):

        #---------------------------------
        # From Bob Bolton (Nov. 3, 2009)
        #---------------------------------
        self.mps_to_mmph = np.float64(3600000)
        self.mmph_to_mps = (np.float64(1) / np.float64(3600000))
        self.forever     = np.float64(999999999)  # [minutes]
        
    #   set_constants()
    #-------------------------------------------------------------------
    def initialize(self, cfg_prefix=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print ' '
            print 'Evaporation component: Initializing...'
        
        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_prefix = cfg_prefix
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()    # (12/3/09)
        self.initialize_config_vars()
        self.read_grid_info()
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars()
        
        #------------------------------------------------------
        # NB! "Sample steps" must be defined before we return
        #     Check all other process modules.
        #------------------------------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print 'Evaporation component: Disabled.'
            self.ET     = self.initialize_scalar(0, dtype='float64')
            self.vol_ET = self.initialize_scalar(0, dtype='float64')
            self.DONE   = True
            self.status = 'initialized'
            return

        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #-----------------------
        # Initialize variables
        #-----------------------
        self.initialize_computed_vars()  # (such as 'ET')
        #---------------------------------------------
        # Must come after initialize_computed_vars()
        # but not needed with new framework. (2/5/13)
        #---------------------------------------------
        ## self.initialize_required_components(mode)
        self.check_input_types()   # (Uses "mp" vars)
        
        self.open_output_files()
        self.status = 'initialized'
        
    #   initialize()
    #-------------------------------------------------------------------
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):
        
        #-------------------------------------------------
        # Note: self.ET already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI)
        
        #-------------------------
        # Update computed values 
        #-------------------------
        self.update_ET_rate()
        self.update_ET_integral()
        self.update_water_balance()

        #---------------------------------------
        # Read next ET vars from input files ?
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
        #-----------------------------------------------
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
        if (self.comp_status == 'Enabled'):
            self.close_input_files()   ##  TopoFlow input "data streams"
            self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Evaporation component')

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
##    def get_cca_port_info(self):
##
##        self.cca_port_names = ['meteorology', 'channels',
##                               'snow', 'infil', 'satzone']
##        self.cca_port_short = ['mp', 'cp', 'sp', 'ip', 'gp']
##        self.cca_port_type  = "IRFPort"
##        self.cca_project    = "edu.csdms.models"
##
##    #   get_cca_port_info()
##    #-------------------------------------------------------------------
##    def embed_child_components(self):
##
##        #------------------------------------------------
##        # Note: Don't call this if (self.CCA == True).
##        #------------------------------------------------
##        # Instantiate and embed "process components"
##        # in the place of the CCA ports.
##        #------------------------------------------------
##        # But how do we choose a given process method
##        # such as "kinematic wave" ??
##        #------------------------------------------------
##        import met_base
##        import channels_kinematic_wave
##        import snow_degree_day
##        import infil_green_ampt
##        import GW_darcy_layers
##
##        self.mp = met_base.met_component()
##        self.cp = channels_kinematic_wave.channels_component()
##        self.sp = snow_degree_day.snow_component()
##        self.ip = infil_green_ampt.infil_component()
##        self.gp = GW_darcy_layers.GW_component()
##
##    #   embed_child_components()
##    #-------------------------------------------------------------------
##    def add_child_ports(self):
##
##        self.add_child_port('cp', 'sp')
##        self.add_child_port('cp', 'ep', SELF=True)
##        self.add_child_port('cp', 'ip')
##        self.add_child_port('cp', 'gp')
####        self.add_child_port('cp', 'dp')   # (NOT READY YET)
##        self.add_child_port('cp', 'mp')
####        self.add_child_port('cp', 'iip')  # (NOT READY YET)
##        #----------------------------------
##        self.add_child_port('sp', 'mp')
##        #----------------------------------
##        self.add_child_port('ip', 'mp')
##        self.add_child_port('ip', 'sp')
##        self.add_child_port('ip', 'ep', SELF=True)
##        self.add_child_port('ip', 'gp')
##        #----------------------------------
##        self.add_child_port('gp', 'cp')
##        self.add_child_port('gp', 'ip')
##        #----------------------------------
##        self.add_child_port('mp', 'sp')
##        
##        ## self.add_child_port('mp', 'pp')   ##### in future ?
##        
##    #   add_child_ports()        
##    #-------------------------------------------------------------------
##    def initialize_ports(self):
##
##        #----------------------------------------------
##        # Initialize the process objects/components
##        # This is also where output files are opened.
##        #----------------------------------------------
##        # Must be before "h_table" is needed (by EVAP) ??
##        if (self.gp.get_status() != 'initialized'):        # groundwater
##            self.gp.initialize( cfg_prefix=self.cfg_prefix )
##            
##        # Must be before snow ??
##        if (self.mp.get_status() != 'initialized'):        # met vars 
##            self.mp.initialize( cfg_prefix=self.cfg_prefix )
##            
##        if (self.cp.get_status() != 'initialized'):        # channel vars
##            self.cp.initialize( cfg_prefix=self.cfg_prefix )
##            
##        if (self.sp.get_status() != 'initialized'):        # snow vars
##            self.sp.initialize( cfg_prefix=self.cfg_prefix ) 
##
##        #-----------------------------------------------------------
##        # (6/28/10) Note that every initialize() calls initialize_
##        # required_components().  But the latter function only
##        # calls initialize_ports() if component's mode is "driver".
##        # So we don't need to add "initializing" test below to
##        # prevent things from trying to initialize each other.
##        #-----------------------------------------------------------
##        if (self.ip.get_status() != 'initialized'):        # infil vars
##            self.ip.initialize( cfg_prefix=self.cfg_prefix )
##            
##    #   initialize_ports()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing ET.
        #----------------------------------------------------
        are_scalars = np.array([
                         # self.is_scalar('d'),
                         #---------------------------------
                         # self.is_scalar('h_table'),
                         #---------------------------------
                         self.is_scalar('T_air'),
                         self.is_scalar('T_surf') ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #******************************************************
        #  Any faster to use np.empty vs. np.zeros ??
        #******************************************************
        self.ET = np.zeros([self.ny, self.nx], dtype='Float64')
        self.vol_ET = self.initialize_scalar(0, dtype='float64')
        
        #------------------------------------------
        # h_table = water table height
        # Assume h_table is always a grid.
        # h_table, dzw and ET must be compatible.
        #------------------------------------------
##        H_IS_GRID = self.gp.is_grid('h_table')
##        if (H_IS_GRID):
##            self.ET = zeros([self.ny, self.nx], dtype='Float64')            
##        else:
##            self.ET = np.float64(0)
##            print '********* WARNING: water table is not a grid'
##            print '                   but it should be.'

    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_Qc(self):
  
        #---------------------------------------------
        # Compute the conductive energy between the
        # surface and subsurface using Fourier's law
        #---------------------------------------------
        # soil_x is converted from [cm] to [m] when
        # it is read from the GUI and then stored
        #---------------------------------------------
        ## T_surf  = self.get_port_data('T_surf', self.mp)
        T_surf  = self.T_surf ## (2/3/13)
        self.Qc = self.K_soil * (self.T_soil_x - T_surf) / self.soil_x
        
    #   update_Qc()
    #-------------------------------------------------------------------
    def update_ET_rate(self):

        #-------------------------------------------------
        # Note: This method should be overridden by some
        #       method for computing ET rate.
        #-------------------------------------------------
        pass
        
        #---------------------------------------------------------------
        # NB!  h_snow is needed by the Bulk_Exchange_Coeff function
        #     (which is called by the Energy_Balance_ET_Rate function)
        #     to adjust reference height, z.  Don't just use h0_snow.
        #---------------------------------------------------------------
##        if   (self.method == 0):    
##            ET = np.float64(0.0)
##        elif (self.method == 1):    
##            ET = Priestley_Taylor_ET_Rate(self.alpha, self.K_soil, \
##                                          self.T_soil_x, self.soil_x, \
##                                          self.mp.Qn_SW, self.mp.Qn_LW, \
##                                          self.mp.T_air, self.mp.T_surf)
##        elif (self.method == 2):    
##            ET = Energy_Balance_ET_Rate(self.K_soil, \
##                                        self.T_soil_x, self.soil_x, \
##                                        self.mp.Qn_SW, self.mp.Qn_LW, \
##                                        self.mp.T_air, self.mp.T_surf, \
##                                        self.mp.uz, self.mp.z, \
##                                        self.mp.z0_air, self.mp.rho_air, \
##                                        self.mp.Cp_air, self.sp.h_snow)
##        else:
##            raise RuntimeError('No match found for expression')
##
##        self.ET = ET
        
    #   update_ET_rate()
    #-------------------------------------------------------------------
    def update_ET_integral(self):

        #------------------------------------------------
        # Update mass total for GW, sum over all pixels
        #------------------------------------------------   
        volume = np.double(self.ET * self.da * self.dt)  # [m^3]
        if (np.size( volume ) == 1):
            self.vol_ET += (volume * self.rti.n_pixels)
        else:
            self.vol_ET += np.sum(volume)
            
    #   update_ET_integral()
    #-------------------------------------------------------------------
    def update_water_balance(self):

        #-------------------------------------------------------
        # Note: Computed ET values are generally taken to be
        #       "potential" values which may not be achieved
        #       if there is not enough water at or near the
        #       surface.  This function first tries to consume
        #       the required water from surface water (depth)
        #       and then goes on to extract water from the
        #       top soil layer (subsurface).
        #-------------------------------------------------------
        # Note: ET = ET rate with units of [m/s].
        #       ev = ET_vars = structure
        #       mv = met_vars = structure
        #       gv = gw_vars = structure
        #       iv = infil_vars = structure
        #        d = depth of surface water [m]
        #        h = water table height above datum
        #        y = thicknesses [m] of all soil layers
        #            when using Darcy subsurface flow
        #-------------------------------------------------------

        #-------------------------------------------
        # If Richards' equation is being used for
        # infiltration, then don't need to remove
        # water from layers as done in remainder
        # and y (wetted thicknesses) is not needed
        #-------------------------------------------
        # But still need to remove surface water
        # first !!  This isn't done yet. ********
        #-------------------------------------------
        RICHARDS_1D = (self.ip.get_scalar_long('RICHARDS'))
##        print 'In evap_base.update_water_balance(),'
##        print '   RICHARDS_1D is set to: ', RICHARDS_1D
##        print ' '
        if (RICHARDS_1D): return
        
        #-------------------------------------
        # Depth of water to be removed by ET
        #-------------------------------------------
        # (8/25/09) Does it make sense to allow ET
        # ET and dzw to be scalars ??
        #-------------------------------------------
        dzw = (self.dt * self.ET)
        ## print 'size(dzw) =', np.size(dzw)
        
        #----------------
        # For debugging
        #----------------
        #if (np.size(dzw) == 1) then begin
        #    msg = [' ','ERROR: dzw is not an array. ', ' ']
        #    result = GUI_Message(msg, /INFO)
        #    STOP
        #endif

        ## depth = self.cp.get_grid_double('d')      # (channel depth)
        ## depth = self.get_port_data('d', self.cp)  # (channel depth)
        depth = self.depth    # (2/3/13, "d@channel")
        UPDATE_DEPTH = False
        
        wL  = np.where( dzw <= depth )
        nwL = np.size( wL[0] )
        wG  = np.where( dzw > depth )
        nwG = np.size( wG[0] )

        if (nwL != 0):    
            #---------------------------------
            # Reduce the surface water depth
            #---------------------------------
            depth[wL]    = (depth[wL] - dzw[wL])
            UPDATE_DEPTH = True
            dzw[wL]      = np.float64(0)
        
        if (nwG != 0):    
            #-----------------------------
            # Save a copy of initial dzw
            #-----------------------------
            dzw0 = dzw.copy()
            
            #-------------------------------------
            # Consume all surface water first
            # This doesn't account for channels.
            #-------------------------------------
            dzw[wG]      = dzw[wG] - depth[wG]
            depth[wG]    = np.float64(0)
            UPDATE_DEPTH = True
            
            #---------------------------------------
            # Try to take remainder from top layer
            # Compute water content of top layer
            #---------------------------------------
            # Used before 7/13/06
            #----------------------
            # p  = gv.soil_P[0]  ;(top layer porosity)
            # y0 = y[*,*,0]
            # content_1 = (y0[wG] * p)
            #---------------------------------------------
            # self.gp.qs is a 1D array of doubles that
            # gives theta_sat for each soil layer.
            # This is taken equal to porosity here.
            #---------------------------------------------
            # self.gp.y[0,:,:] is a grid of doubles that
            # gives the "wetted thickness" of top layer
            #---------------------------------------------
            p0 = self.p0       # (2/3/13, new framework)
            y0 = self.y0       # (2/3/13, new framework)
            h  = self.h_table  # (2/3/13, new framework)
            #---------------------------------------------            
##            p0 = self.get_port_data('qs[0]', self.gp)
##            y0 = self.get_port_data('y[0,:,:]', self.gp)
##            h  = self.get_port_data('h_table', self.gp)
            #---------------------------------------------         
##            p0 = self.gp.get_scalar_double('qs[0]')   #########
##            y0 = self.gp.get_grid_double('y[0,:,:]')  #########
##            h  = self.gp.get_grid_double('h_table')
            #---------------------------------------------             
##            p0 = gv.qs[0]  # (top layer porosity)
##            y0 = gv.y[0,:,:]
##            h  = gv.h   ################ (5/7/09)

            SCALAR_POROSITY = (np.size(p0) == 1)  # (Always True now)
            if (SCALAR_POROSITY):    
                content_1 = (y0[wG] * p0)
            else:    
                content_1 = (y0[wG] * p0[wG])
            
            wwL  = np.where( dzw[wG] <= content_1 )
            nwwL = np.size( wwL[0] )
            wwG  = np.where( dzw[wG] > content_1 )
            nwwG = np.size( wwG[0] )

            #####################################################
            # See Notes at top regarding "nested WHERE calls".
            #####################################################

            #---------------------------------------------
            # Can get all remaining water from top layer
            # Reduce the water table height
            #---------------------------------------------
            if (nwwL != 0):    
                if (SCALAR_POROSITY):
                    dh = dzw.flat[wwL] / p0
                    #### dh = dzw[wG][wwL] / p0
                else:
                    dh = dzw.flat[wwL] / p0.flat[wwL]
                    #### dh = dzw[wG][wwL] / p0[wG][wwL]

                h.flat[wwL]   = h.flat[wwL]  - dh
                y0.flat[wwL]  = y0.flat[wwL] - dh
                dzw.flat[wwL] = np.float64(0)     # (not really needed ?)
                
##                h[wG][wwL]   = h[wG][wwL] - dh
##                y0[wG][wwL]  = y0[wG][wwL] - dh
##                dzw[wG][wwL] = np.float64(0)   # (not really needed ?)
            
            #-----------------------------------------------
            # Can't get all remaining water from top layer
            #-----------------------------------------------
            # Get what is available, and then redefine ET
            # for mass balance consistency
            #-----------------------------------------------
            if (nwwG != 0):
                dh = y0.flat[wwG]
                h.flat[wwG]   = h.flat[wwG] - dh
                y0.flat[wwG]  = np.float64(0)
                dzw.flat[wwG] = dzw.flat[wwG] - content_1[wwG]
                #################################################
                ##### Is there a problem in above line with
                ##### content_1[wwG] part ???
                #------------------------------------------------
                dzw_used    = dzw0.flat[wwG] - dzw.flat[wwG]
                self.ET.flat[wwG] = (dzw_used / self.dt)
                
##                dh = y0[wG][wwG]
##                h[wG][wwG]   = h[wG][wwG] - dh
##                y0[wG][wwG]  = np.float64(0)
##                dzw[wG][wwG] = dzw[wG][wwG] - content_1[wwG]
##                #--------------------------------------------
##                dzw_used    = dzw0[wG][wwG] - dzw[wG][wwG]
##                self.ET[wG][wwG] = (dzw_used / self.dt)
      
            #-------------------------
            # Replace top layer in y
            #-------------------------
            print '    ET component changing "h_table" in GW component.'
##            print '       type(y0) =', type(y0)
##            print '       type(h)  =', type(h)
            self.set_port_data('y[0,:,:]', y0, self.gp)
            self.set_port_data('h_table', h, self.gp)

        #---------------------------------------------------
        # Update the channel depth in channels component ?
        #---------------------------------------------------
        if (UPDATE_DEPTH):
            print '    ET component changing "d" in Channels component.'
            self.set_port_data('d', depth, self.cp)
            print '       type(depth) =', type(depth)
            
    #   update_water_balance()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #----------------------------------------------------
        # Note: Priestley-Taylor method needs alpha but the
        #       energy balance method doesn't, so they each
        #       override "*_input_files" methods.  (2/5/13)
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
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.ET_gs_file = (self.out_directory + self.er_gs_file)
        #---------------------------------------------------------
        self.ET_ts_file = (self.out_directory + self.er_ts_file)

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_nio()
        self.update_outfile_names()

        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_ER_GRIDS):
            model_output.open_new_gs_file( self, self.ET_gs_file, self.rti,
                                           var_name='ET',
                                           long_name='evaporation_rate',
                                           units_name='mm/hr')
                                           ### units_name='m/s')
            
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_ER_PIXELS):
            model_output.open_new_ts_file( self, self.ER_ts_file, IDs,
                                           var_name='ET',
                                           long_name='evaporation_rate',
                                           units_name='mm/hr')

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
##        if (SAVE_ER_GRIDS  == False) and  \
##           (SAVE_ER_PIXELS == False): return
           
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
    #---------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_ER_GRIDS):  model_output.close_gs_file( self, 'ET')
        #-----------------------------------------------------------------
        if (self.SAVE_ER_PIXELS): model_output.close_gs_file( self, 'ET')  

    #   close_output_files()   
    #---------------------------------------------------------------------  
    def save_grids(self):
        
        #-----------------------------------
        # Save grid stack to a netCDF file
        #---------------------------------------------
        # Note that add_grid() methods will convert
        # var from scalar to grid now, if necessary.
        #--------------------------------------------- 
        if (self.SAVE_ER_GRIDS):
            ET_mmph = self.ET * self.mps_to_mmph    # (Bolton 28 Aug)
            model_output.add_grid( self, ET_mmph, 'ET', self.time_min )

    #   save_grids()            
    #---------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min   ########
         
        if (self.SAVE_ER_PIXELS):
            ET_mmph = self.ET * self.mps_to_mmph    # (Bolton 28 Aug)
            model_output.add_values_at_IDs( self, time, ET_mmph, 'ET', IDs )

    #   save_pixel_values()
    #---------------------------------------------------------------------

    
        
