
## Copyright (c) 2001-2010, Scott D. Peckham
## January, May, July 2009
## May 2010 (changes to unit_test() and read_cfg_file()

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "dynamic wave" channel flow component
#          and related functions.  It inherits from the channels
#          "base class" in "channels_base.py".
#-----------------------------------------------------------------------
#
#  class channels_component
#
#      get_attribute()   # (10/26/11)
#      get_input_var_names()     # (defined in channels_base.py)
#      get_output_var_names()    # (defined in channels_base.py)
#      get_var_name()            # (defined in channels_base.py)
#      get_var_units()           # (defined in channels_base.py)
#      ------------------------
#      update_velocity()
#
#-----------------------------------------------------------------------

import numpy as np

from topoflow.components import channels_base

#-----------------------------------------------------------------------
class channels_component(channels_base.channels_component):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'Channels_Dynamic_Wave',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #------------------------------------------------------
        'comp_name':          'ChannelsDynamWave',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Channels_Dynamic_Wave.cfg.in',
        'cfg_extension':      '_channels_dynamic_wave.cfg',
        'cmt_var_prefix':     '/ChannelsDynamWave/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Channels_Dynamic_Wave.xml',
        'dialog_title':       'Channels: Dynamic Wave Parameters',
        'time_units':         'seconds' }
    
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        #-----------------------------------------------------------
        # This is done in channels_base.set_computed_input_vars()
        #-----------------------------------------------------------
        # self.KINEMATIC_WAVE = False 
        # self.DIFFUSIVE_WAVE = False
        # self.DYNAMIC_WAVE   = True 

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

    #   get_attribute()
    #-------------------------------------------------------------------
    def update_velocity(self):

        #--------------------------------------------------------
        # NOTES:  This update to u uses values of u, d, Q, and
        #         S_free from the last time step, then u is
        #         overwritten with new values.
        #--------------------------------------------------------
        
        #-----------------------------
        # Get the free-surface slope
        #-----------------------------
        self.update_free_surface_slope()
        
        #----------------------------------------
        # Compute the wetted bed area (2/12/07)
        #----------------------------------------
        #pb = d / np.cos(angle)
        #;ww = np.where(np.abs(angle) >= (np.pi/2d), nww)  ;(should disallow)
        #;if (nww ne 0) then pb[ww]=(d * 10d)
        #pw = width + (2d * pb)
        #bA = pw * ds_chan        ;(a 2D array)
        #*** bA = da             ;(for testing only)

        #---------------------------------------
        # Always need to know where (d <= 0),
        # so do it here before everything else
        #--------------------------------------- 
        wp = np.where( self.d > 0 )
        np = np.size( wp[0])
        wn = np.where( self.d <= 0 )
        nn = np.size( wn[0] )
        #-----------------------------------------
        # This makes f=0 and du=0 where (d le 0)
        #-----------------------------------------
        dinv = self.d.copy()  # (just to get size)
        #########################
        if (np != 0):
            dinv[wp] = (np.float64(1) / self.d[wp])
        if (nn != 0):    
            dinv[wn] = np.float64(0)
        
        #--------------------------------
        # Compute an effective depth so
        # that f doesn't blow up ??
        #--------------------------------
        #** deff = (d > 1e-9)
        #** deff = (d > 0.001)   ;(1 mm)
        #----------------------------------
        # deff = self.d
        
        #-----------------------------
        # Compute f for Manning case
        #-----------------------------
        if (self.MANNING):    
            nval  = self.nval
            self.f = self.g * (nval * nval * (dinv ** self.one_third))
            
        #---------------------------------
        # Compute f for Law of Wall case
        #---------------------------------
        if (self.LAW_OF_WALL):
            smoothness = (self.aval / self.z0val) * self.d
            #**** smoothness = (0.476d / z0val) * deff
            smoothness = np.maximum(smoothness, float64(1.1))
            #################################################
            self.f = (self.kappa / np.log(smoothness)) ** np.float64(2)
        
        #-----------------
        # Before 2/13/07
        #----------------------------------------------
        # Increment flow velocities: inputs - outputs
        #----------------------------------------------
        #*** grav = 9.81d * (S_free * d)    ;(USE  deff vs. d ??)
        #grav = 9.81d * (S_free * deff)     ;(USE  deff vs. d ??)
        #fric = u*(R + (f * u))
        #acc  = (grav - fric)           ;(positive or negative)
        #u2   = u + (dt * acc / deff)   ;(before next part)
        
        #------------------------------------------
        # Start with the interior (nonflux) terms
        #------------------------------------------
        # When (n eq 0), should have (d gt 0) and
        # (u eq 0), so (grav gt 0) and (acc gt 0)
        # but fric = Atrm = Rtrm = 0.
        #------------------------------------------
        # Multiply all scalars first, then grids
        # Note:  u, d and Q are always grids
        # Note:  (Pw / wtop) = (A2 / Atop)
        #------------------------------------------
        grav = self.g * (self.S_free * self.d)
        fric = self.f * (self.u ** np.float64(2))
        #---------------------------------------
        #####################################################
        angle  = self.angle
        width  = self.width
        d      = self.d
        u      = self.u
        Q      = self.Q
        #---------------------------------------
        p1     = self.d8.p1   ; w1 = self.d8.w1
        p2     = self.d8.p2   ; w2 = self.d8.w2
        p3     = self.d8.p3   ; w3 = self.d8.w3
        p4     = self.d8.p4   ; w4 = self.d8.w4
        p5     = self.d8.p5   ; w5 = self.d8.w5
        p6     = self.d8.p6   ; w6 = self.d8.w6
        p7     = self.d8.p7   ; w7 = self.d8.w7
        p8     = self.d8.p8   ; w8 = self.d8.w8
        #####################################################
        #---------------------------------------
        wtop = width + (2 * d * np.tan(angle))    #(top width)
        Atop = self.d8.ds * wtop                    #(top area)
##        Atop = ds_chan * wtop                     #(top area)
        Rtrm = (u * self.R) * (self.da / Atop)      #(da = pixel area)
        #---------------------------------------
        Pw = width + (np.float64(2) * d / np.cos(angle))  #(wetted perimeter)
        A2 = self.d8.ds * Pw                          #(wetted surf. area)
##        A2 = ds_chan * Pw                           #(wetted surf. area)
        Atrm = (u * Q) * ((np.float64(1) / Atop) - (np.float64(1) / A2))
        #---------------------------------------
        acc = (grav + Atrm - fric - Rtrm)          #(positive or negative)
        u2  = u + (self.dt * dinv * acc)           #(before next part)
        #----------------------------------
        fac = (self.dt / A2) * dinv                #(always grid)
        uu  = u * (Pw / wtop)                      #(always grid)
            
        #-------------------------------------------
        # Add momentum fluxes from D8 child pixels
        #-------------------------------------------
        if (self.d8.p1_OK):    
            u2[p1] += (u[w1] - uu[p1]) * Q[w1] * fac[p1]
        if (self.d8.p2_OK):    
            u2[p2] += (u[w2] - uu[p2]) * Q[w2] * fac[p2]
        if (self.d8.p3_OK):    
            u2[p3] += (u[w3] - uu[p3]) * Q[w3] * fac[p3]
        if (self.d8.p4_OK):    
            u2[p4] += (u[w4] - uu[p4]) * Q[w4] * fac[p4]
        if (self.d8.p5_OK):    
            u2[p5] += (u[w5] - uu[p5]) * Q[w5] * fac[p5]
        if (self.d8.p6_OK):    
            u2[p6] += (u[w6] - uu[p6]) * Q[w6] * fac[p6]
        if (self.d8.p7_OK):    
            u2[p7] += (u[w7] - uu[p7]) * Q[w7] * fac[p7]
        if (self.d8.p8_OK):    
            u2[p8] += (u[w8] - uu[p8]) * Q[w8] * fac[p8]
        
        #--------------------------------
        # Don't allow u2 to be negative
        #----------------------------------------------
        # If uphill flow is allowed, to which pixel ?
        # This worked when d0 grid was used.
        #----------------------------------------------
        u2 = np.maximum(u2, np.float64(0))
        
        #---------------
        # Copy u2 to u
        #---------------
        self.u = u2
        
    #    update_velocity()                       
    #-------------------------------------------------------------------


         
