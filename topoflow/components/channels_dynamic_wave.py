
#  For details on how velocity is computed, see:
#    Peckham, S.D. (2009) Geomorphometry and spatial hydrologic modeling,
#    In: Hengl, T. and Reuter, H.I. (Eds), Geomorphometry: Concepts, Software
#    and Applications, Chapter 25, Developments in Soil Science,
#    vol. 33, Elsevier, 579-602, 
#    http://dx.doi.org/10.1016/S0166-2481(08)00025-1.

#--------------------------------------------------------------------------------
#  Copyright (c) 2001-2017, Scott D. Peckham
#
#  Feb 2017.  Cleaned up and re-tested.  Fixed missing dinv = 1/d.
#  Sep 2014.  New standard names and BMI updates and testing.
#  Nov 2013.  Converted TopoFlow to a Python package.
#  Feb 2013.  Adapted to use EMELI framework.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  May 2010.  Changes to unit_test() and read_cfg_file().
#  Jul 2009.  Updates.
#  May 2009.  Updates.
#  Jan 2009.  Converted from IDL to Python with I2PY.
#
#-----------------------------------------------------------------------
#  NOTES:  This file defines a "dynamic wave" channel flow component
#          and related functions.  It inherits from the channels
#          "base class" in "channels_base.py".
#-----------------------------------------------------------------------
#
#  class channels_component
#
#      get_component_name()
#      get_attribute()           # (10/26/11)
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
    def get_component_name(self):
  
        return 'TopoFlow_Channels_Dynamic_Wave'

    #   get_component_name()      
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
        # Note:  This update to u uses values of u, d, Q, and
        #        S_free from the last time step, then u is
        #        overwritten with new values.
        #--------------------------------------------------------
        # Note:  update_free_surface_slope() is called in
        #        channels_base.update() before calling this.
        #--------------------------------------------------------
        # Note:  For speed, multiply scalars first, then grids.
        #        u, d and Q are always grids
        #--------------------------------------------------------

		#------------------------------------------
		# Start with the interior (nonflux) terms
		#------------------------------------------
		# When (n eq 0), should have (d gt 0) and
		# (u eq 0), so (grav gt 0) and (acc gt 0)
		# but fric = Atrm = Rtrm = 0.
		#------------------------------------------
		angle  = self.angle
		width  = self.width
		vol    = self.vol
		d      = self.d       # (need for Pw)
		u      = self.u
		Q      = self.Q
		f      = self.f

        #-------------------------
        # Make sure u=0 when d=0
        #-------------------------
		## wg = self.d_is_pos
		wb = self.d_is_neg
		u[ wb ] = 0.0

		## d[ self.d8.noflow_IDs ] = 0.0

		#----------------------------------------
		# Update inverse of channel flow volume
		#------------------------------------------------
		# Wherever vol = 0, we should have u = 0, Q = 0.
		#------------------------------------------------
        # Note: Subscripting with Boolean arrays works
        #       even if all values are False.
		#------------------------------------------------
		Vinv = vol.copy()
		wg = (vol > 0)         # (Boolean array)
		### wg = (vol > 1e-5)         # (Boolean array)
		wb = np.invert( wg )   # (opposite Boolean array)
		Vinv[ wg ] = 1 / vol[ wg]
		Vinv[ wb ] = 0.0

		#------------------------------------
		# Compute the areas that are needed
		#------------------------------------------------------------
		# Note:  Ap = "wetted surface area", which is not the same
		#        as "wetted area" used to compute hydraulic radius.
		#------------------------------------------------------------
		# Note:  (Atop / wtop ) = (Ap / Pw) = d8.ds,
		#        so (Pw / wtop) = (Ap / Atop).
		#------------------------------------------------------------
		## wtop = width + (2 * d * np.tan(angle))   # (top width)
		## Atop = self.d8.ds * wtop                 # (top area)
		## A1 = self.d8.ds * width                  # (bottom area)
        ## Acell = self.da                          # (grid cell area)
        #--------------------------------------
        # See update_trapezoid_Rh() for P_wet
        #--------------------------------------
		## Pw = width + (np.float64(2) * d / np.cos(angle))
		Pw = self.P_wet          # (wetted perimeter)
		Ap = self.d8.ds * Pw     # (wetted surf. area)

		#-----------------------------------------------
		# Compute the "Q term" (momentum outflow rate)
		#-----------------------------------------------
		Qtrm = (u * Q)

		#----------------------------------------------
		# Compute the gravitational acceleration term
		#--------------------------------------------------------
		# Only the grav term will contribute at first, when d=0
		#--------------------------------------------------------
		grav = self.g * self.S_free * vol

        #--------------------------------------------------
        # Compute the frictional loss (shear stress) term
        #------------------------------------------------------
        # This can be large enough to make momentum negative,
        # for some reason, whether we use Manning or Log Law
        # of the Wall.  Using a smaller dt helps.
        #------------------------------------------------------
        # For Manning, we have:
        # f = g * (n^2) / d^(1/3)    (or d -> Rh)
        #------------------------------------------------------
		fric = (f * u * u * Ap)

        #------------------------------------------------------
        # Other ideas for how to compute frictional loss term
        #------------------------------------------------------
# 		dinv = d.copy()
# 		wg = (d > 0)         # (Boolean array)
# 		wb = np.invert( wg )   # (opposite Boolean array)
# 		dinv[ wg ] = 1 / d[ wg]
# 		dinv[ wb ] = 0.0
# 		n      = self.nval
# 		S_fric = n * n * u * u * (dinv**(4./3))
# 		fric   = (self.g * S_fric * vol)
        #------------------------------------------------------
		## fric = (self.tau * Ap) / self.rho_H2O  # (doesn't work)
        #------------------------------------------------------
		## fric = (f * u * np.abs(u) * Ap)
        #------------------------------------------------------
		## fric = np.minimum( fric, grav )   ##############

		#----------------------------------------------------
		# Compute the "acceleration" (positive or negative)
		#----------------------------------------------------
		dM = (grav - fric - Qtrm) * self.dt  
		## du = (grav - fric) * self.dt    # (for testing)  ###########

		#---------------------------------------------------------
		# Compute the new flow momentum (divided by rho) due
		# only to contributions from "interior" (nonflux) terms.
		# We'll modify this with flux terms below.
		#---------------------------------------------------------
		M = (u * vol) + dM

		#----------------------------------------
		# Now prepare to compute flux terms
        # which also contribute to the velocity
		#----------------------------------------
		# Note that "compound statements" with 
		# a semicolon are OK.
		#---------------------------------------
		p1 = self.d8.p1   ;  w1 = self.d8.w1
		p2 = self.d8.p2   ;  w2 = self.d8.w2
		p3 = self.d8.p3   ;  w3 = self.d8.w3
		p4 = self.d8.p4   ;  w4 = self.d8.w4
		p5 = self.d8.p5   ;  w5 = self.d8.w5
		p6 = self.d8.p6   ;  w6 = self.d8.w6
		p7 = self.d8.p7   ;  w7 = self.d8.w7
		p8 = self.d8.p8   ;  w8 = self.d8.w8

        #------------------------------------
        # Use this to compute "uQ_in" below
        #------------------------------------
		## M0 = M.copy()  #######

		#-------------------------------------------
		# Add momentum fluxes from D8 child pixels
		#-------------------------------------------
		if (self.d8.p1_OK):
			M[p1] += u[w1] * Q[w1] * self.dt
		if (self.d8.p2_OK):
			M[p2] += u[w2] * Q[w2] * self.dt
		if (self.d8.p3_OK):
			M[p3] += u[w3] * Q[w3] * self.dt  
		if (self.d8.p4_OK):
			M[p4] += u[w4] * Q[w4] * self.dt 
		if (self.d8.p5_OK):
			M[p5] += u[w5] * Q[w5] * self.dt  
		if (self.d8.p6_OK):
			M[p6] += u[w6] * Q[w6] * self.dt   
		if (self.d8.p7_OK):
			M[p7] += u[w7] * Q[w7] * self.dt   
		if (self.d8.p8_OK):
			M[p8] += u[w8] * Q[w8] * self.dt    

        #------------------------------------------
        # How much uQ was added from all D8 kids?
        #------------------------------------------
		## uQ_in = (M - M0) / self.dt  #######

        #--------------
        # For testing
        #--------------
# 		if ((self.time_index % 500) == 0):
# 			print '## M:     min, max =', M.min(),     M.max()
# 			print '## dM:    min, max =', dM.min(),    dM.max()
# 			print '## grav:  min, max =', grav.min(),  grav.max()
# 			print '## fric:  min, max =', fric.min(),  fric.max()
# 			print '## Qtrm:  min, max =', Qtrm.min(),  Qtrm.max()
# 			print '##'
# 			print '## Q:     min, max =', Q.min(),     Q.max()
# 			print '## vol:   min, max =', vol.min(),   vol.max()
# 			print '## d:     min, max =', d.min(),     d.max()
# 			print '## u:     min, max =', u.min(),     u.max()
# 			print '## f:     min, max =', f.min(),     f.max()
# 			print '## angle: min, max =', angle.min(), angle.max()
# 			print '##==================================================='

        #-----------------------------------------
        # Check where momentum is negative
        # Can this be avoided with smaller dt ??
        #--------------------------------------------------------------
        # For Treynor, starting from d=0 everywhere, a large fraction
        # of cells have M < 0 in the beginning.  Later on, there are
        # two cells on ridges (no D8 kids) for which M is slightly < 0,
        # with dt = 0.05 sec.  Both of them lie outside of the main
        # basin -- (col,row) = (4,19) & (22,23).
        #--------------------------------------------------------------
        # Other trouble spots early on include (col, row) =
        #   (8,38), (10,41), (19,25), (6,27), (5,29), (22,24)
        #--------------------------------------------------------------
# 		wb = np.where( M < 0 )
# 		nb = wb[0].size
# 		if (nb > 0):
# 			print 'WARNING: Total momentum < 0 in ' + str(nb) + ' cells.'
# 			if (nb < 5):
# 				print '## cols =', wb[1]
# 				print '## rows =', wb[0]
# 			print '## M:     min, max =', M.min(),     M.max()
# 			print '## grav:  min, max =', grav.min(),  grav.max()
# 			print '## fric:  min, max =', fric.min(),  fric.max()
# 			print '## Qtrm:  min, max =', Qtrm.min(),  Qtrm.max()
# 			print '## uQ_in: min, max =', uQ_in.min(), uQ_in.max()
# 			print '##==================================================='

        #-----------------------------------------------
        # Don't allow losses to make momentum negative
        #-----------------------------------------------
        ###########################  
		M = np.maximum(M, 0.0)
        ###########################  

        #----------------------------------------------
        # Compute new velocity as: momentum / (rho V)
        #----------------------------------------------
		u2 = M * Vinv

		#----------------------
		# Check the u2 values
		#----------------------
# 		w1 = np.where(u2 < 0)
# 		n1 = np.size(w1[0])
# 		print '### Number of negative values in u2 =', n1

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
		self.u[:] = u2

		## print '## u:     min, max =', u.min(),    u.max()
        
    #   update_velocity() 
    #-------------------------------------------------------------------
#     def update_velocity0(self):
# 
#         #--------------------------------------------------------
#         # Note:  This update to u uses values of u, d, Q, and
#         #        S_free from the last time step, then u is
#         #        overwritten with new values.
#         #--------------------------------------------------------
#         # Note:  update_free_surface_slope() is called in
#         #        channels_base.update() before calling this.
#         #--------------------------------------------------------
#         # Note:  For speed, multiply scalars first, then grids.
#         #        u, d and Q are always grids
#         #--------------------------------------------------------
# 
# 		#------------------------------------------
# 		# Start with the interior (nonflux) terms
# 		#------------------------------------------
# 		# When (n eq 0), should have (d gt 0) and
# 		# (u eq 0), so (grav gt 0) and (acc gt 0)
# 		# but fric = Atrm = Rtrm = 0.
# 		#------------------------------------------
# 		angle  = self.angle
# 		width  = self.width
# 		vol    = self.vol
# 		d      = self.d       # (need for Pw)
# 		u      = self.u
# 		Q      = self.Q
# 		f      = self.f
# 
# 		## d[ self.d8.noflow_IDs ] = 0.0
# 
# 		#----------------------------------------
# 		# Update inverse of channel flow volume
# 		#------------------------------------------------
# 		# Wherever vol = 0, we should have u = 0, Q = 0.
# 		#------------------------------------------------
#         # Note: Subscripting with Boolean arrays works
#         #       even if all values are False.
# 		#------------------------------------------------
# 		Vinv = vol.copy()
# 		wg = (vol > 0)         # (Boolean array)
# 		### wg = (vol > 1e-2)         # (Boolean array)
# 		wb = np.invert( wg )   # (opposite Boolean array)
# 		Vinv[ wg ] = 1 / vol[ wg]
# 		Vinv[ wb ] = 0.0
# 
# 		#------------------------------------
# 		# Compute the areas that are needed
# 		#------------------------------------------------------------
# 		# Note:  Ap = "wetted surface area", which is not the same
# 		#        as "wetted area" used to compute hydraulic radius.
# 		#------------------------------------------------------------
# 		# Note:  (Atop / wtop ) = (Ap / Pw) = d8.ds,
# 		#        so (Pw / wtop) = (Ap / Atop).
# 		#------------------------------------------------------------
# 		## wtop = width + (2 * d * np.tan(angle))   # (top width)
# 		## Atop = self.d8.ds * wtop                 # (top area)
# 		## A1 = self.d8.ds * width                  # (bottom area)
#         ## Acell = self.da                          # (grid cell area)
#         #-----------------------------------------------
# 		Pw = width + (np.float64(2) * d / np.cos(angle))  # (wetted perimeter)
# 		Ap = self.d8.ds * Pw                              # (wetted surf. area)
# 
# 		#-----------------------
# 		# Compute the "R term"
# 		#-----------------------
# 		Rtrm = (u * self.R) * (self.da * Vinv)
# 
# 		#-------------------------------------
# 		# Compute gravity and friction terms
# 		#--------------------------------------------------------
# 		# Only the grav term will contribute at first, when d=0
# 		#--------------------------------------------------------
# 		grav = self.g * self.S_free
# 		fric = (f * u * u * Ap * Vinv)
# 		### fric = np.minimum( fric, grav )   ####################
# 
# 		#----------------------------------------------------
# 		# Compute the "acceleration" (positive or negative)
# 		#----------------------------------------------------
# 		## du = (grav - fric - Rtrm) * self.dt  
# 		du = (grav - fric) * self.dt    # (for testing)  ###############################
# 
# 		#-------------------------------------------------
# 		# Compute the new flow velocity just due to
# 		# contributions from "interior" (nonflux) terms.
# 		# We'll modify this with flux terms.
# 		#-------------------------------------------------
# 		u2 = u + du
# 		### u2  = u.copy() + du    
# 
# # 		print '## Q:     min, max =', Q.min(),     Q.max()
# # 		print '## vol:   min, max =', vol.min(),   vol.max()
# # 		print '## d:     min, max =', d.min(),     d.max()
# # 		print '## u:     min, max =', u.min(),     u.max()
# # 		print '## f:     min, max =', f.min(),     f.max()
# # 		print '## grav:  min, max =', grav.min(),  grav.max()
# # 		print '## fric:  min, max =', fric.min(),  fric.max()
# # 		print '## Rtrm:  min, max =', Rtrm.min(),  Rtrm.max()
# # 		print '## angle: min, max =', angle.min(), angle.max()
# # 		print '## du:    min, max =', du.min(),    du.max()
# # 		print '##==================================================='
# 
# 		#----------------------------------------
# 		# Now prepare to compute flux terms
#         # which also contribute to the velocity
# 		#----------------------------------------
# 		fac = Vinv * self.dt
# 
# 		#---------------------------------------
# 		# Note that "compound statements" with 
# 		# a semicolon are OK.
# 		#---------------------------------------
# 		p1 = self.d8.p1   ;  w1 = self.d8.w1
# 		p2 = self.d8.p2   ;  w2 = self.d8.w2
# 		p3 = self.d8.p3   ;  w3 = self.d8.w3
# 		p4 = self.d8.p4   ;  w4 = self.d8.w4
# 		p5 = self.d8.p5   ;  w5 = self.d8.w5
# 		p6 = self.d8.p6   ;  w6 = self.d8.w6
# 		p7 = self.d8.p7   ;  w7 = self.d8.w7
# 		p8 = self.d8.p8   ;  w8 = self.d8.w8
# 
# 		#-------------------------------------------
# 		# Add momentum fluxes from D8 child pixels
# 		#-------------------------------------------
# 		if (self.d8.p1_OK):
# 			u2[p1] += (u[w1] - u[p1]) * Q[w1] * fac[p1]
# 		if (self.d8.p2_OK):
# 			u2[p2] += (u[w2] - u[p2]) * Q[w2] * fac[p2]
# 		if (self.d8.p3_OK):
# 			u2[p3] += (u[w3] - u[p3]) * Q[w3] * fac[p3]   
# 		if (self.d8.p4_OK):
# 			u2[p4] += (u[w4] - u[p4]) * Q[w4] * fac[p4]  
# 		if (self.d8.p5_OK):
# 			u2[p5] += (u[w5] - u[p5]) * Q[w5] * fac[p5]   
# 		if (self.d8.p6_OK):
# 			u2[p6] += (u[w6] - u[p6]) * Q[w6] * fac[p6]    
# 		if (self.d8.p7_OK):
# 			u2[p7] += (u[w7] - u[p7]) * Q[w7] * fac[p7]   
# 		if (self.d8.p8_OK):
# 			u2[p8] += (u[w8] - u[p8]) * Q[w8] * fac[p8]    
# 
# 		#----------------------
# 		# Check the u2 values
# 		#----------------------
# # 		w1 = np.where(u2 < 0)
# # 		n1 = np.size(w1[0])
# # 		print '### Number of negative values in u2 =', n1
# 
# 		#--------------------------------
# 		# Don't allow u2 to be negative
# 		#----------------------------------------------
# 		# If uphill flow is allowed, to which pixel ?
# 		# This worked when d0 grid was used.
# 		#----------------------------------------------
#  		u2 = np.maximum(u2, np.float64(0))
# 
# 		#---------------
# 		# Copy u2 to u
# 		#---------------
# 		self.u[:] = u2
# 
# 		## print '## u:     min, max =', u.min(),    u.max()
#         
#     #   update_velocity0()                       
    #-------------------------------------------------------------------
#     def update_velocity_last(self):
# 
#         #--------------------------------------------------------
#         # Note:  This update to u uses values of u, d, Q, and
#         #        S_free from the last time step, then u is
#         #        overwritten with new values.
#         #--------------------------------------------------------
#         # Note:  update_free_surface_slope() is called in
#         #        channels_base.update() before calling this.
#         #--------------------------------------------------------
#         # Note:  For speed, multiply scalars first, then grids.
#         #        u, d and Q are always grids
#         #--------------------------------------------------------
# 
# 		#------------------------------------------
# 		# Start with the interior (nonflux) terms
# 		#------------------------------------------
# 		# When (n eq 0), should have (d gt 0) and
# 		# (u eq 0), so (grav gt 0) and (acc gt 0)
# 		# but fric = Atrm = Rtrm = 0.
# 		#------------------------------------------
# 		angle  = self.angle
# 		width  = self.width
# 		d      = self.d
# 		u      = self.u
# 		Q      = self.Q
# 		f      = self.f
# 
# 		## d[ self.d8.noflow_IDs ] = 0.0
# 
# 		#-------------------------------
# 		# Update inverse of flow depth
# 		#----------------------------------------------
# 		# Wherever d = 0, we should have u = 0, Q = 0.
# 		#----------------------------------------------
# 		dinv = d.copy()
# 		wg = np.where( d > 0 )
# 		ng = np.size( wg[0] )
# 		if (ng > 0):
# 			dinv[ wg ] = 1 / d[ wg ]
# 		#-----------------------------
# 		wb = np.where( d <= 0 )
# 		nb = np.size( wb[0] )
# 		if (nb > 0):
# 			dinv[ wb ] = 0.0
# 		#-----------------------------
# # 		wz = np.where( d <= 0 )
# # 		nz = np.size( wz[0] )
# # 		dinv = (1 / d)
# # 		if (nz > 0):
# # 			dinv[ wz] = 0
# # 			## u[ wz ]   = 0
# # 			## Q[ wz ]   = 0
# # 		-------------------------------------------
# # 		self.d_is_pos = (d > 0)
# # 		self.d_is_neg = np.invert( self.d_is_pos )
# 
# 		#-----------------------------------------
# 		# Compute the "R term" first (need Atop)
# 		#-----------------------------------------
# 		wtop = width + (2 * d * np.tan(angle))     # (top width)
# 		Atop = self.d8.ds * wtop                   # (top area)
# 		Rtrm = (u * self.R) * (self.da / Atop)     # (da = pixel area)
# 		Rtrm = Rtrm * dinv
# 
# 		#-----------------------
# 		# Compute the "A term"
# 		#------------------------------------------------------------
# 		# Note:  A2 = "wetted surface area", which is not the same
# 		#        as "wetted area" used to compute hydraulic radius.
# 		#------------------------------------------------------------
# 		Pw = width + (np.float64(2) * d / np.cos(angle))  # (wetted perimeter)
# 		A2 = self.d8.ds * Pw                              # (wetted surf. area)
# 		Atrm = (u * Q) * ((np.float64(1) / Atop) - (np.float64(1) / A2))
# 		Atrm = Atrm * dinv
# 
# 		#---------------------------------------------------
# 		# Compute gravity and friction terms (divided by d)
# 		#--------------------------------------------------------
# 		# Only the grav term will contribute at first, when d=0
# 		#--------------------------------------------------------
# 		grav = self.g * self.S_free
# 		fric = (f * u * u * dinv)
# 		## fric = f * (u ** np.float64(2)) * dinv
# 
# 
# 		#----------------------------------------------------
# 		# Compute the "acceleration" (positive or negative)
# 		#----------------------------------------------------
# 		## acc = (grav - fric)   ## for testing
# 		## acc = (grav - fric - Rtrm)
# 		acc = (Atrm - Rtrm + grav - fric)  #################################
# 		### acc = (grav + Atrm - fric - Rtrm)
# 
# # 		print '## d:     min, max =', d.min(),    d.max()
# # 		print '## u:     min, max =', u.min(),    u.max()
# # 		print '## f:     min, max =', f.min(),    f.max()
# # 		print '## grav:  min, max =', grav.min(), grav.max()
# # 		print '## fric:  min, max =', fric.min(), fric.max()
# # 		print '## Rtrm:  min, max =', Rtrm.min(), Rtrm.max()
# # 		print '## acc:   min, max =', acc.min(),  acc.max()
# # 		print '##==================================================='
# 
# 		#-------------------------------------------------
# 		# Compute the new flow velocity just due to
# 		# contributions from "interior" (nonflux) terms.
# 		# We'll modify this with flux terms.
# 		#-------------------------------------------------
# 		## du = (self.dt * acc) 
# 		u2  = u + (self.dt * acc)     # (before next part)
# 
# 		#------------------------------------
# 		# Now prepare to compute flux terms
# 		#------------------------------------
# 		# Compute new velocity, (u + du)
# 		#---------------------------------------------
# 		# Note:  (Atop / wtop ) = (A2 / Pw) = d8.ds,
# 		#        so (Pw / wtop) = (A2 / Atop).
# 		#---------------------------------------------
# 		# Note: fac and uu will always be grids.
# 		#---------------------------------------------
# 		fac = (self.dt / A2) * dinv
# 		uu  = u * (Pw / wtop)
# 		## uu  = u2 * (Pw / wtop)   ######## Use this instead ??
# 
# 		#---------------------------------------
# 		# Note that "compound statements" with 
# 		# a semicolon are OK.
# 		#---------------------------------------
# 		p1 = self.d8.p1   ;  w1 = self.d8.w1
# 		p2 = self.d8.p2   ;  w2 = self.d8.w2
# 		p3 = self.d8.p3   ;  w3 = self.d8.w3
# 		p4 = self.d8.p4   ;  w4 = self.d8.w4
# 		p5 = self.d8.p5   ;  w5 = self.d8.w5
# 		p6 = self.d8.p6   ;  w6 = self.d8.w6
# 		p7 = self.d8.p7   ;  w7 = self.d8.w7
# 		p8 = self.d8.p8   ;  w8 = self.d8.w8
# 
# 		#-------------------------------------------
# 		# Add momentum fluxes from D8 child pixels
# 		#-------------------------------------------
# 		if (self.d8.p1_OK):
# 			u2[p1] += (u[w1] - uu[p1]) * Q[w1] * fac[p1]
# 		if (self.d8.p2_OK):
# 			u2[p2] += (u[w2] - uu[p2]) * Q[w2] * fac[p2]
# 		if (self.d8.p3_OK):    
# 			u2[p3] += (u[w3] - uu[p3]) * Q[w3] * fac[p3]
# 		if (self.d8.p4_OK):    
# 			u2[p4] += (u[w4] - uu[p4]) * Q[w4] * fac[p4]
# 		if (self.d8.p5_OK):    
# 			u2[p5] += (u[w5] - uu[p5]) * Q[w5] * fac[p5]
# 		if (self.d8.p6_OK):    
# 			u2[p6] += (u[w6] - uu[p6]) * Q[w6] * fac[p6]
# 		if (self.d8.p7_OK):    
# 			u2[p7] += (u[w7] - uu[p7]) * Q[w7] * fac[p7]
# 		if (self.d8.p8_OK):    
# 			u2[p8] += (u[w8] - uu[p8]) * Q[w8] * fac[p8]
# 
# #         #-------------------------------------------
# #         # Add momentum fluxes from D8 child pixels
# #         #-------------------------------------------
# #         if (self.d8.p1_OK):
# #             ## print '## u[w1]: min, max =', u[w1].min(), u[w1].max()  #######
# #             u2[p1] += u[w1] * Q[w1] * fac[p1]
# #         if (self.d8.p2_OK):
# #             u2[p2] += u[w2] * Q[w2] * fac[p2]
# #         if (self.d8.p3_OK):    
# #             u2[p3] += u[w3] * Q[w3] * fac[p3]
# #         if (self.d8.p4_OK):    
# #             u2[p4] += u[w4] * Q[w4] * fac[p4]
# #         if (self.d8.p5_OK):    
# #             u2[p5] += u[w5] * Q[w5] * fac[p5]
# #         if (self.d8.p6_OK):    
# #             u2[p6] += u[w6] * Q[w6] * fac[p6]
# #         if (self.d8.p7_OK):    
# #             u2[p7] += u[w7] * Q[w7] * fac[p7]
# #         if (self.d8.p8_OK):    
# #             u2[p8] += u[w8] * Q[w8] * fac[p8]
# 
# 		#--------------------------------
# 		# Don't allow u2 to be negative
# 		#----------------------------------------------
# 		# If uphill flow is allowed, to which pixel ?
# 		# This worked when d0 grid was used.
# 		#----------------------------------------------
# 		u2 = np.maximum(u2, np.float64(0))
# 
# 		#---------------
# 		# Copy u2 to u
# 		#---------------
# 		self.u[:] = u2
# 		### self.u = u2
#         
#     #   update_velocity_last()                       
    #-------------------------------------------------------------------


         
