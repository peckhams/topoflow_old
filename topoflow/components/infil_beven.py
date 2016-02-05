
## Copyright (c) 2009-2013, Scott D. Peckham
## May 2009
## May 2010 (changes to unit_test())

####################################
#   THIS IS NOT READY YET.
####################################

#-----------------------------------------------------------------------
#  NOTES:  This file defines a Beven (Exponential K) infiltration
#          component and related functions.  It inherits from the
#          infiltration "base class" in "infil_base.py".

#-----------------------------------------------------------------------
#
#  class infil_beven_exp_K
#
#      get_attribute()       # (10/26/11)
#      update_infil_rate()

#  Functions:
#      Beven_Exp_K_Infil_Rate_v1  (May 2006)
#      Beven_Exp_K_Infil_Rate_1D  (Not written)
#      Beven_Exp_K_Infil_Rate_3D  (Not written)
#      Beven_Exp_K_Infil_Rate_v2  (Not written)

#-----------------------------------------------------------------------

import numpy as np

# import sys, os, time, glob, platform

from topoflow.components import infil_base

#-----------------------------------------------------------------------
class infil_beven( infil_base.infil_component ):

    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        map = {'comp_name':          'InfilBeven',
               'version':            '3.1',
               'model_name':         'Infiltration_Beven',
               'model_family':       'TopoFlow',
               'cfg_template_file':  'Infil_Beven.cfg.in',
               'cfg_extension':      '_infil_beven.cfg',
               'cmt_var_prefix':     '/InfilBeven/Input/Var/',
               'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Infil_Beven.xml',
               'dialog_title':       'Infiltration: Beven Method Parameters',
               'time_step_type':     'fixed',
               'time_units':         'seconds',
               'mesh_type':          'none',
               'author_name':        'Scott Peckham'}

        try:
            return map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

    #   get_attribute() 
    #-------------------------------------------------------------------
    def update_infil_rate(self, P, SM, ET):

        #-----------------------------------------------------
        # This function is not totally correct but is stable
        # and gives results similar to the correct method.
        #-----------------------------------------------------
        r = (P + SM)
        self.IN = Beven_Exp_K_Infil_Rate_v1(self, r)

        #-----------------------------------------------------
        # Numerically more correct method (??) but unstable
        #-----------------------------------------------------        
        ## self.IN = Beven_Exp_K_Infil_Rate_v2(self, r, r_last, n)
   
    #   update_infil_rate()
    #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def Beven_Exp_K_Infil_Rate_v1(iv, r):

    #------------------------------------------------------------
    #Notes:  IN  = infiltration rate [m/s]
    #        Ks  = saturated hydraulic conductivity [m/s]
    #        Ki  = initial hydraulic conductivity [m/s]
    #        qs  = soil moisture content [dimless]
    #        qi  = soil initial moisture content [dimless]
    #        C   = "storage suction factor" [m]
    #            = (delta_psi * delta_theta) > 0
    #        f   = parameter in K*(z) = K0 * exp(f*z) [1/m]
    #              K* < Ks, is "eff. K behind wetting front"
    #              (f < 0, roughly between -1 and -13)
    #              Is f roughly equal to (-gamma/G) in the
    #              Smith-Parlange method?  (Compare denoms.)
    #              f is sometimes written as (del_theta/m).
    #        K0  = K*(0) = coeff. in previous equation [m/s]
    #         I  = cum. infiltration depth (since reset) [m]
    #         P  = precipitation rate [m/s]
    #        SM  = snowmelt rate [m/s]
    #         r  = (P + SM)  [m/s]

    #        This comes from: Beven, K. (1984) "Infiltration
    #        into a class of vertically non-uniform soils",
    #        Hydrol. Sciences J., 29(4), 425-434.

    #        Note that the infiltration rate has a max possible
    #        value of (P + SM).  It seems that this method does
    #        not asymptote to Ks as I increases.  Can we simply
    #        add Ks?  We should reset I between "events", but
    #        haven't decided how to do this yet.

    #        Total infiltrated depth, I, is incremented in the
    #        calling function, called Infiltration.
    #------------------------------------------------------------

    #--------------------------------------------
    #Note that Ks is used to store K* and c is
    #different than c used for Richards' method.
    #Need to add f to the set of infil vars.
    #Also need to add GUI to collect all vars.
    #Note that t1 and t3 are both < 0 & c > 0.
    #--------------------------------------------
    dq = (iv.qs - iv.qi)
    t1 = (iv.Ks * (iv.f)) / dq
    t2 = (iv.c * (iv.I))
    t3 = (np.float64(1) - np.exp(-float64(1) * (iv.f) * (iv.I) / dq))
    IN = t1 * t2 / t3
    
    #-------------------------------------------------
    #Initially, IN = r, and all of the incoming
    #water infiltrates.  IN cannot exceed r.
    #Ponding time, Tp, is time until (IN lt r).
    #-------------------------------------------------
    IN = np.minimum(IN, r)
    
    #-----------------------------------
    #Is P_total less than Ks anywhere ?
    #If so, set IN = P_total there.
    #-----------------------------------
    iv, IN, r = Check_Low_Rainrate(iv, IN, r)
    
    return IN
    
#   Beven_Exp_K_Infil_Rate_v1
#-----------------------------------------------------------------------



