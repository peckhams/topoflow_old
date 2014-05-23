#! /usr/bin/env python
#
# Copyright (c) 2001-2013, Scott D. Peckham
#
#-----------------------------------------------------------------------

import numpy as np
import os
import tempfile
# See:  http://docs.python.org/2/library/tempfile.html

from topoflow.framework import framework
## from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
#
#  topoflow_test()    # Use framework to run TopoFlow.
#  erode_test()
#
#  ref_test()         # For passing references between components.
#
#  framework_test1()  # Test some basic framework functions.
#  framework_test2()
#
#  bobs_erode_test()  # Use framework to run Erode.
#
#-----------------------------------------------------------------------
def topoflow_test( driver_port_name='hydro_model',
                   cfg_prefix=None, cfg_directory=None,
                   time_interp_method='Linear'):

    #-----------------------------------------------------
    # Note: The "driver_port_name" defaults to using a
    #       component of "hydro_model" type as the driver.
    #       The component of this type specified in the
    #       provider_file will be the driver component.
    #
    #       Any other component in the provider_file can
    #       also be used as the driver.  Examples are:
    #          meteorology
    #          channels
    #          snow
    #          satzone
    #          evap
    #          infil
    #          diversions
    #          ice
    #-----------------------------------------------------
    
    #----------------------------------------------
    # Get path to the current file (framework.py)
    # At top need: "#! /usr/bin/env python" ??
    #----------------------------------------------
    paths = framework.get_package_paths()
    framework_dir = paths['framework']
    examples_dir  = paths['examples']
    
    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'June_20_67'
    if (cfg_directory == None):
        cfg_directory = examples_dir + 'Treynor_Iowa/'

    #-------------------------------------------------------
    # (2/6/13) Since the framework runs the clock now, do
    # we still need to specify a "driver_port" ??
    # Might still be necessary for use in CSDMS framework.
    #-------------------------------------------------------
    f = framework.framework()

    # No longer needed or used.
    ## f.repo_dir = framework_dir  #### repo_dir vs. repo_path ??
    
    #------------------------------
    # Run the full TopoFlow model
    #------------------------------
    f.run_model( driver_port_name=driver_port_name,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

    #----------------
    # Run the model
    #----------------
##    f.run_model_old( driver_name=driver_name,
##                     cfg_prefix=cfg_prefix,
##                     cfg_directory=cfg_directory,
##                     time_interp_method=time_interp_method)

#   topoflow_test()
#-----------------------------------------------------------------------
def erode_test( cfg_prefix=None, cfg_directory=None,
                time_interp_method='Linear'):

    #----------------------------------------------
    # Get path to the current file (framework.py)
    # At top need: "#! /usr/bin/env python" ??
    #----------------------------------------------
    paths = framework.get_package_paths()
    framework_dir = paths['framework']
    examples_dir  = paths['examples']
    
    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'Test'
    if (cfg_directory == None):
        cfg_directory = examples_dir + 'Erode_Test'
            
    f = framework.framework()
    driver_port_name = 'LEM'

    # No longer needed or used.
    ## f.repo_dir = framework_dir  #### repo_dir vs. repo_path ??

    #----------------------
    # Run the Erode model
    #----------------------   
    f.run_model( driver_port_name=driver_port_name,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   erode_test()
#-----------------------------------------------------------------------
def ref_test():

    #---------------------------------------------------------
    # Notes: See "initialize_scalar()" method in BMI.base.py
    #        which is based on this test.
    #---------------------------------------------------------
    import random
    
    class comp1():
        def initialize(self):
            #---------------------------------
            # Case 1. 0-d numpy array.
            # Mutable, if done carefully.
            # Prints without array brackets.
            #---------------------------------
            self.x = np.array(1.0, dtype='float64')

            #--------------------------------------
            # Case 2. 1-d numpy array, 1 element.
            # Mutable, if done carefully.
            # Prints with array brackets.
            #--------------------------------------
            ## self.x = np.array([1.0], dtype='float64')

            #------------------------------------
            # Case 3. numpy scalar;  immutable.
            #------------------------------------           
            ## self.x = np.float64(1)
            
        #-------------------------------------------------
        def update(self):
            #--------------------------
            # These work with Case 1.
            #--------------------------
            # self.x += 1
            # self.x += random.random()
            # self.x *= random.random()
            self.x.fill( random.random() )
            
            #----------------------------
            # Doesn't work with Case 1.
            #----------------------------            
            # self.x = self.x + 1
            # self.x = random.random()
            # self.x[:] = random.random()   #(can't slice 0D)
            # setattr( self, 'x', random.random() )
            #----------------------------------------------------------
            # new_x = np.array( random.random(), dtype='float64' )
            # setattr( self, 'x', new_x )
            #----------------------------------------------------------
            # new_x = np.array( random.random(), dtype='float64' )
            # setattr( self, 'x[:]', new_x )
            
            #----------------------------
            # This works with Case 2.
            #----------------------------            
            # self.x[0] = random.random()
            # self.x[0] += 1  # (works)
            # self.x[:] = random.random()  # (works for 1D, 2D, etc.)
            
            ## self.x = self.x + 1  # (ref is lost; not in place)
        #-------------------------------------------------
        def get_0d_double(self, var_name):
            return getattr( self, var_name )
            #-----------------------------------
            # Using exec like this also works.
            #-----------------------------------
            # exec("result = self." + var_name) in globals(), locals()
            # return result
            
            #------------------------------------------
            # These next three "break" the reference.
            #------------------------------------------
            # return np.float64( result )
            # return result.astype('float64')
            # return np.array( result, dtype='float64' )

    #-----------------------------------------------------
    class comp2():
        def print_x(self):
            print 'After c1.update(), c2.x =', self.x
        #-------------------------------------------------
        def set_0d_double(self, var_name, scalar):
            setattr( self, var_name, scalar )
            #-----------------------------------
            # Using exec like this also works.
            #-----------------------------------
            # exec("self." + var_name + " = scalar") in globals(), locals()
    #-----------------------------------------------------
        
    #---------------------------
    # Instantiate 2 components
    #---------------------------
    c1 = comp1()
    c2 = comp2()
    c1.initialize()
    
    #-------------------------------------
    # Copy reference to x from c1 to c2.
    #-------------------------------------
    vector = c1.get_0d_double('x')
    c2.set_0d_double('x', vector)
    #------------------------------------------    
    # c2.x = c1.x   # (works)
    # exec("c2.x = c1.x") in globals(), locals()  # (works)
    # exec("c2.x = c1.x")       # (problem?)
    # exec("c2." + "x = c1.x")  # (problem?)
    
    for k in xrange(3):
        c1.update()
        c2.print_x()
    
#   ref_test()
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def framework_test1():
    
    #----------------------------------------------
    # Get path to the current file (framework.py)
    # At top need: "#! /usr/bin/env python" ??
    #----------------------------------------------
    paths = framework.get_package_paths()
    framework_dir = paths['framework']
    examples_dir  = paths['examples']

    f = framework.framework()

    # No longer needed or used.
    ## f.repo_dir = framework_dir  #### repo_dir vs. repo_path ??
    ## repo_file = framework_dir + 'component_repository.xml'
    
    ## driver_port_name = 'hydro_model'  # (TopoFlow Driver test)

    f.read_repository( SILENT=False )

    #-----------------------------------------
    # Set the working directory for test run
    #-----------------------------------------
    cfg_directory = examples_dir + 'Treynor_Iowa/'
    os.chdir( cfg_directory )     ############ IS THIS NEEDED?
    cfg_prefix    = 'June_20_67'
    
    print 'Components in repository:'
    for comp_name in f.repo_list:
        print '   ' + comp_name
    print ' '

    tf_comp_info = f.comp_info[ 'topoflow_driver' ]  
    print 'Checking some info for "topoflow_driver":'
    print '    comp_name  =', tf_comp_info.comp_name
    print '    model_name =', tf_comp_info.model_name
    print '    uses_ports =', tf_comp_info.uses_ports
    print ' '

    #--------------------------------------------
    # Instantiate a complete set of components.
    #--------------------------------------------
    # Now the instantiate() method only allows
    # one component of each "type" (port_name).
    #--------------------------------------------
    comp_list = ['tf_meteorology', 
                 'tf_channels_kin_wave',
                 'tf_infil_green_ampt',
                 'tf_evap_priestley_taylor',
                 'tf_snow_degree_day',
                 'tf_ice_gc2d',
                 'tf_satzone_darcy_layers',
                 'tf_diversions_fraction_method',
                 'topoflow_driver' ]
    for comp_name in comp_list:
        f.instantiate( comp_name, SILENT=False )
    print ' '
    print 'ALL_PYTHON =', f.ALL_PYTHON
    print ' '

    #------------------------------------------
    # Check if all uses ports have a provider
    #------------------------------------------
    COMPLETE = f.comp_set_complete( SILENT=False, REPORT=True )
    
    #------------------------------------------------------
    # Try to instantiate another with port_name == snow.
    # It should issue an error message.
    # Also, Treynor_Iowa doesn't have this CFG file yet.
    #------------------------------------------------------
    f.instantiate( 'tf_snow_energy_balance', SILENT=False )

    #------------------------------------------------
    # Remove existing snow component and try again.
    #------------------------------------------------
    f.remove( 'tf_snow_degree_day', SILENT=False )
    f.instantiate( 'tf_snow_energy_balance', SILENT=False )
    
    #-----------------------------------
    # Instantiate all components and
    # store instances in self.comp_set.
    #-----------------------------------
##    skip_list = ['tf_data_his', 'erode_d8_global',
##                 'erode_d8_local', 'd8_global', 'dem_smoother']
##    for comp_name in f.repo_list:
##        if (comp_name not in skip_list):
##            f.instantiate( comp_name )
##            print 'Instantiated component named:', comp_name
##    print ' '
    #-------------------------------------------------------------
##    try:
##        for comp_name in f.repo_list:
##            f.instantiate( comp_name )
##            print 'Instantiated component named:', comp_name
##    except:
##        print 'ERROR: Could not instantiate component.'
##    print ' '
 
    #-----------------------------
    # Initialize some components
    #-------------------------------------------------------------
    # (3/15/12) The TopoFlow components currently call a method:
    # "initialize_required_components()" that in turn calls:
    # "initialize_ports()", inherited from CSDMS_base.py.
    # But this won't be done with the new approach.
    #-------------------------------------------------------------
    f.initialize( 'meteorology', cfg_prefix )
    print 'Initialized component of type: meteorology.'
    f.initialize( 'snow', cfg_prefix )
    print 'Initialized component of type: snow.'   
    print ' '

    #--------------------------------------------
    # Copy references from Meteorology to Snow.
    #--------------------------------------------
    print 'Copying references from Meteorology to Snow component...'
    provider_name = 'meteorology'
    user_name     = 'snow'
    var_name_list = ['air__temperature',
                     'land_surface__net_irradiation_flux',
                     'land_surface__temperature',
                     'water__density' ]

    for long_var_name in var_name_list:
        f.connect( provider_name, user_name,
                   long_var_name )
    
    #--------------------------------------------
    # Copy references from Snow to Meteorology.
    #--------------------------------------------
    print 'Copying references from Snow to Meteorology component...'
    provider_name = 'snow'
    user_name     = 'meteorology'
    var_name_list = [ 'snow__density',
                      'snow__depth',
                      'snow__liquid_equivalent_depth',
                      'snow__melt_rate' ]
    for long_var_name in var_name_list:
        f.connect( provider_name, user_name,
                   long_var_name )

    #---------------------------------------
    # Check whether the references worked.
    #---------------------------------------
    snow_comp = f.comp_set['snow']
    met_comp  = f.comp_set['meteorology']
    print ' '
    print 'From meteorology, snow got: density of water =', \
          snow_comp.rho_H2O
          ### snow_comp.met_vars.rho_H2O
    print 'From snow, meteorology got: density of snow  =', \
          met_comp.rho_snow
          ### met_comp.snow_vars.rho_snow
    
    #----------------
    # Final message
    #----------------
    print 'Finished with framework_test1.'
    print ' '
    
#   framework_test1()
#-----------------------------------------------------------------------
def framework_test2():

    #----------------------------------------------
    # Get path to the current file (framework.py)
    # At top need: "#! /usr/bin/env python" ??
    #----------------------------------------------
    paths = framework.get_package_paths()
    framework_dir = paths['framework']
    examples_dir  = paths['examples']

    f = framework.framework()

    # No longer needed or used.
    ## f.repo_dir = framework_dir  #### repo_dir vs. repo_path ??
    ## repo_file = framework_dir + 'component_repository.xml'
        
    ## driver_port_name = 'hydro_model'  # (TopoFlow Driver test)

    f.read_repository( SILENT=False )

    #-----------------------------------------
    # Set the working directory for test run
    #-----------------------------------------
    cfg_prefix    = 'June_20_67'
    cfg_directory = examples_dir + 'Treynor_Iowa/'
    os.chdir( cfg_directory )     ############ IS THIS NEEDED?
    f.cfg_prefix = cfg_prefix  ## (needed by initialize())
    
    print 'Components in repository:'
    for comp_name in f.repo_list:
        print '   ' + comp_name
    print ' '
 
    #-----------------------------------------------------
    # Set self.comp_set_list and self.provider_list
    # from info in the provider file.
    #-----------------------------------------------------
    f.provider_file = (cfg_prefix + '_providers.txt')
    f.read_provider_file()
    
    #--------------------------------------------
    # Instantiate a complete set of components.
    #--------------------------------------------
    # Now the instantiate() method only allows
    # one component of each "type" (port_name).
    #--------------------------------------------    
    for comp_name in f.comp_set_list:
        f.instantiate( comp_name, SILENT=False )

    #---------------------------------------------
    # Try to automatically connect every user to
    # a corresponding provider in the comp_set.
    #---------------------------------------------------
    # Provider components are initialized in the order
    # of provider_list and then set references in each
    # component that uses one or more of their vars.
    #---------------------------------------------------
    OK = f.initialize_and_connect_comp_set( REPORT=True )
    if not(OK):
        return
    
    #------------------------------------------------
    # Call update() for each provider, in order.
    # Don't worry about reconciling time steps yet.
    #------------------------------------------------
    print ' '
    # f.update_all()  # (not ordered)

##    for k in xrange(10):
##        for port_name in f.provider_list:
##            print 'Calling update() for port_name =', port_name
##            if (port_name != 'hydro_model'):
##                bmi = f.comp_set[ port_name ]
##                bmi.update( -1.0 )
            
    #----------------------------------------------
    # Call the driver component's "update" method
    #----------------------------------------------
##    print ' '
##    print "Calling driver's update() method..."
##    f.update( 'hydro_model' )
    
    #----------------
    # Final message
    #----------------
    print 'Finished with framework_test2.'
    print ' '

#   framework_test2()
#-----------------------------------------------------------------------
def bobs_erode_test( driver_port_name='LEM',
                     cfg_prefix=None, cfg_directory=None,
                     time_interp_method='Linear'):
    
    #----------------------------------------------
    # Get path to the current file (framework.py)
    # At top need: "#! /usr/bin/env python" ??
    #----------------------------------------------
    paths = framework.get_package_paths()
    framework_dir = paths['framework']
    ## examples_dir  = paths['examples']
    
    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'Erode_Test_LCP'
    if (cfg_directory == None):
        cfg_directory = os.getenv("HOME") + '/Erode_Tests/Bob_LCP/'
        ## cfg_directory = '~/Erode_Tests/Bob_LCP/'  # (doesn't work)
    
    #-------------------------------------------------------
    # (2/6/13) Since the framework runs the clock now, do
    # we still need to specify a "driver_port" ??
    # Might still be necessary for use in CSDMS framework.
    #-------------------------------------------------------
    f = framework.framework()

    # No longer needed or used.
    ## f.repo_dir = framework_dir  #### repo_dir vs. repo_path ??
    
    #------------------------------
    # Run the full TopoFlow model
    #------------------------------
    f.run_model( driver_port_name=driver_port_name,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   bobs_erode_test()
#-----------------------------------------------------------------------


