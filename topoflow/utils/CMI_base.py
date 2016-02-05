#
# NOTE:  This is not currently used.  See emeli.py in frameworks folder.
#      
# Copyright (c) 2010-2012, Scott D. Peckham
#
# February 2012  (Complete CMI, starting from CSDMS_base.py.)
#
#-----------------------------------------------------------------------
#
#  Notes:  This file defines a "base class" with a CMI (Component
#          Model Interface) for CSDMS "process" components.  These
#          methods allow a CSDMS component written in Python to
#          communicate with other CSDMS components and to be used
#          within a CSDMS/CCA framework.
#
#          CMI methods often call BMI (Basic Model Interface) methods.
#          The BMI interface is designed to be completely framework
#          independent, and therefore contains no methods related to
#          CCA framework concepts such as ports.
#
#          Some "private" utility methods are defined at the end.
#
#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class CMI_component
#
#      __init__()
#
#      --------------------------------------------
#      CMI methods related to ports and framework
#      --------------------------------------------
#      set_services()         # (not a CMI method)
#      release_services()     # (not a CMI method)
#      get_cca_port_info()
#      get_cca_ports()
#      release_cca_ports()
#
#      -----------------------------------
#      CMI/BMI methods to get model info
#      -----------------------------------
#      get_status()
#      set_status()
#      get_mode()
#      set_mode()
#      get_attribute()   # (component and model attributes)
#
#      --------------------------------------
#      CMI/BMI methods to get variable info
#      --------------------------------------
#      get_values()
#      set_values()
#      get_values_at_indices()
#      set_values_at_indices()
#      -------------------------
#      get_input_var_names()
#      get_output_var_names()
#      -------------------------
#      get_var_rank()
#      get_var_type()
#      get_var_units()
#
#      ----------------------------------
#      CMI/BMI methods to get grid info
#      ----------------------------------
#      get_grid_shape()
#      get_grid_spacing()
#      get_grid_lower_left_corner()
#      get_grid_attribute()                 ### NEW CMI METHOD ??
#
#      ----------------------------------
#      CMI/BMI methods to get time info
#      ----------------------------------
#      get_time_step()
#      get_time_units()
#      get_time()   (current, start, or stop)
#          get_start_time()
#          get_current_time()
#          get_end_time()
#
#      --------------------------------------
#      CMI methods for fine-grained control
#      --------------------------------------
#      initialize()
#      run_until()
#      run_for()
#      finalize()
#      ------------------
#      run_model()          # (not required for CMI)
#      go()
#
#      ------------------------------
#      Private methods or utilities
#      ------------------------------
#      check_finished()     # (called by run_model(), update().)
#      print_traceback()
#
#-----------------------------------------------------------------------

# import BMI_base

from topoflow.components import topoflow_driver

import sys
import traceback

# import numpy
# import os
# import time

## import edu.csdms.models.IRFPort         # (in IMPL file only)
## import edu.csdms.tools.ConfigDialog
## import edu.csdms.tools.IRFPortQueue
## import edu.csdms.tools.PrintQueue
## import edu.csdms.tools.TemplateFiles

#-----------------------------------------------------------------------
def unit_test():

    c = CMI_component()
    print 'Instantiated CMI component.'

#   unit_test()
#-----------------------------------------------------------------------
### class CMI_component( BMI_base.BMI_component ):

class CMI_component:

    #-------------------------------------------------------------------
    def __init__(self):

        self.bmi = topoflow_driver.topoflow_driver()
        
##        ## self.DEBUG    = True
##        self.DEBUG       = False
##        self.SKIP_ERRORS = False
##        self.SILENT      = True        # (new default: 11/16/11) 
##        self.REPORT      = False
##        self.DONE        = False
##        self.status      = 'created'   # (OpenMI 2.0 conventions)
##
##        self.in_directory     = None
##        self.out_directory    = None
##        self.site_prefix      = None
##        self.case_prefix      = None
##        self.comp_status      = 'Enabled'
##        
##    #   __init__()
    #-------------------------------------------------------------------
    def set_services( self, services ):

        #-------------------------------------------------------
        # Notes: This is to be called by IMPL's setServices(),
        #        similar to how boccaSetServices() is called.
        #-------------------------------------------------------
        
        #-----------------------
        # Get model attributes
        #-----------------------
        comp_name    = self.bmi.get_attribute( 'comp_name' )
        gui_xml_file = self.bmi.get_attribute( 'gui_xml_file' )
        dialog_title = self.bmi.get_attribute( 'dialog_title' )
        self.mode = 'nondriver'  # (default)

        try:
            print 'SUCCESS:', comp_name, 'component status =', self.bmi.get_status()
        except:
            print 'FAILURE:', comp_name, 'component has not been instantiated yet.'
            return    # (do not return a value)

        #-------------------------------------------------------
        # Call the function that creates the tabbed-dialog GUI
        #-------------------------------------------------------
        # user_input is an object of type gov.cca.TypeMap
        #-------------------------------------------------------
        button_label = 'Configure'
        dialog_tool  = edu.csdms.tools.ConfigDialog.ConfigDialog()
        dialog_tool.make_dialog( services, gui_xml_file, button_label, dialog_title )
        self.user_input = dialog_tool.get_user_input()

        # self.bmi.get_user_input( services, self.d_services )


    #   set_services()
    #-------------------------------------------------------------------
    def release_services( self, services ):

        #--------------------------------------------------------
        # Notes: This is to be called by IMPL's setServices(),
        #        similar to how boccaReleaseServices is called.
        #        But we aren't doing anything in IMPL file yet.
        #--------------------------------------------------------
        pass
    
    #   release_services()
    #-------------------------------------------------------------------
    def get_cca_port_info(self):

        #-------------------------------------------------------
        # Read required CCA port names from a file?
        # This is done already by the PortQueue service class.
        #-------------------------------------------------------
        self.cca_port_names = ['']
        self.cca_port_short = ['']
        self.cca_port_type  = "IRFPort"
        self.cca_project    = "edu.csdms.models"
        ## self.cca_project    = "edu.csdms.ports"
        
    #   get_cca_port_info()
    #-------------------------------------------------------------------
    def get_cca_ports(self, d_services):

        #-------------------------------------------------------------
        # Example call in initialize method of "TopoFlow_Impl.py":
        #
        # OK = self.tf.get_cca_ports( self.d_services )
        # if not(OK): return -1
        #-------------------------------------------------------------
        # print 'CALLING get_cca_ports()...'

        ############################################        
        self.get_cca_port_info()
        ############################################
        port_names   = self.cca_port_names
        short_names  = self.cca_port_short
        port_type    = self.cca_port_type
        project_name = self.cca_project

        exec("import " + project_name)
        
        #--------------------------------------------
        # Does this component have any uses ports ?
        #--------------------------------------------
        if (port_names[0] == ''):
            if (self.DEBUG):
                print '------------------------------------------'
                print 'NOTE: This component has no uses ports.'
                print '------------------------------------------'
            SUCCESS = True
            return SUCCESS
        
        #--------------------------------------------
        # Note:  Use "d_services" added by Bocca to
        #        an Impl file to get a CCA port
        #--------------------------------------------
        SUCCESS = True
        str2 = project_name + "." + port_type + "." + port_type
        str3 = "( port )"
        
        for k in xrange(len(port_names)):
            #-----------------------------------
            # Try to get a "generic" CCA port.
            #-----------------------------------
            name     = port_names[k]
            name_str = '"' + name + '".'
            try:
                port = d_services.getPort( name )
                print "SUCCESS: Got CCA port: " + name_str
                ## print "*** type(port) =", type(port)
            except:
                print '#####################################################'
                print ' ERROR:  Unable to get CCA port named: ' + name_str
                print ' '
                print ' Please select a component from the palette that'
                print ' provides a port named ' + name_str
                print '#####################################################'
                print ' '
                ## message = "FAILURE:  Unable to get CCA port: " + name
                ## print message
                port    = None
                SUCCESS = False

            if not(SUCCESS):
                break

            #------------------------------------------
            # Try to typecast the port to "port_type"
            # and then store it within "self"
            #------------------------------------------
            ### str1 = "self." + short_names[k]
            str1 = "self.bmi." + short_names[k]   ############## (3/12/12)
            exec( str1 + " = " + str2 + str3 )
 
            exec( "UNABLE = (" + str1 + " == None)" )
            if (UNABLE):
                print 'FAILURE: Unable to cast CCA port: ' + name
                exec( "d_services.releasePort( "  + name + " )")
                SUCCESS = False
           
        return SUCCESS
        
    #   get_cca_ports()
    #-------------------------------------------------------------------
    def release_cca_ports(self, d_services):
   
        #----------------------
        #  Release all ports
        #----------------------
        for name in self.cca_port_names:
            d_services.releasePort( name )
                 
    #   release_cca_ports()
    #-------------------------------------------------------------------
    def get_status( self ):

        #-----------------------------------------------------
        # Notes: Return component status as a string.  The
        #        possible return values are from OpenMI 2.0:
        #
        #            created, initializing, initialized,
        #            updating, updated, finalizing, finalized,
        #            failed (could add "stopped").
        #-----------------------------------------------------
        return self.bmi.get_status()

    #   get_status()
    #-------------------------------------------------------------------
    def set_status( self, status ):

        #-----------------------------------------------------
        # Note: Valid values for the status string are:
        #
        #           created, initializing, initialized,
        #           updating, updated, finalizing, finalized,
        #           failed (could add "stopped").
        #       Could check if valid before setting.
        #-----------------------------------------------------
        self.status = status

    #   set_status()
    #-------------------------------------------------------------------
    def get_mode( self ):

        #-------------------------------------------------------
        # Maybe get mode string from get_attribute() instead ?
        #-------------------------------------------------------
        return self.mode    # ('driver' or 'nondriver')

    #   get_mode()
    #-------------------------------------------------------------------
    def set_mode( self, mode ):

        mode_list = ['driver', 'nondriver']
        if (mode.lower() in mode_list):
            self.mode = mode.lower()   # (store directly into self ??)
        else:
            print '################################################'
            print ' ERROR: In CMI_base.set_mode():'
            print '        Mode must be "driver" or "nondriver".'
            print '################################################'
            print ' '

         #############################################
         # Should we also set mode at BMI level ??
         #############################################
         # self.bmi.set_mode( 'driver' )
         
    #   set_mode()
    #-------------------------------------------------------------------
    def get_attribute( self, att_name ):

        #------------------------------------------------------------
        # Note: There is a BMI.get_attribute() method, but perhaps
        #       this one should read component-specific attributes
        #       from a config file (vs. model-specific).
        #------------------------------------------------------------
        try:
            #----------------------------------------------
            # This will work for model attributes, like:
            #     model_name, version, time_step_type,
            #     time_units, mesh_type, author_name, etc.
            #----------------------------------------------
            return self.bmi.get_attribute( att_name )
        except:
            pass

        #--------------------------------------------------
        # Read key-value pairs from a file and then build
        # a dictionary like this ??
        #--------------------------------------------------
        amap = {'comp_name':          'ErodeGlobal',
                'model_family':       'Erode', 
                 ## 'cfg_extension':      '_erode_global.cfg',   # (In BMI-level ??)
                'cmt_var_prefix':     '/ErodeGlobal/Input/Var/',
                'cfg_template_file':  'Erode_Global.cfg.in',
                'gui_xml_file':       '/home/csdms/cca/erode/0.5/src/share/cmt/gui/Erode_Global.xml',
                'dialog_title':       'LEM: Erode Global-Timestep Parameters' }

        try:
            return amap[ att_name.lower() ]
        except:
            print '#####################################################'
            print ' ERROR: Could not find CMI attribute: ' + att_name
            print '#####################################################'
            print ' '       
                    
    #   get_attribute()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_values( self, long_var_name ):

        #---------------------------------------------------------
        # Note: This method returns a NumPy "ndarray" object
        #       that Babel is able to pass to other components
        #       as a SIDL generic array.
        #---------------------------------------------------------
        #       The functions called by this one below are "port
        #       functions" (or methods) and therefore have
        #       static return types.
        #---------------------------------------------------------
        #       This component (self) often won't have a port
        #       connection to the caller and so can't call the
        #       methods of its caller.  Regridding and unit
        #       conversion is done by PortQueue just prior to
        #       calling CMT.set_values().
        #---------------------------------------------------------
        return self.bmi.get_values( long_var_name )
        
#         try:
#             #----------------------------------------------
#             # Get data type and rank for long_var_name.
#             # Assume that NumPy dtype string is returned.
#             #----------------------------------------------
#             type_name = self.bmi.get_var_type( long_var_name )
#             rank      = self.bmi.get_var_rank( long_var_name )
# 
#             #-----------------------------------------------------------
#             # Call the appropriate BMI method.  We could do this other
#             # ways in Python, but this is the general approach.
#             #-----------------------------------------------------------
#             if (type_name == 'float64'):
#                 if (rank == 0):
#                     values = self.bmi.get_0d_double( long_var_name )
#                 elif (rank == 1):
#                     values = self.bmi.get_1d_double( long_var_name )
#                 elif (rank == 2):
#                     values = self.bmi.get_2d_double( long_var_name )
#                 elif (rank == 3):
#                     values = self.bmi.get_3d_double( long_var_name )
#     ##            elif (rank == 4):
#     ##                values = self.bmi.get_4d_double( long_var_name )
#                     #-------------------------------------------------------------
#             elif (type_name == 'int32'):
#                 if (rank == 0):
#                     values = self.bmi.get_0d_int( long_var_name )
#                 elif (rank == 1):
#                     values = self.bmi.get_1d_int( long_var_name )
#                 elif (rank == 2):
#                     values = self.bmi.get_2d_int( long_var_name )
#                 elif (rank == 3):
#                     values = self.bmi.get_3d_int( long_var_name )
#     ##            elif (rank == 4):
#     ##                values = self.bmi.get_4d_int( long_var_name )
#                     #-------------------------------------------------------------
#             else:
#                 print '############################################'
#                 print ' ERROR: In CMI_base.get_values():'
#                 print '        Type_name of "' + type_name + '"'
#                 print '        is not yet supported.'
#                 print '############################################'
#                 print ' '
#                 
#             #-----------------------------------------------
#             # Return values on this component's own grid
#             # and with its own units for variables.
#             # Set values will do conversions when needed ?
#             #-----------------------------------------------
#             return values
# 
#         except:
#             print '######################################'
#             print ' ERROR: In CMI_base.get_values():'
#             print '        A BMI method call failed.'
#             print '        Returning a value of zero.'
#             print '######################################'
#             print ' '
#             return numpy.float64( 0 )

    #   get_values()
    #-------------------------------------------------------------------
    def set_values( self, long_var_name, values ):

        self.bmi.set_values( long_var_name, values )
        
#         #----------------------------------------------
#         # Get data type and rank for long_var_name.
#         # Assume that NumPy dtype string is returned.
#         #----------------------------------------------
#         dtype = self.bmi.get_var_type( long_var_name )
#         rank  = self.bmi.get_var_rank( long_var_name )
# 
#         #---------------------------------------------
#         # Check data type and rank for input values.
#         #---------------------------------------------
#         val_dtype = str( values.dtype )
#         val_rank  = numpy.rank( values )
#         
#         if (val_rank != rank):
#             print '#############################################'
#             print ' ERROR: In CMI_base.set_values():'
#             print '        Rank of values does not match the'
#             print '        rank of target variable named:'
#             print '          ' + long_var_name
#             print '#############################################'
#             print ' '
#             return
#         
#         if (val_dtype != dtype):
#             #------------------------------------
#             # Attempt to cast values to dtype ?
#             #------------------------------------
#             new_vals = values.astype( dtype )
#         else:
#             new_vals = values
#             
# 	#------------------------------------------------------------------
# 	# First check the type of the data that is stored in the
# 	# generic array.  Cast the generic array to a SIDL array of
# 	# that type.  Call the appropriate BMI set_value function and
# 	# pass it a pointer to the start of the data. This is pseudo
# 	# code that looks more or less like what it would be in C.
# 	# The underlying data is contiguous in memory and of unit stride.
# 	#------------------------------------------------------------------
#         if (dtype == 'float64'):
#             if (rank == 0):
#                 self.bmi.set_0d_double( long_var_name, new_vals )
#             elif (rank == 1):
#                 self.bmi.set_1d_double( long_var_name, new_vals )
#             elif (rank == 2):
#                 self.bmi.set_2d_double( long_var_name, new_vals )
#             elif (rank == 3):
#                 self.bmi.set_3d_double( long_var_name, new_vals )
# 	#------------------------------------------------------------------
#         elif (dtype == 'int32'):
#             if (rank == 0):
#                 self.bmi.set_0d_int( long_var_name, new_vals )
#             elif (rank == 1):
#                 self.bmi.set_1d_int( long_var_name, new_vals )
#             elif (rank == 2):
#                 self.bmi.set_2d_int( long_var_name, new_vals )
#             elif (rank == 3):
#                 self.bmi.set_3d_int( long_var_name, new_vals )
#         #-------------------------------------------------------------
#         else:
#             print '############################################'
#             print ' ERROR: In CMI_base.set_values():'
#             print '        Type_name of "' + dtype + '"'
#             print '        is not yet supported.'
#             print '############################################'
#             print ' '
                
    #   set_values()
    #-------------------------------------------------------------------
    def get_values_at_indices( self, long_var_name, indices ):

	#-------------------------------------------------------------
	# Notes:  For the specified variable, get the values at the
	#         specified indices and return them.  If the wrapped
	#         model is raster, then each raster cell has a long
	#         integer, calendar-style index and these are used
	#         for indices.  If the wrapped model is ugrid, then
	#         each cell in the grid has a unique, long-integer
	#         ID and these are used for indices.
	#-------------------------------------------------------------
        self.bmi.get_values_at_indices( long_var_name, indices )
    
# 	#----------------------------------------------
# 	# Get data type and rank for long_var_name.
# 	# Assume that NumPy dtype string is returned.
# 	#----------------------------------------------
# 	dtype = self.bmi.get_var_type( long_var_name )
# 	rank  = self.bmi.get_var_rank( long_var_name )
# 
# 	if (dtype == 'float64'):
# 		if (rank == 2):
# 			return self.bmi.get_2d_double_at_indices( long_var_name, indices )
# 	elif (dtype == 'int32'):
# 		if (rank == 2):
# 			return self.bmi.get_2d_int_at_indices( long_var_name, indices )
            
    #   get_values_at_indices()
    #-------------------------------------------------------------------
    def set_values_at_indices( self, long_var_name, indices, values ):

        self.bmi.set_values_at_indices( long_var_name, indices, values )
        
#         #----------------------------------------------
#         # Get data type and rank for long_var_name.
#         # Assume that NumPy dtype string is returned.
#         #----------------------------------------------
#         dtype = self.bmi.get_var_type( long_var_name )
#         rank  = self.bmi.get_var_rank( long_var_name )
# 
#         if (dtype == 'float64'):
#             if (rank == 2):
#                 self.bmi.set_2d_double_at_indices( long_var_name, indices, values )
#                 return
#         elif (dtype == 'int32'):
#             if (rank == 2):
#                 self.bmi.set_2d_int_at_indices( long_var_name, indices, values )
#                 return
# 
#         print '#####################################################'
#         print ' ERROR: CMI_base.set_values_at_indices() currently'
#         print '        requires target variable to be a 2D array.'
#         print '#####################################################'
#         print ' '
                
    #   set_values_at_indices()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #------------------------------------------------------
        # Note:  NumPy dtype of returned array will be "|Sn",
        #        where n is the length of the longest string.
        #------------------------------------------------------
        return self.bmi.get_input_var_names()
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):

        #------------------------------------------------------
        # Note:  NumPy dtype of returned array will be "|Sn",
        #        where n is the length of the longest string.
        #------------------------------------------------------
        return self.bmi.get_output_var_names()
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_type( self, long_var_name ):

        #-----------------------------------------------------
        # Notes:  This may need to be in CMI as well as BMI.
        #-----------------------------------------------------
        return self.bmi.get_var_type( long_var_name )
    
    #   get_var_type()
    #-------------------------------------------------------------------
    def get_var_units( self, long_var_name ):

        #-------------------------------------------------------------
        # This needs to be in CMI as well as BMI to be called by a
        # service component that converts units.  See CMI.initialize.
        #-------------------------------------------------------------
        return self.bmi.get_var_units( long_var_name )
    
    #   get_var_units()
    #-------------------------------------------------------------------
    def get_var_rank( self, long_var_name ):

        #-----------------------------------------------------
        # Notes:  This may need to be in CMI as well as BMI.
        #-----------------------------------------------------
        return self.bmi.get_var_rank( long_var_name )
    
    #   get_var_rank()
    #-------------------------------------------------------------------   
    #-------------------------------------------------------------------
    def get_grid_shape(self, long_var_name):

        #-------------------------------------------------------
        # Note: Return ndarray with [nx, ny, nz] as int32.
        #       We may want to reverse the order of the
        #       returned values in this method and the next 2.
        #-------------------------------------------------------
        return self.bmi.get_grid_shape()
    
    #   get_grid_shape()
    #-------------------------------------------------------------------
    def get_grid_spacing(self, long_var_name):

        #--------------------------------------------------------
        # Note: Return ndarray with [dx, dy, dz] as float64.
        #       xres and yres could have angular units like
        #       arcseconds, if (info.pixel_geom == 0).  In
        #       that case, xres (or dx) will vary with latitude
        #       and this method isn't really appropriate.
        #-------------------------------------------------------- 
        return self.bmi.get_grid_spacing()

    #   get_grid_spacing()
    #-------------------------------------------------------------------
    def get_grid_lower_left_corner(self, long_var_name):

        #--------------------------------------------------------
        # Note: Return ndarray with [x0, y0, z0] as float64.

        #       What about units for [x0,y0,z0] and [dx,dy,dz]?
        #--------------------------------------------------------        
        return self.bmi.get_grid_lower_left_corner( long_var_name )

    #   get_grid_lower_left_corner()
    #-------------------------------------------------------------------
    def get_grid_attribute(self, long_var_name, att_name):

        return self.bmi.get_grid_attribute( long_var_name, att_name ) 

    #   get_grid_attribute()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # BMI methods to get time-related info
    #-------------------------------------------------------------------
    def get_time_step(self):

        return self.bmi.get_time_step()
    
    #   get_time_step()
    #-------------------------------------------------------------------        
    def get_time_units(self):

        #-----------------------------------------------------------
        # Notes: The get_attribute() method should always provide
        #        a model's time units.  However, some components
        #        have an initialize_time_vars() method that saves
        #        a units string in self.time_units.  These strings
        #        are lower-case, e.g. 'seconds', 'minutes',
        #        or 'years'.  They should conform to UDUNITS.
        #-----------------------------------------------------------
        unit_str = self.get_attribute( 'time_units' )
        return unit_str
    
    #   get_time_units()
    #-------------------------------------------------------------------
    def get_time(self, when='current'):

        #---------------------------------------------------------
        # Note: The "when" argument can be: current, start, end.
        #---------------------------------------------------------
        # return self.bmi.get_time( when )

        if (when == 'current'):
            return self.bmi.get_current_time()
        elif (when == 'start'):
            return self.bmi.get_start_time()
        elif (when == 'stop'):
            return self.bmi.get_end_time()
        elif (when == 'end'):
            return self.bmi.get_end_time()
        else:
            return self.bmi.get_current_time()  # (e.g. 'now')
        
    #   get_time()
    #-------------------------------------------------------------------   
    #-------------------------------------------------------------------
    def initialize( self, cfg_prefix ):

        #-----------------------------------------------------------
        # NB!  We currently ignore cfg_prefix argument and assume
        #      that a string variable called "case_prefix is
        #      available from the GUI.
        #-----------------------------------------------------------
        #      We could use "cfg_dir = os.getcwd()" instead of "."
        #-----------------------------------------------------------
        self.set_status( 'initializing' )
        DEBUG = False
        
        #----------------------------------------------
        # Get CCA ports for all required components.
        # Later use the PortQueue service here.
        #----------------------------------------------
        OK = self.get_cca_ports( self.d_services )
        if not(OK):
            comp_name = self.get_attribute( 'comp_name' )
##            comp_name = self.bmi.get_attribute( 'comp_name' )
            print 'ERROR during', comp_name + '.get_cca_ports().'
            return
            ## return -1  # (not in initialize())

        #----------------------
        # Set up a Port Queue (something like this)
        #----------------------
        # port_names = ???
        # self.pq = edu.csdms.tools.IRFPortQueue.PortQueue( port_names )
        
        #--------------------------------------------
        # Build name of new model config file (CFG)
        #--------------------------------------------
        var_prefix   = self.get_attribute( 'cmt_var_prefix' )
        cfg_ext      = self.get_attribute( 'cfg_extension' )
        cfg_template = self.get_attribute( 'cfg_template_file')
##        var_prefix   = self.bmi.get_attribute( 'cmt_var_prefix' )
##        cfg_ext      = self.bmi.get_attribute( 'cfg_extension' )
##        cfg_template = self.bmi.get_attribute( 'cfg_template_file')
        #---------------------------------------------------------------
        cfg_dir  = "."
        var_name   = var_prefix + "case_prefix"  # (need full name)
        cfg_prefix = self.user_input.getString( var_name, "DEFAULT" )
        if (DEBUG):
            print '#### Read from GUI: cfg_prefix =', cfg_prefix
        cfg_file = cfg_prefix + cfg_ext

        #--------------------------------------------------------
        # Read values from tabbed-dialog GUI and then use them
        # to create a valid input file from a "template file".
        # Template files have extension ".in", e.g. ".cfg.in".
        # All info read from GUI is in self.user_input.
        #--------------------------------------------------------
        template = edu.csdms.tools.TemplateFiles.TemplateFiles()
        template.add_file( cfg_template, cfg_file )
        template.substitute( self.user_input, var_prefix, cfg_dir )

        #------------------------------------
        # Initialize the wrapped model.
        # Currently need to set the "mode".
        #------------------------------------
        self.bmi.initialize( cfg_prefix, mode=self.mode)

        self.set_status( 'initialized' )
        
    #   initialize()
    #-------------------------------------------------------------------
    def run_until( self ):

        self.set_status( 'updating' )


        self.set_status( 'updated' )

    #   run_until()
    #-------------------------------------------------------------------
    def run_for( self ):

        self.set_status( 'updating' )


        self.set_status( 'updated' )

    #   run_for()
    #-------------------------------------------------------------------
    def finalize( self ):

        self.set_status( 'finalizing' )
        self.bmi.finalize()
        self.bmi.release_cca_ports( self.d_services )
        self.set_status( 'finalized' )
        
    #   finalize()
    #-------------------------------------------------------------------
    def run_model(self):
        
        #-------------------------------------------------------------
        # Use CMI methods vs. BMI.run_model. (10/26/11)
        #-------------------------------------------------------------
        # Note: CMI.initialize ignores arg and sets cfg_prefix.
        #-------------------------------------------------------------
        ### self.mode = 'driver'
        self.set_mode( 'driver' )
        comp_name = self.get_attribute( 'comp_name' ) 
        self.initialize( '' )

        ##########################################
        # Use something like this instead soon
        ##########################################
##        try:
##            self.run_until()  ########
##        except:
##            print 'ERROR in CMI.run_model() at:'
##            print '    time_index =', self.bmi.time_index
##            self.print_traceback( comp_name + '.run_model()' )
##            return -1   # (component can still be used)
            
        while not(self.bmi.DONE):  #################
            try: 
                self.update( self.bmi.time_sec )    # (check argument; maybe use None?)
            except:
                print 'ERROR in CMI.run_model() at:'
                print '    time_index =', self.bmi.time_index
                self.print_traceback( comp_name + '.run_model()' )
                return -1   # (component can still be used)

        #--------------------------------------------
        # Later do report with a service component.
        # NB!  BMI.finalize() closes the log_file.
        #--------------------------------------------
        self.bmi.print_final_report(comp_name=comp_name + ' component', mode='driver')
        self.set_mode( 'nondriver' )  # (restore default)
        ## self.mode = 'nondriver'  # (restore default)
        self.finalize()          # (very last thing)
        return 0

    #   run_model() 
    #-------------------------------------------------------------------
    def go(self):

        self.run_model()

    #   go()
    #-------------------------------------------------------------------
    # Private methods or utilities; not part of BMI interface.
    #-------------------------------------------------------------------
    def check_finished(self):

        #------------------------------------------------------
        # Note: If self.DONE has already been set to True by
        #       another function or component, this function
        #       preserves that setting (see below).
        #------------------------------------------------------
        if (self.stop_code == 0):
            #---------------------
            # Stop after n_steps
            #---------------------
            TIMES_UP = (self.time_index >= self.n_steps)
        elif (self.stop_code == 1):
            #-----------------------
            # Stop after stop_time 
            #-----------------------
            TIMES_UP = (self.time >= self.stop_time)
        elif (self.stop_code == 2):
            #-----------------------------------------
            # Stop if "steady-state", but only works
            # as written here for global timesteps. 
            #-----------------------------------------
            TIMES_UP = (self.dz_max < self.dz_tolerance)   #######

        self.DONE = (self.DONE or TIMES_UP)

    #   check_finished()
    #-------------------------------------------------------------------
    def print_traceback(self, comp_name='None'):

        if (comp_name == 'None'):
            comp_name = self.get_attribute( 'comp_name' )

        print '################################################'
        print ' ERROR encountered in ' + comp_name
        print '       Please check your input parameters.'
        print '################################################'

        #-----------------------------------------------
        # Print the Python exception traceback message
        #-----------------------------------------------
        traceback.print_exc()
    
    #   print_traceback()                    
    #-------------------------------------------------------------------
