
## Copyright (c) 2001-2010, Scott D. Peckham
## January, August 2009
## May 2010 (changes to unit_test(), initialize(), etc.)

################################################################
# NB!  The update_volume_in() method ONLY tracks precip now.
#      Add get_ports() method and "precip port" ???
#      "channels_base.py" now has update_volume_out() also.

################################################################

#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class basins_component  (inherits from CSDMS_base.py)
#
#      initialize()
#      update()          # (non-OpenMI arguments)
#      finalize()
#      get_cfg_extension()     ############### OBSOLETE SOON ?
#      read_cfg_file()         ############### OBSOLETE SOON ?
#      -----------------------
#      get_cca_port_info()      # (5/11/10)
#      write_outlet_file()      # (11/7/11)
#      read_outlet_data()
#      read_main_basin_IDs()
#      check_outlet_IDs()
#      -----------------------
#      update_volume_in()    # (commented out)
#      update_volume_out()   # (commented out)
#
#-----------------------------------------------------------------------

from numpy import *
import numpy
import os
import os.path

## import cfg_files as cfg
import BMI_base
import tf_utils

#-----------------------------------------------------------------------
def unit_test():

    b = basins_component()
    b.CCA   = False
    b.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    os.chdir( cfg_directory )
    cfg_prefix      = tf_utils.TF_Test_case_prefix()
    b.site_prefix = 'Treynor'  #############       

    b.initialize( cfg_prefix=cfg_prefix, mode='driver' )

    print 'outlet_ID    =', b.outlet_ID
    print 'basin_area   =', b.basin_area
    print 'basin_relief =', b.basin_relief
    print ' '
    print 'n_outlets    =', b.n_outlets
    print 'outlet_cols  =', b.outlet_cols
    print 'outlet_rows  =', b.outlet_rows
    print 'reliefs      =', b.basin_reliefs
    print 'areas        =', b.basin_areas
    print ' '
    print 'nx           =', b.nx
    print 'ny           =', b.ny
    print ' '
    print "get_status()            = ", b.get_status()
    print "is_scalar('n_outlets')  = ", b.is_scalar('n_outlets')
    print "is_grid('n_outlets')    = ", b.is_grid('n_outlets')
    # Next one has double size, since its really a tuple.
    print "get_vector_long('outlet_IDs') = ", \
          b.get_vector_long('outlet_IDs')
    print "get_scalar_double('basin_area') = ", \
          b.get_scalar_double('basin_area')
##    print "get_vector_double('basin_areas') = ", \
##          b.get_vector_double('basin_areas')
    print 'Finished with unit_test().'
    print ' '
    
#   unit_test()   
#-----------------------------------------------------------------------
class basins_component( BMI_base.BMI_component):
    
    def initialize(self, cfg_prefix=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print 'Basins component: Initializing...'

        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_prefix = cfg_prefix
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        ## self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()
 
        #---------------------------------------------
        # Read outlet IDs (IDs of monitored pixels)
        # and their attributes like area and relief.
        # Then read IDs of all cells in the first
        # (or main) basin, i.e. above first outlet.
        #---------------------------------------------
        self.read_outlet_data()   # (uses nx and ny)
        
        # self.read_main_basin_IDs()  # (see below)

        #-------------------------------------------
        # Prepare to track total water in and out
        # of the main basin (using basin RTM file)
        #-------------------------------------------
        self.get_pvolume = False   ####
        TRACK_VOLUME = (self.get_pvolume and (self.basin_RTM_file != ''))
        self.TRACK_VOLUME = TRACK_VOLUME
        if (TRACK_VOLUME):
            #-------------------------------------------
            # Prepare to track total water in and out
            # of the main basin (using basin RTM file)
            #----------------------------------------------
            # This requires knowing the IDs of the pixels
            # that lie within the basin (basin_IDs).
            #----------------------------------------------
            self.read_main_basin_IDs()
            self.volume_in  = self.initialize_scalar( 0, dtype='float64')
            self.volume_out = self.initialize_scalar( 0, dtype='float64')

        self.status = 'initialized'
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, Q, time, dt, da, pv):

        self.status = 'updating'  # (OpenMI)
        
        if (self.TRACK_VOLUME):
            self.update_volume_out(Q, dt)
            self.update_volume_in(time, dt, da, pv)

        #------------------------
        # Update internal clock
        #------------------------
        # self.update_time()

        self.status = 'updated'  # (OpenMI)
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalized'  # (OpenMI)

    #   finalize()
    #-------------------------------------------------------------------
    def get_cfg_extension(self):

        return '_basins.cfg'
    
    #   get_cfg_extension()
    #-------------------------------------------------------------------
    def read_config_file(self):

        #-----------------------
        # 6/28/10.  Need this.
        #-----------------------
        pass
    
    #   read_config_file()
    #-------------------------------------------------------------------
    def read_cfg_file(self):

        #---------------------------------
        # (5/11/10) Needs to be here for
        # initialize_input_vars()
        #---------------------------------
        # Need to read:
        #     site_prefix        ###############
        #     outlet_file
        #     basin_RTM_file
        #---------------------------------
        pass

    #   read_cfg_file()
    #-------------------------------------------------------------------
    def get_cca_port_info(self):

        self.cca_port_names = ['None']   # (has no uses ports)
        self.cca_port_short = ['np']
        self.cca_port_type  = "IRFPort"
        self.cca_project    = "edu.csdms.models"

    #   get_cca_port_info()  
    #-------------------------------------------------------------------
    def write_outlet_file(self, pixel_file):
 
        #---------------------------------------------------
        # Written to provide default if missing. (11/7/11)
        #---------------------------------------------------
        dash_line = ''.rjust(60, "-")  #########
        unit = open( pixel_file, 'w' )
        unit.write('\n')
        unit.write( dash_line + '\n')
        unit.write(' Monitored Grid Cell (Outlet) Information\n')
        unit.write( dash_line + '\n')
        format1 = '%10s%10s%16s%16s'
        format2 = '%10i%10i%16.4f%16.4f'
        header = format1 % ('Column', 'Row', 'Area [km^2]', 'Relief [m]')
        unit.write( header + '\n' )
        unit.write( dash_line + '\n')
        info_str = format2 % (self.nx/2, self.ny/2, 0.0, 0.0)
        unit.write( info_str + '\n')
        unit.close()
 
    #   write_outlet_file()
    #-------------------------------------------------------------------
    def read_outlet_data(self):

        #----------------------------------------------------------
        # Note: outlet_IDs and basin_IDs are 1D arrays or vectors
        #       but have very different sizes.  Recall that they
        #       are (row,col) tuples and not just long ints.
        #----------------------------------------------------------
        # Note:  "basin_IDs" gives IDs of all grid cells that lie
        #        within the basin that drains to outlet_ID.
        #----------------------------------------------------------
        # Notes: outlet data is stored in a multi-column
        #        textfile as:
        #            Column, Row, Area [km^2], Relief [m]
        #---------------------------------------------------
        if (hasattr(self, 'pixel_file')):
            self.outlet_file = self.pixel_file    # (10/25/11)
        else:
            self.outlet_file = (self.in_directory +
                                self.case_prefix + '_outlets.txt')

        #--------------------------------------
        # Does outlet_file exist ? (10/25/11)
        #--------------------------------------
        if not(os.path.exists( self.outlet_file )):
            hash_line = ''.rjust(60, "#")
            print hash_line 
            print ' ERROR: Could not find monitored pixel file:'
            print ' ' + self.outlet_file
            print ' Creating default version of file. '
            print hash_line 
            print ' '
            self.write_outlet_file( self.outlet_file ) 
            # return
 
##        file_unit = open(self.outlet_file, 'r')
##        cfg.skip_header(file_unit, n_lines=6)
        
        file_unit = open(self.outlet_file, 'r')    
        lines     = file_unit.readlines()
        file_unit.close()
        lines = lines[6:]   # (skip over 6 header lines)
        n_lines = len(lines)
        
        self.outlet_cols   = numpy.zeros(n_lines, dtype='Int64')  ###
        self.outlet_rows   = numpy.zeros(n_lines, dtype='Int64')
        self.basin_areas   = numpy.zeros(n_lines, dtype='Float64')
        self.basin_reliefs = numpy.zeros(n_lines, dtype='Float64')
        
        n = 0
        for line in lines:
            words = line.split()    
            if (len(words) >= 4):
                self.outlet_cols[n]   = int64(float64(words[0]))
                self.outlet_rows[n]   = int64(float64(words[1]))
                self.basin_areas[n]   = float64(words[2])
                self.basin_reliefs[n] = float64(words[3])
                n += 1
        self.n_outlets = n
        
        #-------------------------------------------
        # Save IDs as a tuple of row indices and
        # calendar indices, "numpy.where" style
        #------------------------------------------- 
##        self.outlet_IDs    = (self.outlet_rows,
##                              self.outlet_cols)
##        self.outlet_ID     = (self.outlet_rows[0],
##                              self.outlet_cols[0])
        
        #------------------------------------------
        # Save IDs as a 1D array of long-integer,
        # calendar-style indices
        #----------------------------------------------------
        # NB! Must return as Int32 vs. Int64 currently !!!!     SEEMS OBSOLETE.
        #----------------------------------------------------
        # (11/20/11) If outlet_file only has one outlet,
        # then using int32() here converts vector with one
        # element to a scalar and produces an error.
        #-----------------------------------------------------
        ## self.outlet_IDs = int32(self.outlet_rows * self.nx) + int32(self.outlet_cols)
        self.outlet_IDs = (self.outlet_rows * self.nx) + self.outlet_cols
        self.outlet_ID  = self.outlet_IDs[0]
        ## print 'outlet_ID =', self.outlet_ID
        
        self.basin_area   = self.basin_areas[0]
        self.basin_relief = self.basin_reliefs[0]

        #------------------------------------------
        # Are all the outlet IDs inside the DEM ?
        #------------------------------------------
        OK = self.check_outlet_IDs()
        if not(OK):
            print 'ERROR: Some outlet_IDs lie outside of DEM.'
            print ' '
            return
        
    #   read_outlet_data()
    #-------------------------------------------------------------------
    def read_main_basin_IDs(self):

        #----------------------------------------------------------
        # Note: outlet_IDs and basin_IDs are 1D arrays or vectors
        #       but have very different sizes.  Recall that they
        #       are (row,col) tuples and not just long ints.
        #----------------------------------------------------------
        # Note:  "basin_IDs" gives IDs of all grid cells that lie
        #        within the basin that drains to outlet_ID.
        #----------------------------------------------------------
        self.basin_RTM_file = (self.in_directory +
                               self.site_prefix + '_basin.rtm')
        
        #----------------------------------
        # Read basin pixels from RTM_file
        #----------------------------------
        RTM_unit = open(self.basin_RTM_file, 'rb')
        RTM_filesize = os.path.getsize(RTM_unit.name)
        n_IDs = RTM_filesize / int32(4)

        basin_IDs = fromfile(RTM_unit, count=n_IDs, dtype='Int32')
        if (self.rti.SWAP_ENDIAN):
            basin_IDs.byteswap(True)
        RTM_unit.close()

        #-----------------------------------------
        # NB! basin_IDs is now 1D array of longs
        #-----------------------------------------
        wb  = where(basin_IDs >= 0)
        nwb = size(wb)   # (see note)
        ### nwb = size(wb[0])
        if (nwb != 0):    
            basin_IDs = basin_IDs[wb]
        else:    
            basin_IDs = self.outlet_ID  ###

        #------------------------------------------
        # Save IDs as a 1D array of long-integer,
        # calendar-style indices
        #------------------------------------------
        self.basin_IDs = basin_IDs    # (use Q.flat[basin_IDs])

        #-------------------------------------------
        # Save IDs as a tuple of row indices and
        # calendar indices, "numpy.where" style
        #-------------------------------------------        
        ### nx = self.nx
        ### self.basin_IDs = (basin_IDs / nx, basin_IDs % nx)

    #   read_main_basin_IDs()
    #-------------------------------------------------------------------
    def check_outlet_IDs(self):

        OK = True
        ## if not(self.PIXEL_OUTPUT): return OK
     
        wbad = where(logical_or((self.outlet_IDs < 0), \
                               (self.outlet_IDs > (self.rti.n_pixels - 1))))
        n_bad = size(wbad[0])
        
        if (n_bad != 0):    
            msg = array(['SORRY, ', ' ', \
                         'One or more of the monitored pixel IDs ', \
                         'are not in the range of valid values. ', \
                         ' ', \
                         'You can use hydrologic GIS software to get ', \
                         'the outlet ID for the main basin. ', ' '])
            result = GUI_Message(msg, INFO=True, TITLE='Missing Input')
            OK = False
            
        return OK

    #   check_outlet_IDs()
    #-------------------------------------------------------------------
    # NOTE:  May be better to move this into "precip_base.py".
    #-------------------------------------------------------------------
##    def update_volume_in(self, time, dt, da, pv):
##       
##        #----------------------------------------------------------
##        # Notes:  This procedure integrates precip. over the main
##        #         basin using the model vs. sampling timestep.
##
##        #         Recall that da is a grid [km^2].
##        #----------------------------------------------------------
##        if (pv.method == 0): return
##        
##        if (pv.method == 1):
##            #------------------------------------------------
##            # In this case, pv.rates and pv.durations are
##            # 1D vectors (vs. scalar or grid), which does
##            # not conform to the general approach now used
##            # throughout TopoFlow and by PRECIP_METHOD 2.
##            #------------------------------------------------
##            wd  = where(time < pv.duration_sums)
##            nwd = size(wd[0])
##            if (nwd != 0):
##                # rate = pv.rates[wd[0]]     ########
##                rate = pv.rates[wd[0][0]]     ######################
##                dvol = dt * rate * pv.basin_area * float64(1000000)
##                self.volume_in += dvol
##        else:    
##            #----------------------------------------------------
##            # If pv.durations is a scalar, then duration_sums
##            # is equal to the same scalar, but this still works
##            # as written (3/20/07)
##            #----------------------------------------------------
##            n_rates = size(pv.rates)
##            if (n_rates == 1):    
##                P_rates = pv.rates
##            else:    
##                P_rates = pv.rates[self.basin_IDs]
##            #-------------------------------------------------------
##            n_durs = size(pv.duration_sums)
##            if (time <= pv.duration_sums[n_durs - 1]):    
##                if (size(da) == 1):    
##                    nb   = size(self.basin_IDs[0]) ### BUG FIX.
##                    dvol = dt * sum(double(P_rates * da * nb))
##                else:    
##                    dvol = dt * sum(double(P_rates * da[self.basin_IDs]))
##                self.volume_in += dvol
##        
##    #   update_volume_in()
##    #-------------------------------------------------------------------
##    def update_volume_out(self, Q, dt):
##        
##        self.volume_out += (Q[self.outlet_ID] * dt)
##          
##    #   update_volume_out()               
    #-------------------------------------------------------------------
        


