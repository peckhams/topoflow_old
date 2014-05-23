"""
This module reads and writes a series of data cubes to a NetCDF file.


=== Set up some parameters ===

Define the shape of the data cube.
>>> shape = (4, 5, 3)
>>> res = (1., 100., 200.)

Set the values of the data cube.
>>> var = numpy.random.random (shape)

Write cubes at these times.
>>> times = range (5)

This is the name of the variable that we will read/write.
>>> var_name = 'grain_size'


=== Write a cube ===

>>> import nccs_files
>>> f = nccs_files.nccs_file ()

Open the file for writing.
>>> f.open_new_file ('test_file.nc', shape=shape, res=res,
...                  dtype='float64',
...                  var_name=var_name,
...                  long_name="Sediment grain size",
...                  units_name="phe",
...                  comment="Cube of sediment grain size")
True


The name of the file may not be the same as what was specified in the call
to open_new_file.  If the file already exists, a number will be appended to
the file name (before the extension).

>>> file_name = f.file_name
>>> print file_name # doctest: +ELLIPSIS
test_file....nc

Write the variable
>>> f.add_cube (var, var_name)

Close the file.
>>> f.close_file ()

=== Read a cube ===

>>> import nccs_files
>>> f = nccs_files.nccs_file ()
>>> f.open_file (file_name)
True

Read variable from file and compare it to what we wrote.
>>> var_from_file = f.get_cube (var_name, 0)
>>> (var_from_file == var).all ()
True
>>> f.close_file ()

=== Write a series of cubes ===

Open the file for writing.
>>> f.open_new_file ('test_file.nc', shape=shape, res=res,
...                  dtype='float64',
...                  var_name=var_name,
...                  long_name="Sediment grain size",
...                  units_name="phe",
...                  comment="Cube of sediment grain size")
True

>>> file_name = f.file_name
>>> print file_name # doctest: +ELLIPSIS
test_file....nc

Write the variable
>>> for time in times:
...   f.add_cube (var+time, var_name)

Close the file.
>>> f.close_file ()

=== Read a series of cubes ===

>>> import nccs_files
>>> f = nccs_files.nccs_file ()
>>> f.open_file (file_name)
True

Read variable from file and compare it to what we wrote.
>>> values_match = True
>>> for time in times:
...   var_from_file = f.get_cube (var_name, time)
...   values_match &= (var_from_file == var+time).all ()
>>> values_match
True
>>> f.close_file ()

"""

import os
import sys
import time

import numpy
import bov_files
import file_utils
import rti_files

#---------------------------------------------------------------------
# This class is for I/O of time-indexed 3D arrays to netCDF files.
#---------------------------------------------------------------------
#
#   unit_test()
#   save_nccs_cube()   # (not ready yet)
#
#   class nccs_file():
#
#       import_nio()
#       open_file()
#       get_nio_type_map()
#       open_new_file()
#       update_time_index()
#----------------------------
#       add_cube()
#       get_cube()
#----------------------------
#       close_file()
#       close()
#
#-------------------------------------------------------------------
def save_nccs_cube(nccs_file_name=None, rtg_file_name=None):

    nccs = nccs_file ()
    OK = nccs.open_file(nccs_file_name)
    if not OK:
        return

    var_name  = 'H'
    time_index = 200
    cube = nccs.get_cube(var_name, time_index)
    nccs.close()
    
    cube = numpy.array(cube)
    print 'min(cube), max(cube) =', cube.min(), cube.max()

    rtg_unit = open(rtg_file_name, 'wb')
    cube.tofile(unit)
    rtg_unit.cube()

#   save_nccs_cube()
#-------------------------------------------------------------------
class nccs_file():

    #----------------------------------------------------------
    # Note:  ncgs = NetCDF Grid Stack (used by CSDMS)
    #----------------------------------------------------------
    def import_nio(self):

        try:
            import Nio  # (a module in the PyNIO package) 
            # print 'Imported Nio version: ' + Nio.__version__
            return Nio
        except:
##            python_version = sys.version[:3]
##            print ' '
##            print 'SORRY, Cannot write netCDF files because'
##            print 'the "Nio" package cannot be imported.'
##            print ' '
##            if (python_version != '2.6'):
##                print 'Note that "PyNIO" is only installed for'
##                print 'Python version 2.6 on "beach".'
##                print 'The current Python version is:', python_version
##                print ' '
            return False
        
    #   import_nio()
    #----------------------------------------------------------
    def open_file(self, file_name):

        #--------------------------------------------------
        # Try to import the Nio module from PyNIO package
        #--------------------------------------------------
        Nio = self.import_nio ()
        if not Nio:
          return
        
        #-------------------------
        # Open file to read only
        #-------------------------
        try:
            nccs_unit = Nio.open_file(file_name, mode="r")
            self.nccs_unit = nccs_unit
            return True
        except:
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def get_nio_type_map(self):

        #----------------------------------------
        # Possible settings for "nio_type_code"
        #-------------------------------------------
        # nio_type_code = "d"  # (double, Float64)
        # nio_type_code = "f"  # (float,  Float32)
        # nio_type_code = "l"  # (long,   Int64)
        # nio_type_code = "i"  # (int,    Int32)
        # nio_type_code = "h"  # (short,  Int16)
        # nio_type_code = "b"  # (byte,   Int8)
        # nio_type_code = "S1" # (char)
        #-------------------------------------------
        nio_type_map = {'float64':'d', 'float32':'f',
                        'int64':'l', 'int32':'i',
                        'int16':'s', 'int8':'b',
                        'S|100':'S1'}  # (check last entry)                      

        return nio_type_map
    
    #   get_nio_type_map()
    #----------------------------------------------------------
    def open_new_file(self, file_name, info=None,
                      var_name='X',
                      long_name=None,
                      units_name='None',
                      dtype='float64',
                      ### dtype='float64'
                      time_units='minutes',
                      comment='',
                      shape=(1,1,1),
                      res=(1.,1.,1.),
                      MAKE_RTI=True, MAKE_BOV=False):
            
        #--------------------------------------------------
        # Try to import the Nio module from PyNIO package
        #--------------------------------------------------
        Nio = self.import_nio ()
        if not Nio:
            return False

        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        self.file_name = file_name
        
        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.format     = 'nccs'
        self.file_name  = file_name
        self.time_index = 0
        self.var_name   = var_name
        self.shape      = shape
        self.res        = res
        
        if (long_name is None):
            long_name = var_name
        self.long_name  = long_name
        self.units_name = units_name
        self.dtype      = dtype

        #-----------------------------------
        # Get Nio type code for this dtype
        #------------------------------------
        nio_type_map  = self.get_nio_type_map()        
        nio_type_code = nio_type_map[ dtype.lower() ]        
        self.nio_type_code = nio_type_code
        
        #-------------------------------------
        # Open a new netCDF file for writing
        #-------------------------------------
        # Sample output from time.asctime():
        #     "Thu Oct  8 17:10:18 2009"
        #-------------------------------------
        opt = Nio.options()
        opt.PreFill = False            # (for efficiency)
        opt.HeaderReserveSpace = 4000  # (4000 bytes, for efficiency)
        history = "Created using PyNIO " + Nio.__version__ + " on "
        history = history + time.asctime() + ". " 
        history = history + comment
        # print 'MADE IT PAST history BLOCK'
        
        try:
            nccs_unit = Nio.open_file (file_name, mode="w",
                                       options=opt, history=history)
            OK = True
        except:
            OK = False
            return OK

        #----------------------------------------------
        # Create grid dimensions nx and ny, plus time
        #----------------------------------------------
        # Without using "int()" here, we get this:
        #     TypeError: size must be None or integer
        #----------------------------------------------
        nccs_unit.create_dimension("nz", self.shape[0])
        nccs_unit.create_dimension("ny", self.shape[1])
        nccs_unit.create_dimension("nx", self.shape[2])
        nccs_unit.create_dimension("time", None)   # (unlimited dimension)
        # print 'MADE IT PAST create_dimension CALLS.'
        
        #-------------------------
        # Create a time variable
        #------------------------------------------
        #('d' = float64; must match in add_cube()
        #------------------------------------------
        tvar = nccs_unit.create_variable ('time', 'd', ("time",))
        nccs_unit.variables['time'].units = time_units
        
        #--------------------------------
        # Create a variable in the file
        #----------------------------------
        # Returns "var" as a PyNIO object
        #----------------------------------
        var = nccs_unit.create_variable (var_name, nio_type_code,
                                         ("time", "nz", "ny", "nx"))

        #----------------------------------
        # Specify a "nodata" fill value ?
        #----------------------------------
        var._FillValue = -9999.0    ## Does this jive with Prefill above ??
        
        #------------------------------------
        # Create attributes of the variable
        #------------------------------------
        nccs_unit.variables[var_name].long_name = long_name
        nccs_unit.variables[var_name].units = units_name
        nccs_unit.variables[var_name].dz = self.res[0]
        nccs_unit.variables[var_name].dy = self.res[1]
        nccs_unit.variables[var_name].dx = self.res[2]
        nccs_unit.variables[var_name].y_south_edge = 0.
        nccs_unit.variables[var_name].y_north_edge = self.res[1]*self.shape[1]
        nccs_unit.variables[var_name].x_west_edge = 0.
        nccs_unit.variables[var_name].x_east_edge = self.res[2]*self.shape[2]
        nccs_unit.variables[var_name].z_bottom_edge = 0.
        nccs_unit.variables[var_name].z_top_edge = self.res[0]*self.shape[0]
        
        self.nccs_unit = nccs_unit
        return OK
    
    #   open_new_file()
    #----------------------------------------------------------
    def add_cube(self, grid, var_name, time=None,
                 time_index=-1):

        #---------------------------------
        # Assign a value to the variable
        #-------------------------------------------
        # This syntax works for scalars and grids
        #-------------------------------------------
        # nc_unit.variables[var_name].assign_value( grid )
 
        #-------------------------------------
        # Can use time_index to overwrite an
        # existing grid vs. simple append.
        #-------------------------------------
        if time_index == -1:
            time_index = self.time_index
        if time == None:
            time = numpy.float64(time_index)
            
        #---------------------------------------
        # Write a time to existing netCDF file
        #---------------------------------------
        times = self.nccs_unit.variables['time']
        times[time_index] = time
        
        #---------------------------------------
        # Write a grid to existing netCDF file
        #---------------------------------------
        var = self.nccs_unit.variables[var_name]
        if numpy.rank(grid) == 0:
            #-----------------------------------------------
            # "grid" is actually a scalar (dynamic typing)
            # so convert it to a grid before saving
            #-----------------------------------------------
            grid2 = grid + numpy.zeros([self.nz, self.ny, self.nx],
                                        dtype=self.dtype)
            var[time_index] = grid2.astype(self.dtype)
        else:
            var[time_index] = grid.astype(self.dtype)

        #---------------------------
        # Increment the time index
        #---------------------------            
        self.time_index += 1

    #   add_cube()
    #----------------------------------------------------------
    def get_cube(self, var_name, time_index):

        var = self.nccs_unit.variables[var_name]
        return var[time_index]
        
    #   get_cube()
    #-------------------------------------------------------------------
    def close_file(self):

        self.nccs_unit.close ()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        self.nccs_unit.close ()

    #   close()
    #-------------------------------------------------------------------
   
if __name__ == "__main__":
    import doctest
    doctest.testmod()
 
