#!/usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup

## from distutils.core import setup

setup(name='topoflow',
      version='3.4',
      description='d8-based, spatial hydrologic model',
      author='Scott D. Peckham',
      author_email='Scott.Peckham@colorado.edu',
      license='MIT',
      url='http://csdms.colorado.edu/wiki/Model:TopoFlow',
      packages=['topoflow',
                'topoflow.components',
                'topoflow.components.tests',
                'topoflow.examples',
                'topoflow.framework',   # (later in REQUIRES)
                'topoflow.framework.tests',
                'topoflow.gui',         # (11/8/13)
                'topoflow.utils',       # (later in REQUIRES)
                'topoflow.utils.tests'], 
      install_requires=['numpy', 'scipy', 'h5py', 'netCDF4'],
      entry_points={
          'console_scripts': [
              'topoflow = topoflow.components.topoflow:main',
          ]
       },
      #-------------------------------------------------------
      # There is debate online about the pros and cons of
      # using "install_requires" since it can interfere with
      # a user's existing installed packages.
      #-------------------------------------------------------
      # PyNIO allows reading and writing of NetCDF files.
      # scimath.units (from Canopy) allows unit conversion.
      #-------------------------------------------------------      
      # Right now, the topoflow package includes subpackages
      # called "utils" and "framework" that may later be
      # distributed as separate Python packages.
      #-------------------------------------------------------      
      # install_requires=["numpy","scimath","PyNIO"],
      # provides=[],
      # obsoletes=[]
      # package_data={'':['']},
      # data_files=[],   # (include topoflow.examples here instead?)
      # test_suite='topoflow.tests',
     )
