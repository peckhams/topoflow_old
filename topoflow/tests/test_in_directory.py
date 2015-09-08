#!/usr/bin/env python
#
# Tests using absolute and relative paths for the TopoFlow .cfg file.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_meteorology.cfg
# - Treynor.rti
# - June_20_67_rain_rates.txt
# - Treynor_slope.bin
# - Treynor_aspect.bin

from shutil import rmtree
from os import makedirs, listdir, chdir, getcwd
from os.path import join, dirname, exists
from nose.tools import assert_is_instance, assert_is_not_none
from numpy.testing import assert_almost_equal
from topoflow.components.met_base import met_component as Model
from . import input_dir, output_dir, time_factor


def setup_module():
    global comp
    comp = Model()
    if exists(output_dir) is False:
        makedirs(output_dir)


def teardown_module():
    parent_dir = dirname(output_dir)
    if exists(parent_dir):
        rmtree(parent_dir)


def test_initialize_from_elsewhere():
    cfg_file = join(input_dir, 'June_20_67_meteorology.cfg')
    print '\n***cfg_file:', cfg_file
    comp.initialize(cfg_file)


def test_initialize_from_current_directory():
    current_dir = getcwd()
    chdir(input_dir)
    cfg_file = 'June_20_67_meteorology.cfg'
    print '\n***cfg_file:', cfg_file
    comp.initialize(cfg_file)
    chdir(current_dir)
