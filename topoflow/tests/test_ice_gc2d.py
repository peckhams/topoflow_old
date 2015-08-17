#!/usr/bin/env python
#
# Tests the GC2D TopoFlow component.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_ice_valley_glacier.cfg
# - Treynor.rti
# - Treynor_DEM.rtg

from shutil import rmtree
from os import makedirs, listdir
from os.path import join, dirname, exists
from nose.tools import raises, assert_is_instance, assert_is_not_none
from topoflow.components.ice_base import ice_component as Model
from . import input_dir, output_dir


def setup_module():
    global comp
    comp = Model()
    if exists(output_dir) is False:
        makedirs(output_dir)


def teardown_module():
    parent_dir = dirname(output_dir)
    if exists(parent_dir):
        rmtree(parent_dir)


def test_is_instance():
    assert_is_instance(comp, Model)


# Error using filter2d in gc2d.py, line 985.
@raises(ValueError)
def test_irf():
    cfg_file = join(input_dir, 'June_20_67_ice_valley_glacier.cfg')

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
