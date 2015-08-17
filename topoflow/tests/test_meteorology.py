#!/usr/bin/env python
#
# Tests the Meteorology TopoFlow component.
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
from os import makedirs, listdir
from os.path import join, dirname, exists
from nose.tools import assert_is_instance, assert_is_not_none
from topoflow.components.met_base import met_component as Model
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


def test_irf():
    cfg_file = join(input_dir, 'June_20_67_meteorology.cfg')

    # This attributes is needed but not provided.
    comp.h_snow = 1

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
