#!/usr/bin/env python
#
# Tests the InfilGreenAmpt TopoFlow component.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_infil_green_ampt.cfg
# - Treynor.rti

from shutil import rmtree
from os import makedirs, listdir
from os.path import join, dirname, exists
from nose.tools import raises, assert_is_instance, assert_is_not_none
import numpy as np
from topoflow.components.infil_green_ampt import infil_component
from . import input_dir, output_dir


def setup_module():
    global comp
    comp = infil_component()
    if exists(output_dir) is False:
        makedirs(output_dir)


def teardown_module():
    parent_dir = dirname(output_dir)
    if exists(parent_dir):
        rmtree(parent_dir)


def test_is_instance():
    assert_is_instance(comp, infil_component)


def test_irf():
    cfg_file = join(input_dir, 'June_20_67_infil_green_ampt.cfg')

    # @mperignon found these attributes are needed but not provided.
    comp.P_rain = 1
    comp.SM = 0
    comp.h_table = 1
    comp.elev = 10

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
