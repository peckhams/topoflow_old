#!/usr/bin/env python
#
# Tests the EvapPriestleyTaylor TopoFlow component.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_evap_priestley_taylor.cfg
# - Treynor.rti

from shutil import rmtree
from os import makedirs, listdir
from os.path import join, dirname, exists
from nose.tools import assert_is_instance, assert_is_not_none
from topoflow.components.evap_priestley_taylor import evap_component
from . import input_dir, output_dir


def setup_module():
    global comp
    comp = evap_component()
    if exists(output_dir) is False:
        makedirs(output_dir)


def teardown_module():
    parent_dir = dirname(output_dir)
    if exists(parent_dir):
        rmtree(parent_dir)


def test_is_instance():
    assert_is_instance(comp, evap_component)


def test_irf():
    cfg_file = join(input_dir, 'June_20_67_evap_priestley_taylor.cfg')

    # @mperignon found these attributes are needed but not provided.
    comp.Qn_SW = 1
    comp.Qn_LW = 1
    comp.T_air = 10
    comp.T_surf = 20

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
