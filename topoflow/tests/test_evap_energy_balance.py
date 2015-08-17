#!/usr/bin/env python
#
# Tests the EvapEnergyBalance TopoFlow component.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_evap_energy_balance.cfg
# - Treynor.rti

from shutil import rmtree
from os import makedirs, listdir
from os.path import join, dirname, exists
from nose.tools import assert_is_instance, assert_is_not_none
import numpy as np
from topoflow.components.evap_energy_balance import evap_component as Model
from . import input_dir, output_dir


def setup_module():
    global comp
    comp = Model()
    if not exists(output_dir):
        makedirs(output_dir)


def teardown_module():
    parent_dir = dirname(output_dir)
    if exists(parent_dir):
        rmtree(parent_dir)


def test_is_instance():
    assert_is_instance(comp, Model)


def test_irf():
    cfg_file = join(input_dir, 'June_20_67_evap_energy_balance.cfg')

    # @mperignon found these attributes are needed but not provided.
    comp.h_snow = 1
    comp.Q_sum = 1
    comp.Qe = 1
    comp.T_air = 10
    comp.T_surf = 20

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
