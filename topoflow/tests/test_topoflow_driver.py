#!/usr/bin/env python
#
# Tests the TopoFlow driver component.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_topoflow.cfg
# - Treynor.rti

from shutil import rmtree
from os import makedirs, listdir
from os.path import join, dirname, exists
import numpy as np
from nose.tools import assert_is_instance, assert_is_not_none
from topoflow.components.topoflow_driver import topoflow_driver as Model
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
    cfg_file = join(input_dir, 'June_20_67_topoflow.cfg')

    # @mperignon and @mdpiper found these attributes are needed but not provided.
    comp.Q_outlet = np.array(1)
    comp.Q_peak = 1
    comp.T_peak = 1
    comp.u_peak = 1
    comp.Tu_peak = 1
    comp.d_peak = 1
    comp.Td_peak = 1
    comp.P_max = 1
    comp.vol_P = 1
    comp.vol_Q = 1
    comp.vol_SM = 1
    comp.vol_MR = 1
    comp.vol_ET = 1
    comp.vol_IN = 1
    comp.vol_Rg = 1
    comp.vol_GW = 1
    comp.vol_R = 1
    comp.Q_min = 0
    comp.Q_max = 1
    comp.u_min = 0
    comp.u_max = 1
    comp.d_min = 0
    comp.d_max = 1

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
