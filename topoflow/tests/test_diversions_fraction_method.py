#!/usr/bin/env python
#
# Tests the DiversionsFraction TopoFlow component.
#
# site_prefix = 'Treynor'
# case_prefix = 'June_20_67'
#
# Required input files:
# - June_20_67_diversions_fraction_method.cfg
# - Treynor.rti

from shutil import rmtree
from os import makedirs, listdir
from os.path import join, dirname, exists
from nose.tools import raises, assert_is_instance, assert_is_not_none
from topoflow.components.diversions_fraction_method import diversions_component
from . import input_dir, output_dir


def setup_module():
    global comp
    comp = diversions_component()
    if not exists(output_dir):
        makedirs(output_dir)


def teardown_module():
    parent_dir = dirname(output_dir)
    if exists(parent_dir):
        rmtree(parent_dir)


def test_is_instance():
    assert_is_instance(comp, diversions_component)


def test_irf():
    cfg_file = join(input_dir, 'June_20_67_diversions_fraction_method.cfg')

    comp.initialize(cfg_file)
    comp.update()
    comp.finalize()


def test_has_output():
    assert_is_not_none(listdir(output_dir))
