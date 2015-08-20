"""Nosetests for TopoFlow components."""

from os.path import join, dirname, realpath, expanduser


input_dir = join(dirname(realpath(__file__)), '..', 'examples', 'Treynor_Iowa')
output_dir = join(expanduser('~'), 'TopoFlow_Tests', 'Test1')
time_factor = 5.0
