
S.D. Peckham
September 25, 2014

This folder contains files that were used to test the Diversions component when it was first developed.  The CFG files are of the "old style" and will need to be converted to the "new style" before these tests can be used again.

Some functionality of the Diversions component has been moved into channels_base.py as the update_diversions() method, but this has not been tested yet.  The original code is in: topoflow/components/diversions_fraction_method_LAST.py.

The plan is to add a test to topoflow/framework/tests/test_framework.py that uses the Test_Plane_Canal data set to test the Diversions component.
