
SILENT = False
if not(SILENT):
    print 'Importing TopoFlow packages:'
    print '   topoflow.utils'
    print '   topoflow.utils.tests'
    print '   topoflow.components'
    print '   topoflow.components.tests'
    print '   topoflow.framework'
    print '   topoflow.framework.tests'
    print '   topoflow.gui (unfinished)'
    print ' '

import topoflow.utils
import topoflow.utils.tests
#-----------------------------------
import topoflow.components
import topoflow.components.tests
#-----------------------------------
import topoflow.framework
import topoflow.framework.tests
import topoflow.framework.tests.test_framework
#-----------------------------------
import topoflow.gui

#---------------------------------------
# Import the topoflow_test() function,
# so it can be executed with:
#    >>> import topoflow
#    >>> topoflow.topoflow_test().
#---------------------------------------
# from topoflow.framework.tests.test_framework import topoflow_test

#------------------------------------------------
# Another idea. This can also be executed with:
# >>> topoflow.topoflow_test()
#-------------------------------------------------------------
##def topoflow_test():
##
##    topoflow.framework.tests.test_framework.topoflow_test()
##
###   topoflow_test() 
#-------------------------------------------------------------
