
# S.D. Peckham
# June 10, 2009
#-------------------

import time
import numpy

def loop_test(n=1000000):

    #---------------------------------
    # For loop using Python's range
    #---------------------------------
    start  = time.time()
    my_sum = 0
    
    for k in range(n):
        my_sum += 1

    run_time = (time.time() - start)
    print 'run time with range =', run_time
    
    #---------------------------------
    # For loop using Python's xrange
    #---------------------------------
    start  = time.time()
    my_sum = 0
    
    for k in xrange(n):
        my_sum += 1

    run_time = (time.time() - start)
    print 'run time with xrange =', run_time

    #--------------------------------
    # For loop using NumPy's arange
    #---------------------------------
    start  = time.time()
    my_sum = 0
    
    for k in numpy.arange(n):
        my_sum += 1

    run_time = (time.time() - start)
    print 'run time with arange =', run_time
