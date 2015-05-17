import cython
import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def testarray(np.ndarray[DTYPE_t, ndim=1] a):
   cdef int i,j,k,l
   cdef DTYPE_t sum=0.0

   for i in range(100):
    for j in range(100):
        for k in range(100):
            for l in range(1000):
                sum+=1.0
            #a[0]+=sum #this is the trouble line that makes the code slow.