# cython: profile=False
cimport cython
import numpy as np
cimport numpy as np
from math import sqrt
import sys

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.profile(False)

## Regularisers for NNLSL1
# Proximal
cpdef prox_nnl1(np.ndarray[np.float64_t] x, double lam, int compartment_start, int compartment_size) :
    """
    Author: Matteo Frigo - athena @ Inria
    """
    # NB: this preserves non-negativity
    cdef:
        np.ndarray[np.float64_t] v
        size_t i
    v = x.copy()
    for i in range(compartment_start, compartment_size):
        if v[i] <= lam:
            v[i] = 0.0
        else:
            v[i] -= lam
    return v


## Regularisers for HNNLS
# Penalty term
cpdef omega_hierarchical(np.ndarray[np.float64_t] v, np.ndarray[object] subtree, np.ndarray[np.float64_t] weight, double lam, double n) :
    """
    Author: Matteo Frigo - athena @ Inria
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    cdef:
        int nG = weight.size
        size_t k, i
        double xn, tmp = 0.0

    if lam != 0:
        if n == 2:
            for k in range(nG):
                idx = subtree[k]
                xn = 0.0
                for i in idx:
                    xn += v[i]*v[i]
                    tmp += weight[k] * sqrt( xn )
        elif n == np.Inf:
            for k in range(nG):
                idx = subtree[k]
                tmp += weight[k] * max( v[idx] )
    return lam*tmp

# Proximal operator of the penalty term
cpdef prox_hierarchical( np.ndarray[np.float64_t] x, np.ndarray[object] subtree, np.ndarray[np.float64_t] weight, double lam, double n ) :
    """
    Author: Matteo Frigo - athena @ Inria
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    cdef:
        np.ndarray[np.float64_t] v
        int nG = weight.size, N, rho
        size_t k, i
        double r, xn, theta

    v = x.copy()
    v[v<0] = 0.0

    if lam != 0:
        if n == np.inf :
            for k in range(nG) :
                idx = subtree[k]
                # xn = max( v[idx] )
                r = weight[k] * lam
                for i in idx :
                    if v[i] <= r:
                        v[i] = 0.0
                    else :
                        v[i] -= r
        if n == 2:
            for k in range(nG):
                idx = subtree[k]
                xn = 0.0
                for i in idx:
                    xn += v[i]*v[i]
                    xn = sqrt(xn)
                r = weight[k] * lam
                if xn > r:
                    r = (xn-r)/xn
                    for i in idx :
                        v[i] *= r
                else:
                    for i in idx:
                        v[i] = 0.0
    return v
