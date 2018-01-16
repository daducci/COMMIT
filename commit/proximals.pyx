cimport cython
import numpy as np
cimport numpy as np
from math import sqrt
import sys

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.profile(False)

#####  #####   ####  #    # # #    #   ##   #       ####
#    # #    # #    #  #  #  # ##  ##  #  #  #      #
#    # #    # #    #   ##   # # ## # #    # #       ####
#####  #####  #    #   ##   # #    # ###### #           #
#      #   #  #    #  #  #  # #    # #    # #      #    #
#      #    #  ####  #    # # #    # #    # ######  ####
cpdef non_negativity(np.ndarray[np.float64_t] x, int compartment_start, int compartment_size):
    """
    POCS for the first orthant (non-negativity)

    Author: Matteo Frigo - athena @ Inria
    """
    cdef:
        np.ndarray[np.float64_t] v
        size_t i
    v = x.copy()
    for i in range(compartment_start, compartment_size):
        if v[i] < 0.0:
            v[i] = 0.0
    return v

cpdef soft_thresholding(np.ndarray[np.float64_t] x, double lam, int compartment_start, int compartment_size) :
    """
    Proximal of L1 norm

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

cpdef projection_onto_l2_ball(np.ndarray[np.float64_t] x, double lam, int compartment_start, int compartment_size) :
    """
    Proximal of L2 norm

    Author: Matteo Frigo - athena @ Inria
    """
    # NB: this preserves non-negativity
    cdef:
        np.float64_t xn
        np.ndarray[np.float64_t] v
        size_t i
    v = x.copy()
    xn = sqrt(sum(v[compartment_start:compartment_size]**2))
    if xn > lam
        for i in range(compartment_start, compartment_size):
            v[i] = v[i]/xn*lam
    return v

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
