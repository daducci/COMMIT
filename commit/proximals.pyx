#!python
#cython: boundscheck=False, wraparound=False
"""
Author: Matteo Frigo - lts5 @ EPFL and Dep. of CS @ Univ. of Verona

This structure is based on the previous work of Rafael Carrillo and was
supported by the LTS5 laboratory at EPFL, Lausanne.
"""
cimport cython
import numpy as np
cimport numpy as np
from math import sqrt
import sys


cpdef non_negativity(np.ndarray[np.float64_t] x, int compartment_start, int compartment_size):
    """
    POCS for the first orthant (non-negativity)
    """
    cdef:
        np.ndarray[np.float64_t] v
        size_t i
    v = x.copy()
    for i in range(compartment_start, compartment_start+compartment_size):
        if v[i] < 0.0:
            v[i] = 0.0
    return v


cpdef soft_thresholding(np.ndarray[np.float64_t] x, double lam, int compartment_start, int compartment_size) :
    """
    Proximal of L1 norm
    """
    # NB: this preserves non-negativity
    cdef:
        np.ndarray[np.float64_t] v
        size_t i
    v = x.copy()
    for i in range(compartment_start, compartment_start+compartment_size):
        if v[i] <= lam:
            v[i] = 0.0
        else:
            v[i] -= lam
    return v


cpdef projection_onto_l2_ball(np.ndarray[np.float64_t] x, double lam, int compartment_start, int compartment_size) :
    """
    Proximal of L2 norm
    """
    # NB: this preserves non-negativity
    cdef:
        np.float64_t xn
        np.ndarray[np.float64_t] v
        size_t i
    v = x.copy()
    xn = sqrt(sum(v[compartment_start:compartment_start+compartment_size]**2))
    if xn > lam:
        for i in range(compartment_start, compartment_start+compartment_size):
            v[i] = v[i]/xn*lam
    return v


cpdef omega_group_sparsity(np.ndarray[np.float64_t] v, np.ndarray[object] subtree, np.ndarray[np.float64_t] weight, double lam, double n) :
    """
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    cdef:
        int nG = weight.size
        size_t k
        double tmp = 0.0

    if lam != 0:
        if n == 2:
            for k in range(nG):
                idx = subtree[k]
                tmp += weight[k] * sqrt( sum(v[idx]**2) )
        elif n == np.Inf:
            for k in range(nG):
                idx = subtree[k]
                tmp += weight[k] * max( v[idx] )
    return lam*tmp


cpdef prox_group_sparsity( np.ndarray[np.float64_t] x, np.ndarray[object] subtree, np.ndarray[np.float64_t] weight, double lam, double n ) :
    """
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    cdef:
        np.ndarray[np.float64_t] v
        int nG = weight.size
        size_t k, i
        double r, xn

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
                xn = sqrt( sum(v[idx]**2) )
                r = weight[k] * lam
                if xn > r:
                    r = (xn-r)/xn
                    for i in idx :
                        v[i] *= r
                else:
                    for i in idx:
                        v[i] = 0.0
    return v
