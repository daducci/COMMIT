#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False
"""
Author: Matteo Frigo - lts5 @ EPFL and Dep. of CS @ Univ. of Verona

This structure is based on the previous work of Rafael Carrillo and was
supported by the LTS5 laboratory at EPFL, Lausanne.
"""
cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt


cpdef non_negativity(double [::1] x, int compartment_start, int compartment_size):
    """
    POCS for the first orthant (non-negativity)
    """
    cdef:
        int i
    for i in xrange(compartment_start, compartment_start+compartment_size):
        if x[i] <= 0.0 :
            x[i] = 0.0
    return np.asarray( x )


cpdef soft_thresholding(double [::1] x, double lam, int compartment_start, int compartment_size) :
    """
    Proximal of L1 norm
    """
    # NB: this preserves non-negativity
    cdef:
        int i
    for i in xrange(compartment_start, compartment_start+compartment_size):
        if x[i] <= lam:
            x[i] = 0.0
        else:
            x[i] = x[i] - lam
    return np.asarray( x )


cpdef projection_onto_l2_ball(double [::1] x, double lam, int compartment_start, int compartment_size) :
    """
    Proximal of L2 norm
    """
    # NB: this preserves non-negativity
    cdef:
        double xn = 0.0, k
        int i
    for i in xrange(compartment_start, compartment_start+compartment_size):
        xn += x[i]*x[i]
    xn = sqrt(xn)
    if xn > lam :
        k = 1. - lam/xn
        for i in xrange(compartment_start, compartment_start+compartment_size):
            x[i] = x[i]*k
    else :
        for i in xrange(compartment_start, compartment_start+compartment_size):
            x[i] = 0
    return np.asarray( x )


cpdef omega_group_sparsity(double [::1] x, int [::1] group_idx, int [::1] group_size, double [::1] group_weight, double lam, double n) :
    """
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    cdef:
        int nG = group_size.size, N
        int k, i, j = 0
        double omega = 0.0, gNorm, x_i

    if lam != 0:
        if n == 2:
            for k in xrange(nG):
                N = group_size[k]
                gNorm = 0.0
                for i in xrange(j,j+N) :
                    x_i = x[group_idx[i]]
                    gNorm += x_i*x_i
                omega += group_weight[k] * sqrt( gNorm )
                j += N
        elif n == np.inf:
            for k in xrange(nG):
                N = group_size[k]
                gNorm = x[group_idx[j]]
                for i in xrange(j+1,j+N) :
                    x_i = x[group_idx[i]]
                    if x_i > gNorm :
                        gNorm = x_i
                omega += group_weight[k] * gNorm
                j += N
    return lam*omega


cpdef prox_group_sparsity( double [::1] x, int [::1] group_idx, int [::1] group_size, double [::1] group_weight, double lam, double n ) :
    """
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    cdef:
        int nG = group_size.size, N
        int k, i, j = 0
        double wl, gNorm, x_i

    k = x.size
    for i in xrange(k):
        if x[i] <= 0.0:
            x[i] = 0.0

    if lam != 0:
        if n == 2 :
            for k in xrange(nG) :
                N = group_size[k]
                gNorm = 0.0
                for i in xrange(j,j+N) :
                    x_i = x[group_idx[i]]
                    gNorm += x_i*x_i
                gNorm = sqrt( gNorm )

                wl = group_weight[k] * lam
                if gNorm <= wl :
                    for i in xrange(j,j+N) :
                        x[ group_idx[i] ] = 0.0
                else :
                    wl = (gNorm-wl)/gNorm
                    for i in xrange(j,j+N) :
                        x[ group_idx[i] ] *= wl
                j += N
        # elif n == np.inf :
        # [TODO] TO be correctly implemented
        #     for k in range(nG) :
        #         idx = subtree[k]
        #         # xn = max( v[idx] )
        #         r = weight[k] * lam
        #         for i in idx :
        #             if v[i] <= r:
        #                 v[i] = 0.0
        #             else :
        #                 v[i] -= r
    return np.asarray( x )