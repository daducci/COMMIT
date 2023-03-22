#!python
# cython: boundscheck=False, wraparound=False, profile=False
import cython
from cython.operator import dereference
from libcpp.vector cimport vector
cimport numpy as np

import numpy as np
import nibabel

from dipy.io.stateful_tractogram import StatefulTractogram as sft
from dipy.io.streamline import save_tractogram


# Interface to actual C code
cdef extern from "spline_smoothing_c.cpp":
    int do_spline_smoothing(
        float* ptr_npaFiberI, int nP, float* ptr_npaFiberO, float ratio, float segment_len
        ) nogil

# cdef extern from "spline_smoothing_c.cpp":
#     int compute_control_p(
#         float* ptr_npaFiberI, int nP, float* ptr_npaFiberO, float ratio, float segment_len
#         ) nogil
    

# cpdef compute_control_points(Sft_IN, segment_len = 1, control_point_ratio = 0.25, verbose = False, tol = 1):
#     if control_point_ratio < 0 or control_point_ratio > 1 :
#         raise ValueError( "'control_point_ratio' parameter must be in [0..1]" )

#     # create the structure for the input and output polyline
#     cdef float [:, ::1] origFib
#     cdef float* ptr_origFib
#     cdef float [:, ::1] knotsFib = np.ascontiguousarray( np.zeros( (3*10000,1) ).astype(np.float32) )
#     cdef float* ptr_knotsFib = &knotsFib[0,0]

    
#     control_p= []
#     original_fib = []
    
#     count_not_smoothed = 0
#     for f in Sft_IN.streamlines :
#         streamline = np.ascontiguousarray( f.copy())
#         origFib = streamline
#         ptr_origFib = &origFib[0,0]

#         n = compute_control_p( ptr_origFib, f.shape[0], ptr_knotsFib, control_point_ratio, segment_len, tol )
#         if n != 0 :
#             streamline = np.reshape( knotsFib[:3*n].copy(), (n,3) )
#             control_p.append( streamline )
            
#             original_fib.append(f)
            
#         else :
#             # control_p.append( f )
#             count_not_smoothed = count_not_smoothed + 1

#     Reduced_sft = sft(control_p, Sft_IN, Sft_IN.space)
#     save_tractogram(Reduced_sft, "Reduced_sft_configuration.trk", bbox_valid_check=True)

#     return Reduced_sft

cpdef perform_smoothing(Sft_IN, REF, segment_len = 1, control_point_ratio = 1, verbose = False, tol = 1) :

    """Perform the conversion of a tractoram to the sparse data-structure internally
    used by COMMIT to perform the matrix-vector multiplications with the operator A
    during the inversion of the linear system.

    Parameters
    ----------
    trk_filename_in : string
        Path to the .trk file containing the tractogram to process.

    trk_filename_out : string
        Path to the .trk where to store the filtered tractogram.

    control_point_ratio : float
        Percent of control points to use in the interpolating spline (default : 0.25).

    segment_len : float
        Sampling resolution of the final streamline after interpolation (default : 1.0).

    verbose : boolean
        Print information and progess (default : False).
    """

    if control_point_ratio < 0 or control_point_ratio > 1 :
        raise ValueError( "'control_point_ratio' parameter must be in [0..1]" )

    # create the structure for the input and output polyline
    cdef float [:, ::1] npaFiberI
    cdef float* ptr_npaFiberI
    cdef float [:, ::1] npaFiberO = np.ascontiguousarray( np.zeros( (3*10000,1) ).astype(np.float32) )
    cdef float* ptr_npaFiberO = &npaFiberO[0,0]

    
    trk_fiber_out = []
    old_coord = []
    
    count_not_smoothed = 0
    for f in Sft_IN.streamlines :
        streamline = np.ascontiguousarray( f.copy())
        npaFiberI = streamline
        ptr_npaFiberI = &npaFiberI[0,0]

        n = do_spline_smoothing( ptr_npaFiberI, f.shape[0], ptr_npaFiberO, control_point_ratio, segment_len )

        if n != 0 :
            streamline = np.reshape( npaFiberO[:3*n].copy(), (n,3) )
            trk_fiber_out.append( streamline )
            
            old_coord.append(f)
            
        else :
            trk_fiber_out.append( f )
            count_not_smoothed = count_not_smoothed + 1
        

    if count_not_smoothed >0 :
        print '%d streamlines not smoothed due to \"segment_len\" > length of streamline' % count_not_smoothed

    
    Start_sft = sft.from_sft(trk_fiber_out, Sft_IN)
    # save_tractogram(Start_sft, "Smoothed_start_configuration.trk", bbox_valid_check=True)
    
    # Old_sft = sft.from_sft(old_coord, REF)
    # save_tractogram(Sft_IN, "NON_Smoothed_start_configuration.trk", bbox_valid_check=True)

    # return trk_fiber_out, old_coord
    return Start_sft


    