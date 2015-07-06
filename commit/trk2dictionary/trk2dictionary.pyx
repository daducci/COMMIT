#!python
# cython: boundscheck=False, wraparound=False, profile=False

import cython
import numpy as np
cimport numpy as np
import nibabel
from os.path import join, exists
from os import makedirs
import time

# Interface to actual C code
cdef extern from "trk2dictionary_c.cpp":
    int trk2dictionary(
        char* strTRKfilename, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars, int n_properties, float fiber_shift, int points_to_skip,
        float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
        float* _ptrMASK, float* ptrTDI, char* path_out, int c
    ) nogil


cpdef run( filename_trk, path_out, filename_peaks = None, filename_mask = None, do_intersect = True, fiber_shift = 0, points_to_skip = 0, vf_THR = 0.1, flip_peaks = [False,False,False] ):
    """Perform the conversion of a tractoram to the sparse data-structure internally
    used by COMMIT to perform the matrix-vector multiplications with the operator A
    during the inversion of the linear system.

    Parameters
    ----------
    filename_trk : string
        Path to the .trk file containing the tractogram to convert.

    path_out : string
        Path to the folder where to store the sparse data structure.

    filename_peaks : string
        Path to the NIFTI file containing the peaks to use as extra-cellular contributions.
        The data matrix should be 4D with last dimension 3*N, where N is the number
        of peaks in each voxel. (default : no extra-cellular contributions)

    filename_mask : string
        Path to a binary mask to use to restrict the analysis to specific areas.

    do_intersect : boolean
        If True then fiber segments that intersect voxel boundaries are splitted (default).
        If False then the centroid of the segment is used as its voxel position.

    fiber_shift : float
        If necessary, apply a translation to fiber coordinates (default : 0) to account
        for differences between the reference system of the tracking algorithm and COMMIT.
        The value is specified in voxel units, eg 0.5 translates by half voxel.

    points_to_skip : integer
        If necessary, discard first points at beginning/end of a fiber (default : 0).

    vf_THR : float
        Discard peaks smaller than vf_THR * max peak (default : 0.1).

    flip_peaks : list of three boolean
        If necessary, flips peak orientations along each axis (default : no flipping).
    """
    tic = time.time()
    print '\n-> Creating the dictionary from tractogram:'
    print '\t* Segment position = %s' % ( 'COMPUTE INTERSECTIONS' if do_intersect else 'CENTROID' )
    print '\t* Fiber shift      = %.3f (voxel-size units)' % fiber_shift
    print '\t* Points to skip   = %d' % points_to_skip

    print '\t* Loading data:'

    # fiber-tracts from .trk
    print '\t\t* tractogram'
    try :
        _, trk_hdr = nibabel.trackvis.read( filename_trk, as_generator=True )
    except :
        raise IOError( 'Track file not found' )
    Nx = trk_hdr['dim'][0]
    Ny = trk_hdr['dim'][1]
    Nz = trk_hdr['dim'][2]
    Px = trk_hdr['voxel_size'][0]
    Py = trk_hdr['voxel_size'][1]
    Pz = trk_hdr['voxel_size'][2]
    print '\t\t\t- %d x %d x %d' % ( Nx, Ny, Nz )
    print '\t\t\t- %.4f x %.4f x %.4f' % ( Px, Py, Pz )
    print '\t\t\t- %d fibers' % trk_hdr['n_count']

    # white-matter mask
    cdef float* ptrMASK
    cdef float [:, :, ::1] niiMASK_img
    if filename_mask is not None :
        print '\t\t* filtering mask'
        niiMASK = nibabel.load( filename_mask )
        niiMASK_hdr = niiMASK.header if nibabel.__version__ >= '2.0.0' else niiMASK.get_header()
        print '\t\t\t- %d x %d x %d' % ( niiMASK.shape[0], niiMASK.shape[1], niiMASK.shape[2] )
        print '\t\t\t- %.4f x %.4f x %.4f' % ( niiMASK_hdr['pixdim'][1], niiMASK_hdr['pixdim'][2], niiMASK_hdr['pixdim'][3] )
        if ( Nx!=niiMASK.shape[0] or Ny!=niiMASK.shape[1] or Nz!=niiMASK.shape[2] or
             abs(Px-niiMASK_hdr['pixdim'][1])>1e-3 or abs(Py-niiMASK_hdr['pixdim'][2])>1e-3 or abs(Pz-niiMASK_hdr['pixdim'][3])>1e-3 ) :
            print '\t\t  [WARNING] WM dataset does not have the same geometry as the tractogram'
        niiMASK_img = np.ascontiguousarray( niiMASK.get_data().astype(np.float32) )
        ptrMASK  = &niiMASK_img[0,0,0]
    else :
        print '\t\t* no mask specified to filter IC compartments'
        ptrMASK = NULL

    # peaks file for EC contributions
    cdef float* ptrPEAKS
    cdef float [:, :, :, ::1] niiPEAKS_img
    cdef int Np
    if filename_peaks is not None :
        print '\t\t* EC orientations'
        niiPEAKS = nibabel.load( filename_peaks )
        niiPEAKS_hdr = niiPEAKS.header if nibabel.__version__ >= '2.0.0' else niiPEAKS.get_header()
        print '\t\t\t- %d x %d x %d x %d' % ( niiPEAKS.shape[0], niiPEAKS.shape[1], niiPEAKS.shape[2], niiPEAKS.shape[3] )
        print '\t\t\t- %.4f x %.4f x %.4f' % ( niiPEAKS_hdr['pixdim'][1], niiPEAKS_hdr['pixdim'][2], niiPEAKS_hdr['pixdim'][3] )
        print '\t\t\t- ignoring peaks < %.2f * MaxPeak' % vf_THR
        print '\t\t\t- flipping axes : [ x=%s, y=%s, z=%s ]' % ( flip_peaks[0], flip_peaks[1], flip_peaks[2] )
        if ( Nx!=niiPEAKS.shape[0] or Ny!=niiPEAKS.shape[1] or Nz!=niiPEAKS.shape[2] or
             abs(Px-niiPEAKS_hdr['pixdim'][1])>1e-3 or abs(Py-niiPEAKS_hdr['pixdim'][2])>1e-3 or abs(Pz-niiPEAKS_hdr['pixdim'][3])>1e-3 ) :
            print "\t\t  [WARNING] PEAKS dataset does not have the same geometry as the tractogram"
        if niiPEAKS.shape[3] % 3 :
            raise RuntimeError( 'PEAKS dataset must have 3*k volumes' )
        if vf_THR < 0 or vf_THR > 1 :
            raise RuntimeError( 'vf_THR must be between 0 and 1' )
        niiPEAKS_img = np.ascontiguousarray( niiPEAKS.get_data().astype(np.float32) )
        ptrPEAKS = &niiPEAKS_img[0,0,0,0]
        Np = niiPEAKS.shape[3]/3
    else :
        print '\t\t* no dataset specified for EC compartments'
        Np = 0
        ptrPEAKS = NULL

    # output path
    print '\t\t* output written to "%s"' % path_out
    if not exists( path_out ):
        makedirs( path_out )


    # create TDI mask
    cdef float [:, :, ::1] niiTDI_img = np.ascontiguousarray( np.zeros((Nx,Ny,Nz),dtype=np.float32) )
    cdef float* ptrTDI  = &niiTDI_img[0,0,0]

    # calling actual C code
    ret = trk2dictionary( filename_trk,
        trk_hdr['dim'][0], trk_hdr['dim'][1], trk_hdr['dim'][2],
        trk_hdr['voxel_size'][0], trk_hdr['voxel_size'][1], trk_hdr['voxel_size'][2],
        trk_hdr['n_count'], trk_hdr['n_scalars'], trk_hdr['n_properties'], fiber_shift, points_to_skip,
        ptrPEAKS, Np, vf_THR, -1 if flip_peaks[0] else 1, -1 if flip_peaks[1] else 1, -1 if flip_peaks[2] else 1,
        ptrMASK, ptrTDI, path_out, 1 if do_intersect else 0 );
    if ret == 0 :
        print '   [ DICTIONARY not generated ]'
        return None
    print '   [ %.1f seconds ]' % ( time.time() - tic )

    # save TDI map
    affine = None
    if ptrPEAKS != NULL :
        affine = niiPEAKS.affine if nibabel.__version__ >= '2.0.0' else niiPEAKS.get_affine()
    niiTDI = nibabel.Nifti1Image( niiTDI_img, affine )
    nibabel.save( niiTDI, join(path_out,'dictionary_tdi.nii.gz') )
