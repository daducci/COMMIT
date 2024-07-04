#!python
# cython: language_level=3, c_string_type=str, c_string_encoding=ascii, boundscheck=False, wraparound=False, profile=False

from libc.stdlib cimport malloc, free
from libcpp cimport bool
cimport numpy as np

import nibabel
import numpy as np

import os
from os.path import join, exists, splitext, dirname, isdir, isfile
from os import makedirs, remove

import amico

from dicelib.clustering import run_clustering
from dicelib.ui import _in_notebook
from dicelib.ui import ProgressBar, setup_logger
from dicelib import ui
from dicelib.utils import format_time

import pickle

from pkg_resources import get_distribution

import shutil

import time

 
logger = setup_logger('trk2dictionary')

# Interface to actual C code
cdef extern from "trk2dictionary_c.cpp":
    int trk2dictionary(
        char* filename_tractogram, int data_offset, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars,
        int n_properties, float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len,  float max_fiber_len,
        float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
        float* _ptrMASK, float* _ptrISO, double** ptrTDI, char* path_out, int c, double* ptrPeaksAffine,
        int nReplicas, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, bool* ptrBlurApplyTo,
        float* ptrTractsAffine, short* prtHashTable, int threads_count, int verbose
    ) nogil

def _get_header( niiFILE ):
    return niiFILE.header if nibabel.__version__ >= '2.0.0' else niiFILE.get_header()

def _get_affine( niiFILE ):
    return niiFILE.affine if nibabel.__version__ >= '2.0.0' else niiFILE.get_affine()

cpdef cat_function( infilename, outfilename ):
    """ Concatenate binary file """
    with open( outfilename, "wb" ) as outfile:
        for fname in infilename:
            with open( fname, 'rb' ) as inFile:
                shutil.copyfileobj( inFile, outfile )
            remove( fname )

cpdef compute_tdi( np.uint32_t[::1] v, np.float32_t[::1] l, int nx, int ny, int nz, int verbose):
    cdef np.float32_t [::1] tdi = np.zeros( nx*ny*nz, dtype=np.float32 )
    cdef unsigned long long i=0
    with ProgressBar(total=v.size, disable=verbose<3, hide_on_exit=True) as pbar:

        for i in range(v.size):
            tdi[ v[i] ] += l[i]
            pbar.update()

    return tdi


cpdef run( filename_tractogram=None, path_out=None, filename_peaks=None, filename_mask=None, filename_ISO=None,
            do_intersect=True, fiber_shift=0, min_seg_len=1e-3, min_fiber_len=0.0, max_fiber_len=250.0,
            vf_THR=0.1, peaks_use_affine=False, flip_peaks=[False,False,False], blur_clust_groupby=None,
            blur_clust_thr=0, blur_spacing=0.25, blur_core_extent=0.0, blur_gauss_extent=0.0,
            blur_gauss_min=0.1, blur_apply_to=None, TCK_ref_image=None, ndirs=500, n_threads=-1,
            keep_temp=False, verbose=3
            ):
    """Perform the conversion of a tractoram to the sparse data-structure internally
    used by COMMIT to perform the matrix-vector multiplications with the operator A
    during the inversion of the linear system.

    Parameters
    ----------
    filename_tractogram : string
        Path to the tractogram (.trk or .tck) containing the streamlines to load.

    TCK_ref_image: string
        When loading a .tck tractogram, path to the NIFTI file containing the information about
        the geometry to be used for the tractogram to load. If not specified, it will try to use
        the information from filename_peaks or filename_mask.

    path_out : string
        Path to the folder for storing the sparse data structure. If not specified (default),
        a folder name "COMMIT" will be created in the same folder of the tractogram.

    filename_mask : string
        Path to a binary mask for restricting the analysis to specific areas.
        Segments outside this mask are discarded. If not specified (default),
        the mask is created from all voxels intersected by the tracts.

    filename_ISO : string
        Path to a binary mask for WIP models.

    do_intersect : boolean
        If True then streamline segments that intersect voxel boundaries are splitted (default).
        If False then the centroid of the segment is used as its voxel position.

    fiber_shift : float or list of three float
        If necessary, apply a translation to streamline coordinates (default : 0) to account
        for differences between the reference system of the tracking algorithm and COMMIT.
        The value is specified in voxel units, eg 0.5 translates by half voxel.

    min_seg_len : float
        Discard segments <= than this length in mm (default : 1e-3).

    min_fiber_len : float
        Discard streamlines <= than this length in mm (default : 0.0).

    max_fiber_len : float
        Discard streamlines >= than this length in mm (default : 250.0).

    filename_peaks : string
        Path to the NIFTI file containing the peaks to use as extra-cellular contributions.
        The data matrix should be 4D with last dimension 3*N, where N is the number
        of peaks in each voxel. (default : no extra-cellular contributions).

    peaks_use_affine : boolean
        Whether to rotate the peaks according to the affine matrix (default : False).

    vf_THR : float
        Discard peaks smaller than vf_THR * max peak (default : 0.1).

    flip_peaks : list of three boolean
        If necessary, flips peak orientations along each axis (default : no flipping).

    blur_spacing : float
        To obtain the blur effect, streamlines are duplicated and organized in a cartesian grid;
        this parameter controls the spacing of the grid in mm (defaut : 0.25).

    blur_core_extent: float
        Extent of the core inside which the segments have equal contribution to the central one (default : 0.0).

    blur_gauss_extent: float
        Extent of the gaussian damping at the border (default : 0.0).

    blur_gauss_min: float
        Minimum value of the Gaussian to consider when computing the sigma (default : 0.1).

    blur_apply_to: array of bool
        For each input streamline, decide whether blur is applied or not to it (default : None, meaning apply to all).

    blur_clust_thr: float
        Clustering threshold used to remove redundant streamlines from the input tractogram

    ndirs : int
        Number of orientations on the sphere used to discretize the orientation of each
        each segment in a streamline (default : 500).

    n_threads: int
        Number of threads. If nothing is specified, the value used is the number of CPUs available
        in the system (default : -1).

    keep_temp: boolean
        If True, the temporary files are not deleted (default : False).

    verbose: int
        The verbosity level (0: only errors, 1: errors and warnings, 2: errors, warnings and info, 3: errors, warnings, info and progress bars, 4: errors, warnings, info, progress bars and debug)
    """

    ui.set_verbose('trk2dictionary' ,verbose)

    # check the value of ndirs
    if not amico.lut.is_valid(ndirs):
        logger.error( 'Unsupported value for ndirs.\nNote: Supported values for ndirs are [500 (default), 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 32761]' )

    # check conflicts of fiber_shift
    if np.isscalar(fiber_shift) :
        fiber_shiftX = fiber_shift
        fiber_shiftY = fiber_shift
        fiber_shiftZ = fiber_shift
    elif len(fiber_shift) == 3 :
        fiber_shiftX = fiber_shift[0]
        fiber_shiftY = fiber_shift[1]
        fiber_shiftZ = fiber_shift[2]
    else :
        logger.error( '"fiber_shift" must be a scalar or a vector with 3 elements' )

    # check for invalid parameters in the blur
    if blur_core_extent < 0 :
        logger.error( 'The extent of the core must be non-negative' )

    if blur_gauss_extent < 0 :
        logger.error( 'The extent of the blur must be non-negative' )

    if blur_gauss_extent > 0 or blur_core_extent > 0:
        if blur_spacing <= 0 :
            logger.error( 'The grid spacing of the blur must be positive' )

    tic = time.time()
    logger.info( 'Creating data structure from tractogram' )

    logger.subinfo( 'Configuration:', indent_char='*', indent_lvl=1 )
    logger.subinfo( f'Segment position = {"COMPUTE INTERSECTIONS" if do_intersect else "CENTROID"}', indent_lvl=2, indent_char='-' )
    logger.subinfo( f'Coordinates shift in X = {fiber_shiftX:.3f} (voxel-size units)', indent_lvl=2, indent_char='-' )
    logger.subinfo( f'Coordinates shift in Y = {fiber_shiftY:.3f} (voxel-size units)', indent_lvl=2, indent_char='-' )
    logger.subinfo( f'Coordinates shift in Z = {fiber_shiftZ:.3f} (voxel-size units)', indent_lvl=2, indent_char='-' )
    if min_seg_len >= 1e-3:
        logger.subinfo( f'Min segment len    = {min_seg_len:.3f} mm', indent_lvl=2, indent_char='-' )
    else:
        logger.subinfo( f'Min segment len    = {min_seg_len:.2e} mm', indent_lvl=2, indent_char='-' )
    logger.subinfo( f'Min streamline len = {min_fiber_len:.2f} mm', indent_lvl=2, indent_char='-' )
    logger.subinfo( f'Max streamline len = {max_fiber_len:.2f} mm', indent_lvl=2, indent_char='-' )

    # check blur params
    cdef :
        double [:] blurRho
        double [:] blurAngle
        double [:] blurWeights
        bool [:] blurApplyTo
        int nReplicas
        float blur_sigma
        int i = 0

    if (blur_gauss_extent==0 and blur_core_extent==0) or (blur_spacing==0) :
        nReplicas = 1
        blurRho = np.array( [0.0], np.double )
        blurAngle = np.array( [0.0], np.double )
        blurWeights = np.array( [1], np.double )
    else:
        tmp = np.arange(0,blur_core_extent+blur_gauss_extent+1e-6,blur_spacing)
        tmp = np.concatenate( (tmp,-tmp[1:][::-1]) )
        x, y = np.meshgrid( tmp, tmp )
        r = np.sqrt( x*x + y*y )
        idx = (r <= blur_core_extent+blur_gauss_extent)
        blurRho = r[idx]
        blurAngle = np.arctan2(y,x)[idx]
        nReplicas = blurRho.size

        blurWeights = np.empty( nReplicas, np.double  )
        if blur_gauss_extent == 0 :
            blurWeights[:] = 1.0
        else:
            blur_sigma = blur_gauss_extent / np.sqrt( -2.0 * np.log( blur_gauss_min ) )
            for i in xrange(nReplicas):
                if blurRho[i] <= blur_core_extent :
                    blurWeights[i] = 1.0
                else:
                    blurWeights[i] = np.exp( -(blurRho[i] - blur_core_extent)**2 / (2.0*blur_sigma**2) )

    if nReplicas == 1 :
        logger.subinfo( 'Do not blur streamlines', indent_lvl=2, indent_char='-' )
    else :
        logger.subinfo( 'Blur streamlines:', indent_lvl=2, indent_char='-' )
        logger.subinfo( f'core extent  = {blur_core_extent:.3f}', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'gauss extent = {blur_gauss_extent:.3f} (sigma = {blur_sigma:.3f})', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'grid spacing = {blur_spacing:.3f}' , indent_lvl=3, indent_char='-' )
        logger.subinfo( f'weights = [ {np.min(blurWeights):.3f} ... {np.max(blurWeights):.3f} ]', indent_lvl=3, indent_char='-' )

    if min_seg_len < 0 :
        logger.error( '"min_seg_len" must be >= 0' )
    if min_fiber_len < 0 :
        logger.error( '"min_fiber_len" must be >= 0' )
    if max_fiber_len < min_fiber_len :
        logger.error( '"max_fiber_len" must be >= "min_fiber_len"' )

    if filename_tractogram is None:
        logger.error( '"filename_tractogram" not defined' )

    if path_out is None:
        path_out = dirname(filename_tractogram)
        if path_out == '':
            path_out = '.'
        if not isdir(path_out):
            logger.error( '"path_out" cannot be inferred from "filename_tractogram"' )
        path_out = join(path_out,'COMMIT')

    # create output path
    logger.subinfo( f'Output written to "{path_out}"', indent_char='-', indent_lvl=2 )
    if not exists( path_out ):
        makedirs( path_out )

    path_temp = join(path_out, 'temp')
    logger.subinfo( f'Temporary files written to "{path_temp}"', indent_char='-', indent_lvl=2 )

    if exists(path_temp):
        shutil.rmtree(path_temp, ignore_errors=True)
    makedirs(path_temp, exist_ok=True)

    if n_threads == 0 or n_threads > 255 :
        logger.error( 'Number of n_threads must be between 1 and 255' )

    if n_threads == -1 :
        # Set to the number of CPUs in the system
        try :
            n_threads = os.cpu_count()
        except :
            n_threads = 1

    logger.subinfo( f'Using parallel computation with {n_threads} threads', indent_char='-', indent_lvl=2 )

    if np.isscalar(blur_clust_thr):
        blur_clust_thr = np.array( [blur_clust_thr] )

    if blur_clust_thr[0]> 0:
        logger.subinfo( 'Reducing streamlines redundancy:', indent_char='*', indent_lvl=1)
        logger.subinfo( f'Input tractogram "{filename_tractogram}"', indent_lvl=2, indent_char='-' )
        input_hdr = nibabel.streamlines.load( filename_tractogram, lazy_load=True ).header
        input_n_count = int(input_hdr['count'])

        extension = splitext(filename_tractogram)[1]

        if filename_mask is None and TCK_ref_image is None:
            if extension == ".tck":
                logger.error( 'TCK files do not contain information about the geometry. Use "TCK_ref_image" for that' )
            else:
                logger.error( 'Unknown file extension. Use "filename_mask" or "TCK_ref_image" for that' )

        input_tractogram = os.path.basename(filename_tractogram)[:len(os.path.basename(filename_tractogram))-4]
        filename_out = join( path_out, f'{input_tractogram}_clustered_thr_{float(blur_clust_thr[0])}.tck' )

        if blur_clust_groupby:
            hdr = nibabel.streamlines.load( filename_tractogram, lazy_load=True ).header
            temp_idx = np.arange(int(hdr['count']))
            log_list = []
            ret_subinfo = logger.subinfo(f'Clustering with threshold = {blur_clust_thr[0]}', indent_lvl=2, indent_char='-', with_progress=verbose>2)
            with ProgressBar(disable=verbose<3, hide_on_exit=True, subinfo=ret_subinfo, log_list=log_list):
                idx_centroids = run_clustering(tractogram_in=filename_tractogram, tractogram_out=filename_out,
                                            temp_folder=path_temp, atlas=blur_clust_groupby, clust_thr=blur_clust_thr[0],
                                            n_threads=n_threads, keep_temp_files=True, force=True, verbose=1, log_list=log_list)
        else:
            log_list = []
            ret_subinfo = logger.subinfo(f'Clustering with threshold = {blur_clust_thr[0]}', indent_lvl=2, indent_char='-', with_progress=verbose>2)
            with ProgressBar(disable=verbose<3, hide_on_exit=True, subinfo=ret_subinfo, log_list=log_list):
                idx_centroids = run_clustering(tractogram_in=filename_tractogram, tractogram_out=filename_out,
                                            temp_folder=path_temp, clust_thr=blur_clust_thr[0],
                                            keep_temp_files=True, force=True, verbose=1)
        filename_tractogram = filename_out


    # Load data from files
    logger.subinfo( 'Loading data:', indent_char='*', indent_lvl=1 )
    cdef short [:] htable = amico.lut.load_precomputed_hash_table(ndirs)
    cdef short* ptrHashTable = &htable[0]

    # Streamlines from tractogram
    logger.subinfo( 'Tractogram:', indent_lvl=2, indent_char='-' )

    if not exists(filename_tractogram):
        logger.error( f'Tractogram file not found: {filename_tractogram}' )
    extension = splitext(filename_tractogram)[1]
    if extension != ".trk" and extension != ".tck":
        logger.error( 'Invalid input file: only .trk and .tck are supported' )

    hdr = nibabel.streamlines.load( filename_tractogram, lazy_load=True ).header


    if extension == ".trk":
        logger.subinfo ( f'geometry taken from "{filename_tractogram}"', indent_lvl=3, indent_char='-' )
        Nx = int(hdr['dimensions'][0])
        Ny = int(hdr['dimensions'][1])
        Nz = int(hdr['dimensions'][2])
        Px = hdr['voxel_sizes'][0]
        Py = hdr['voxel_sizes'][1]
        Pz = hdr['voxel_sizes'][2]

        data_offset = 1000
        n_count = hdr['nb_streamlines']
        n_scalars = hdr['nb_scalars_per_point']
        n_properties = hdr['nb_properties_per_streamline']

    else:
        if TCK_ref_image is None:
            if filename_peaks is not None:
                TCK_ref_image = filename_peaks
            elif filename_mask is not None:
                TCK_ref_image = filename_mask
            else:
                logger.error( 'TCK files do not contain information about the geometry. Use "TCK_ref_image" for that' )
        logger.subinfo ( f'geometry taken from "{TCK_ref_image}"', indent_lvl=3, indent_char='-' )

        niiREF = nibabel.load( TCK_ref_image )
        niiREF_hdr = _get_header( niiREF )
        Nx = niiREF.shape[0]
        Ny = niiREF.shape[1]
        Nz = niiREF.shape[2]
        Px = niiREF_hdr['pixdim'][1]
        Py = niiREF_hdr['pixdim'][2]
        Pz = niiREF_hdr['pixdim'][3]
        data_offset = int(hdr['_offset_data'])
        n_count = int(hdr['count'])
        n_scalars = 0
        n_properties = 0

    if n_threads > n_count:
        logger.warning( 'Reducing the number of threads to the number of streamlines' )
        n_threads = n_count

    logger.subinfo( f'{Nx} x {Ny} x {Nz}', indent_lvl=3, indent_char='-' )
    logger.subinfo( f'{Px:.4f} x {Py:.4f} x {Pz:.4f}', indent_lvl=3, indent_char='-' )
    if blur_clust_thr[0]> 0:
        logger.subinfo( f'{input_n_count} streamlines', indent_lvl=3, indent_char='-' )
    else:
        logger.subinfo( f'{n_count} streamlines', indent_lvl=3, indent_char='-' )
    if Nx >= 2**16 or Nz >= 2**16 or Nz >= 2**16 :
        logger.error( 'The max dim size is 2^16 voxels' )

    # check copmpatibility between blurApplyTo and number of streamlines
    if blur_apply_to is None:
        blur_apply_to = np.repeat([True], n_count)
    else :
        if blur_apply_to.size != n_count :
            logger.error( '"blur_apply_to" must have one value per streamline' )
        logger.subinfo( f'{sum(blur_apply_to)} blurred streamlines', indent_lvl=3, indent_char='-' )
    blurApplyTo = blur_apply_to

    # get toVOXMM matrix (remove voxel scaling from affine) in case of TCK
    cdef float [:] toVOXMM
    cdef float* ptrToVOXMM
    if extension == ".tck":
        M = _get_affine( niiREF ).copy()
        M[:3, :3] = M[:3, :3].dot( np.diag([1./Px,1./Py,1./Pz]) )
        toVOXMM = np.ravel(np.linalg.inv(M)).astype('<f4')
        ptrToVOXMM = &toVOXMM[0]

    # white-matter mask
    cdef float* ptrMASK
    cdef float [:, :, ::1] niiMASK_img
    if filename_mask is not None :
        logger.subinfo( 'Filtering mask:', indent_lvl=2, indent_char='-' )
        niiMASK = nibabel.load( filename_mask )
        niiMASK_hdr = _get_header( niiMASK )
        logger.subinfo( f'{niiMASK.shape[0]} x {niiMASK.shape[1]} x {niiMASK.shape[2]}', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'{niiMASK_hdr["pixdim"][1]:.4f} x {niiMASK_hdr["pixdim"][2]:.4f} x {niiMASK_hdr["pixdim"][3]:.4f}', indent_lvl=3, indent_char='-' )
        if ( Nx!=niiMASK.shape[0] or Ny!=niiMASK.shape[1] or Nz!=niiMASK.shape[2] or 
            abs(Px-niiMASK_hdr['pixdim'][1])>1e-3 or abs(Py-niiMASK_hdr['pixdim'][2])>1e-3 or abs(Pz-niiMASK_hdr['pixdim'][3])>1e-3 ) :
            logger.warning( 'Dataset does not have the same geometry as the tractogram' )
        niiMASK_img = np.ascontiguousarray( np.asanyarray( niiMASK.dataobj ).astype(np.float32) )
        ptrMASK  = &niiMASK_img[0,0,0]
    else :
        logger.subinfo( 'No mask specified to filter IC compartments', indent_lvl=2, indent_char='-' )
        ptrMASK = NULL

    # peaks file for EC contributions
    cdef float* ptrPEAKS
    cdef float [:, :, :, ::1] niiPEAKS_img
    cdef int Np
    cdef double [:,:, :,::1] niiTDI_img = np.ascontiguousarray( np.zeros((n_threads,Nx,Ny,Nz),dtype=np.float64) )
    cdef double** ptrTDI = <double**>malloc( n_threads * sizeof(double*) )
    for i in range(n_threads):
        ptrTDI[i] = &niiTDI_img[i,0,0,0]

    cdef double [:, ::1] peaksAffine
    cdef double* ptrPeaksAffine

    if filename_peaks is not None :
        logger.subinfo( 'EC orientations:', indent_lvl=2, indent_char='-' )
        niiPEAKS = nibabel.load( filename_peaks )
        niiPEAKS_hdr = _get_header( niiPEAKS )
        logger.subinfo( f'{niiPEAKS.shape[0]} x {niiPEAKS.shape[1]} x {niiPEAKS.shape[2]} x {niiPEAKS.shape[3]}', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'{niiPEAKS_hdr["pixdim"][1]:.4f} x {niiPEAKS_hdr["pixdim"][2]:.4f} x {niiPEAKS_hdr["pixdim"][3]:.4f}', indent_lvl=3, indent_char='-' )
        
        logger.subinfo( f'ignoring peaks < {vf_THR:.2f} * MaxPeak', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'{"" if peaks_use_affine else "not "}using affine matrix', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'flipping axes : [ x={flip_peaks[0]}, y={flip_peaks[1]}, z={flip_peaks[2]} ]', indent_lvl=3, indent_char='-' )
        if ( Nx!=niiPEAKS.shape[0] or Ny!=niiPEAKS.shape[1] or Nz!=niiPEAKS.shape[2] or
            abs(Px-niiPEAKS_hdr['pixdim'][1])>1e-3 or abs(Py-niiPEAKS_hdr['pixdim'][2])>1e-3 or abs(Pz-niiPEAKS_hdr['pixdim'][3])>1e-3 ) :
            logger.warning( "Dataset does not have the same geometry as the tractogram" )
        if niiPEAKS.shape[3] % 3 :
            logger.error( 'PEAKS dataset must have 3*k volumes' )
        if vf_THR < 0 or vf_THR > 1 :
            logger.error( '"vf_THR" must be between 0 and 1' )
        niiPEAKS_img = np.ascontiguousarray( np.asanyarray( niiPEAKS.dataobj ).astype(np.float32) )
        ptrPEAKS = &niiPEAKS_img[0,0,0,0]
        Np = niiPEAKS.shape[3]/3

        # affine matrix to rotate gradien directions (if required)
        if peaks_use_affine :
            peaksAffine = np.ascontiguousarray( niiPEAKS.affine[:3,:3].T )
        else :
            peaksAffine = np.ascontiguousarray( np.eye(3) )
        ptrPeaksAffine = &peaksAffine[0,0]
    else :
        logger.subinfo( 'No dataset specified for EC compartments', indent_lvl=2, indent_char='-' )
        Np = 0
        ptrPEAKS = NULL
        ptrPeaksAffine = NULL
    
    # ISO map for isotropic compartment
    cdef float* ptrISO
    cdef float [:, :, ::1] niiISO_img
    if filename_ISO is not None :
        logger.subinfo( 'Restricted ISO map', indent_lvl=2, indent_char='-' )
        niiISO = nibabel.load( filename_ISO )
        niiISO_hdr = _get_header( niiISO )
        logger.subinfo( f'{niiISO.shape[0]} x {niiISO.shape[1]} x {niiISO.shape[2]}', indent_lvl=3, indent_char='-' )
        logger.subinfo( f'{niiISO_hdr["pixdim"][1]:.4f} x {niiISO_hdr["pixdim"][2]:.4f} x {niiISO_hdr["pixdim"][3]:.4f}', indent_lvl=3, indent_char='-' )
        if ( Nx!=niiISO.shape[0] or Ny!=niiISO.shape[1] or Nz!=niiISO.shape[2] or
            abs(Px-niiISO_hdr['pixdim'][1])>1e-3 or abs(Py-niiISO_hdr['pixdim'][2])>1e-3 or abs(Pz-niiISO_hdr['pixdim'][3])>1e-3 ) :
            logger.warning( 'Dataset does not have the same geometry as the tractogram' )
        niiISO_img = np.ascontiguousarray( np.asanyarray( niiISO.dataobj ).astype(np.float32) )
        ptrISO  = &niiISO_img[0,0,0]
    else :
        logger.subinfo( 'No ISO map specified, using tdi', indent_lvl=2, indent_char='-' )
        ptrISO = NULL

    # write dictionary information info file
    dictionary_info = {}
    dictionary_info['filename_tractogram'] = filename_tractogram
    if blur_clust_thr[0]> 0:
        idx_centroids = np.array(idx_centroids, dtype=np.uint32)
        dictionary_info['n_count'] = input_n_count
        dictionary_info['tractogram_centr_idx'] = idx_centroids
    dictionary_info['TCK_ref_image'] = TCK_ref_image
    dictionary_info['path_out'] = path_out
    dictionary_info['filename_peaks'] = filename_peaks
    dictionary_info['filename_mask'] = filename_mask
    dictionary_info['do_intersect'] = do_intersect
    dictionary_info['fiber_shift'] = fiber_shift
    dictionary_info['min_seg_len'] = min_seg_len
    dictionary_info['min_fiber_len'] = min_fiber_len
    dictionary_info['max_fiber_len'] = max_fiber_len
    dictionary_info['vf_THR'] = vf_THR
    dictionary_info['peaks_use_affine'] = peaks_use_affine
    dictionary_info['flip_peaks'] = flip_peaks
    dictionary_info['blur_core_extent'] = blur_core_extent
    dictionary_info['blur_gauss_extent'] = blur_gauss_extent
    dictionary_info['blur_gauss_min'] = blur_gauss_min
    dictionary_info['blur_spacing'] = blur_spacing
    dictionary_info['blur_sigma'] = blur_sigma
    dictionary_info['blur_apply_to'] = blur_apply_to
    dictionary_info['ndirs'] = ndirs
    dictionary_info['n_threads'] = n_threads
    with open( join(path_out,'dictionary_info.pickle'), 'wb+' ) as dictionary_info_file:
        pickle.dump(dictionary_info, dictionary_info_file, protocol=2)

    # calling actual C code
    ret = trk2dictionary( filename_tractogram, data_offset,
        Nx, Ny, Nz, Px, Py, Pz, n_count, n_scalars, n_properties,
        fiber_shiftX, fiber_shiftY, fiber_shiftZ, min_seg_len, min_fiber_len, max_fiber_len,
        ptrPEAKS, Np, vf_THR, -1 if flip_peaks[0] else 1, -1 if flip_peaks[1] else 1, -1 if flip_peaks[2] else 1,
        ptrMASK, ptrISO, ptrTDI, path_temp, 1 if do_intersect else 0, ptrPeaksAffine,
        nReplicas, &blurRho[0], &blurAngle[0], &blurWeights[0], &blurApplyTo[0], ptrToVOXMM, ptrHashTable, n_threads, verbose if not _in_notebook() else 0 );
    if ret == 0 :
        logger.warning( 'DICTIONARY not generated' )
        return None

    # NOTE: this is to ensure flushing all the output from the cpp code
    print(end='', flush=True)
    time.sleep(0.5)

    # Concatenate files together
    log_list = []
    ret_subinfo = logger.subinfo( 'Saving data structure', indent_char='*', indent_lvl=1, with_progress=verbose>2 )
    cdef int discarded = 0
    with ProgressBar(disable=verbose<3, hide_on_exit=True, subinfo=ret_subinfo, log_list=log_list):
        for j in range(n_threads-1):
            path_IC_f = join(path_temp, f'dictionary_IC_f_{j+1}.dict')
            kept = np.fromfile( join(path_temp, f'dictionary_TRK_kept_{j}.dict'), dtype=np.uint8 )
            IC_f = np.fromfile( join(path_temp, f'dictionary_IC_f_{j+1}.dict'), dtype=np.uint32 )
            discarded += np.count_nonzero(kept==0)
            IC_f -= discarded
            IC_f_save = np.memmap( path_IC_f, dtype="uint32", mode='w+', shape=IC_f.shape )
            IC_f_save[:] = IC_f[:]
            IC_f_save.flush()
            del IC_f_save
            # np.save( path_out + f'/dictionary_IC_f_{j+1}.dict', IC_f, allow_pickle=True)

        fileout = join(path_out, 'dictionary_TRK_kept.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_TRK_kept_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_TRK_norm.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_TRK_norm_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_TRK_len.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_TRK_len_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_TRK_lenTot.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_TRK_lenTot_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_IC_f.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_IC_f_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_IC_v.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_IC_v_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_IC_o.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_IC_o_{j}.dict') ]
        cat_function( dict_list, fileout )

        fileout = join(path_out, 'dictionary_IC_len.dict')
        dict_list = []
        for j in range(n_threads):
            dict_list += [ join(path_temp, f'dictionary_IC_len_{j}.dict') ]
        cat_function( dict_list, fileout )


        # save TDI and MASK maps
        if TCK_ref_image is not None:
            TDI_affine = _get_affine( niiREF )
        elif filename_mask is not None :
            TDI_affine = _get_affine( niiMASK )
        elif filename_peaks is not None :
            TDI_affine = _get_affine( niiPEAKS )
        else :
            TDI_affine = np.diag( [Px, Py, Pz, 1] )

    # save TDI and MASK maps
    cdef np.float32_t [::1] niiTDI_mem = np.zeros( Nx*Ny*Nz, dtype=np.float32 )
    log_list = []
    ret_subinfo = logger.subinfo( 'Saving TDI and MASK maps', indent_char='*', indent_lvl=1, with_progress=verbose>2)
    with ProgressBar(disable=verbose<3, hide_on_exit=True, subinfo=ret_subinfo, log_list=log_list):
        v = np.fromfile( join(path_out, 'dictionary_IC_v.dict'),   dtype=np.uint32 )
        l = np.fromfile( join(path_out, 'dictionary_IC_len.dict'), dtype=np.float32 )
        
        niiTDI_mem = compute_tdi( v, l, Nx, Ny, Nz, verbose )
        niiTDI_img_save = np.reshape( niiTDI_mem, (Nx,Ny,Nz), order='F' )

        niiTDI = nibabel.Nifti1Image( niiTDI_img_save, TDI_affine )
        niiTDI_hdr = _get_header( niiTDI )
        niiTDI_hdr['descrip'] = f'Created with COMMIT {get_distribution("dmri-commit").version}'
        nibabel.save( niiTDI, join(path_out,'dictionary_tdi.nii.gz') )

        if filename_mask is not None :
            niiMASK = nibabel.Nifti1Image( niiMASK_img, TDI_affine )
        else :
            niiMASK = nibabel.Nifti1Image( (np.asarray(niiTDI_img_save)>0).astype(np.float32), TDI_affine )

        niiMASK_hdr = _get_header( niiMASK )
        niiMASK_hdr['descrip'] = niiTDI_hdr['descrip']
        nibabel.save( niiMASK, join(path_out,'dictionary_mask.nii.gz') )

        if not keep_temp:
            shutil.rmtree(path_temp)

    logger.info( f'[ {format_time(time.time() - tic)} ]' )