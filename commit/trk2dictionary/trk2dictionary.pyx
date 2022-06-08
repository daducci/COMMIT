#!python
# cython: language_level=3, c_string_type=str, c_string_encoding=ascii, boundscheck=False, wraparound=False, profile=False
from __future__ import print_function
import cython
import numpy as np
cimport numpy as np
import nibabel
from os.path import join, exists, splitext, dirname, isdir
from os import makedirs, remove
import time
import amico
import pickle
from amico.util import LOG, NOTE, WARNING, ERROR
from pkg_resources import get_distribution

from libcpp cimport bool

# Interface to actual C code
cdef extern from "trk2dictionary_c.cpp":
    int trk2dictionary(
        char* filename_tractogram, int data_offset, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars,
        int n_properties, float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len,  float max_fiber_len,
        float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
        float* _ptrMASK, float* ptrTDI, char* path_out, int c, double* ptrPeaksAffine,
        int nReplicas, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, bool* ptrBlurApplyTo,
        float* ptrTractsAffine, unsigned short ndirs, short* prtHashTable
    ) nogil

def _get_header( niiFILE ):
    return niiFILE.header if nibabel.__version__ >= '2.0.0' else niiFILE.get_header()

def _get_affine( niiFILE ):
    return niiFILE.affine if nibabel.__version__ >= '2.0.0' else niiFILE.get_affine()

cpdef run( filename_tractogram=None, path_out=None, filename_peaks=None, filename_mask=None, do_intersect=True,
    fiber_shift=0, min_seg_len=1e-3, min_fiber_len=0.0, max_fiber_len=250.0,
    vf_THR=0.1, peaks_use_affine=False, flip_peaks=[False,False,False],
    blur_spacing=0.25, blur_core_extent=0.0, blur_gauss_extent=0.0, blur_gauss_min=0.1, blur_apply_to=None,
    filename_trk=None, gen_trk=None, TCK_ref_image=None, ndirs=32761
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

    ndirs : int
        Number of orientations on the sphere used to discretize the orientation of each
        each segment in a streamline (default : 32761).

    filename_trk : string
        DEPRECATED. Use filename_tractogram instead.

    gen_trk : string
        DEPRECATED. No tractogram will be saved any more, but the returned coefficients will account
        for the streamlines that were pre-filtered in this function.
    """

    # check the value of ndirs
    if not amico.lut.is_valid(ndirs):
        ERROR( 'Unsupported value for ndirs.\nNote: Supported values for ndirs are [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 32761 (default)]' )

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
        ERROR( '"fiber_shift" must be a scalar or a vector with 3 elements' )

    # check for invalid parameters in the blur
    if blur_core_extent < 0 :
        ERROR( 'The extent of the core must be non-negative' )

    if blur_gauss_extent < 0 :
        ERROR( 'The extent of the blur must be non-negative' )

    if blur_gauss_extent > 0 or blur_core_extent > 0:
        if blur_spacing <= 0 :
            ERROR( 'The grid spacing of the blur must be positive' )

    tic = time.time()
    LOG( '\n-> Creating the dictionary from tractogram:' )

    LOG( '\n   * Configuration:' )
    print( f'\t- Segment position = {"COMPUTE INTERSECTIONS" if do_intersect else "CENTROID"}' )
    print( f'\t- Coordinates shift in X = {fiber_shiftX:.3f} (voxel-size units)' )
    print( f'\t- Coordinates shift in Y = {fiber_shiftY:.3f} (voxel-size units)' )
    print( f'\t- Coordinates shift in Z = {fiber_shiftZ:.3f} (voxel-size units)' )
    if min_seg_len >= 1e-3:
        print( f'\t- Min segment len  = {min_seg_len:.3f} mm' )
    else:
        print( f'\t- Min segment len  = {min_seg_len:.2e} mm' )
    print( f'\t- Min streamline len    = {min_fiber_len:.2f} mm' )
    print( f'\t- Max streamline len    = {max_fiber_len:.2f} mm' )

    # check blur params
    cdef :
        double [:] blurRho
        double [:] blurAngle
        double [:] blurWeights
        bool [:] blurApplyTo
        int nReplicas
        float blur_sigma

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
        print( '\t- Do not blur streamlines' )
    else :
        print( '\t- Blur streamlines:' )
        print( f'\t\t- core extent  = {blur_core_extent:.3f}' )
        print( f'\t\t- gauss extent = {blur_gauss_extent:.3f} (sigma = {blur_sigma:.3f})' )
        print( f'\t\t- grid spacing = {blur_spacing:.3f}' )
        print( f'\t\t- weights = [ {np.min(blurWeights):.3f} ... {np.max(blurWeights):.3f} ]' )

    if min_seg_len < 0 :
        ERROR( '"min_seg_len" must be >= 0' )
    if min_fiber_len < 0 :
        ERROR( '"min_fiber_len" must be >= 0' )
    if max_fiber_len < min_fiber_len :
        ERROR( '"max_fiber_len" must be >= "min_fiber_len"' )

    if filename_trk is None and filename_tractogram is None:
        ERROR( '"filename_tractogram" not defined' )

    if filename_trk is not None and filename_tractogram is not None:
        WARNING('"filename_trk" will not be considered, "filename_tractogram" will be used')

    if filename_trk is not None and filename_tractogram is None:
        filename_tractogram = filename_trk
        WARNING('"filename_trk" parameter is deprecated, use "filename_tractogram" instead')

    if path_out is None:
        path_out = dirname(filename_tractogram)
        if path_out == '':
            path_out = '.'
        if not isdir(path_out):
            ERROR( '"path_out" cannot be inferred from "filename_tractogram"' )
        path_out = join(path_out,'COMMIT')

    if gen_trk is not None:
        WARNING('"gen_trk" parameter is deprecated')

    # create output path
    print( f'\t- Output written to "{path_out}"' )
    if not exists( path_out ):
        makedirs( path_out )

    # Load data from files
    LOG( '\n   * Loading data:' )
    cdef short [:] htable = amico.lut.load_precomputed_hash_table(ndirs)
    cdef short* ptrHashTable = &htable[0]

    # Streamlines from tractogram
    print( '\t- Tractogram' )

    if not exists(filename_tractogram):
        ERROR( f'Tractogram file not found: {filename_tractogram}' )
    extension = splitext(filename_tractogram)[1]
    if extension != ".trk" and extension != ".tck":
        ERROR( 'Invalid input file: only .trk and .tck are supported' )

    hdr = nibabel.streamlines.load( filename_tractogram, lazy_load=True ).header

    if extension == ".trk":
        print ( f'\t\t- geometry taken from "{filename_tractogram}"' )
        Nx = hdr['dimensions'][0]
        Ny = hdr['dimensions'][1]
        Nz = hdr['dimensions'][2]
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
                ERROR( 'TCK files do not contain information about the geometry. Use "TCK_ref_image" for that' )
        print ( f'\t\t- geometry taken from "{TCK_ref_image}"' )

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

    print( f'\t\t- {Nx} x {Ny} x {Nz}' )
    print( f'\t\t- {Px:.4f} x {Py:.4f} x {Pz:.4f}' )
    print( f'\t\t- {n_count} streamlines' )
    if Nx >= 2**16 or Nz >= 2**16 or Nz >= 2**16 :
        ERROR( 'The max dim size is 2^16 voxels' )

    # check copmpatibility between blurApplyTo and number of streamlines
    if blur_apply_to is None:
        blur_apply_to = np.repeat([True], n_count)
    else :
        if blur_apply_to.size != n_count :
            ERROR( '"blur_apply_to" must have one value per streamline' )
        print( f'\t\t\t- {sum(blur_apply_to)} blurred streamlines' )
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
        print( '\t- Filtering mask' )
        niiMASK = nibabel.load( filename_mask )
        niiMASK_hdr = _get_header( niiMASK )
        print( f'\t\t- {niiMASK.shape[0]} x {niiMASK.shape[1]} x {niiMASK.shape[2]}' )
        print( f'\t\t- {niiMASK_hdr["pixdim"][1]:.4f} x {niiMASK_hdr["pixdim"][2]:.4f} x {niiMASK_hdr["pixdim"][3]:.4f}' )
        if ( Nx!=niiMASK.shape[0] or Ny!=niiMASK.shape[1] or Nz!=niiMASK.shape[2] or
            abs(Px-niiMASK_hdr['pixdim'][1])>1e-3 or abs(Py-niiMASK_hdr['pixdim'][2])>1e-3 or abs(Pz-niiMASK_hdr['pixdim'][3])>1e-3 ) :
            WARNING( 'Dataset does not have the same geometry as the tractogram' )
        niiMASK_img = np.ascontiguousarray( np.asanyarray( niiMASK.dataobj ).astype(np.float32) )
        ptrMASK  = &niiMASK_img[0,0,0]
    else :
        print( '\t- No mask specified to filter IC compartments' )
        ptrMASK = NULL

    # peaks file for EC contributions
    cdef float* ptrPEAKS
    cdef float [:, :, :, ::1] niiPEAKS_img
    cdef int Np
    cdef float [:, :, ::1] niiTDI_img = np.ascontiguousarray( np.zeros((Nx,Ny,Nz),dtype=np.float32) )
    cdef float* ptrTDI  = &niiTDI_img[0,0,0]
    cdef double [:, ::1] peaksAffine
    cdef double* ptrPeaksAffine
    if filename_peaks is not None :
        print( '\t- EC orientations' )
        niiPEAKS = nibabel.load( filename_peaks )
        niiPEAKS_hdr = _get_header( niiPEAKS )
        print( f'\t\t- {niiPEAKS.shape[0]} x {niiPEAKS.shape[1]} x {niiPEAKS.shape[2]} x {niiPEAKS.shape[3]}' )
        print( f'\t\t- {niiPEAKS_hdr["pixdim"][1]:.4f} x {niiPEAKS_hdr["pixdim"][2]:.4f} x {niiPEAKS_hdr["pixdim"][3]:.4f}' )

        print( f'\t\t- ignoring peaks < {vf_THR:.2f} * MaxPeak' )
        print( f'\t\t- {"" if peaks_use_affine else "not "}using affine matrix' )
        print( f'\t\t- flipping axes : [ x={flip_peaks[0]}, y={flip_peaks[1]}, z={flip_peaks[2]} ]' )
        if ( Nx!=niiPEAKS.shape[0] or Ny!=niiPEAKS.shape[1] or Nz!=niiPEAKS.shape[2] or
            abs(Px-niiPEAKS_hdr['pixdim'][1])>1e-3 or abs(Py-niiPEAKS_hdr['pixdim'][2])>1e-3 or abs(Pz-niiPEAKS_hdr['pixdim'][3])>1e-3 ) :
            WARNING( "Dataset does not have the same geometry as the tractogram" )
        if niiPEAKS.shape[3] % 3 :
            ERROR( 'PEAKS dataset must have 3*k volumes' )
        if vf_THR < 0 or vf_THR > 1 :
            ERROR( '"vf_THR" must be between 0 and 1' )
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
        print( '\t- No dataset specified for EC compartments' )
        Np = 0
        ptrPEAKS = NULL
        ptrPeaksAffine = NULL

    # write dictionary information info file
    dictionary_info = {}
    dictionary_info['filename_tractogram'] = filename_tractogram
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
    with open( join(path_out,'dictionary_info.pickle'), 'wb+' ) as dictionary_info_file:
        pickle.dump(dictionary_info, dictionary_info_file, protocol=2)

    # calling actual C code
    ret = trk2dictionary( filename_tractogram, data_offset,
        Nx, Ny, Nz, Px, Py, Pz, n_count, n_scalars, n_properties,
        fiber_shiftX, fiber_shiftY, fiber_shiftZ, min_seg_len, min_fiber_len, max_fiber_len,
        ptrPEAKS, Np, vf_THR, -1 if flip_peaks[0] else 1, -1 if flip_peaks[1] else 1, -1 if flip_peaks[2] else 1,
        ptrMASK, ptrTDI, path_out, 1 if do_intersect else 0, ptrPeaksAffine,
        nReplicas, &blurRho[0], &blurAngle[0], &blurWeights[0], &blurApplyTo[0], ptrToVOXMM, ndirs, ptrHashTable  );
    if ret == 0 :
        WARNING( 'DICTIONARY not generated' )
        return None

    # save TDI and MASK maps
    if TCK_ref_image is not None:
        TDI_affine = _get_affine( niiREF )
    elif filename_mask is not None :
        TDI_affine = _get_affine( niiMASK )
    elif filename_peaks is not None :
        TDI_affine = _get_affine( niiPEAKS )
    else :
        TDI_affine = np.diag( [Px, Py, Pz, 1] )

    niiTDI = nibabel.Nifti1Image( niiTDI_img, TDI_affine )
    niiTDI_hdr = _get_header( niiTDI )
    niiTDI_hdr['descrip'] = f'Created with COMMIT {get_distribution("dmri-commit").version}'
    nibabel.save( niiTDI, join(path_out,'dictionary_tdi.nii.gz') )

    if filename_mask is not None :
        niiMASK = nibabel.Nifti1Image( niiMASK_img, TDI_affine )
    else :
        niiMASK = nibabel.Nifti1Image( (np.asarray(niiTDI_img)>0).astype(np.float32), TDI_affine )
    niiMASK_hdr = _get_header( niiMASK )
    niiMASK_hdr['descrip'] = niiTDI_hdr['descrip']
    nibabel.save( niiMASK, join(path_out,'dictionary_mask.nii.gz') )

    LOG( f'\n   [ {time.time() - tic:.1f} seconds ]' )