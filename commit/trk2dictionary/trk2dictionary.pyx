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


# Interface to actual C code
cdef extern from "trk2dictionary_c.cpp":
    int trk2dictionary(
        char* filename_tractogram, int data_offset, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars, 
        int n_properties, float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len,  float max_fiber_len,
        float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
        float* _ptrMASK, float* ptrTDI, char* path_out, int c, double* ptrPeaksAffine,
        int nBlurRadii, double blurSigma, double* ptrBlurRadii, int* ptrBlurSamples, double* ptrBlurWeights,  float* ptrTractsAffine, unsigned short ndirs, short* prtHashTable
    ) nogil


cpdef run( filename_tractogram=None, path_out=None, filename_peaks=None, filename_mask=None, do_intersect=True,
    fiber_shift=0, min_seg_len=1e-3, min_fiber_len=0.0, max_fiber_len=250.0,
    vf_THR=0.1, peaks_use_affine=False, flip_peaks=[False,False,False], 
    blur_radii=[], blur_samples=[], blur_sigma=0.0,
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
        If True then fiber segments that intersect voxel boundaries are splitted (default).
        If False then the centroid of the segment is used as its voxel position.

    fiber_shift : float or list of three float
        If necessary, apply a translation to fiber coordinates (default : 0) to account
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

    blur_radii : list of float
        Translate each segment to given radii to assign a broader fiber contribution (default : []).
    
    blur_samples : list of integer
        Segments are duplicated along a circle at a given radius; this parameter controls the
        number of samples to take over a given circle (defaut : []).

    blur_sigma: float
        The contributions of the segments at different radii are damped as a Gaussian (default : 0.0).
    
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
    if type(blur_radii)==list:
        blur_radii = np.array(blur_radii, np.double)
    elif type(blur_radii)!=np.ndarray:
        ERROR( '"blur_radii" must be a list of floats' )
    if type(blur_samples)==list:
        blur_samples = np.array(blur_samples, np.int32)
    elif type(blur_samples)!=np.ndarray:
        ERROR( '"blur_samples" must be a list of integers' )

    if blur_sigma > 0 :
        if blur_radii.size != blur_samples.size :
            ERROR( 'The number of blur radii and blur samples must match' )

        if np.count_nonzero( blur_radii<=0 ):
            ERROR( 'A blur radius was <= 0; only positive radii can be used' )

        if np.count_nonzero( blur_samples<1 ):
            ERROR( 'Please specify at least 1 sample per blur radius' )

    tic = time.time()
    LOG( '\n-> Creating the dictionary from tractogram:' )
    
    LOG( '\n   * Configuration:' )
    print( '\t- Segment position = %s' % ( 'COMPUTE INTERSECTIONS' if do_intersect else 'CENTROID' ) )
    print( '\t- Fiber shift X    = %.3f (voxel-size units)' % fiber_shiftX )
    print( '\t- Fiber shift Y    = %.3f (voxel-size units)' % fiber_shiftY )
    print( '\t- Fiber shift Z    = %.3f (voxel-size units)' % fiber_shiftZ )
    if min_seg_len >= 1e-3:
        print( '\t- Min segment len  = %.3f mm' % min_seg_len )
    else:
        print( '\t- Min segment len  = %.2e mm' % min_seg_len )
    print( '\t- Min fiber len    = %.2f mm' % min_fiber_len )
    print( '\t- Max fiber len    = %.2f mm' % max_fiber_len )

    # check blur params
    cdef :
        double [:] blurRadii
        int [:] blurSamples
        double [:] blurWeights
        double* ptrBlurRadii
        int* ptrBlurSamples
        double* ptrBlurWeights
        int nBlurRadii
        float [:] ArrayInvM
        float* ptrArrayInvM
    
    # add a fake radius for original segment
    if blur_sigma == 0:
        nBlurRadii = 1
        blurRadii = np.array( [0.0], np.double )
        blurSamples = np.array( [1], np.int32 )
        blurWeights = np.array( [1], np.double )
    else:
        nBlurRadii = len(blur_radii)+1
        blurRadii = np.insert( blur_radii, 0, 0.0 ).astype(np.double)
        blurSamples = np.insert( blur_samples, 0, 1 ).astype(np.int32)

        # compute weights for gaussian damping
        blurWeights = np.empty_like( blurRadii )
        for i in xrange(nBlurRadii):
            blurWeights[i] = np.exp( -blurRadii[i]**2 / (2.0*blur_sigma**2) )

    if nBlurRadii == 1 :
        print( '\t- Do not blur fibers' )
    else :
        print( '\t- Blur fibers:' )
        print( '\t\t- sigma = %.3f' % blur_sigma )
        print( '\t\t- radii =   [ ', end="" )
        for i in xrange( 1, blurRadii.size ) :
            print( '%.3f ' % blurRadii[i], end="" )
        print( ']' )
        print( '\t\t- weights = [ ', end="" )
        for i in xrange( 1, blurWeights.size ) :
            print( '%.3f ' % blurWeights[i], end="" )
        print( ']' )
        print( '\t\t- samples = [ ', end="" )
        for i in xrange( 1, blurSamples.size ) :
            print( '%5d ' % blurSamples[i], end="" )
        print( ']' )

    ptrBlurRadii   = &blurRadii[0]
    ptrBlurSamples = &blurSamples[0]
    ptrBlurWeights = &blurWeights[0]

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
    print( '\t- Output written to "%s"' % path_out )
    if not exists( path_out ):
        makedirs( path_out )

    # Load data from files
    LOG( '\n   * Loading data:' )
    cdef short [:] htable = amico.lut.load_precomputed_hash_table(ndirs)
    cdef short* ptrHashTable = &htable[0]

    # Streamlines from tractogram
    print( '\t- Tractogram' )
    
    if not exists(filename_tractogram):
        ERROR( 'Tractogram file not found: %s' % filename_tractogram )        
    extension = splitext(filename_tractogram)[1]
    if extension != ".trk" and extension != ".tck":
        ERROR( 'Invalid input file: only .trk and .tck are supported' )
    
    hdr = nibabel.streamlines.load( filename_tractogram, lazy_load=True ).header
            
    if extension == ".trk":
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

    if extension == ".tck":
        if TCK_ref_image is None:
            if filename_peaks is not None:
                TCK_ref_image = filename_peaks
            elif filename_mask is not None:
                TCK_ref_image = filename_mask
            else:
                ERROR( 'TCK files do not contain information about the geometry. Use "TCK_ref_image" for that' )

        print ('\t\t- geometry taken from "%s"' %TCK_ref_image)

        nii_image = nibabel.load(TCK_ref_image)
        nii_hdr = nii_image.header if nibabel.__version__ >= '2.0.0' else nii_image.get_header()
        Nx = nii_image.shape[0]
        Ny = nii_image.shape[1]
        Nz = nii_image.shape[2]
        Px = nii_hdr['pixdim'][1]
        Py = nii_hdr['pixdim'][2]
        Pz = nii_hdr['pixdim'][3]
        data_offset = int(hdr['_offset_data'])  #set offset
        n_count = int(hdr['count'])  #set number of fibers
        n_scalars = 0
        n_properties = 0
        
    print( '\t\t- %d x %d x %d' % ( Nx, Ny, Nz ) )
    print( '\t\t- %.4f x %.4f x %.4f' % ( Px, Py, Pz ) )
    print( '\t\t- %d fibers' % n_count )
    if Nx >= 2**16 or Nz >= 2**16 or Nz >= 2**16 :
        ERROR( 'The max dim size is 2^16 voxels' )
    
    # get the affine matrix
    if extension == ".tck":
        scaleMat = np.diag(np.divide(1.0, [Px,Py,Pz]))
        M = nii_hdr.get_best_affine()

        # Affine matrix without scaling, i.e. diagonal is 1
        M[:3, :3] = np.dot(scaleMat, M[:3, :3])
        M = M.astype('<f4') # affine matrix in float value
        invM = np.linalg.inv(M) # inverse affine matrix
        #create a vector of inverse matrix M
        ArrayInvM = np.ravel(invM)
        ptrArrayInvM = &ArrayInvM[0]

    # white-matter mask
    cdef float* ptrMASK
    cdef float [:, :, ::1] niiMASK_img
    if filename_mask is not None :
        print( '\t- Filtering mask' )
        niiMASK = nibabel.load( filename_mask )
        niiMASK_hdr = niiMASK.header if nibabel.__version__ >= '2.0.0' else niiMASK.get_header()
        print( '\t\t- %d x %d x %d' % ( niiMASK.shape[0], niiMASK.shape[1], niiMASK.shape[2] ) )
        print( '\t\t- %.4f x %.4f x %.4f' % ( niiMASK_hdr['pixdim'][1], niiMASK_hdr['pixdim'][2], niiMASK_hdr['pixdim'][3] ) )
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
    cdef double [:, ::1] affine
    cdef double* ptrAFFINE
    if filename_peaks is not None :
        print( '\t- EC orientations' )
        niiPEAKS = nibabel.load( filename_peaks )
        niiPEAKS_hdr = niiPEAKS.header if nibabel.__version__ >= '2.0.0' else niiPEAKS.get_header()
        print( '\t\t- %d x %d x %d x %d' % ( niiPEAKS.shape[0], niiPEAKS.shape[1], niiPEAKS.shape[2], niiPEAKS.shape[3] ) )
        print( '\t\t- %.4f x %.4f x %.4f' % ( niiPEAKS_hdr['pixdim'][1], niiPEAKS_hdr['pixdim'][2], niiPEAKS_hdr['pixdim'][3] ) )
        print( '\t\t- ignoring peaks < %.2f * MaxPeak' % vf_THR )
        print( '\t\t- %susing affine matrix' % ( "" if peaks_use_affine else "not " ) )
        print( '\t\t- flipping axes : [ x=%s, y=%s, z=%s ]' % ( flip_peaks[0], flip_peaks[1], flip_peaks[2] ) )
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
            affine = np.ascontiguousarray( niiPEAKS.affine[:3,:3].T )
        else :
            affine = np.ascontiguousarray( np.eye(3) )
        ptrAFFINE = &affine[0,0]
    else :
        print( '\t- No dataset specified for EC compartments' )
        Np = 0
        ptrPEAKS = NULL
        ptrAFFINE = NULL

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
    dictionary_info['blur_radii'] = blur_radii
    dictionary_info['blur_samples'] = blur_samples
    dictionary_info['blur_sigma'] = blur_sigma    
    dictionary_info['ndirs'] = ndirs
    with open( join(path_out,'dictionary_info.pickle'), 'wb+' ) as dictionary_info_file:
        pickle.dump(dictionary_info, dictionary_info_file, protocol=2)

    # calling actual C code
    ret = trk2dictionary( filename_tractogram, data_offset,
        Nx, Ny, Nz, Px, Py, Pz, n_count, n_scalars, n_properties,
        fiber_shiftX, fiber_shiftY, fiber_shiftZ, min_seg_len, min_fiber_len, max_fiber_len,
        ptrPEAKS, Np, vf_THR, -1 if flip_peaks[0] else 1, -1 if flip_peaks[1] else 1, -1 if flip_peaks[2] else 1,
        ptrMASK, ptrTDI, path_out, 1 if do_intersect else 0, ptrAFFINE,
        nBlurRadii, blur_sigma, ptrBlurRadii, ptrBlurSamples, ptrBlurWeights, ptrArrayInvM, ndirs, ptrHashTable  );
    if ret == 0 :
        WARNING( 'DICTIONARY not generated' )
        return None

    # save TDI and MASK maps
    if filename_mask is not None :
        affine = niiMASK.affine if nibabel.__version__ >= '2.0.0' else niiMASK.get_affine()
    elif filename_peaks is not None :
        affine = niiPEAKS.affine if nibabel.__version__ >= '2.0.0' else niiPEAKS.get_affine()
    else :
        affine = np.diag( [Px, Py, Pz, 1] )

    niiTDI = nibabel.Nifti1Image( niiTDI_img, affine )
    nii_hdr = niiTDI.header if nibabel.__version__ >= '2.0.0' else niiTDI.get_header()
    nii_hdr['descrip'] = 'Created with COMMIT %s'%get_distribution('dmri-commit').version
    nibabel.save( niiTDI, join(path_out,'dictionary_tdi.nii.gz') )

    if filename_mask is not None :
        niiMASK = nibabel.Nifti1Image( niiMASK_img, affine )
    else :
        niiMASK = nibabel.Nifti1Image( (np.asarray(niiTDI_img)>0).astype(np.float32), affine )
    nii_hdr = niiMASK.header if nibabel.__version__ >= '2.0.0' else niiMASK.get_header()
    nii_hdr['descrip'] = 'Created with COMMIT %s'%get_distribution('dmri-commit').version
    nibabel.save( niiMASK, join(path_out,'dictionary_mask.nii.gz') )

    LOG( '\n   [ %.1f seconds ]' % ( time.time() - tic ) )