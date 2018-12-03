a#!python
# cython: boundscheck=False, wraparound=False, profile=False

import cython
import numpy as np
cimport numpy as np
import nibabel
from os.path import join, exists
from os import makedirs, remove
import time

# Interface to actual C code
cdef extern from "trk2dictionary_c.cpp":
    int trk2dictionary(
        char* filename_tractogram, int data_offset, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars, int n_properties, float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, int points_to_skip, float min_seg_len,
        float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
        float* _ptrMASK, float* ptrTDI, char* path_out, int c, double* ptrAFFINE,
        int nBlurRadii, double blurSigma, double* ptrBlurRadii, int* ptrBlurSamples, double* ptrBlurWeights
    ) nogil


cpdef run( filename_tractogram, path_out, TCK_ref_image = None, filename_peaks = None, filename_mask = None, do_intersect = True,
    fiber_shift = 0, points_to_skip = 0, vf_THR = 0.1, peaks_use_affine = False,
    flip_peaks = [False,False,False], min_seg_len = 1e-3, gen_trk = True,
    blur_radii = [], blur_samples = [], blur_sigma = 1.0
    ):
    """Perform the conversion of a tractoram to the sparse data-structure internally
    used by COMMIT to perform the matrix-vector multiplications with the operator A
    during the inversion of the linear system.

    Parameters
    ----------
    filename_trk : string
        Path to the .trk or .tck file containing the tractogram to convert.

    path_out : string
        Path to the folder where to store the sparse data structure.

    filename_peaks : string
        Path to the NIFTI file containing the peaks to use as extra-cellular contributions.
        The data matrix should be 4D with last dimension 3*N, where N is the number
        of peaks in each voxel. (default : no extra-cellular contributions)

    filename_mask : string
        Path to a binary mask to restrict the analysis to specific areas. Segments
        outside this mask are discarded. If not specified (default), the mask is created from
        all voxels intersected by the tracts.

    do_intersect : boolean
        If True then fiber segments that intersect voxel boundaries are splitted (default).
        If False then the centroid of the segment is used as its voxel position.

    fiber_shift : float or list of three float
        If necessary, apply a translation to fiber coordinates (default : 0) to account
        for differences between the reference system of the tracking algorithm and COMMIT.
        The value is specified in voxel units, eg 0.5 translates by half voxel.
        Do noth use if you are using fiber_shiftX or fiber_shiftY or fiber_shiftZ.

    points_to_skip : integer
        If necessary, discard first points at beginning/end of a fiber (default : 0).

    vf_THR : float
        Discard peaks smaller than vf_THR * max peak (default : 0.1).

    peaks_use_affine : boolean
        Whether to rotate the peaks according to the affine matrix (default : False).

    flip_peaks : list of three boolean
        If necessary, flips peak orientations along each axis (default : no flipping).

    min_seg_len : float
        Discard segments <= than this length in mm (default : 1e-3)

    gen_trk : boolean
        If True then generate a .trk file in the 'path_out' containing the fibers used in the dictionary (default : True)
    blur_radii : list of float
        Translate each segment to given radii to assign a broader fiber contribution (default : [])
    blur_samples : list of integer
        Segments are duplicated along a circle at a given radius; this parameter controls the number of samples to take over a given circle (defaut : [])
    blur_sigma
        The contributions of the segments at different radii are damped as a Gaussian (default : 1.0)
    """

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
        raise RuntimeError( 'fiber_shift must be a scalar or a vector with 3 elements' )

    tic = time.time()
    print '\n-> Creating the dictionary from tractogram:'
    print '\t* Segment position = %s' % ( 'COMPUTE INTERSECTIONS' if do_intersect else 'CENTROID' )
    print '\t* Fiber shift X    = %.3f (voxel-size units)' % fiber_shiftX
    print '\t* Fiber shift Y    = %.3f (voxel-size units)' % fiber_shiftY
    print '\t* Fiber shift Z    = %.3f (voxel-size units)' % fiber_shiftZ
    print '\t* Points to skip   = %d' % points_to_skip
    print '\t* Min segment len  = %.2e' % min_seg_len

    # check blur params
    cdef :
        double [:] blurRadii
        int [:] blurSamples
        double [:] blurWeights
        double* ptrBlurRadii
        int* ptrBlurSamples
        double* ptrBlurWeights
        int nBlurRadii

    if len(blur_radii) != len(blur_samples) :
        raise RuntimeError( 'number of radii and samples must match' )

    # convert to numpy arrays (add fake radius for original segment)
    nBlurRadii = len(blur_radii)+1
    blurRadii = np.array( [0.0]+blur_radii, np.double )
    blurSamples = np.array( [1]+blur_samples, np.int32 )

    # compute weights for gaussian damping
    blurWeights = np.empty_like( blurRadii )
    for i in xrange(nBlurRadii):
        blurWeights[i] = np.exp( -blurRadii[i]**2 / (2.0*blur_sigma**2) )

    if nBlurRadii == 1 :
        print '\t* Do not blur fibers'
    else :
        print '\t* Blur fibers :'
        print '\t\t- sigma = %.3f' % blur_sigma
        print '\t\t- radii =   [',
        for i in xrange( 1, blurRadii.size ) :
            print '%.3f' % blurRadii[i],
        print ']'
        print '\t\t- samples = [',
        for i in xrange( 1, blurSamples.size ) :
            print '%5d' % blurSamples[i],
        print ']'
        print '\t\t- weights = [',
        for i in xrange( 1, blurWeights.size ) :
            print '%.3f' % blurWeights[i],
        print ']'

    ptrBlurRadii   = &blurRadii[0]
    ptrBlurSamples = &blurSamples[0]
    ptrBlurWeights = &blurWeights[0]

    # minimum segment length
    if min_seg_len < 0 :
        raise RuntimeError( 'min_seg_len must be >= 0' )


    print '\t* Loading data:'
    
    # fiber-tracts from tractogram
    print '\t\t* tractogram'

    # fiber-tracts from tractogram
    print '\t\t* tractogram'
    extension = os.path.splitext(filename_tractogram)[1]  #take extension of file
    
    if extension != ".trk" and extension != ".tck" :
        raise IOError( 'Invalid input file. Please enter tractogram file .trk or .tck' )
    try :
        if (extension == ".trk") : # if file tractogram is .trk
            _, trk_hdr = nibabel.trackvis.read( filename_tractogram )
        else : # if file tractogram is .tck
            tck_hdr = nibabel.streamlines.tck.TckFile._read_header(filename_tractogram)
    except :
        raise IOError( 'Tractogram file not found' )
        
    if (extension == ".trk"):
        Nx = trk_hdr['dim'][0]
        Ny = trk_hdr['dim'][1]
        Nz = trk_hdr['dim'][2]
        Px = trk_hdr['voxel_size'][0]
        Py = trk_hdr['voxel_size'][1]
        Pz = trk_hdr['voxel_size'][2]

        data_offset = 1000
        n_count = trk_hdr['n_count']
        n_scalars = trk_hdr['n_scalars']
        n_properties = trk_hdr['n_properties']
        
    if (extension == ".tck"):
        #open file .nii and get header of this to get info on the structure
        if TCK_ref_image is None:
            if filename_peaks is not None:
                TCK_ref_image = filename_peaks
            elif filename_mask is not None:
                TCK_ref_image = filename_mask
            else:
                raise RuntimeError( 'TCK files do not contain info on the geometry. Use "TCK_ref_image" for that.' )

        #load the TCK_ref_image( .nii file) with nibabel
        nii_image = nibabel.load(TCK_ref_image)
        #read the header of nii file
        nii_hdr = nii_image.header if nibabel.__version__ >= '2.0.0' else nii_image.get_header()

        #set shape's of tractogram
        Nx = nii_image.shape[0]
        Ny = nii_image.shape[1]
        Nz = nii_image.shape[2]

        #set distance's of control points
        Px = nii_hdr['pixdim'][1]
        Py = nii_hdr['pixdim'][2]
        Pz = nii_hdr['pixdim'][3]

        #set offset and number of streamlines
        data_offset = int(tck_hdr['_offset_data'])  #set offset
        n_count = int(tck_hdr['count'])  #set number of fibers

        #set number of proprieties and numebr of scalar to zero, because there are not present in .tck file
        n_scalars = 0
        n_properties = 0

        
    print '\t\t\t- %d x %d x %d' % ( Nx, Ny, Nz )
    print '\t\t\t- %.4f x %.4f x %.4f' % ( Px, Py, Pz )
    print '\t\t\t- %d fibers' % trk_hdr['n_count']
    if Nx >= 2**16 or Nz >= 2**16 or Nz >= 2**16 :
        raise RuntimeError( 'The max dim size is 2^16 voxels' )

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
            print '\t\t  [WARNING] dataset does not have the same geometry as the tractogram'
        niiMASK_img = np.ascontiguousarray( niiMASK.get_data().astype(np.float32) )
        ptrMASK  = &niiMASK_img[0,0,0]
    else :
        print '\t\t* no mask specified to filter IC compartments'
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
        print '\t\t* EC orientations'
        niiPEAKS = nibabel.load( filename_peaks )
        niiPEAKS_hdr = niiPEAKS.header if nibabel.__version__ >= '2.0.0' else niiPEAKS.get_header()
        print '\t\t\t- %d x %d x %d x %d' % ( niiPEAKS.shape[0], niiPEAKS.shape[1], niiPEAKS.shape[2], niiPEAKS.shape[3] )
        print '\t\t\t- %.4f x %.4f x %.4f' % ( niiPEAKS_hdr['pixdim'][1], niiPEAKS_hdr['pixdim'][2], niiPEAKS_hdr['pixdim'][3] )
        print '\t\t\t- ignoring peaks < %.2f * MaxPeak' % vf_THR
        print '\t\t\t- %susing affine matrix' % ( "" if peaks_use_affine else "not " )
        print '\t\t\t- flipping axes : [ x=%s, y=%s, z=%s ]' % ( flip_peaks[0], flip_peaks[1], flip_peaks[2] )
        if ( Nx!=niiPEAKS.shape[0] or Ny!=niiPEAKS.shape[1] or Nz!=niiPEAKS.shape[2] or
             abs(Px-niiPEAKS_hdr['pixdim'][1])>1e-3 or abs(Py-niiPEAKS_hdr['pixdim'][2])>1e-3 or abs(Pz-niiPEAKS_hdr['pixdim'][3])>1e-3 ) :
            print "\t\t  [WARNING] dataset does not have the same geometry as the tractogram"
        if niiPEAKS.shape[3] % 3 :
            raise RuntimeError( 'PEAKS dataset must have 3*k volumes' )
        if vf_THR < 0 or vf_THR > 1 :
            raise RuntimeError( 'vf_THR must be between 0 and 1' )
        niiPEAKS_img = np.ascontiguousarray( niiPEAKS.get_data().astype(np.float32) )
        ptrPEAKS = &niiPEAKS_img[0,0,0,0]
        Np = niiPEAKS.shape[3]/3

        # affine matrix to rotate gradien directions (if required)
        if peaks_use_affine :
            affine = np.ascontiguousarray( niiPEAKS.affine[:3,:3].T )
        else :
            affine = np.ascontiguousarray( np.eye(3) )
        ptrAFFINE = &affine[0,0]
    else :
        print '\t\t* no dataset specified for EC compartments'
        Np = 0
        ptrPEAKS = NULL
        ptrAFFINE = NULL

    # output path
    print '\t\t* output written to "%s"' % path_out
    if not exists( path_out ):
        makedirs( path_out )

    # calling actual C code
    ret = trk2dictionary( filename_trk,
        trk_hdr['dim'][0], trk_hdr['dim'][1], trk_hdr['dim'][2],
        trk_hdr['voxel_size'][0], trk_hdr['voxel_size'][1], trk_hdr['voxel_size'][2],
        trk_hdr['n_count'], trk_hdr['n_scalars'], trk_hdr['n_properties'], fiber_shiftX, fiber_shiftY, fiber_shiftZ, points_to_skip, min_seg_len,
        ptrPEAKS, Np, vf_THR, -1 if flip_peaks[0] else 1, -1 if flip_peaks[1] else 1, -1 if flip_peaks[2] else 1,
        ptrMASK, ptrTDI, path_out, 1 if do_intersect else 0, ptrAFFINE,
        nBlurRadii, blur_sigma, ptrBlurRadii, ptrBlurSamples, ptrBlurWeights );
    if ret == 0 :
        print '   [ DICTIONARY not generated ]'
        return None

    # create new TRK with only fibers in the WM mask
    if gen_trk :
        print '\t* Generate tractogram matching the dictionary: '
        fib, _ = nibabel.trackvis.read( filename_trk, as_generator=True )
        fibKept = []
        file_kept = np.fromfile( join(path_out,'dictionary_TRK_kept.dict'), dtype=np.bool_ )
        ind = 0
        for f in fib:
            if file_kept[ind]:
                fibKept.append( (f[0],None, None) )
            ind = ind+1
        nibabel.trackvis.write( join(path_out,'dictionary_TRK_fibers.trk'), fibKept, trk_hdr )
        print '\t  [ %d fibers kept ]' % np.count_nonzero( file_kept )
    print '   [ %.1f seconds ]' % ( time.time() - tic )

    # save TDI and MASK maps
    if filename_mask is not None :
        affine = niiMASK.affine if nibabel.__version__ >= '2.0.0' else niiMASK.get_affine()
    elif filename_peaks is not None :
        affine = niiPEAKS.affine if nibabel.__version__ >= '2.0.0' else niiPEAKS.get_affine()
    else :
        affine = np.diag( [Px, Py, Pz, 1] )

    niiTDI = nibabel.Nifti1Image( niiTDI_img, affine )
    nibabel.save( niiTDI, join(path_out,'dictionary_tdi.nii.gz') )

    if filename_mask is not None :
        niiMASK = nibabel.Nifti1Image( niiMASK_img, affine )
    else :
        niiMASK = nibabel.Nifti1Image( (np.asarray(niiTDI_img)>0).astype(np.float32), affine )
    nibabel.save( niiMASK, join(path_out,'dictionary_mask.nii.gz') )


cpdef convert_old_dictionary( path ):
    """Perform the conversion of the files representing a dictionary, i.e. dictionary_*.dict,
    from the old format to the new one, where the files *_{vx,vy,vz}.dict are replaced
    by a single file *_v.dict (same for the files *_{ox,oy}.dict).

    Parameters
    ----------
    path : string
        Path to the folder containing the dictionary_*.dict files.
    """
    if not exists( join(path,'dictionary_IC_vx.dict') ):
        raise RuntimeError( 'Folder does not contain dictionary files in the old format' )

    niiTDI = nibabel.load( join(path,'dictionary_tdi.nii.gz') )
    Nx, Ny, Nz = niiTDI.shape[:3]
    x = np.fromfile( join(path,'dictionary_IC_vx.dict'), dtype=np.uint16 ).astype(np.uint32)
    y = np.fromfile( join(path,'dictionary_IC_vy.dict'), dtype=np.uint16 ).astype(np.uint32)
    z = np.fromfile( join(path,'dictionary_IC_vz.dict'), dtype=np.uint16 ).astype(np.uint32)
    v = x + Nx * ( y + Ny * z )
    v.tofile( join(path,'dictionary_IC_v.dict') )
    remove( join(path,'dictionary_IC_vx.dict') )
    remove( join(path,'dictionary_IC_vy.dict') )
    remove( join(path,'dictionary_IC_vz.dict') )

    x = np.fromfile( join(path,'dictionary_EC_vx.dict'), dtype=np.uint8 ).astype(np.uint32)
    y = np.fromfile( join(path,'dictionary_EC_vy.dict'), dtype=np.uint8 ).astype(np.uint32)
    z = np.fromfile( join(path,'dictionary_EC_vz.dict'), dtype=np.uint8 ).astype(np.uint32)
    v = x + Nx * ( y + Ny * z )
    v.tofile( join(path,'dictionary_EC_v.dict') )
    remove( join(path,'dictionary_EC_vx.dict') )
    remove( join(path,'dictionary_EC_vy.dict') )
    remove( join(path,'dictionary_EC_vz.dict') )

    x = np.fromfile( join(path,'dictionary_IC_ox.dict'), dtype=np.uint8 ).astype(np.uint16)
    y = np.fromfile( join(path,'dictionary_IC_oy.dict'), dtype=np.uint8 ).astype(np.uint16)
    v = y + 181 * x
    v.tofile( join(path,'dictionary_IC_o.dict') )
    remove( join(path,'dictionary_IC_ox.dict') )
    remove( join(path,'dictionary_IC_oy.dict') )

    x = np.fromfile( join(path,'dictionary_EC_ox.dict'), dtype=np.uint8 ).astype(np.uint16)
    y = np.fromfile( join(path,'dictionary_EC_oy.dict'), dtype=np.uint8 ).astype(np.uint16)
    v = y + 181 * x
    v.tofile( join(path,'dictionary_EC_o.dict') )
    remove( join(path,'dictionary_EC_ox.dict') )
    remove( join(path,'dictionary_EC_oy.dict') )
