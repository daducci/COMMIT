#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False

import numpy as np

from dicelib.ui import setup_logger

# Interfaces to actual C code performing the multiplications
cdef extern void COMMIT_A(
    int _nF, int _nE, int _nV, int _nS, int _ndirs,
    double *_v_in, double *_v_out,
    unsigned int *_ICf, unsigned int *_ICv, unsigned short *_ICo, float *_ICl,
    unsigned int *_ECv, unsigned short *_ECo,
    unsigned int *_ISOv,
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
    unsigned int* _ICthreads, unsigned int* _ECthreads, unsigned int* _ISOthreads,
    unsigned int _nIC, unsigned int _nEC, unsigned int _nISO, unsigned int _nThreads
) nogil

cdef extern void COMMIT_At(
    int _nF, int _n, int _nE, int _nV, int _nS, int _ndirs,
    double *_v_in, double *_v_out,
    unsigned int *_ICf, unsigned int *_ICv, unsigned short *_ICo, float *_ICl,
    unsigned int *_ECv, unsigned short *_ECo,
    unsigned int *_ISOv,
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
    unsigned char* _ICthreadsT, unsigned int* _ECthreadsT, unsigned int* _ISOthreadsT,
    unsigned int _nIC, unsigned int _nEC, unsigned int _nISO, unsigned int _nThreads
) nogil

cdef extern void Tikhonov(
    int _nF,
    double *_v_in, double *_v_out,
    unsigned int *_ISOv,
    double _lambda,
    unsigned int *_ISOthreads, unsigned int _nISO, unsigned int _nThreads,
    unsigned int *_neighbors, unsigned int *_indptr
) nogil

cdef extern void Tikhonov_t(
    int _nF,
    double *_v_in, double *_v_out,
    unsigned int *_ISOv,
    double _lambda,
    unsigned int *_ISOthreadsT, unsigned int _nISO, unsigned int _nThreads,
    unsigned int *_neighbors, unsigned int *_indptr
) nogil

cdef extern void COMMIT_A_nolut(
    int _nF,
    double *_v_in, double *_v_out,
    unsigned int *_ICf, unsigned int *_ICv, float *_ICl,
    unsigned int *_ISOv,
    unsigned int* _ICthreads, unsigned int* _ISOthreads,
    unsigned int _nISO, unsigned int _nThreads
) nogil

cdef extern void COMMIT_At_nolut(
    int _nF, int _n,
    double *_v_in, double *_v_out,
    unsigned int *_ICf, unsigned int *_ICv, float *_ICl,
    unsigned int *_ISOv,
    unsigned char* _ICthreadsT, unsigned int* _ISOthreadsT,
    unsigned int _nISO, unsigned int _nThreads
) nogil

logger = setup_logger('operator')


def precompute_neighbours(mask, mask_ix, mask_iy, mask_iz):
    idx = mask.ravel(order='F').nonzero()[0]
    N = idx.size  # Number of voxels in the mask

    # Create the lookup table (lut)
    lut = np.full(mask.size, -1, dtype=np.int32)
    lut[idx] = np.arange(N)

    # Define neighbour offsets (6-connectivity)
    neighbor_offsets = np.array([
        [-1, 0, 0],  # Left
        [1, 0, 0],   # Right
        [0, -1, 0],  # Down
        [0, 1, 0],   # Up
        [0, 0, -1],  # Back
        [0, 0, 1]    # Front
    ], dtype=np.int32)

    # Prepare voxel coordinates
    x_i = mask_ix[:, np.newaxis]
    y_i = mask_iy[:, np.newaxis]
    z_i = mask_iz[:, np.newaxis]

    # Compute neighbour coordinates
    x_n = x_i + neighbor_offsets[:, 0]
    y_n = y_i + neighbor_offsets[:, 1]
    z_n = z_i + neighbor_offsets[:, 2]

    # Check boundaries
    valid = (
        (x_n >= 0) & (x_n < mask.shape[0]) &
        (y_n >= 0) & (y_n < mask.shape[1]) &
        (z_n >= 0) & (z_n < mask.shape[2])
    )

    # Flatten the valid mask and coordinates
    valid_flat = valid.ravel()
    x_n_flat = x_n.ravel()[valid_flat]
    y_n_flat = y_n.ravel()[valid_flat]
    z_n_flat = z_n.ravel()[valid_flat]

    # Original voxel indices repeated for each neighbour
    voxel_indices = np.repeat(np.arange(N), neighbor_offsets.shape[0])[valid_flat]

    # Convert neighbour coordinates to flat indices
    idx_n = np.ravel_multi_index(
        (x_n_flat, y_n_flat, z_n_flat),
        mask.shape,
        order='F'
    )

    # Map to local indices
    neighbor_local_indices = lut[idx_n]

    # Filter out neighbours not in the mask
    valid_neighbors = neighbor_local_indices >= 0
    voxel_indices = voxel_indices[valid_neighbors]
    neighbor_local_indices = neighbor_local_indices[valid_neighbors]

    # Sort voxel_indices and neighbor_local_indices based on voxel_indices
    sort_order = np.argsort(voxel_indices)
    voxel_indices = voxel_indices[sort_order]
    neighbor_local_indices = neighbor_local_indices[sort_order]

    # Build the indptr array
    counts = np.bincount(voxel_indices, minlength=N)
    indptr = np.zeros(N + 1, dtype=np.uint32)
    indptr[1:] = np.cumsum(counts)

    # The neighbour indices array
    neighbours = neighbor_local_indices.astype(np.uint32)


    return neighbours, indptr




cdef class LinearOperator :
    """This class is a wrapper to the C code for performing marix-vector multiplications
    with the COMMIT linear operator A. The multiplications are done using C code
    that uses information from the DICTIONARY, KERNELS and THREADS data structures.
    """
    cdef int nS, nF, nR, nE, nT, nV, nI, n, ndirs
    cdef public int adjoint, n1, n2
    cdef public double tikhonov_lambda

    cdef DICTIONARY
    cdef KERNELS
    cdef THREADS
    cdef nolut

    cdef unsigned int*   ICf
    cdef float*          ICl
    cdef unsigned int*   ICv
    cdef unsigned short* ICo
    cdef unsigned int*   ECv
    cdef unsigned short* ECo
    cdef unsigned int*   ISOv

    cdef float* LUT_IC
    cdef float* LUT_EC
    cdef float* LUT_ISO

    cdef unsigned int*  ICthreads
    cdef unsigned int*  ECthreads
    cdef unsigned int*  ISOthreads

    cdef unsigned char* ICthreadsT
    cdef unsigned int*  ECthreadsT
    cdef unsigned int*  ISOthreadsT

    cdef unsigned int*  neighbours
    cdef unsigned int*  indptr

    cdef public         neigh
    cdef public         neigh_nptr

    def __init__( self, DICTIONARY, KERNELS, THREADS, tikhonov_lambda=0, nolut=False ) :
        """Set the pointers to the data structures used by the C code."""
        self.DICTIONARY         = DICTIONARY
        self.KERNELS            = KERNELS
        self.THREADS            = THREADS
        self.nolut              = nolut

        self.nF                 = DICTIONARY['IC']['nF']    # number of FIBERS
        self.nR                 = KERNELS['wmr'].shape[0]   # number of FIBER RADII
        self.nE                 = DICTIONARY['EC']['nE']    # number of EC segments
        self.nT                 = KERNELS['wmh'].shape[0]   # number of EC TORTUOSITY values
        self.nV                 = DICTIONARY['nV']          # number of VOXELS
        self.nI                 = KERNELS['iso'].shape[0]   # number of ISO contributions
        self.n                  = DICTIONARY['IC']['n']     # numbner of IC segments
        self.ndirs              = KERNELS['wmr'].shape[1]   # number of directions
        self.tikhonov_lambda    = tikhonov_lambda           # equalizer parameter of the Tikhonov regularization term

        if KERNELS['wmr'].size > 0 :
            self.nS = KERNELS['wmr'].shape[2]       # number of SAMPLES
        elif KERNELS['wmh'].size > 0 :
            self.nS = KERNELS['wmh'].shape[2]
        else :
            self.nS = KERNELS['wmr'].shape[1]

        self.adjoint    = 0                         # direct of inverse product

        self.n1 = self.nV*self.nS
        self.n2 = self.nR*self.nF + self.nT*self.nE + self.nI*self.nV

        # get C pointers to arrays in DICTIONARY
        cdef unsigned int [::1]   ICf  = DICTIONARY['IC']['fiber']
        self.ICf = &ICf[0]
        cdef float [::1]          ICl  = DICTIONARY['IC']['len']
        self.ICl = &ICl[0]
        cdef unsigned int [::1]   ICv  = DICTIONARY['IC']['v']
        self.ICv = &ICv[0]
        cdef unsigned short [::1] ICo  = DICTIONARY['IC']['o']
        self.ICo = &ICo[0]
        cdef unsigned int [::1]   ECv  = DICTIONARY['EC']['v']
        self.ECv = &ECv[0]
        cdef unsigned short [::1] ECo  = DICTIONARY['EC']['o']
        self.ECo = &ECo[0]
        cdef unsigned int [::1]   ISOv = DICTIONARY['ISO']['v']
        self.ISOv = &ISOv[0]

        # get C pointers to arrays in KERNELS
        cdef float [:, :, ::1] wmrSFP = KERNELS['wmr']
        self.LUT_IC  = &wmrSFP[0,0,0]

        cdef float [:, :, ::1] wmhSFP = KERNELS['wmh']
        self.LUT_EC  = &wmhSFP[0,0,0]
        cdef float [:, ::1] isoSFP = KERNELS['iso']
        self.LUT_ISO = &isoSFP[0,0]

        # get C pointers to arrays in THREADS
        cdef unsigned int [::1] ICthreads = THREADS['IC']
        self.ICthreads  = &ICthreads[0]
        cdef unsigned int [::1] ECthreads = THREADS['EC']
        self.ECthreads  = &ECthreads[0]
        cdef unsigned int [::1] ISOthreads = THREADS['ISO']
        self.ISOthreads = &ISOthreads[0]

        cdef unsigned char [::1] ICthreadsT = THREADS['ICt']
        self.ICthreadsT  = &ICthreadsT[0]
        cdef unsigned int  [::1] ECthreadsT = THREADS['ECt']
        self.ECthreadsT  = &ECthreadsT[0]
        cdef unsigned int  [::1] ISOthreadsT = THREADS['ISOt']
        self.ISOthreadsT = &ISOthreadsT[0]

        cdef unsigned int [::1] neighbours
        cdef unsigned int [::1] indptr

        if tikhonov_lambda > 0:
            # precompute neighbours
            self.neigh, self.neigh_nptr = precompute_neighbours(DICTIONARY['MASK'], DICTIONARY['MASK_ix'], DICTIONARY['MASK_iy'], DICTIONARY['MASK_iz'])

            neighbours = self.neigh
            self.neighbours = &neighbours[0]

            indptr = self.neigh_nptr
            self.indptr = &indptr[0]


    @property
    def T( self ) :
        """Transpose of the explicit matrix."""
        C = LinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS, self.tikhonov_lambda, self.nolut )
        C.adjoint = 1 - C.adjoint
        return C


    @property
    def shape( self ) :
        """Size of the explicit matrix."""
        if not self.adjoint :
            return ( self.n1, self.n2 )
        else :
            return ( self.n2, self.n1 )


    def dot( self, double [::1] v_in  ):
        """Wrapper to C code for efficiently performing the matrix-vector multiplications.

        Parameters
        ----------
        v_in : 1D numpy.array of double
            Input vector for the matrix-vector multiplication

        Returns
        -------
        v_out : 1D numpy.array of double
            Results of the multiplication
        """

        # Permit only matrix-vector multiplications
        if v_in.size != self.shape[1] :
            raise RuntimeError( "A.dot(): dimensions do not match" )

        # Create output array
        cdef double [::1] v_out = np.zeros( self.shape[0], dtype=np.float64 )

        cdef unsigned int nthreads = self.THREADS['n']
        cdef unsigned int nIC = self.KERNELS['wmr'].shape[0]
        cdef unsigned int nEC = self.KERNELS['wmh'].shape[0]
        cdef unsigned int nISO = self.KERNELS['iso'].shape[0]

        # Call the cython function to read the memory pointers
        if not self.adjoint :
            # DIRECT PRODUCT A*x
            if self.nolut:
                with nogil:
                    COMMIT_A_nolut(
                        self.nF,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICv, self.ICl,
                        self.ISOv,
                        self.ICthreads, self.ISOthreads,
                        nISO, nthreads
                    )
            else:
                with nogil:
                    COMMIT_A(
                        self.nF, self.nE, self.nV, self.nS, self.ndirs,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICv, self.ICo, self.ICl,
                        self.ECv, self.ECo,
                        self.ISOv,
                        self.LUT_IC, self.LUT_EC, self.LUT_ISO,
                        self.ICthreads, self.ECthreads, self.ISOthreads,
                        nIC, nEC, nISO, nthreads
                    )
                if self.tikhonov_lambda > 0:
                    with nogil:
                        Tikhonov(
                            self.nF,
                            &v_in[0], &v_out[0],
                            self.ISOv,
                            self.tikhonov_lambda,
                            self.ISOthreads, nISO, nthreads,
                            self.neighbours, self.indptr
                        )
        else :
            # INVERSE PRODUCT A'*y
            if self.nolut:
                with nogil:
                    COMMIT_At_nolut(
                        self.nF, self.n,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICv, self.ICl,
                        self.ISOv,
                        self.ICthreadsT, self.ISOthreadsT,
                        nISO, nthreads
                    )
            else:
                with nogil:
                    COMMIT_At(
                        self.nF, self.n, self.nE, self.nV, self.nS, self.ndirs,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICv, self.ICo, self.ICl,
                        self.ECv, self.ECo,
                        self.ISOv,
                        self.LUT_IC, self.LUT_EC, self.LUT_ISO,
                        self.ICthreadsT, self.ECthreadsT, self.ISOthreadsT,
                        nIC, nEC, nISO, nthreads
                    )
                if self.tikhonov_lambda > 0:
                    with nogil:
                        Tikhonov_t(
                            self.nF,
                            &v_in[0], &v_out[0],
                            self.ISOv,
                            self.tikhonov_lambda,
                            self.ISOthreads, nISO, nthreads,
                            self.neighbours, self.indptr
                        )

        return v_out
