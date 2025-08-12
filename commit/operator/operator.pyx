#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False

import numpy as np

from dicelib.ui import setup_logger

# Interfaces to actual C code performing the multiplications
cdef extern void COMMIT_A(
    int _nF, int _nE, int _nV, int _nS, int _ndirs,
    double *_v_in, double *_v_out,
    unsigned int *_ICf, unsigned int *_ICeval, unsigned int *_ICv, unsigned short *_ICo, float *_ICl,
    unsigned int *_ECv, unsigned short *_ECo,
    unsigned int *_ISOv,
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
    unsigned int* _ICthreads, unsigned int* _ECthreads, unsigned int* _ISOthreads,
    unsigned int _nIC, unsigned int _nEC, unsigned int _nISO, unsigned int _nThreads
) nogil

cdef extern void COMMIT_At(
    int _nF, int _n, int _nE, int _nV, int _nS, int _ndirs,
    double *_v_in, double *_v_out,
    unsigned int *_ICf, unsigned int *_ICeval, unsigned int *_ICv, unsigned short *_ICo, float *_ICl,
    unsigned int *_ECv, unsigned short *_ECo,
    unsigned int *_ISOv,
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
    unsigned char* _ICthreadsT, unsigned int* _ECthreadsT, unsigned int* _ISOthreadsT,
    unsigned int _nIC, unsigned int _nEC, unsigned int _nISO, unsigned int _nThreads
) nogil

cdef extern void COMMIT_A_nolut(
    int _nF,
    double *_v_in, double *_v_out,
    unsigned int *_ICf,  unsigned int *_ICeval, unsigned int *_ICv, float *_ICl,
    unsigned int *_ISOv,
    unsigned int* _ICthreads, unsigned int* _ISOthreads,
    unsigned int _nISO, unsigned int _nThreads
) nogil

cdef extern void COMMIT_At_nolut(
    int _nF, int _n,
    double *_v_in, double *_v_out,
    unsigned int *_ICf,  unsigned int *_ICeval, unsigned int *_ICv, float *_ICl,
    unsigned int *_ISOv,
    unsigned char* _ICthreadsT, unsigned int* _ISOthreadsT,
    unsigned int _nISO, unsigned int _nThreads
) nogil

logger = setup_logger('operator')

cdef class LinearOperator :
    """This class is a wrapper to the C code for performing matrix-vector multiplications
    with the COMMIT linear operator A. The calculations are done using C code
    that uses information from the DICTIONARY, KERNELS and THREADS data structures.
    """
    cdef public int adjoint, nROWS, nCOLS

    cdef:
        DICTIONARY, KERNELS, THREADS
        int nSAMPLES, ndirs, nolut
        int ICn, ICnVOX, ICnSTR, ICnRF, ECn, ECnRF, ISOn, ISOnRF
        unsigned int *ICf, *ICv, *ECv, *ISOv, *ICeval
        float *ICl, *LUT_IC, *LUT_EC, *LUT_ISO
        unsigned short *ICo, *ECo
        unsigned int *ICthreads, *ECthreads, *ISOthreads, *ECthreadsT, *ISOthreadsT
        unsigned char *ICthreadsT


    def __init__( self, DICTIONARY, KERNELS, THREADS, nolut=False ) :
        """Set the pointers to the data structures used by the C code."""
        self.DICTIONARY = DICTIONARY
        self.KERNELS    = KERNELS
        self.THREADS    = THREADS
        self.nolut      = nolut

        self.ICn        = DICTIONARY['IC']['n']     # number of IC contributions (i.e. segments)
        self.ICnSTR     = DICTIONARY['IC']['nSTR']  # number of IC streamlines
        self.ICnRF      = KERNELS['wmr'].shape[0]   # number of IC response functions
        self.ICnVOX     = DICTIONARY['IC']['nVOX']  # number of IC voxels
        self.ECn        = DICTIONARY['EC']['n']     # number of EC contributions (i.e. peaks)
        self.ECnRF      = KERNELS['wmh'].shape[0]   # number of EC response functions
        self.ISOn       = DICTIONARY['ISO']['n']    # number of ISO contributions (i.e. voxels)
        self.ISOnRF     = KERNELS['iso'].shape[0]   # number of ISO response functions

        # number of SAMPLES and SAMPLING DIRECTIONS
        if KERNELS['wmr'].size > 0 :
            self.nSAMPLES = KERNELS['wmr'].shape[2]
        elif KERNELS['wmh'].size > 0 :
            self.nSAMPLES = KERNELS['wmh'].shape[2]
        else :
            self.nSAMPLES = KERNELS['wmr'].shape[1]
        self.ndirs = KERNELS['wmr'].shape[1]   # number of directions

        self.adjoint    = 0 # direct of inverse product
        self.nROWS = self.ICnVOX*self.nSAMPLES
        self.nCOLS = self.ICnSTR*self.ICnRF + self.ECn*self.ECnRF + self.ISOn*self.ISOnRF

        # get C pointers to arrays in DICTIONARY
        cdef unsigned int [::1]   ICf  = DICTIONARY['IC']['str']
        self.ICf = &ICf[0]
        cdef unsigned int [::1]   ICeval = DICTIONARY['IC']['eval']
        self.ICeval = &ICeval[0]
        cdef float [::1]          ICl  = DICTIONARY['IC']['len']
        self.ICl = &ICl[0]
        cdef unsigned int [::1]   ICv  = DICTIONARY['IC']['vox']
        self.ICv = &ICv[0]
        cdef unsigned short [::1] ICo  = DICTIONARY['IC']['dir']
        self.ICo = &ICo[0]
        cdef unsigned int [::1]   ECv  = DICTIONARY['EC']['vox']
        self.ECv = &ECv[0]
        cdef unsigned short [::1] ECo  = DICTIONARY['EC']['dir']
        self.ECo = &ECo[0]
        cdef unsigned int [::1]   ISOv = DICTIONARY['ISO']['vox']
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


    @property
    def T( self ) :
        """Transpose of the explicit matrix."""
        C = LinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS, self.nolut )
        C.adjoint = 1 - C.adjoint
        return C


    @property
    def shape( self ) :
        """Size of the explicit matrix."""
        if not self.adjoint :
            return ( self.nROWS, self.nCOLS )
        else :
            return ( self.nCOLS, self.nROWS )


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
                        self.ICnSTR,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICeval, self.ICv, self.ICl,
                        self.ISOv,
                        self.ICthreads, self.ISOthreads,
                        nISO, nthreads
                    )
            else:
                with nogil:
                    COMMIT_A(
                        self.ICnSTR, self.ECn, self.ISOn, self.nSAMPLES, self.ndirs,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICeval, self.ICv, self.ICo, self.ICl,
                        self.ECv, self.ECo,
                        self.ISOv,
                        self.LUT_IC, self.LUT_EC, self.LUT_ISO,
                        self.ICthreads, self.ECthreads, self.ISOthreads,
                        nIC, nEC, nISO, nthreads
                    )
        else :
            # INVERSE PRODUCT A'*y
            if self.nolut:
                with nogil:
                    COMMIT_At_nolut(
                        self.ICnSTR, self.ICn,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICeval, self.ICv, self.ICl,
                        self.ISOv,
                        self.ICthreadsT, self.ISOthreadsT,
                        nISO, nthreads
                    )
            else:
                with nogil:
                    COMMIT_At(
                        self.ICnSTR, self.ICn, self.ECn, self.ISOn, self.nSAMPLES, self.ndirs,
                        &v_in[0], &v_out[0],
                        self.ICf, self.ICeval, self.ICv, self.ICo, self.ICl,
                        self.ECv, self.ECo,
                        self.ISOv,
                        self.LUT_IC, self.LUT_EC, self.LUT_ISO,
                        self.ICthreadsT, self.ECthreadsT, self.ISOthreadsT,
                        nIC, nEC, nISO, nthreads
                    )

        return v_out
