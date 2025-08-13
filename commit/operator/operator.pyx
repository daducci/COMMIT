#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False

import numpy as np

from dicelib.ui import setup_logger

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
        unsigned int [::1] ICf
        unsigned int [::1] ICv
        unsigned int [::1] ECv
        unsigned int [::1] ISOv
        unsigned int [::1] ICeval
        float [::1]        ICl
        float [:, :, ::1] LUT_IC
        float [:, :, ::1] LUT_EC
        float [:, ::1]    LUT_ISO
        unsigned short [::1] ICo
        unsigned short [::1] ECo
        unsigned int [::1]  ICthreads
        unsigned int [::1]  ECthreads
        unsigned int [::1]  ISOthreads
        unsigned int [::1]  ECthreadsT
        unsigned int [::1]  ISOthreadsT
        unsigned char [::1] ICthreadsT
        unsigned int nThreads


    def __init__( self, DICTIONARY, KERNELS, THREADS, nolut=False ) :
        """Set the pointers to the data structures used by the C code."""
        self.DICTIONARY = DICTIONARY
        self.KERNELS    = KERNELS
        self.THREADS    = THREADS

        self.ICn    = DICTIONARY['IC']['n']     # number of IC contributions (i.e. segments)
        self.ICnSTR = DICTIONARY['IC']['nSTR']  # number of IC streamlines
        self.ICnRF  = KERNELS['wmr'].shape[0]   # number of IC response functions
        self.ICnVOX = DICTIONARY['IC']['nVOX']  # number of IC voxels
        self.ECn    = DICTIONARY['EC']['n']     # number of EC contributions (i.e. peaks)
        self.ECnRF  = KERNELS['wmh'].shape[0]   # number of EC response functions
        self.ISOn   = DICTIONARY['ISO']['n']    # number of ISO contributions (i.e. voxels)
        self.ISOnRF = KERNELS['iso'].shape[0]   # number of ISO response functions

        # number of SAMPLES
        if KERNELS['wmr'].size > 0 :
            self.nSAMPLES = KERNELS['wmr'].shape[2]
        elif KERNELS['wmh'].size > 0 :
            self.nSAMPLES = KERNELS['wmh'].shape[2]
        else :
            self.nSAMPLES = KERNELS['wmr'].shape[1]
        # number of SAMPLING DIRECTIONS
        self.ndirs = KERNELS['wmr'].shape[1]
        # use lut or not
        self.nolut = nolut
        # direct of inverse product
        self.adjoint = 0
        # actual size of matrix A
        self.nROWS   = self.ICnVOX*self.nSAMPLES
        self.nCOLS   = self.ICnSTR*self.ICnRF + self.ECn*self.ECnRF + self.ISOn*self.ISOnRF

        # get pointers to arrays in DICTIONARY
        self.ICf    = DICTIONARY['IC']['str']
        self.ICl    = DICTIONARY['IC']['len']
        self.ICv    = DICTIONARY['IC']['vox']
        self.ICo    = DICTIONARY['IC']['dir']
        self.ICeval = DICTIONARY['IC']['eval']
        self.ECv    = DICTIONARY['EC']['vox']
        self.ECo    = DICTIONARY['EC']['dir']
        self.ISOv   = DICTIONARY['ISO']['vox']

        # get pointers to arrays in KERNELS
        self.LUT_IC  = KERNELS['wmr']
        self.LUT_EC  = KERNELS['wmh']
        self.LUT_ISO = KERNELS['iso']

        # get pointers to arrays in THREADS
        self.nThreads    = THREADS['n']
        self.ICthreads   = THREADS['IC']
        self.ECthreads   = THREADS['EC']
        self.ISOthreads  = THREADS['ISO']
        self.ICthreadsT  = THREADS['ICt']
        self.ECthreadsT  = THREADS['ECt']
        self.ISOthreadsT = THREADS['ISOt']


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

        # Call the c++ function to perform the multiplications
        with nogil:
            if not self.adjoint :
                # DIRECT PRODUCT A*x
                if self.nolut:
                    COMMIT_A_nolut(
                        &v_in[0], &v_out[0],
                        self.ICnSTR,
                        &self.ICf[0], &self.ICeval[0], &self.ICv[0], &self.ICl[0],
                        &self.ISOv[0],
                        &self.ICthreads[0], &self.ISOthreads[0],
                        self.ISOnRF, self.nThreads
                    )
                else:
                    COMMIT_A(
                        &v_in[0], &v_out[0],
                        self.nSAMPLES, self.ndirs,
                        self.ICnSTR, self.ECn, self.ISOn,
                        &self.ICf[0], &self.ICeval[0], &self.ICv[0], &self.ICo[0], &self.ICl[0],
                        &self.ECv[0], &self.ECo[0],
                        &self.ISOv[0],
                        &self.LUT_IC[0,0,0], &self.LUT_EC[0,0,0], &self.LUT_ISO[0,0],
                        &self.ICthreads[0], &self.ECthreads[0], &self.ISOthreads[0],
                        self.ICnRF, self.ECnRF, self.ISOnRF, self.nThreads
                    )
            else :
                # INVERSE PRODUCT A'*y
                if self.nolut:
                    COMMIT_At_nolut(
                        &v_in[0], &v_out[0],
                        self.ICnSTR, self.ICn,
                        &self.ICf[0], &self.ICeval[0], &self.ICv[0], &self.ICl[0],
                        &self.ISOv[0],
                        &self.ICthreadsT[0], &self.ISOthreadsT[0],
                        self.ISOnRF, self.nThreads
                    )
                else:
                    COMMIT_At(
                        &v_in[0], &v_out[0],
                        self.nSAMPLES, self.ndirs,
                        self.ICnSTR, self.ICn, self.ECn, self.ISOn,
                        &self.ICf[0], &self.ICeval[0], &self.ICv[0], &self.ICo[0], &self.ICl[0],
                        &self.ECv[0], &self.ECo[0],
                        &self.ISOv[0],
                        &self.LUT_IC[0,0,0], &self.LUT_EC[0,0,0], &self.LUT_ISO[0,0],
                        &self.ICthreadsT[0], &self.ECthreadsT[0], &self.ISOthreadsT[0],
                        self.ICnRF, self.ECnRF, self.ISOnRF, self.nThreads
                    )

        return v_out


# Interfaces of the external C functions performing the multiplications
# =====================================================================
cdef extern void COMMIT_A(
    double *, double *,
    int, int, int, int, int,
    unsigned int *, unsigned int *, unsigned int *, unsigned short *, float *,
    unsigned int *, unsigned short *,
    unsigned int *,
    float *, float *, float *,
    unsigned int*, unsigned int*, unsigned int*,
    unsigned int, unsigned int, unsigned int, unsigned int
) nogil

cdef extern void COMMIT_At(
    double *, double *,
    int, int, int, int, int, int,
    unsigned int *, unsigned int *, unsigned int *, unsigned short *, float *,
    unsigned int *, unsigned short *,
    unsigned int *,
    float *, float *, float *,
    unsigned char*, unsigned int*, unsigned int*,
    unsigned int, unsigned int, unsigned int, unsigned int
) nogil

cdef extern void COMMIT_A_nolut(
    double *, double *,
    int,
    unsigned int *,  unsigned int *, unsigned int *, float *,
    unsigned int *,
    unsigned int *, unsigned int *,
    unsigned int, unsigned int
) nogil

cdef extern void COMMIT_At_nolut(
    double *, double *,
    int, int,
    unsigned int *,  unsigned int *, unsigned int *, float *,
    unsigned int *,
    unsigned char*, unsigned int*,
    unsigned int, unsigned int
) nogil
