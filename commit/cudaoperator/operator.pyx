#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False

import cython
import numpy as np
cimport numpy as np
from amico.util import ERROR, LOG

cdef extern from "operator_withCUDA.cuh":
    int checkCompatibility(np.uint64_t, int)

cdef extern from "operator_withCUDA.cuh":
    cdef cppclass C_CudaLinearOperator "CudaLinearOperator":
        C_CudaLinearOperator(
            np.uint32_t*,
            np.uint32_t*,
            np.uint16_t*,
            np.float32_t*,
            np.float32_t*,

            np.uint32_t*,
            np.uint16_t*,
            np.float32_t*,

            np.float32_t*,

            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            
            int,
            int)

        int   getCudaStatus()
        void  setTransposeData(np.uint32_t*, np.uint32_t*, np.uint16_t*, np.float32_t*)
        void  destroy()
        void  dot(np.float64_t*, np.float64_t*)
        void Tdot(np.float64_t*, np.float64_t*)

cdef class CudaLinearOperator :
    """This class is a wrapper to the CUDA C++ code for performing marix-vector multiplications
    with the COMMIT linear operator A in a CUDA GPU. The multiplications are done using CUDA C++ code
    that uses information from the DICTIONARY and KERNELS data structures.
    """
    cdef int nS, nF, nR, nE, nT, nV, nI, n, ndirs, gpu_id
    cdef public int adjoint, n1, n2

    cdef DICTIONARY
    cdef KERNELS
    cdef THREADS

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

    # pointer to this operator in GPU memory
    cdef C_CudaLinearOperator* thisptr

    # these should be always None, they remain for compatibility
    cdef unsigned int*   ICthreads
    cdef unsigned int*   ECthreads
    cdef unsigned int*   ISOthreads
    cdef unsigned char*  ICthreadsT
    cdef unsigned int*   ECthreadsT
    cdef unsigned int*   ISOthreadsT


    def __init__( self, DICTIONARY, KERNELS, THREADS, fcall = 0 ) :
        """Set the pointers to the data structures used by the C code."""
        self.DICTIONARY = DICTIONARY
        self.KERNELS    = KERNELS
        self.THREADS    = THREADS

        self.nF         = DICTIONARY['IC']['nF']    # number of FIBERS
        self.nR         = KERNELS['wmr'].shape[0]   # number of FIBER RADII
        self.nE         = DICTIONARY['EC']['nE']    # number of EC segments
        self.nT         = KERNELS['wmh'].shape[0]   # number of EC TORTUOSITY values
        self.nV         = DICTIONARY['nV']          # number of VOXELS
        self.nI         = KERNELS['iso'].shape[0]   # number of ISO contributions
        self.n          = DICTIONARY['IC']['n']     # numbner of IC segments
        self.ndirs      = KERNELS['wmr'].shape[1]   # number of directions
        self.gpu_id     = THREADS['GPUID']          # id of the CUDA GPU

        if KERNELS['wmr'].size > 0 :
            self.nS = KERNELS['wmr'].shape[2]       # number of SAMPLES
        elif KERNELS['wmh'].size > 0 :
            self.nS = KERNELS['wmh'].shape[2]
        else :
            self.nS = KERNELS['wmr'].shape[1]

        self.adjoint = 0                            # direct of inverse product

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

        LOG( '\n-> Checking CUDA GPU:' )
        #cdef unsigned long long required_mem = 28*self.n + 6*self.nzeppelins + 8.0*(size_t)nfibers + 16.0*(size_t)nvoxels + 4.0*((size_t)size_lutic + (size_t)size_lutec + (size_t)size_lutiso + (size_t)this->nrows + (size_t)this->ncols)
        cdef int ans = checkCompatibility(0, self.gpu_id)
        if ans == 1:
            ERROR( 'The selected GPU is not detected; check "gpu_id" in "set_threads()"' )
        elif ans == 2:
            ERROR( 'Impossible to set GPU with ID=%d' % gpu_id )
        elif ans == 3:
            ERROR( 'Impossible to get properties from GPU with ID=%d' % gpu_id )
        elif ans == 4:
            ERROR( 'Compute capability must be at least 5.0' )

        # create the operator in GPU memory
        self.thisptr = new C_CudaLinearOperator(
            &ICv[0],
            &ICf[0],
            &ICo[0],
            &ICl[0],
            &wmrSFP[0,0,0],

            &ECv[0],
            &ECo[0],
            &wmhSFP[0,0,0],

            &isoSFP[0,0],

            self.n,
            self.nV,
            self.nF,
            self.nE,
            self.ndirs,
            self.nS,
            self.nR,
            self.nT,
            self.nI,
            
            fcall,
            self.gpu_id)

        # create the transpose of the operator in GPU memory
        if fcall == 1:
            idx = np.lexsort( [np.array(self.DICTIONARY['IC']['o']), np.array(self.DICTIONARY['IC']['fiber'])] )

            self.DICTIONARY['IC']['v']     = self.DICTIONARY['IC']['v'][ idx ]
            self.DICTIONARY['IC']['o']     = self.DICTIONARY['IC']['o'][ idx ]
            self.DICTIONARY['IC']['fiber'] = self.DICTIONARY['IC']['fiber'][ idx ]
            self.DICTIONARY['IC']['len']   = self.DICTIONARY['IC']['len'][ idx ]

            ICf = self.DICTIONARY['IC']['fiber']
            ICl = self.DICTIONARY['IC']['len']
            ICv = self.DICTIONARY['IC']['v']
            ICo = self.DICTIONARY['IC']['o']

            self.ICf = &ICf[0]
            self.ICl = &ICl[0]
            self.ICv = &ICv[0]
            self.ICo = &ICo[0]

            self.thisptr.setTransposeData(&self.ICv[0], &self.ICf[0], &self.ICo[0], &self.ICl[0])

    @property
    def T( self ) :
        """Transpose of the explicit matrix."""
        C = CudaLinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS )
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
            ERROR( "A.dot(): dimensions do not match" )

        # Create output array
        cdef double [::1] v_out = np.zeros( self.shape[0], dtype=np.float64 )

        # Call the cython function to read the memory pointers
        if not self.adjoint :
            # DIRECT PRODUCT A*x
            self.thisptr.dot(&v_in[0], &v_out[0])
        else :
            # INVERSE PRODUCT A'*y
            self.thisptr.Tdot(&v_in[0], &v_out[0])

        return v_out

    def destroy( self ):
        """Free all memory of the CUDA GPU"""
        self.thisptr.destroy()
