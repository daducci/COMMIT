#!python
#cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False, binding=False
cimport cython
import numpy as np
cimport numpy as np

import time
import glob
import sys
from os import makedirs, remove
from os.path import exists, join as pjoin
import nibabel
import cPickle
import commit.models
import commit.solvers
import amico.scheme
import amico.lut
import pyximport
pyximport.install( reload_support=True )


def setup( lmax = 12 ) :
    """General setup/initialization of the COMMIT framework."""
    amico.lut.precompute_rotation_matrices( lmax )


cdef class Evaluation :
    """Class to hold all the information (data and parameters) when performing an
    evaluation with the COMMIT framework.
    """
    cdef public niiDWI
    cdef public niiDWI_img
    cdef public scheme
    cdef public model
    cdef public KERNELS
    cdef public DICTIONARY
    cdef public THREADS
    cdef public A
    cdef public x
    cdef public CONFIG

    def __init__( self, study_path, subject ) :
        """Setup the data structures with default values.

        Parameters
        ----------
        study_path : string
            The path to the folder containing all the subjects from one study
        subject : string
            The path (relative to previous folder) to the subject folder
        """
        self.niiDWI     = None # set by "load_data" method
        self.scheme     = None # set by "load_data" method
        self.model      = None # set by "set_model" method
        self.KERNELS    = None # set by "load_kernels" method
        self.DICTIONARY = None # set by "load_dictionary" method
        self.THREADS    = None # set by "set_threads" method
        self.A          = None # set by "build_operator" method
        self.x          = None # set by "fit" method

        # store all the parameters of an evaluation with COMMIT
        self.CONFIG = {}
        self.set_config('study_path', study_path)
        self.set_config('subject', subject)
        self.set_config('DATA_path', pjoin( study_path, subject ))

        self.set_config('doNormalizeSignal', True)
        self.set_config('doMergeB0', False)
        self.set_config('doNormalizeKernels', True)
        self.set_config('doDemean', False)
        self.set_config('doNormalizeMaps', False)



    def set_config( self, key, value ) :
        self.CONFIG[ key ] = value

    def get_config( self, key ) :
        return self.CONFIG.get( key )


    def load_data( self, dwi_filename = 'DWI.nii', scheme_filename = 'DWI.scheme', b0_thr = 0 ) :
        """Load the diffusion signal and its corresponding acquisition scheme.

        Parameters
        ----------
        dwi_filename : string
            The file name of the DWI data, relative to the subject folder (default : 'DWI.nii')
        scheme_filename : string
            The file name of the corresponding acquisition scheme (default : 'DWI.scheme')
        b0_thr : float
            The threshold below which a b-value is considered a b0 (default : 0)
        """

        # Loading data and acquisition scheme
        tic = time.time()
        print '\n-> Loading data:'

        print '\t* DWI signal...'
        self.set_config('dwi_filename', dwi_filename)
        self.niiDWI  = nibabel.load( pjoin( self.get_config('DATA_path'), dwi_filename) )
        self.niiDWI_img = self.niiDWI.get_data().astype(np.float32)
        if self.niiDWI_img.ndim ==3 :
            self.niiDWI_img = np.expand_dims( self.niiDWI_img, axis=3 )
        hdr = self.niiDWI.header if nibabel.__version__ >= '2.0.0' else self.niiDWI.get_header()
        self.set_config('dim', self.niiDWI_img.shape[0:3])
        self.set_config('pixdim', tuple( hdr.get_zooms()[:3] ))
        print '\t\t- dim    = %d x %d x %d x %d' % self.niiDWI_img.shape
        print '\t\t- pixdim = %.3f x %.3f x %.3f' % self.get_config('pixdim')

        print '\t* Acquisition scheme...'
        self.set_config('scheme_filename', scheme_filename)
        self.set_config('b0_thr', b0_thr)
        self.scheme = amico.scheme.Scheme( pjoin( self.get_config('DATA_path'), scheme_filename), b0_thr )
        print '\t\t- %d samples, %d shells' % ( self.scheme.nS, len(self.scheme.shells) )
        print '\t\t- %d @ b=0' % ( self.scheme.b0_count ),
        for i in xrange(len(self.scheme.shells)) :
            print ', %d @ b=%.1f' % ( len(self.scheme.shells[i]['idx']), self.scheme.shells[i]['b'] ),
        print

        if self.scheme.nS != self.niiDWI_img.shape[3] :
            raise ValueError( 'Scheme does not match with DWI data' )

        if self.scheme.dwi_count == 0 :
            raise ValueError( 'There are no DWI volumes in the data' )

        print '   [ %.1f seconds ]' % ( time.time() - tic )

        # Preprocessing
        tic = time.time()
        print '\n-> Preprocessing:'

        if self.get_config('doNormalizeSignal') :
            if self.scheme.b0_count > 0 :
                print '\t* Normalizing to b0...',
                sys.stdout.flush()
                mean = np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 )
                idx = mean <= 0
                mean[ idx ] = 1
                mean = 1 / mean
                mean[ idx ] = 0
                for i in xrange(self.scheme.nS) :
                    self.niiDWI_img[:,:,:,i] *= mean
            else :
                print '\t* There are no b0 volume(s) for normalization...',
            print '[ min=%.2f,  mean=%.2f, max=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.mean(), self.niiDWI_img.max() )

        if self.scheme.b0_count > 1 :
            if self.get_config('doMergeB0') :
                print '\t* Merging multiple b0 volume(s)...',
                mean = np.expand_dims( np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 ), axis=3 )
                self.niiDWI_img = np.concatenate( (mean, self.niiDWI_img[:,:,:,self.scheme.dwi_idx]), axis=3 )
            else :
                print '\t* Keeping all b0 volume(s)...',
            print '[ %d x %d x %d x %d ]' % self.niiDWI_img.shape

        if self.get_config('doDemean') :
            print '\t* Demeaning signal...',
            sys.stdout.flush()
            mean = np.repeat( np.expand_dims(np.mean(self.niiDWI_img,axis=3),axis=3), self.niiDWI_img.shape[3], axis=3 )
            self.niiDWI_img = self.niiDWI_img - mean
            print '[ min=%.2f,  mean=%.2f, max=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.mean(), self.niiDWI_img.max() )

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def set_model( self, model_name ) :
        """Set the model to use to describe the signal contributions in each voxel.

        Parameters
        ----------
        model_name : string
            The name of the model (must match a class name in "commit.models" module)
        """
        # Call the specific model constructor
        if hasattr(commit.models, model_name ) :
            self.model = getattr(commit.models,model_name)()
        else :
            raise ValueError( 'Model "%s" not recognized' % model_name )

        self.set_config('ATOMS_path', pjoin( self.get_config('study_path'), 'kernels', self.model.id ))


    def generate_kernels( self, regenerate = False, lmax = 12 ) :
        """Generate the high-resolution response functions for each compartment.
        Dispatch to the proper function, depending on the model.

        Parameters
        ----------
        regenerate : boolean
            Regenerate kernels if they already exist (default : False)
        lmax : int
            Maximum SH order to use for the rotation procedure (default : 12)
        """
        if self.scheme is None :
            raise RuntimeError( 'Scheme not loaded; call "load_data()" first.' )
        if self.model is None :
            raise RuntimeError( 'Model not set; call "set_model()" method first.' )

        # store some values for later use
        self.set_config('lmax', lmax)
        self.model.scheme = self.scheme

        print '\n-> Simulating with "%s" model:' % self.model.name

        # check if kernels were already generated
        tmp = glob.glob( pjoin(self.get_config('ATOMS_path'),'A_*.npy') )
        if len(tmp)>0 and not regenerate :
            print '   [ Kernels already computed. Call "generate_kernels( regenerate=True )" to force regeneration. ]'
            return

        # create folder or delete existing files (if any)
        if not exists( self.get_config('ATOMS_path') ) :
            makedirs( self.get_config('ATOMS_path') )
        else :
            for f in glob.glob( pjoin(self.get_config('ATOMS_path'),'*') ) :
                remove( f )

        # auxiliary data structures
        aux = amico.lut.load_precomputed_rotation_matrices( lmax )
        idx_IN, idx_OUT = amico.lut.aux_structures_generate( self.scheme, lmax )

        # Dispatch to the right handler for each model
        tic = time.time()
        self.model.generate( self.get_config('ATOMS_path'), aux, idx_IN, idx_OUT )
        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def load_kernels( self ) :
        """Load rotated kernels and project to the specific gradient scheme of this subject.
        Dispatch to the proper function, depending on the model.
        """
        if self.model is None :
            raise RuntimeError( 'Model not set; call "set_model()" method first.' )
        if self.scheme is None :
            raise RuntimeError( 'Scheme not loaded; call "load_data()" first.' )

        tic = time.time()
        print '\n-> Resampling LUT for subject "%s":' % self.get_config('subject')

        # auxiliary data structures
        idx_OUT, Ylm_OUT = amico.lut.aux_structures_resample( self.scheme, self.get_config('lmax') )

        # Dispatch to the right handler for each model
        if self.get_config('doMergeB0') :
            print '\t* Merging multiple b0 volume(s)...',
        else :
            print '\t* Keeping all b0 volume(s)...',
        self.KERNELS = self.model.resample( self.get_config('ATOMS_path'), idx_OUT, Ylm_OUT, self.get_config('doMergeB0') )
        nIC  = self.KERNELS['wmr'].shape[0]
        nEC  = self.KERNELS['wmh'].shape[0]
        nISO = self.KERNELS['iso'].shape[0]
        print '[ OK ]'


        # ensure contiguous arrays for C part
        self.KERNELS['wmr'] = np.ascontiguousarray( self.KERNELS['wmr'] )
        self.KERNELS['wmh'] = np.ascontiguousarray( self.KERNELS['wmh'] )
        self.KERNELS['iso'] = np.ascontiguousarray( self.KERNELS['iso'] )

        # De-mean kernels
        if self.get_config('doDemean') :
            print '\t* Demeaning signal...',
            for j in xrange(181) :
                for k in xrange(181) :
                    for i in xrange(nIC) :
                        self.KERNELS['wmr'][i,j,k,:] -= self.KERNELS['wmr'][i,j,k,:].mean()
                    for i in xrange(nEC) :
                        self.KERNELS['wmh'][i,j,k,:] -= self.KERNELS['wmh'][i,j,k,:].mean()
            for i in xrange(nISO) :
                self.KERNELS['iso'][i] -= self.KERNELS['iso'][i].mean()
            print '[ OK ]'

        # Normalize atoms
        if self.get_config('doNormalizeKernels') :
            print '\t* Normalizing...',

            self.KERNELS['wmr_norm'] = np.zeros( nIC )
            for i in xrange(nIC) :
                self.KERNELS['wmr_norm'][i] = np.linalg.norm( self.KERNELS['wmr'][i,0,0,:] )
                for j in xrange(181) :
                    for k in xrange(181) :
                        self.KERNELS['wmr'][i,j,k,:] /= self.KERNELS['wmr_norm'][i]

            self.KERNELS['wmh_norm'] = np.zeros( nEC )
            for i in xrange(nEC) :
                self.KERNELS['wmh_norm'][i] = np.linalg.norm( self.KERNELS['wmh'][i,0,0,:] )
                for j in xrange(181) :
                    for k in xrange(181) :
                        self.KERNELS['wmh'][i,j,k,:] /= self.KERNELS['wmh_norm'][i]

            self.KERNELS['iso_norm'] = np.zeros( nISO )
            for i in xrange(nISO) :
                self.KERNELS['iso_norm'][i] = np.linalg.norm( self.KERNELS['iso'][i,:] )
                self.KERNELS['iso'][i,:] /= self.KERNELS['iso_norm'][i]

            print '[ OK ]'

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    cpdef load_dictionary( self, path, use_mask = False ) :
        """Load the sparse structure previously created with "trk2dictionary" script.

        Parameters
        ----------
        path : string
            Folder containing the output of the trk2dictionary script (relative to subject path)
        use_mask : boolean
            If False (default) the optimization will be conducted only on the voxels actually
            traversed by tracts. If True, the mask specified in trk2dictionary
            (i.e. "filename_mask" paramater) will be used instead.
            NB: if no mask was specified in trk2dictionary, the "tdi" and
            "mask" masks are equivalent and this parameter is not influent.
        """
        if self.niiDWI is None :
            raise RuntimeError( 'Data not loaded; call "load_data()" first.' )

        tic = time.time()
        print '\n-> Loading the dictionary:'
        self.DICTIONARY = {}
        self.set_config('TRACKING_path', pjoin(self.get_config('DATA_path'),path))

        # load mask
        self.set_config('dictionary_mask', 'mask' if use_mask else 'tdi' )
        mask_filename = pjoin(self.get_config('TRACKING_path'),'dictionary_%s.nii'%self.get_config('dictionary_mask'))
        if not exists( mask_filename ) :
            mask_filename += '.gz'
            if not exists( mask_filename ) :
                raise RuntimeError( 'Dictionary not found. Execute ''trk2dictionary'' script first.' );
        niiMASK = nibabel.load( mask_filename )
        self.DICTIONARY['MASK'] = (niiMASK.get_data() > 0).astype(np.uint8)

        # segments from the tracts
        # ------------------------
        print '\t* segments from the tracts...',
        sys.stdout.flush()

        self.DICTIONARY['TRK'] = {}
        self.DICTIONARY['TRK']['norm'] = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_TRK_norm.dict'), dtype=np.float32 )
        self.DICTIONARY['TRK']['len']  = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_TRK_len.dict'), dtype=np.float32 )

        self.DICTIONARY['IC'] = {}
        self.DICTIONARY['IC']['fiber'] = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_IC_f.dict'), dtype=np.uint32 )
        self.DICTIONARY['IC']['v']     = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_IC_v.dict'), dtype=np.uint32 )
        self.DICTIONARY['IC']['o']     = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_IC_o.dict'), dtype=np.uint16 )
        self.DICTIONARY['IC']['len']   = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_IC_len.dict'), dtype=np.float32 )
        self.DICTIONARY['IC']['n']     = self.DICTIONARY['IC']['fiber'].size
        self.DICTIONARY['IC']['nF']    = self.DICTIONARY['TRK']['norm'].size

        # reorder the segments based on the "v" field
        idx = np.argsort( self.DICTIONARY['IC']['v'], kind='mergesort' )
        self.DICTIONARY['IC']['v']     = self.DICTIONARY['IC']['v'][ idx ]
        self.DICTIONARY['IC']['o']     = self.DICTIONARY['IC']['o'][ idx ]
        self.DICTIONARY['IC']['fiber'] = self.DICTIONARY['IC']['fiber'][ idx ]
        self.DICTIONARY['IC']['len']   = self.DICTIONARY['IC']['len'][ idx ]
        del idx

        # divide the length of each segment by the fiber length so that all the columns of the libear operator will have same length
        # NB: it works in conjunction with the normalization of the kernels
        cdef :
            np.float32_t [:] sl = self.DICTIONARY['IC']['len']
            np.float32_t [:] tl = self.DICTIONARY['TRK']['norm']
            np.uint32_t  [:] f  = self.DICTIONARY['IC']['fiber']
            int s
        if self.get_config('doNormalizeKernels') :
            for s in xrange(self.DICTIONARY['IC']['n']) :
                sl[s] /= tl[ f[s] ]

        print '[ %d fibers and %d segments ]' % ( self.DICTIONARY['IC']['nF'], self.DICTIONARY['IC']['n'] )

        # segments from the peaks
        # -----------------------
        print '\t* segments from the peaks...',
        sys.stdout.flush()

        self.DICTIONARY['EC'] = {}
        self.DICTIONARY['EC']['v']  = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_EC_v.dict'), dtype=np.uint32 )
        self.DICTIONARY['EC']['o']  = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_EC_o.dict'), dtype=np.uint16 )
        self.DICTIONARY['EC']['nE'] = self.DICTIONARY['EC']['v'].size

        # reorder the segments based on the "v" field
        idx = np.argsort( self.DICTIONARY['EC']['v'], kind='mergesort' )
        self.DICTIONARY['EC']['v'] = self.DICTIONARY['EC']['v'][ idx ]
        self.DICTIONARY['EC']['o'] = self.DICTIONARY['EC']['o'][ idx ]
        del idx

        print ' [ %d segments ]' % self.DICTIONARY['EC']['nE']

        # isotropic compartments
        # ----------------------
        print '\t* isotropic contributions...',
        sys.stdout.flush()

        self.DICTIONARY['ISO'] = {}

        self.DICTIONARY['nV'] = self.DICTIONARY['MASK'].sum()

        vx, vy, vz = ( self.DICTIONARY['MASK'] > 0 ).nonzero() # [TODO] find a way to avoid using int64 (not necessary and waste of memory)
        vx = vx.astype(np.int32)
        vy = vy.astype(np.int32)
        vz = vz.astype(np.int32)
        self.DICTIONARY['ISO']['v'] = vx + self.get_config('dim')[0] * ( vy + self.get_config('dim')[1] * vz )
        del vx, vy, vz

        # reorder the segments based on the "v" field
        idx = np.argsort( self.DICTIONARY['ISO']['v'], kind='mergesort' )
        self.DICTIONARY['ISO']['v'] = self.DICTIONARY['ISO']['v'][ idx ]
        del idx

        print ' [ %d voxels ]' % self.DICTIONARY['nV']

        # post-processing
        # ---------------
        print '\t* post-processing...',
        sys.stdout.flush()

        # get the indices to extract the VOI as in MATLAB (in place of DICTIONARY.MASKidx)
        idx = self.DICTIONARY['MASK'].ravel(order='F').nonzero()[0]
        self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] = np.unravel_index( idx, self.DICTIONARY['MASK'].shape, order='F' )

        lut = np.zeros( self.get_config('dim'), dtype=np.uint32 ).ravel()
        for i in xrange(idx.size) :
            lut[ idx[i] ] = i
        self.DICTIONARY['IC'][ 'v'] = lut[ self.DICTIONARY['IC'][ 'v'] ]
        self.DICTIONARY['EC'][ 'v'] = lut[ self.DICTIONARY['EC'][ 'v'] ]
        self.DICTIONARY['ISO']['v'] = lut[ self.DICTIONARY['ISO']['v'] ]

        print '         [ OK ]'

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def set_threads( self, n = None ) :
        """Set the number of threads to use for the matrix-vector operations with A and A'.

        Parameters
        ----------
        n : integer
            Number of threads to use (default : number of CPUs in the system)
        """
        if n is None :
            # Set to the number of CPUs in the system
            try :
                import multiprocessing
                n = multiprocessing.cpu_count()
            except :
                n = 1

        if n < 1 or n > 255 :
            raise RuntimeError( 'Number of threads must be between 1 and 255' )
        if self.DICTIONARY is None :
            raise RuntimeError( 'Dictionary not loaded; call "load_dictionary()" first.' )
        if self.KERNELS is None :
            raise RuntimeError( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first.' )

        self.THREADS = {}
        self.THREADS['n'] = n

        cdef :
            long [:] C
            long t, tot, i1, i2, N, c
            int i

        tic = time.time()
        print '\n-> Distributing workload to different threads:'
        print '\t* number of threads : %d' % n

        # Distribute load for the computation of A*x product
        print '\t* A operator...',
        sys.stdout.flush()

        if self.DICTIONARY['IC']['n'] > 0 :
            self.THREADS['IC'] = np.zeros( n+1, dtype=np.uint32 )
            if n > 1 :
                N = np.floor( self.DICTIONARY['IC']['n']/n )
                t = 1
                tot = 0
                C = np.bincount( self.DICTIONARY['IC']['v'] )
                for c in C :
                    tot += c
                    if tot >= N :
                        self.THREADS['IC'][t] = self.THREADS['IC'][t-1] + tot
                        t += 1
                        tot = 0
            self.THREADS['IC'][n] = self.DICTIONARY['IC']['n']

            # check if some threads are not assigned any segment
            if np.count_nonzero( np.diff( self.THREADS['IC'].astype(np.int32) ) <= 0 ) :
                self.THREADS = None
                raise RuntimeError( 'Too many threads for the IC compartments to evaluate; try decreasing the number.' )
        else :
            self.THREADS['IC'] = None

        if self.DICTIONARY['EC']['nE'] > 0 :
            self.THREADS['EC'] = np.zeros( n+1, dtype=np.uint32 )
            for i in xrange(n) :
                self.THREADS['EC'][i] = np.searchsorted( self.DICTIONARY['EC']['v'], self.DICTIONARY['IC']['v'][ self.THREADS['IC'][i] ] )
            self.THREADS['EC'][n] = self.DICTIONARY['EC']['nE']

            # check if some threads are not assigned any segment
            if np.count_nonzero( np.diff( self.THREADS['EC'].astype(np.int32) ) <= 0 ) :
                self.THREADS = None
                raise RuntimeError( 'Too many threads for the EC compartments to evaluate; try decreasing the number.' )
        else :
            self.THREADS['EC'] = None

        if self.DICTIONARY['nV'] > 0 :
            self.THREADS['ISO'] = np.zeros( n+1, dtype=np.uint32 )
            for i in xrange(n) :
                self.THREADS['ISO'][i] = np.searchsorted( self.DICTIONARY['ISO']['v'], self.DICTIONARY['IC']['v'][ self.THREADS['IC'][i] ] )
            self.THREADS['ISO'][n] = self.DICTIONARY['nV']

            # check if some threads are not assigned any segment
            if np.count_nonzero( np.diff( self.THREADS['ISO'].astype(np.int32) ) <= 0 ) :
                self.THREADS = None
                raise RuntimeError( 'Too many threads for the ISO compartments to evaluate; try decreasing the number.' )
        else :
            self.THREADS['ISO'] = None

        print ' [ OK ]'

        # Distribute load for the computation of At*y product
        print '\t* A\' operator...',
        sys.stdout.flush()

        if self.DICTIONARY['IC']['n'] > 0 :
            self.THREADS['ICt'] = np.full( self.DICTIONARY['IC']['n'], n-1, dtype=np.uint8 )
            if n > 1 :
                idx = np.argsort( self.DICTIONARY['IC']['fiber'], kind='mergesort' )
                C = np.bincount( self.DICTIONARY['IC']['fiber'] )
                t = tot = i1 = i2 = 0
                N = np.floor(self.DICTIONARY['IC']['n']/n)
                for c in C :
                    i2 += c
                    tot += c
                    if tot >= N :
                        self.THREADS['ICt'][ i1:i2 ] = t
                        t += 1
                        if t==n-1 :
                            break
                        i1 = i2
                        tot = c
                self.THREADS['ICt'][idx] = self.THREADS['ICt'].copy()

        else :
            self.THREADS['ICt'] = None

        if self.DICTIONARY['EC']['nE'] > 0 :
            self.THREADS['ECt'] = np.zeros( n+1, dtype=np.uint32 )
            N = np.floor( self.DICTIONARY['EC']['nE']/n )
            for i in xrange(1,n) :
                self.THREADS['ECt'][i] = self.THREADS['ECt'][i-1] + N
            self.THREADS['ECt'][n] = self.DICTIONARY['EC']['nE']

            # check if some threads are not assigned any segment
            if np.count_nonzero( np.diff( self.THREADS['ECt'].astype(np.int32) ) <= 0 ) :
                self.THREADS = None
                raise RuntimeError( 'Too many threads for the EC compartments to evaluate; try decreasing the number.' )
        else :
            self.THREADS['ECt'] = None

        if self.DICTIONARY['nV'] > 0 :
            self.THREADS['ISOt'] = np.zeros( n+1, dtype=np.uint32 )
            N = np.floor( self.DICTIONARY['nV']/n )
            for i in xrange(1,n) :
                self.THREADS['ISOt'][i] = self.THREADS['ISOt'][i-1] + N
            self.THREADS['ISOt'][n] = self.DICTIONARY['nV']

            # check if some threads are not assigned any segment
            if np.count_nonzero( np.diff( self.THREADS['ISOt'].astype(np.int32) ) <= 0 ) :
                self.THREADS = None
                raise RuntimeError( 'Too many threads for the ISO compartments to evaluate; try decreasing the number.' )
        else :
            self.THREADS['ISOt'] = None

        print '[ OK ]'

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def build_operator( self ) :
        """Compile/build the operator for computing the matrix-vector multiplications by A and A'
        using the informations from self.DICTIONARY, self.KERNELS and self.THREADS.
        NB: needs to call this function to update pointers to data structures in case
            the data is changed in self.DICTIONARY, self.KERNELS or self.THREADS.
        """
        if self.DICTIONARY is None :
            raise RuntimeError( 'Dictionary not loaded; call "load_dictionary()" first.' )
        if self.KERNELS is None :
            raise RuntimeError( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first.' )
        if self.THREADS is None :
            raise RuntimeError( 'Threads not set; call "set_threads()" first.' )

        tic = time.time()
        print '\n-> Building linear operator A:'

        # need to pass these parameters at runtime for compiling the C code
        from commit.operator import config
        config.nTHREADS = self.THREADS['n']
        config.model    = self.model.id
        config.nIC      = self.KERNELS['wmr'].shape[0]
        config.nEC      = self.KERNELS['wmh'].shape[0]
        config.nISO     = self.KERNELS['iso'].shape[0]
        if not 'commit.operator.operator' in sys.modules :
            import commit.operator.operator
        else :
            reload( sys.modules['commit.operator.operator'] )
        self.A = sys.modules['commit.operator.operator'].LinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS )

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def fit( self, tol_fun = 1e-3, max_iter = 100, verbose = 1, x0 = None ) :
        """Fit the model to the data.

        Parameters
        ----------
        tol_fun : float
            Tolerance on the objective function (default : 1e-3)
        max_iter : integer
            Maximum number of iterations (default : 100)
        verbose : integer
            Level of verbosity: 0=no print, 1=print progress (default : 1)
        """
        if self.niiDWI is None :
            raise RuntimeError( 'Data not loaded; call "load_data()" first.' )
        if self.DICTIONARY is None :
            raise RuntimeError( 'Dictionary not loaded; call "load_dictionary()" first.' )
        if self.KERNELS is None :
            raise RuntimeError( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first.' )
        if self.THREADS is None :
            raise RuntimeError( 'Threads not set; call "set_threads()" first.' )
        if self.A is None :
            raise RuntimeError( 'Operator not built; call "build_operator()" first.' )
        if x0 is not None :
            if x0.shape[0] != self.A.shape[1] :
                raise RuntimeError( 'x0: dimension do not match' )

        self.CONFIG['optimization'] = {}
        self.CONFIG['optimization']['tol_fun']  = tol_fun
        self.CONFIG['optimization']['max_iter'] = max_iter
        self.CONFIG['optimization']['verbose']  = verbose

        # run solver
        t = time.time()
        print '\n-> Fit model using "nnls":'
        Y = self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float64)
        self.x, OPT_det = commit.solvers.nnls( Y, self.A, tol_fun=tol_fun, max_iter=max_iter, verbose=verbose, x0=x0 )
        self.CONFIG['optimization']['fit_time'] = round(time.time()-t, 3)
        self.CONFIG['optimization']['||Ax-y||'] = OPT_det['||Ax-y||']
        self.CONFIG['optimization']['Cost function'] = OPT_det['Cost function']
        self.CONFIG['optimization']['Error abs'] = OPT_det['Abs error']
        self.CONFIG['optimization']['Error rel'] = OPT_det['Rel error']
        self.CONFIG['optimization']['x abs'] = OPT_det['Abs x']
        self.CONFIG['optimization']['x rel'] = OPT_det['Rel x']
        self.CONFIG['optimization']['iteration'] = OPT_det['iteration']
        print '   [ %s ]' % ( time.strftime("%Hh %Mm %Ss", time.gmtime(self.CONFIG['optimization']['fit_time']) ) )


    def save_results( self, path_suffix = None ) :
        """Save the output (coefficients, errors, maps etc).

        Parameters
        ----------
        path_suffix : string
            Text to be appended to "Results" to create the output path (default : None)
        """
        if self.x is None :
            raise RuntimeError( 'Model not fitted to the data; call "fit()" first.' )

        RESULTS_path = 'Results_' + self.model.id
        if path_suffix :
            self.set_config('path_suffix', path_suffix)
            RESULTS_path = RESULTS_path + path_suffix

        print '\n-> Saving results to "%s/*":' % RESULTS_path
        tic = time.time()

        # create folder or delete existing files (if any)
        RESULTS_path = pjoin( self.get_config('TRACKING_path'), RESULTS_path )
        if not exists( RESULTS_path ) :
            makedirs( RESULTS_path )
        else :
            for f in glob.glob( pjoin(RESULTS_path,'*') ) :
                remove( f )
        self.set_config('RESULTS_path', RESULTS_path)

        # Configuration and results
        print '\t* configuration and results...',
        sys.stdout.flush()
        nF = self.DICTIONARY['IC']['nF']
        nE = self.DICTIONARY['EC']['nE']
        nV = self.DICTIONARY['nV']
        norm_fib = np.ones( nF )
        # x is the x of the original problem
        # self.x is the x preconditioned
        if self.get_config('doNormalizeKernels') :
            # renormalize the coefficients
            norm1 = np.repeat(self.KERNELS['wmr_norm'],nF)
            norm2 = np.repeat(self.KERNELS['wmh_norm'],nE)
            norm3 = np.repeat(self.KERNELS['iso_norm'],nV)
            norm_fib = np.kron(np.ones(self.KERNELS['wmr'].shape[0]), self.DICTIONARY['TRK']['norm'])
            x = self.x / np.hstack( (norm1*norm_fib,norm2,norm3) )
        else :
            x = self.x

        np.savetxt( pjoin(RESULTS_path,'weights.txt'), x, fmt='%.4f' )
        print '[ OK ]'

        # Map of wovelwise errors
        print '\t* fitting errors:'

        not_NaN = np.ones( self.get_config('dim'), dtype=np.float32 ) * 1e-16 # avoid division by 0

        niiMAP_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        affine = self.niiDWI.affine if nibabel.__version__ >= '2.0.0' else self.niiDWI.get_affine()
        niiMAP     = nibabel.Nifti1Image( niiMAP_img, affine )
        niiMAP_hdr = niiMAP.header if nibabel.__version__ >= '2.0.0' else niiMAP.get_header()

        y_mea = np.reshape( self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32), (nV,-1) )
        y_est = np.reshape( self.A.dot(self.x), (nV,-1) ).astype(np.float32)

        print '\t\t- RMSE...',
        sys.stdout.flush()
        tmp = np.sqrt( np.mean((y_mea-y_est)**2,axis=1) )
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP_hdr['cal_min'] = 0
        niiMAP_hdr['cal_max'] = tmp.max()
        nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_RMSE.nii.gz') )
        print ' [ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() )
        self.CONFIG['map'] = {}
        self.CONFIG['map']['RMSE mean'] = round(tmp.mean(), 3)
        self.CONFIG['map']['RMSE std'] = round(tmp.std(), 3)

        print '\t\t- NRMSE...',
        sys.stdout.flush()
        tmp = np.sum(y_mea**2,axis=1)
        idx = np.where( tmp < 1E-12 )
        tmp[ idx ] = 1
        tmp = np.sqrt( np.sum((y_mea-y_est)**2,axis=1) / tmp )
        tmp[ idx ] = 0
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP_hdr['cal_min'] = 0
        niiMAP_hdr['cal_max'] = 1
        nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_NRMSE.nii.gz') )
        print '[ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() )
        self.CONFIG['map']['NRMSE mean'] = round(tmp.mean(), 3)
        self.CONFIG['map']['NRMSE std'] = round(tmp.std(), 3)

        # Map of compartment contributions
        print '\t* voxelwise contributions:'

        print '\t\t- intra-axonal',
        sys.stdout.flush()
        niiIC_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        if len(self.KERNELS['wmr']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0]
            tmp = ( x[:offset].reshape( (-1,nF) ) * norm_fib.reshape( (-1,nF) ) ).sum( axis=0 )
            xv = np.bincount( self.DICTIONARY['IC']['v'], minlength=nV,
                weights=tmp[ self.DICTIONARY['IC']['fiber'] ] * self.DICTIONARY['IC']['len']
            ).astype(np.float32)
            niiIC_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        print '[ OK ]'

        print '\t\t- extra-axonal',
        sys.stdout.flush()
        niiEC_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        if len(self.KERNELS['wmh']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0]
            tmp = x[offset:offset+nE*len(self.KERNELS['wmh'])].reshape( (-1,nE) ).sum( axis=0 )
            xv = np.bincount( self.DICTIONARY['EC']['v'], weights=tmp, minlength=nV ).astype(np.float32)
            niiEC_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        print '[ OK ]'

        print '\t\t- isotropic',
        sys.stdout.flush()
        niiISO_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        if len(self.KERNELS['iso']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0] + nE * self.KERNELS['wmh'].shape[0]
            xv = x[offset:].reshape( (-1,nV) ).sum( axis=0 )
            niiISO_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        print '   [ OK ]'

        if self.get_config('doNormalizeMaps') :
                niiIC = nibabel.Nifti1Image( niiIC_img / ( niiIC_img + niiEC_img + niiISO_img + not_NaN), affine )
                niiEC = nibabel.Nifti1Image( niiEC_img / ( niiIC_img + niiEC_img + niiISO_img + not_NaN), affine )
                niiISO = nibabel.Nifti1Image( niiISO_img / ( niiIC_img + niiEC_img + niiISO_img + not_NaN), affine )
        else:
                niiIC = nibabel.Nifti1Image( niiIC_img, affine )
                niiEC = nibabel.Nifti1Image( niiEC_img, affine )
                niiISO = nibabel.Nifti1Image( niiISO_img, affine )

        nibabel.save( niiIC , pjoin(RESULTS_path,'compartment_IC.nii.gz') )
        nibabel.save( niiEC , pjoin(RESULTS_path,'compartment_EC.nii.gz') )
        nibabel.save( niiISO , pjoin(RESULTS_path,'compartment_ISO.nii.gz') )

        with open( pjoin(RESULTS_path,'results.pickle'), 'wb+' ) as fid :
            cPickle.dump( [self.CONFIG, self.x, x], fid, protocol=2 )

        print '   [ %.1f seconds ]' % ( time.time() - tic )
