#!python
#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False, binding=False
from __future__ import print_function
cimport cython
import numpy as np
cimport numpy as np

import time
import glob
import sys
from os import makedirs, remove, getcwd, listdir
from os.path import exists, join as pjoin, isfile, isdir
import nibabel
import pickle
import commit.models
import commit.solvers
import amico.scheme
import amico.lut
import pyximport
from pkg_resources import get_distribution

from amico.util import LOG, NOTE, WARNING, ERROR


def setup( lmax=12, ndirs=32761 ) :
    """General setup/initialization of the COMMIT framework.

    Parameters
    ----------
    lmax : int
        Maximum SH order to use for the rotation phase (default : 12)
    ndirs : int
        Number of directions on the half of the sphere representing the possible orientations of the response functions (default : 32761)
    """

    if not amico.lut.is_valid(ndirs):
        ERROR( 'Unsupported value for ndirs.\nNote: Supported values for ndirs are [1, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 32761 (default)]' )

    amico.lut.precompute_rotation_matrices( lmax, ndirs )


def load_dictionary_info( filename ):
    """Function to load dictionary info file

    Parameters
    ----------
    filename : string
        This value is always COMMIT_PATH + dictionary_info.pickle
    """
    if not isfile( filename ):
        ERROR( 'Dictionary is outdated or not found. Execute "trk2dictionary" script first' )
    with open( filename, 'rb' ) as dictionary_info_file:
        if sys.version_info.major == 3:
            aux = pickle.load( dictionary_info_file, fix_imports=True, encoding='bytes' )
            # Pickle files written by Python 2 are loaded with byte
            # keys, whereas those written by Python 3 are loaded with
            # str keys, even when both are written using protocol=2
            result_aux = {(k.decode() if hasattr(k,"decode") else k): v for k, v in aux.items()}
            return result_aux
        else:
            return pickle.load( dictionary_info_file )


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
    cdef public confidence_map_img
    
    def __init__( self, study_path, subject ) :
        """Setup the data structures with default values.

        Parameters
        ----------
        study_path : string
            The path to the folder containing all the subjects from one study
        subject : string
            The path (relative to previous folder) to the subject folder
        """
        self.niiDWI             = None # set by "load_data" method
        self.scheme             = None # set by "load_data" method
        self.model              = None # set by "set_model" method
        self.KERNELS            = None # set by "load_kernels" method
        self.DICTIONARY         = None # set by "load_dictionary" method
        self.THREADS            = None # set by "set_threads" method
        self.A                  = None # set by "build_operator" method
        self.x                  = None # set by "fit" method
        self.confidence_map_img = None # set by "fit" method

        # store all the parameters of an evaluation with COMMIT
        self.CONFIG = {}
        self.set_config('version', get_distribution('dmri-commit').version)
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


    def load_data( self, dwi_filename='DWI.nii', scheme_filename='DWI.scheme', b0_thr=0, b0_min_signal=0 ) :
        """Load the diffusion signal and its corresponding acquisition scheme.

        Parameters
        ----------
        dwi_filename : string
            The file name of the DWI data, relative to the subject folder (default : 'DWI.nii')
        scheme_filename : string
            The file name of the corresponding acquisition scheme (default : 'DWI.scheme')
        b0_thr : float
            The threshold below which a b-value is considered a b0 (default : 0)
        b0_min_signal : float
            Crop to zero the signal in voxels where the b0 <= b0_min_signal * mean(b0[b0>0]) (default : 0)
        """

        # Loading data and acquisition scheme
        tic = time.time()
        LOG( '\n-> Loading data:' )

        print( '\t* DWI signal:' )
        self.set_config('dwi_filename', dwi_filename)
        self.niiDWI  = nibabel.load( pjoin( self.get_config('DATA_path'), dwi_filename) )
        self.niiDWI_img = np.asanyarray( self.niiDWI.dataobj ).astype(np.float32)
        if self.niiDWI_img.ndim ==3 :
            self.niiDWI_img = np.expand_dims( self.niiDWI_img, axis=3 )
        hdr = self.niiDWI.header if nibabel.__version__ >= '2.0.0' else self.niiDWI.get_header()
        self.set_config('dim', self.niiDWI_img.shape[0:3])
        self.set_config('pixdim', tuple( hdr.get_zooms()[:3] ))
        print( '\t\t- dim    : %d x %d x %d x %d' % self.niiDWI_img.shape )
        print( '\t\t- pixdim : %.3f x %.3f x %.3f' % self.get_config('pixdim') )
        print( '\t\t- values : min=%.2f, max=%.2f, mean=%.2f' % ( self.niiDWI_img.min(), self.niiDWI_img.max(), self.niiDWI_img.mean() ) )

        print( '\t* Acquisition scheme:' )
        self.set_config('scheme_filename', scheme_filename)
        self.set_config('b0_thr', b0_thr)
        self.scheme = amico.scheme.Scheme( pjoin( self.get_config('DATA_path'), scheme_filename), b0_thr )
        print( '\t\t- %d samples, %d shells' % ( self.scheme.nS, len(self.scheme.shells) ) )
        print( '\t\t- %d @ b=0' % ( self.scheme.b0_count ), end='' )
        for i in xrange(len(self.scheme.shells)) :
            print( ', %d @ b=%.1f' % ( len(self.scheme.shells[i]['idx']), self.scheme.shells[i]['b'] ), end='' )
        print()

        if self.scheme.nS != self.niiDWI_img.shape[3] :
            ERROR( 'Scheme does not match with DWI data' )

        if self.scheme.dwi_count == 0 :
            ERROR( 'There are no DWI volumes in the data' )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )

        # Preprocessing
        tic = time.time()
        LOG( '\n-> Preprocessing:' )

        if self.get_config('doNormalizeSignal') :
            if self.scheme.b0_count > 0 :
                print( '\t* Normalizing to b0... ', end='' )
                sys.stdout.flush()
                b0 = np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 )
                idx = b0 <= b0_min_signal * b0[b0>0].mean()
                b0[ idx ] = 1
                b0 = 1.0 / b0
                b0[ idx ] = 0
                for i in xrange(self.scheme.nS) :
                    self.niiDWI_img[:,:,:,i] *= b0
                print( '[ min=%.2f, max=%.2f, mean=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.max(), self.niiDWI_img.mean() ) )
                del idx, b0
            else :
                WARNING( 'There are no b0 volumes for normalization' )

        if self.scheme.b0_count > 1 :
            if self.get_config('doMergeB0') :
                print( '\t* Merging multiple b0 volume(s)... ', end='' )
                mean = np.expand_dims( np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 ), axis=3 )
                self.niiDWI_img = np.concatenate( (mean, self.niiDWI_img[:,:,:,self.scheme.dwi_idx]), axis=3 )
                del mean
            else :
                print( '\t* Keeping all b0 volume(s)... ', end='' )
            print( '[ %d x %d x %d x %d ]' % self.niiDWI_img.shape )

        if self.get_config('doDemean') :
            print( '\t* Demeaning signal... ', end='' )
            sys.stdout.flush()
            mean = np.repeat( np.expand_dims(np.mean(self.niiDWI_img,axis=3),axis=3), self.niiDWI_img.shape[3], axis=3 )
            self.niiDWI_img = self.niiDWI_img - mean
            print( '[ min=%.2f, max=%.2f, mean=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.max(), self.niiDWI_img.mean() ) )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )


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
            ERROR( 'Model "%s" not recognized' % model_name )

        self.set_config('ATOMS_path', pjoin( self.get_config('study_path'), 'kernels', self.model.id ))


    def generate_kernels( self, regenerate=False, lmax=12, ndirs=32761 ) :
        """Generate the high-resolution response functions for each compartment.
        Dispatch to the proper function, depending on the model.

        Parameters
        ----------
        regenerate : boolean
            Regenerate kernels if they already exist (default : False)
        lmax : int
            Maximum SH order to use for the rotation procedure (default : 12)
        ndirs : int
            Number of directions on the half of the sphere representing the possible orientations of the response functions (default : 32761)
        """
        if not amico.lut.is_valid(ndirs):
            ERROR( 'Unsupported value for ndirs.\nNote: Supported values for ndirs are [1, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 32761 (default)]' )
        if self.scheme is None :
            ERROR( 'Scheme not loaded; call "load_data()" first' )
        if self.model is None :
            ERROR( 'Model not set; call "set_model()" method first' )

        # store some values for later use
        self.set_config('lmax', lmax)
        self.set_config('ndirs', ndirs)
        self.set_config('model', self.model.get_params())
        self.model.scheme = self.scheme

        LOG( '\n-> Simulating with "%s" model:' % self.model.name )

        # check if kernels were already generated
        tmp = glob.glob( pjoin(self.get_config('ATOMS_path'),'A_*.npy') )
        if len(tmp)>0 and not regenerate :
            LOG( '   [ Kernels already computed. Use option "regenerate=True" to force regeneration ]' )
            return

        # create folder or delete existing files (if any)
        if not exists( self.get_config('ATOMS_path') ) :
            makedirs( self.get_config('ATOMS_path') )
        else :
            for f in glob.glob( pjoin(self.get_config('ATOMS_path'),'*') ) :
                remove( f )

        # auxiliary data structures
        aux = amico.lut.load_precomputed_rotation_matrices( lmax, ndirs )
        idx_IN, idx_OUT = amico.lut.aux_structures_generate( self.scheme, lmax )

        # Dispatch to the right handler for each model
        tic = time.time()
        self.model.generate( self.get_config('ATOMS_path'), aux, idx_IN, idx_OUT, ndirs )
        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )


    def load_kernels( self ) :
        """Load rotated kernels and project to the specific gradient scheme of this subject.
        Dispatch to the proper function, depending on the model.
        """
        if self.model is None :
            ERROR( 'Model not set; call "set_model()" method first' )
        if self.scheme is None :
            ERROR( 'Scheme not loaded; call "load_data()" first' )

        tic = time.time()
        LOG( '\n-> Resampling LUT for subject "%s":' % self.get_config('subject') )

        # auxiliary data structures
        idx_OUT, Ylm_OUT = amico.lut.aux_structures_resample( self.scheme, self.get_config('lmax') )

        # Dispatch to the right handler for each model
        if self.get_config('doMergeB0') :
            print( '\t* Merging multiple b0 volume(s)...' )
        else :
            print( '\t* Keeping all b0 volume(s)...' )
        self.KERNELS = self.model.resample( self.get_config('ATOMS_path'), idx_OUT, Ylm_OUT, self.get_config('doMergeB0'), self.get_config('ndirs') )
        nIC  = self.KERNELS['wmr'].shape[0]
        nEC  = self.KERNELS['wmh'].shape[0]
        nISO = self.KERNELS['iso'].shape[0]
        print( '\t  [ OK ]' )

        # ensure contiguous arrays for C part
        self.KERNELS['wmr'] = np.ascontiguousarray( self.KERNELS['wmr'] )
        self.KERNELS['wmh'] = np.ascontiguousarray( self.KERNELS['wmh'] )
        self.KERNELS['iso'] = np.ascontiguousarray( self.KERNELS['iso'] )

        # De-mean kernels
        if self.get_config('doDemean') :
            print( '\t* Demeaning signal...', end='' )
            for j in xrange(self.get_config('ndirs')) :
                for i in xrange(nIC) :
                    self.KERNELS['wmr'][i,j,:] -= self.KERNELS['wmr'][i,j,:].mean()
                for i in xrange(nEC) :
                    self.KERNELS['wmh'][i,j,:] -= self.KERNELS['wmh'][i,j,:].mean()
            for i in xrange(nISO) :
                self.KERNELS['iso'][i] -= self.KERNELS['iso'][i].mean()
            print( '[ OK ]' )

        # Normalize atoms
        if self.get_config('doNormalizeKernels') :
            print( '\t* Normalizing... ', end='' )

            self.KERNELS['wmr_norm'] = np.zeros( nIC )
            for i in xrange(nIC) :
                self.KERNELS['wmr_norm'][i] = np.linalg.norm( self.KERNELS['wmr'][i,0,:] )
                for j in xrange(self.get_config('ndirs')) :
                    self.KERNELS['wmr'][i,j,:] /= self.KERNELS['wmr_norm'][i]

            self.KERNELS['wmh_norm'] = np.zeros( nEC )
            for i in xrange(nEC) :
                self.KERNELS['wmh_norm'][i] = np.linalg.norm( self.KERNELS['wmh'][i,0,:] )
                for j in xrange(self.get_config('ndirs')) :
                    self.KERNELS['wmh'][i,j,:] /= self.KERNELS['wmh_norm'][i]

            self.KERNELS['iso_norm'] = np.zeros( nISO )
            for i in xrange(nISO) :
                self.KERNELS['iso_norm'][i] = np.linalg.norm( self.KERNELS['iso'][i,:] )
                self.KERNELS['iso'][i,:] /= self.KERNELS['iso_norm'][i]

            print( '[ OK ]' )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )


    cpdef load_dictionary( self, path, use_all_voxels_in_mask=False ) :
        """Load the sparse structure previously created with "trk2dictionary" script.

        Parameters
        ----------
        path : string
            Folder containing the output of the trk2dictionary script (relative to subject path)
        use_all_voxels_in_mask : boolean
            If False (default) the optimization will be conducted only on the voxels actually
            traversed by tracts. If True, then all voxels present in the mask specified in
            trk2dictionary.run(), i.e. "filename_mask" parameter, will be used instead.
            NB: if no mask was specified in trk2dictionary, this parameter is irrelevant.
        """
        if self.niiDWI is None :
            ERROR( 'Data not loaded; call "load_data()" first' )

        tic = time.time()
        LOG( '\n-> Loading the dictionary:' )
        self.DICTIONARY = {}
        self.set_config('TRACKING_path', pjoin(self.get_config('DATA_path'),path))

        # check that ndirs of dictionary matches with that of the kernels
        dictionary_info = load_dictionary_info( pjoin(self.get_config('TRACKING_path'), 'dictionary_info.pickle') )
        if dictionary_info['ndirs'] != self.get_config('ndirs'):
            ERROR( '"ndirs" of the dictionary (%d) does not match with the kernels (%d)' % (dictionary_info['ndirs'], self.get_config('ndirs')) )
        self.DICTIONARY['ndirs'] = dictionary_info['ndirs']

        # load mask
        self.set_config('dictionary_mask', 'mask' if use_all_voxels_in_mask else 'tdi' )
        mask_filename = pjoin(self.get_config('TRACKING_path'),'dictionary_%s.nii'%self.get_config('dictionary_mask'))
        if not exists( mask_filename ) :
            mask_filename += '.gz'
            if not exists( mask_filename ) :
                ERROR( 'Dictionary not found. Execute "trk2dictionary" script first' );
        niiMASK = nibabel.load( mask_filename )
        niiMASK_hdr = niiMASK.header if nibabel.__version__ >= '2.0.0' else niiMASK.get_header()
        if ( self.get_config('dim')[0]!=niiMASK.shape[0] or
             self.get_config('dim')[1]!=niiMASK.shape[1] or
             self.get_config('dim')[2]!=niiMASK.shape[2] or
             abs(self.get_config('pixdim')[0]-niiMASK_hdr['pixdim'][1])>1e-3 or
             abs(self.get_config('pixdim')[1]-niiMASK_hdr['pixdim'][2])>1e-3 or
             abs(self.get_config('pixdim')[2]-niiMASK_hdr['pixdim'][3])>1e-3 ) :
            WARNING( 'Dictionary does not have the same geometry as the dataset' )
        self.DICTIONARY['MASK'] = ( np.asanyarray(niiMASK.dataobj ) > 0).astype(np.uint8)

        # segments from the tracts
        # ------------------------
        print( '\t* Segments from the tracts... ', end='' )
        sys.stdout.flush()

        self.DICTIONARY['TRK'] = {}
        self.DICTIONARY['TRK']['kept']   = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_TRK_kept.dict'), dtype=np.uint8 )
        self.DICTIONARY['TRK']['norm']   = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_TRK_norm.dict'), dtype=np.float32 )
        self.DICTIONARY['TRK']['len']    = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_TRK_len.dict'), dtype=np.float32 )
        self.DICTIONARY['TRK']['lenTot'] = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_TRK_lenTot.dict'), dtype=np.float32 )


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

        # divide the length of each segment by the fiber length so that all the columns of the linear operator will have same length
        # NB: it works in conjunction with the normalization of the kernels
        cdef :
            np.float32_t [:] sl = self.DICTIONARY['IC']['len']
            np.float32_t [:] tl = self.DICTIONARY['TRK']['norm']
            np.uint32_t  [:] f  = self.DICTIONARY['IC']['fiber']
            int s
        if self.get_config('doNormalizeKernels') :
            for s in xrange(self.DICTIONARY['IC']['n']) :
                sl[s] /= tl[ f[s] ]

        print( '[ %d fibers and %d segments ]' % ( self.DICTIONARY['IC']['nF'], self.DICTIONARY['IC']['n'] ) )

        # segments from the peaks
        # -----------------------
        print( '\t* Segments from the peaks...  ', end='' )
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

        print( '[ %d segments ]' % self.DICTIONARY['EC']['nE'] )

        # isotropic compartments
        # ----------------------
        print( '\t* Isotropic contributions...  ', end='' )
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

        print( '[ %d voxels ]' % self.DICTIONARY['nV'] )

        # post-processing
        # ---------------
        print( '\t* Post-processing...          ', end='' )
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

        print( '[ OK ]' )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )


    def set_threads( self, n=None ) :
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
            ERROR( 'Number of threads must be between 1 and 255' )
        if self.DICTIONARY is None :
            ERROR( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.KERNELS is None :
            ERROR( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first' )

        self.THREADS = {}
        self.THREADS['n'] = n

        cdef :
            long [:] C
            long t, tot, i1, i2, N, c
            int i

        tic = time.time()
        LOG( '\n-> Distributing workload to different threads:' )
        print( '\t* Number of threads : %d' % n )

        # Distribute load for the computation of A*x product
        print( '\t* A operator...  ', end='' )
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
                ERROR( 'Too many threads for the IC compartments to evaluate; try decreasing the number', prefix='\n' )
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
                ERROR( 'Too many threads for the EC compartments to evaluate; try decreasing the number', prefix='\n' )
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
                ERROR( 'Too many threads for the ISO compartments to evaluate; try decreasing the number', prefix='\n' )
        else :
            self.THREADS['ISO'] = None

        print( '[ OK ]' )

        # Distribute load for the computation of At*y product
        print( '\t* A\' operator... ', end='' )
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
                ERROR( 'Too many threads for the EC compartments to evaluate; try decreasing the number', prefix='\n' )
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
                ERROR( 'Too many threads for the ISO compartments to evaluate; try decreasing the number', prefix='\n' )
        else :
            self.THREADS['ISOt'] = None

        print( '[ OK ]' )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )


    def build_operator( self, build_dir=None ) :
        """Compile/build the operator for computing the matrix-vector multiplications by A and A'
        using the informations from self.DICTIONARY, self.KERNELS and self.THREADS.
        NB: needs to call this function to update pointers to data structures in case
            the data is changed in self.DICTIONARY, self.KERNELS or self.THREADS.

        Parameters
        ----------
        build_dir : string
            The folder in which to store the compiled files.
            If None (default), they will end up in the .pyxbld directory in the user’s home directory.
            If using this option, it is recommended to use a temporary directory, quit your python
                console between each build, and delete the content of the temporary directory.
        """
        if self.DICTIONARY is None :
            ERROR( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.KERNELS is None :
            ERROR( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first' )
        if self.THREADS is None :
            ERROR( 'Threads not set; call "set_threads()" first' )

        if self.DICTIONARY['IC']['nF'] <= 0 :
            ERROR( 'No streamline found in the dictionary; check your data' )
        if self.DICTIONARY['EC']['nE'] <= 0 and self.KERNELS['wmh'].shape[0] > 0 :
            ERROR( 'The selected model has EC compartments, but no peaks have been provided; check your data' )

        tic = time.time()
        LOG( '\n-> Building linear operator A:' )

        # need to pass these parameters at runtime for compiling the C code
        from commit.operator import config

        compilation_is_needed = False

        if config.nTHREADS is None or config.nTHREADS != self.THREADS['n']:
            compilation_is_needed = True
        if config.nIC is None or config.nIC != self.KERNELS['wmr'].shape[0]:
            compilation_is_needed = True
        if config.model is None or config.model != self.model.id:
            compilation_is_needed = True
        if config.nEC is None or config.nEC != self.KERNELS['wmh'].shape[0]:
            compilation_is_needed = True
        if config.nISO is None or config.nISO != self.KERNELS['iso'].shape[0]:
            compilation_is_needed = True
        if config.build_dir != build_dir:
            compilation_is_needed = True

        if compilation_is_needed or not 'commit.operator.operator' in sys.modules :

            if build_dir is not None:
                if isdir(build_dir) and not len(listdir(build_dir)) == 0:
                    ERROR( '\nbuild_dir is not empty, unsafe build option.' )
                elif config.nTHREADS is not None:
                    ERROR( '\nThe parameter build_dir has changed, unsafe build option.' )
                else:
                    WARNING( '\nUsing build_dir, always quit your python console between COMMIT Evaluation.' )

            config.nTHREADS   = self.THREADS['n']
            config.model      = self.model.id
            config.nIC        = self.KERNELS['wmr'].shape[0]
            config.nEC        = self.KERNELS['wmh'].shape[0]
            config.nISO       = self.KERNELS['iso'].shape[0]
            config.build_dir  = build_dir

            pyximport.install( reload_support=True, language_level=3, build_dir=build_dir, build_in_temp=True, inplace=False )

            if not 'commit.operator.operator' in sys.modules :
                import commit.operator.operator
            else :
                reload( sys.modules['commit.operator.operator'] )

        self.A = sys.modules['commit.operator.operator'].LinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )


    def get_y( self ):
        """
        Returns a numpy array that corresponds to the 'y' vector of the optimisation problem.
        NB: this can be run only after having loaded the dictionary and the data.
        """
        if self.DICTIONARY is None :
            ERROR( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.niiDWI is None :
            ERROR( 'Data not loaded; call "load_data()" first' )
        return self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float64)


    def fit( self, tol_fun=1e-3, tol_x=1e-6, max_iter=100, verbose=True, x0=None, regularisation=None, confidence_map_filename=None, confidence_map_rescale=False ) :
        """Fit the model to the data.

        Parameters
        ----------
        tol_fun : float
            Tolerance on the objective function (default : 1e-3)
        max_iter : integer
            Maximum number of iterations (default : 100)
        verbose : boolean
            Level of verbosity: 0=no print, 1=print progress (default : True)
        x0 : np.array
            Initial guess for the solution of the problem (default : None)
        regularisation : commit.solvers.init_regularisation object
            Python dictionary that describes the wanted regularisation term.
            Check the documentation of commit.solvers.init_regularisation to see
            how to properly define the wanted mathematical formulation
            (default : None)
        confidence_map_filename : string
            Path to the NIFTI file containing a confidence map on the data, 
            relative to the subject folder. The file can be 3D or 4D in 
            the same space as the dwi_filename used (dim and voxel size).
            It should contain float values. 
            (default : None)
        confidence_map_rescale : boolean
            If true, the values of the confidence map will be rescaled to the 
            range [0.0,1.0]. Only the voxels considered in the mask will be affected.
            (default : False)
        """
        if self.niiDWI is None :
            ERROR( 'Data not loaded; call "load_data()" first' )
        if self.DICTIONARY is None :
            ERROR( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.KERNELS is None :
            ERROR( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first' )
        if self.THREADS is None :
            ERROR( 'Threads not set; call "set_threads()" first' )
        if self.A is None :
            ERROR( 'Operator not built; call "build_operator()" first' )
        
        # Confidence map
        self.confidence_map_img = None
        self.set_config('confidence_map_filename', None)
        confidence_array = None
        confidence_array_changed = False

        if confidence_map_filename is not None:
            # Loading confidence map
            tic = time.time()
            LOG( '\n-> Loading confidence map:' )

            if not exists( pjoin( self.get_config('DATA_path'), confidence_map_filename)  ) :            
                ERROR( 'Confidence map not found' )
            
            self.set_config('confidence_map_filename', confidence_map_filename)
            confidence_map  = nibabel.load( pjoin( self.get_config('DATA_path'), confidence_map_filename) )
            self.confidence_map_img = np.asanyarray( confidence_map.dataobj ).astype(np.float32)

            if self.confidence_map_img.ndim not in [3,4]:
                ERROR( 'Confidence map must be 3D or 4D dataset' )

            if self.confidence_map_img.ndim == 3:
                print( '\t* Extending the confidence map volume to match the DWI signal volume(s)... ' )                
                self.confidence_map_img = np.repeat(self.confidence_map_img[:, :, :, np.newaxis], self.niiDWI_img.shape[3], axis=3)
            hdr = confidence_map.header if nibabel.__version__ >= '2.0.0' else confidence_map.get_header()
            confidence_map_dim = self.confidence_map_img.shape[0:3]
            confidence_map_pixdim = tuple( hdr.get_zooms()[:3] )
            print( '\t\t- dim    : %d x %d x %d x %d' % self.confidence_map_img.shape )
            print( '\t\t- pixdim : %.3f x %.3f x %.3f' % confidence_map_pixdim )

            LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )
            
            if ( self.get_config('dim') != confidence_map_dim ):
                ERROR( 'Dataset does not have the same geometry (number of voxels) as the DWI signal' )

            if (self.get_config('pixdim') != confidence_map_pixdim ):
                ERROR( 'Dataset does not have the same geometry (voxel size) as the DWI signal' )
            
            if (self.confidence_map_img.shape != self.niiDWI_img.shape):
                ERROR( 'Dataset does not have the same geometry as the DWI signal' )

            confidence_array = self.confidence_map_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32)

            if(np.isnan(confidence_array).any()): 
                confidence_array[np.isnan(confidence_array)] = 0.
                confidence_array_changed = True
                WARNING('Confidence map contains NaNs. Those values were changed to 0.')
            
            cMAX = np.max(confidence_array)
            cMIN = np.min(confidence_array)  
            if(cMIN == cMAX): 
                self.confidence_map_img = None
                self.set_config('confidence_map_filename', None)
                confidence_array = None
                confidence_array_changed = False
                WARNING('All voxels in the confidence map have the same value. The confidence map will not be used')

            elif(confidence_map_rescale):
                confidence_array = ( confidence_array - cMIN ) / ( cMAX - cMIN )
                LOG ( '\n   [Confidence map interval was scaled from the original [%.1f, %.1f] to the intended [%.1f, %.1f] linearly]' % ( cMIN, cMAX, np.min(confidence_array), np.max(confidence_array) ) ) 
                confidence_array_changed = True
            
            if(confidence_array_changed):
                nV = self.DICTIONARY['nV']
                self.confidence_map_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ] = np.reshape( confidence_array, (nV,-1) ).astype(np.float32)

        if x0 is not None :
            if x0.shape[0] != self.A.shape[1] :
                ERROR( 'x0 dimension does not match the number of columns of the dictionary' )
        if regularisation is None :
            regularisation = commit.solvers.init_regularisation(self)

        self.CONFIG['optimization']                   = {}
        self.CONFIG['optimization']['tol_fun']        = tol_fun
        self.CONFIG['optimization']['tol_x']          = tol_x
        self.CONFIG['optimization']['max_iter']       = max_iter
        self.CONFIG['optimization']['verbose']        = verbose
        self.CONFIG['optimization']['regularisation'] = regularisation

        # run solver
        t = time.time()
        LOG( '\n-> Fit model:' )

        self.x, opt_details = commit.solvers.solve(self.get_y(), self.A, self.A.T, tol_fun = tol_fun, tol_x = tol_x, max_iter = max_iter, verbose = verbose, x0 = x0, regularisation = regularisation, confidence_array = confidence_array)

        self.CONFIG['optimization']['fit_details'] = opt_details
        self.CONFIG['optimization']['fit_time'] = time.time()-t

        LOG( '\n   [ %s ]' % ( time.strftime("%Hh %Mm %Ss", time.gmtime(self.CONFIG['optimization']['fit_time']) ) ) )


    def get_coeffs( self, get_normalized=True ):
        """
        Returns the coefficients estimated from the tractogram passed to trk2dictionary.run(). These coefficients are divided in three
        classes (ic, ec, iso) and correspond to the 'original optimization problem', so if preconditioning was applied via the option
        'doNormalizeKernels', then they are re-scaled accordingly. This behaviour can be overridden using the 'get_normalized' parameter.

        Parameters
        ----------
        get_normalized : boolean
            If True (default), the returned coefficients correspond to the 'original' optimization problem.
            If False, the returned coefficients correspond to the 'preconditioned' optimization problem.
        """
        if self.x is None :
            ERROR( 'Model not fitted to the data; call "fit()" first' )

        nF = self.DICTIONARY['IC']['nF']
        nE = self.DICTIONARY['EC']['nE']
        nV = self.DICTIONARY['nV']

        if get_normalized and self.get_config('doNormalizeKernels') :
            # renormalize the coefficients
            norm1 = np.repeat(self.KERNELS['wmr_norm'],nF)
            norm2 = np.repeat(self.KERNELS['wmh_norm'],nE)
            norm3 = np.repeat(self.KERNELS['iso_norm'],nV)
            norm_fib = np.kron(np.ones(self.KERNELS['wmr'].shape[0]), self.DICTIONARY['TRK']['norm'])
            x = self.x / np.hstack( (norm1*norm_fib,norm2,norm3) )
        else :
            x = self.x

        offset1 = nF * self.KERNELS['wmr'].shape[0]
        offset2 = offset1 + nE * self.KERNELS['wmh'].shape[0]
        kept = np.tile( self.DICTIONARY['TRK']['kept'], self.KERNELS['wmr'].shape[0] )
        xic = np.zeros( kept.size )
        xic[kept==1] = x[:offset1]
        xec = x[offset1:offset2]
        xiso = x[offset2:]

        return xic, xec, xiso


    def save_results( self, path_suffix=None, stat_coeffs='sum', save_est_dwi=False, save_coeff=None, save_opt_details=None ) :
        """Save the output (coefficients, errors, maps etc).

        Parameters
        ----------
        path_suffix : string
            Text to be appended to "Results" to create the output path (default : None)
        stat_coeffs : string
            Stat to be used if more coefficients are estimated for each streamline.
            Options: 'sum', 'mean', 'median', 'min', 'max', 'all' (default : 'sum')
        save_est_dwi : boolean
            Save the estimated DW-MRI signal (default : False)
        save_opt_details : boolean
            DEPRECATED. The details of the optimization and the coefficients are always saved.
        save_coeff : boolean
            DEPRECATED. The estimated weights for the streamlines are always saved.
        """
        RESULTS_path = 'Results_' + self.model.id
        if path_suffix :
            self.set_config('path_suffix', path_suffix)
            RESULTS_path = RESULTS_path + path_suffix

        LOG( '\n-> Saving results to "%s/*":' % RESULTS_path )
        tic = time.time()

        if self.x is None :
            ERROR( 'Model not fitted to the data; call "fit()" first' )

        if save_coeff is not None :
            WARNING('"save_coeff" parameter is deprecated')

        if save_opt_details is not None :
            WARNING('"save_opt_details" parameter is deprecated')

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

        # create folder or delete existing files (if any)
        RESULTS_path = pjoin( self.get_config('TRACKING_path'), RESULTS_path )
        if not exists( RESULTS_path ) :
            makedirs( RESULTS_path )
        else :
            for f in glob.glob( pjoin(RESULTS_path,'*') ) :
                remove( f )
        self.set_config('RESULTS_path', RESULTS_path)

        # Map of voxelwise errors
        print( '\t* Fitting errors:' )

        niiMAP_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        affine = self.niiDWI.affine if nibabel.__version__ >= '2.0.0' else self.niiDWI.get_affine()
        niiMAP     = nibabel.Nifti1Image( niiMAP_img, affine )
        niiMAP_hdr = niiMAP.header if nibabel.__version__ >= '2.0.0' else niiMAP.get_header()
        niiMAP_hdr['descrip'] = 'Created with COMMIT %s'%self.get_config('version')

        y_mea = np.reshape( self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32), (nV,-1) )
        y_est = np.reshape( self.A.dot(self.x), (nV,-1) ).astype(np.float32)

        print( '\t\t- RMSE...  ', end='' )
        sys.stdout.flush()
        tmp = np.sqrt( np.mean((y_mea-y_est)**2,axis=1) )
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP_hdr['cal_min'] = 0
        niiMAP_hdr['cal_max'] = tmp.max()
        nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_RMSE.nii.gz') )
        print( '[ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() ) )

        print( '\t\t- NRMSE... ', end='' )
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
        print( '[ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() ) )
        
        if self.confidence_map_img is not None:
            confidence_array = np.reshape( self.confidence_map_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32), (nV,-1) )
            
            print( '\t\t- RMSE considering the confidence map...  ', end='' )        
            sys.stdout.flush()
            tmp = np.sum(confidence_array,axis=1)
            idx = np.where( tmp < 1E-12 )
            tmp[ idx ] = 1 
            tmp = np.sqrt( np.sum(confidence_array*(y_mea-y_est)**2,axis=1) / tmp )
            tmp[ idx ] = 0
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
            niiMAP_hdr['cal_min'] = 0
            niiMAP_hdr['cal_max'] = tmp.max()
            nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_RMSE_adjusted.nii.gz') )
            print( '[ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() ) )

            print( '\t\t- NRMSE considering the confidence map... ', end='' )
            sys.stdout.flush()
            tmp = np.sum(confidence_array*y_mea**2,axis=1)
            idx = np.where( tmp < 1E-12 )
            tmp[ idx ] = 1
            tmp = np.sqrt( np.sum(confidence_array*(y_mea-y_est)**2,axis=1) / tmp )
            tmp[ idx ] = 0
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
            niiMAP_hdr['cal_min'] = 0
            niiMAP_hdr['cal_max'] = 1
            nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_NRMSE_adjusted.nii.gz') )
            print( '[ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() ) )
            confidence_array = None

        # Map of compartment contributions
        print( '\t* Voxelwise contributions:' )

        print( '\t\t- Intra-axonal... ', end='' )
        sys.stdout.flush()
        niiIC_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        if len(self.KERNELS['wmr']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0]
            tmp = ( x[:offset].reshape( (-1,nF) ) * norm_fib.reshape( (-1,nF) ) ).sum( axis=0 )
            xv = np.bincount( self.DICTIONARY['IC']['v'], minlength=nV,
                weights=tmp[ self.DICTIONARY['IC']['fiber'] ] * self.DICTIONARY['IC']['len']
            ).astype(np.float32)
            niiIC_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        print( '[ OK ]' )

        print( '\t\t- Extra-axonal... ', end='' )
        sys.stdout.flush()
        niiEC_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        if len(self.KERNELS['wmh']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0]
            tmp = x[offset:offset+nE*len(self.KERNELS['wmh'])].reshape( (-1,nE) ).sum( axis=0 )
            xv = np.bincount( self.DICTIONARY['EC']['v'], weights=tmp, minlength=nV ).astype(np.float32)
            niiEC_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        print( '[ OK ]' )

        print( '\t\t- Isotropic... ', end='' )
        sys.stdout.flush()
        niiISO_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        if len(self.KERNELS['iso']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0] + nE * self.KERNELS['wmh'].shape[0]
            xv = x[offset:].reshape( (-1,nV) ).sum( axis=0 )
            niiISO_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        print( '   [ OK ]' )

        if self.get_config('doNormalizeMaps') :
            niiIC = nibabel.Nifti1Image(  niiIC_img  / ( niiIC_img + niiEC_img + niiISO_img + 1e-16), affine, header=niiMAP_hdr )
            niiEC = nibabel.Nifti1Image(  niiEC_img /  ( niiIC_img + niiEC_img + niiISO_img + 1E-16), affine, header=niiMAP_hdr )
            niiISO = nibabel.Nifti1Image( niiISO_img / ( niiIC_img + niiEC_img + niiISO_img + 1E-16), affine, header=niiMAP_hdr )
        else:
            niiIC = nibabel.Nifti1Image(  niiIC_img,  affine, header=niiMAP_hdr )
            niiEC = nibabel.Nifti1Image(  niiEC_img,  affine, header=niiMAP_hdr )
            niiISO = nibabel.Nifti1Image( niiISO_img, affine, header=niiMAP_hdr )

        nibabel.save( niiIC , pjoin(RESULTS_path,'compartment_IC.nii.gz') )
        nibabel.save( niiEC , pjoin(RESULTS_path,'compartment_EC.nii.gz') )
        nibabel.save( niiISO , pjoin(RESULTS_path,'compartment_ISO.nii.gz') )

        # Configuration and results
        print( '\t* Configuration and results:' )

        print( '\t\t- streamline_weights.txt... ', end='' )
        sys.stdout.flush()
        xic, _, _ = self.get_coeffs()
        if stat_coeffs != 'all' and xic.size > 0 :
            xic = np.reshape( xic, (-1,self.DICTIONARY['TRK']['kept'].size) )
            if stat_coeffs == 'sum' :
                xic = np.sum( xic, axis=0 )
            elif stat_coeffs == 'mean' :
                xic = np.mean( xic, axis=0 )
            elif stat_coeffs == 'median' :
                xic = np.median( xic, axis=0 )
            elif stat_coeffs == 'min' :
                xic = np.min( xic, axis=0 )
            elif stat_coeffs == 'max' :
                xic = np.max( xic, axis=0 )
            else :
                ERROR( 'Stat not allowed. Possible values: sum, mean, median, min, max, all', prefix='\n' )

        # scale output weights if blur was used
        dictionary_info = load_dictionary_info( pjoin(self.get_config('TRACKING_path'), 'dictionary_info.pickle') )
        if dictionary_info['blur_sigma'] > 0 :
            if stat_coeffs == 'all' :
                ERROR( 'Not yet implemented. Unable to account for blur in case of multiple streamline constributions.' )
            xic[ self.DICTIONARY['TRK']['kept']==1 ] *= self.DICTIONARY['TRK']['lenTot'] / self.DICTIONARY['TRK']['len']
            
        np.savetxt( pjoin(RESULTS_path,'streamline_weights.txt'), xic, fmt='%.5e' )
        self.set_config('stat_coeffs', stat_coeffs)
        print( '[ OK ]' )

        # Save to a pickle file the following items:
        #   item 0: dictionary with all the configuration details
        #   item 1: np.array obtained through the optimisation process with the normalised kernels
        #   item 2: np.array renormalisation of coeffs in item 1
        print( '\t\t- results.pickle... ', end='' )
        sys.stdout.flush()
        with open( pjoin(RESULTS_path,'results.pickle'), 'wb+' ) as fid :
            pickle.dump( [self.CONFIG, self.x, x], fid, protocol=2 )
        print( '        [ OK ]' )

        if save_est_dwi :
            print( '\t\t- Estimated signal... ', end='' )
            sys.stdout.flush()
            self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ] = y_est
            nibabel.save( nibabel.Nifti1Image( self.niiDWI_img , affine ), pjoin(RESULTS_path,'fit_signal_estimated.nii.gz') )
            self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ] = y_mea
            print( '[ OK ]' )

        LOG( '   [ %.1f seconds ]' % ( time.time() - tic ) )
