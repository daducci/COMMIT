#!python
#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False, binding=False
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
from dicelib.ui import __logger__ as logger
from dicelib.ui import ProgressBar
from dicelib import ui
from importlib import reload, invalidate_caches
import pyximport
from pkg_resources import get_distribution



def setup( lmax=12 ) :
    """General setup/initialization of the COMMIT framework.

    Parameters
    ----------
    lmax : int
        Maximum SH order to use for the rotation phase (default : 12)
    """
    amico.setup( lmax )


def load_dictionary_info( filename ):
    """Function to load dictionary info file

    Parameters
    ----------
    filename : string
        This value is always COMMIT_PATH + dictionary_info.pickle
    """
    if not isfile( filename ):
        logger.error( 'Dictionary is outdated or not found. Execute "trk2dictionary" script first' )
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
    cdef public regularisation_params
    cdef public x
    cdef public CONFIG
    cdef public temp_data
    cdef public confidence_map_img
    cdef public verbose

    def __init__( self, study_path='.', subject='.' ) :
        """Setup the data structures with default values.

        Parameters
        ----------
        study_path : string
            The path to the folder containing all the subjects from one study (default : '.')
        subject : string
            The path (relative to previous folder) to the subject folder (default : '.')
        """
        self.niiDWI                 = None # set by "load_data" method
        self.scheme                 = None # set by "load_data" method
        self.model                  = None # set by "set_model" method
        self.KERNELS                = None # set by "load_kernels" method
        self.DICTIONARY             = None # set by "load_dictionary" method
        self.THREADS                = None # set by "set_threads" method
        self.A                      = None # set by "build_operator" method
        self.regularisation_params  = None # set by "set_regularisation" method
        self.x                      = None # set by "fit" method
        self.confidence_map_img     = None # set by "fit" method
        self.verbose                = 3

        # store all the parameters of an evaluation with COMMIT
        self.CONFIG = {}
        self.temp_data = {}
        self.set_config('version', get_distribution('dmri-commit').version)
        self.set_config('study_path', study_path)
        self.set_config('subject', subject)
        self.set_config('DATA_path', pjoin( study_path, subject ))

        self.set_config('doNormalizeSignal', True)
        self.set_config('doMergeB0', False)
        self.set_config('doNormalizeKernels', True)
        self.set_config('doDemean', False)
        self.set_config('doNormalizeMaps', False)

        ui.set_verbose( self.verbose )


    def set_verbose( self, verbose ) :
        """Set the verbosity level of the logger.

        Parameters
        ----------
        verbose : int
            The verbosity level (0: only errors, 1: errors and warnings, 2: errors, warnings and info, 3: errors, warnings, info and progress bars, 4: errors, warnings, info, progress bars and debug)
        """
        self.verbose = verbose
        ui.set_verbose( verbose )

    def set_config( self, key, value ) :
        self.CONFIG[ key ] = value
        self.temp_data[ key ] = value


    def get_config( self, key ) :
        return self.CONFIG.get( key )


    def get_temp_data( self, key ) :
        return self.temp_data.get( key )


    def load_data( self, dwi_filename, scheme_filename=None, b0_thr=0, b0_min_signal=0, replace_bad_voxels=None ) :
        """Load the diffusion signal and its corresponding acquisition scheme.

        Parameters
        ----------
        dwi_filename : string
            The filename of the DWI data, relative to the subject folder.
        scheme_filename : string
            The file name of the corresponding acquisition scheme.
            If None, assumes the data is scalar, i.e. 1 value per voxel.
        b0_thr : float
            The threshold below which a b-value is considered a b0 (default : 0)
        b0_min_signal : float
            Crop to zero the signal in voxels where the b0 <= b0_min_signal * mean(b0[b0>0]) (default : 0)
        replace_bad_voxels : float, optional
            Value to be used to fill NaN and Inf values in the signal. (default : do nothing)
        """

        # Loading data and acquisition scheme
        tic = time.time()
        logger.subinfo('')
        logger.info( 'Loading data:' )

        logger.subinfo('Acquisition scheme:', indent_lvl=1, indent_char='*' )
        if scheme_filename is not None:
            self.set_config('scheme_filename', scheme_filename)
            self.set_config('b0_thr', b0_thr)
            logger.subinfo('diffusion-weighted signal', indent_char='-', indent_lvl=2)
            self.scheme = amico.scheme.Scheme( pjoin( self.get_config('DATA_path'), scheme_filename), b0_thr )
            logger.subinfo('%d samples, %d shells' % ( self.scheme.nS, len(self.scheme.shells) ) , indent_lvl=2, indent_char='-' )
            scheme_string = f'{self.scheme.b0_count} @ b=0'
            for i in xrange(len(self.scheme.shells)) :
                scheme_string += f', {len(self.scheme.shells[i]["idx"])} @ b={self.scheme.shells[i]["b"]:.1f}'
            logger.subinfo( scheme_string, indent_lvl=2, indent_char='-' )
        else:
            # if no scheme is passed, assume data is scalar
            self.scheme = amico.scheme.Scheme( np.array( [[0,0,0,1000]] ), 0 )
            logger.subinfo('scalar map', indent_char='-', indent_lvl=2)

        logger.subinfo('Signal dataset:', indent_lvl=1, indent_char='*' )
        self.set_config('dwi_filename', dwi_filename)
        self.set_config('b0_min_signal', b0_min_signal)
        self.set_config('replace_bad_voxels', replace_bad_voxels)
        self.niiDWI  = nibabel.load( pjoin( self.get_config('DATA_path'), dwi_filename) )
        self.niiDWI_img = np.asanyarray( self.niiDWI.dataobj ).astype(np.float32)
        if self.niiDWI_img.ndim ==3 :
            self.niiDWI_img = np.expand_dims( self.niiDWI_img, axis=3 )
        hdr = self.niiDWI.header if nibabel.__version__ >= '2.0.0' else self.niiDWI.get_header()
        self.set_config('dim', self.niiDWI_img.shape[0:3])
        self.set_config('pixdim', tuple( hdr.get_zooms()[:3] ))
        logger.subinfo('dim    : %d x %d x %d x %d' % self.niiDWI_img.shape, indent_lvl=2, indent_char='-' )
        logger.subinfo('pixdim : %.3f x %.3f x %.3f' % self.get_config('pixdim'), indent_lvl=2, indent_char='-' )
        logger.subinfo('values : min=%.2f, max=%.2f, mean=%.2f' % ( self.niiDWI_img.min(), self.niiDWI_img.max(), self.niiDWI_img.mean() ), indent_lvl=2, indent_char='-' )

        if self.scheme.nS != self.niiDWI_img.shape[3] :
            logger.error( 'Scheme does not match with input data' )
        if self.scheme.dwi_count == 0 :
            logger.error( 'There are no DWI volumes in the data' )

        # Check for Nan or Inf values in raw data
        if np.isnan(self.niiDWI_img).any() or np.isinf(self.niiDWI_img).any():
            if replace_bad_voxels is not None:
                logger.warning(f'Nan or Inf values in the raw signal. They will be replaced with: {replace_bad_voxels}')
                np.nan_to_num(self.niiDWI_img, copy=False, nan=replace_bad_voxels, posinf=replace_bad_voxels, neginf=replace_bad_voxels)
            else:
                logger.error('Nan or Inf values in the raw signal. Try using the "replace_bad_voxels" or "b0_min_signal" parameters when calling "load_data()"')

        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )

        # Preprocessing
        if self.get_config('scheme_filename') is not None:
            tic = time.time()
            logger.subinfo('')
            logger.info( 'Preprocessing:' )

            if self.get_config('doNormalizeSignal') :
                if self.scheme.b0_count > 0:
                    logger.subinfo(' Normalizing to b0', with_progress=True, indent_lvl=1, indent_char='*')
                    with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
                        b0 = np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 )
                        idx = b0 <= b0_min_signal * b0[b0>0].mean()
                        b0[ idx ] = 1
                        b0 = 1.0 / b0
                        b0[ idx ] = 0
                        for i in xrange(self.scheme.nS) :
                            self.niiDWI_img[:,:,:,i] *= b0
                    logger.subinfo( '[ min=%.2f, max=%.2f, mean=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.max(), self.niiDWI_img.mean() ), indent_lvl=1 )
                    del idx, b0
                else :
                    logger.warning( 'There are no b0 volumes for normalization' )

            if self.scheme.b0_count > 1:
                if self.get_config('doMergeB0') :
                    logger.subinfo('Merging multiple b0 volume(s)', indent_char='*', indent_lvl=1)
                    mean = np.expand_dims( np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 ), axis=3 )
                    self.niiDWI_img = np.concatenate( (mean, self.niiDWI_img[:,:,:,self.scheme.dwi_idx]), axis=3 )
                    del mean
                else :
                    logger.subinfo('Keeping all b0 volume(s)', indent_char='*', indent_lvl=1)
                logger.subinfo('[ %d x %d x %d x %d ]' % self.niiDWI_img.shape )

            if self.get_config('doDemean'):
                mean = np.repeat( np.expand_dims(np.mean(self.niiDWI_img,axis=3),axis=3), self.niiDWI_img.shape[3], axis=3 )
                self.niiDWI_img = self.niiDWI_img - mean
                logger.subinfo('Demeaning signal [ min=%.2f, max=%.2f, mean=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.max(), self.niiDWI_img.mean() ), indent_lvl=1, indent_char='*' )

            # Check for Nan or Inf values in pre-processed data
            if np.isnan(self.niiDWI_img).any() or np.isinf(self.niiDWI_img).any():
                if replace_bad_voxels is not None:
                    logger.warning(f'Nan or Inf values in the signal after the pre-processing. They will be replaced with: {replace_bad_voxels}')
                    np.nan_to_num(self.niiDWI_img, copy=False, nan=replace_bad_voxels, posinf=replace_bad_voxels, neginf=replace_bad_voxels)
                else:
                    logger.error('Nan or Inf values in the signal after the pre-processing. Try using the "replace_bad_voxels" or "b0_min_signal" parameters when calling "load_data()"')

            logger.subinfo( f'[ {( time.time() - tic ):.1f} seconds ]' )


    def set_model( self, model_name ) :
        """Set the model to use to describe the signal contributions in each voxel.

        Parameters
        ----------
        model_name : string
            The name of the model (must match a class name in "commit.models" module)
        """
        # Call the specific model constructor
        if hasattr(commit.models, model_name ) :
            self.model = getattr(commit.models, model_name)()
        else :
            logger.error( 'Model "%s" not recognized' % model_name )

        self.set_config('ATOMS_path', pjoin( self.get_config('study_path'), 'kernels', self.model.id ))


    def generate_kernels( self, regenerate=False, lmax=12, ndirs=500 ) :
        """Generate the high-resolution response functions for each compartment.
        Dispatch to the proper function, depending on the model.

        Parameters
        ----------
        regenerate : boolean
            Regenerate kernels if they already exist (default : False)
        lmax : int
            Maximum SH order to use for the rotation procedure (default : 12)
        ndirs : int
            Number of directions on the half of the sphere representing the possible orientations of the response functions (default : 500)
        """
        logger.subinfo('')
        logger.info( 'Simulating with "%s" model:' % self.model.name )

        if not amico.lut.is_valid( ndirs ):
            logger.error( 'Unsupported value for ndirs.\nNote: Supported values for ndirs are [1, 500 (default), 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 32761]' )
        if self.scheme is None :
            logger.error( 'Scheme not loaded; call "load_data()" first' )
        if self.model is None :
            logger.error( 'Model not set; call "set_model()" method first' )
        if self.model.id=='VolumeFractions' and ndirs!=1:
            ndirs = 1
            logger.subinfo('Forcing "ndirs" to 1 because model is isotropic', indent_char='*', indent_lvl=1)
        if 'commitwipmodels' in sys.modules :
            if self.model.restrictedISO is not None and ndirs!=1:
                ndirs = 1
                logger.subinfo('Forcing "ndirs" to 1 because model is isotropic', indent_char='*', indent_lvl=1)
 
        # store some values for later use
        self.set_config('lmax', lmax)
        self.set_config('ndirs', ndirs)
        self.set_config('model', self.model.get_params())
        self.model.scheme = self.scheme

        # check if kernels were already generated
        tmp = glob.glob( pjoin(self.get_config('ATOMS_path'),'A_*.npy') )
        if len(tmp)>0 and not regenerate :
            logger.subinfo( '[ Kernels already computed. Use option "regenerate=True" to force regeneration ]' )
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
        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )


    def load_kernels( self ) :
        """Load rotated kernels and project to the specific gradient scheme of this subject.
        Dispatch to the proper function, depending on the model.
        """
        if self.model is None :
            logger.error( 'Model not set; call "set_model()" method first' )
        if self.scheme is None :
            logger.error( 'Scheme not loaded; call "load_data()" first' )

        tic = time.time()
        logger.subinfo('')
        logger.info( 'Loading the kernels:' )
        logger.subinfo( 'Resampling LUT for subject "%s":' % self.get_config('subject'), indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            # auxiliary data structures
            idx_OUT, Ylm_OUT = amico.lut.aux_structures_resample( self.scheme, self.get_config('lmax') )

            self.KERNELS = self.model.resample( self.get_config('ATOMS_path'), idx_OUT, Ylm_OUT, self.get_config('doMergeB0'), self.get_config('ndirs') )
            nIC  = self.KERNELS['wmr'].shape[0]
            nEC  = self.KERNELS['wmh'].shape[0]
            nISO = self.KERNELS['iso'].shape[0]

        # Dispatch to the right handler for each model
        if self.get_config('doMergeB0') :
            logger.subinfo( 'Merging multiple b0 volume(s) [-OK-]', indent_char='-', indent_lvl=1 )
        else :
            logger.subinfo( 'Keeping all b0 volume(s)      [-OK-]', indent_char='-', indent_lvl=1 )

        # ensure contiguous arrays for C part
        self.KERNELS['wmr'] = np.ascontiguousarray( self.KERNELS['wmr'] )
        self.KERNELS['wmh'] = np.ascontiguousarray( self.KERNELS['wmh'] )
        self.KERNELS['iso'] = np.ascontiguousarray( self.KERNELS['iso'] )

        # De-mean kernels
        if self.get_config('doDemean') :
            logger.subinfo('Demeaning signal            ', with_progress=True, indent_lvl=1, indent_char='-' )
            with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
                for j in xrange(self.get_config('ndirs')) :
                    for i in xrange(nIC) :
                        self.KERNELS['wmr'][i,j,:] -= self.KERNELS['wmr'][i,j,:].mean()
                    for i in xrange(nEC) :
                        self.KERNELS['wmh'][i,j,:] -= self.KERNELS['wmh'][i,j,:].mean()
                for i in xrange(nISO) :
                    self.KERNELS['iso'][i] -= self.KERNELS['iso'][i].mean()

        # Normalize atoms
        if self.get_config('doNormalizeKernels') :
            logger.subinfo('Normalizing kernels         ', with_progress=True, indent_lvl=1, indent_char='-' )
            with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
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

        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )


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
            logger.error( 'Data not loaded; call "load_data()" first' )

        tic = time.time()
        logger.subinfo('')
        logger.info( 'Loading the dictionary:' )
        self.DICTIONARY = {}
        self.set_config('TRACKING_path', pjoin(self.get_config('DATA_path'),path))

        # check that ndirs of dictionary matches with that of the kernels
        dictionary_info = load_dictionary_info( pjoin(self.get_config('TRACKING_path'), 'dictionary_info.pickle') )
        if dictionary_info['ndirs'] != self.get_config('ndirs'):
            logger.error( '"ndirs" of the dictionary (%d) does not match with the kernels (%d)' % (dictionary_info['ndirs'], self.get_config('ndirs')) )
        self.DICTIONARY['ndirs'] = dictionary_info['ndirs']
        self.DICTIONARY['n_threads'] = dictionary_info['n_threads']

        # load mask
        self.set_config('dictionary_mask', 'mask' if use_all_voxels_in_mask else 'tdi' )
        mask_filename = pjoin(self.get_config('TRACKING_path'),'dictionary_%s.nii'%self.get_config('dictionary_mask'))
        if not exists( mask_filename ) :
            mask_filename += '.gz'
            if not exists( mask_filename ) :
                logger.error( 'Dictionary not found. Execute "trk2dictionary" script first' );
        niiMASK = nibabel.load( mask_filename )
        niiMASK_hdr = niiMASK.header if nibabel.__version__ >= '2.0.0' else niiMASK.get_header()
        if ( self.get_config('dim')[0]!=niiMASK.shape[0] or
             self.get_config('dim')[1]!=niiMASK.shape[1] or
             self.get_config('dim')[2]!=niiMASK.shape[2] or
             abs(self.get_config('pixdim')[0]-niiMASK_hdr['pixdim'][1])>1e-3 or
             abs(self.get_config('pixdim')[1]-niiMASK_hdr['pixdim'][2])>1e-3 or
             abs(self.get_config('pixdim')[2]-niiMASK_hdr['pixdim'][3])>1e-3 ) :
            logger.warning( 'Dictionary does not have the same geometry as the dataset' )
        self.DICTIONARY['MASK'] = ( np.asanyarray(niiMASK.dataobj ) > 0).astype(np.uint8)

        # segments from the tracts
        # ------------------------
        logger.subinfo('Segments from the tracts', indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
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

        logger.subinfo(f"{self.DICTIONARY['IC']['nF']} fibers and {self.DICTIONARY['IC']['n']} segments", indent_char='-', indent_lvl=1 )

        # segments from the peaks
        # -----------------------
        logger.subinfo('Segments from the peaks ', indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            self.DICTIONARY['EC'] = {}
            self.DICTIONARY['EC']['v']  = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_EC_v.dict'), dtype=np.uint32 )
            self.DICTIONARY['EC']['o']  = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_EC_o.dict'), dtype=np.uint16 )
            self.DICTIONARY['EC']['nE'] = self.DICTIONARY['EC']['v'].size

            # reorder the segments based on the "v" field
            idx = np.argsort( self.DICTIONARY['EC']['v'], kind='mergesort' )
            self.DICTIONARY['EC']['v'] = self.DICTIONARY['EC']['v'][ idx ]
            self.DICTIONARY['EC']['o'] = self.DICTIONARY['EC']['o'][ idx ]
            del idx

        logger.subinfo( f"{self.DICTIONARY['EC']['nE']} segments", indent_char='-', indent_lvl=1)

        # isotropic compartments
        # ----------------------
        logger.subinfo('Isotropic contributions ', indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            self.DICTIONARY['ISO'] = {}

            self.DICTIONARY['ISO']['v'] = np.fromfile( pjoin(self.get_config('TRACKING_path'),'dictionary_ISO_v.dict'), dtype=np.uint32 )
            self.DICTIONARY['ISO']['nV'] = self.DICTIONARY['ISO']['v'].size
                
            self.DICTIONARY['nV'] = self.DICTIONARY['MASK'].sum()

            # reorder the segments based on the "v" field
            idx = np.argsort( self.DICTIONARY['ISO']['v'], kind='mergesort' )
            self.DICTIONARY['ISO']['v'] = self.DICTIONARY['ISO']['v'][ idx ]
            del idx

        logger.subinfo( f"{self.DICTIONARY['ISO']['nV']} voxels", indent_char='-', indent_lvl=1 )
        
        # post-processing
        # ---------------
        logger.subinfo('Post-processing         ', indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            # get the indices to extract the VOI as in MATLAB (in place of DICTIONARY.MASKidx)
            idx = self.DICTIONARY['MASK'].ravel(order='F').nonzero()[0]
            self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] = np.unravel_index( idx, self.DICTIONARY['MASK'].shape, order='F' )

            lut = np.zeros( self.get_config('dim'), dtype=np.uint32 ).ravel()
            for i in xrange(idx.size) :
                lut[ idx[i] ] = i
            self.DICTIONARY['IC'][ 'v'] = lut[ self.DICTIONARY['IC'][ 'v'] ]
            self.DICTIONARY['EC'][ 'v'] = lut[ self.DICTIONARY['EC'][ 'v'] ]
            self.DICTIONARY['ISO']['v'] = lut[ self.DICTIONARY['ISO']['v'] ]

        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )



    def set_threads( self, n=None ) :
        """Set the number of threads to use for the matrix-vector operations with A and A'.

        Parameters
        ----------
        n : integer
            Number of threads to use (default : number of CPUs in the system)
        """
        if n is None :
            # Use the same number of threads used in trk2dictionary
            n = self.DICTIONARY['n_threads']


        if n < 1 or n > 255 :
            logger.error( 'Number of threads must be between 1 and 255' )
        if self.DICTIONARY is None :
            logger.error( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.KERNELS is None :
            logger.error( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first' )

        self.THREADS = {}
        self.THREADS['n'] = n

        cdef :
            long [:] C
            long t, tot, i1, i2, N, c
            int i

        tic = time.time()
        logger.subinfo('')
        logger.info( 'Distributing workload to different threads:' )
        logger.subinfo('Number of threads : %d' % n , indent_char='*' )

        # Distribute load for the computation of A*x product
        logger.subinfo('A operator ', indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            if self.DICTIONARY['IC']['n'] > 0 :
                self.THREADS['IC'] = np.zeros( n+1, dtype=np.uint32 )
                if n > 1 :
                    N = np.floor( self.DICTIONARY['IC']['n']/n )
                    t = 1
                    tot = 0
                    C = np.bincount( self.DICTIONARY['IC']['v'] )
                    for c in C :
                        tot += c
                        if tot >= N and t <= n :
                            self.THREADS['IC'][t] = self.THREADS['IC'][t-1] + tot
                            t += 1
                            tot = 0
                self.THREADS['IC'][n] = self.DICTIONARY['IC']['n']

                # check if some threads are not assigned any segment
                if np.count_nonzero( np.diff( self.THREADS['IC'].astype(np.int32) ) <= 0 ) :
                    self.THREADS = None
                    logger.error( 'Too many threads for the IC compartments to evaluate; try decreasing the number', prefix='\n' )
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
                    logger.error( 'Too many threads for the EC compartments to evaluate; try decreasing the number', prefix='\n' )
            else :
                self.THREADS['EC'] = None

            if self.DICTIONARY['nV'] > 0 :
                self.THREADS['ISO'] = np.zeros( n+1, dtype=np.uint32 )
                for i in xrange(n) :
                    self.THREADS['ISO'][i] = np.searchsorted( self.DICTIONARY['ISO']['v'], self.DICTIONARY['IC']['v'][ self.THREADS['IC'][i] ] )
                self.THREADS['ISO'][n] = self.DICTIONARY['ISO']['nV']

                # check if some threads are not assigned any segment
                if np.count_nonzero( np.diff( self.THREADS['ISO'].astype(np.int32) ) <= 0 ) :
                    self.THREADS = None
                    logger.error( 'Too many threads for the ISO compartments to evaluate; try decreasing the number', prefix='\n' )
            else :
                self.THREADS['ISO'] = None


        # Distribute load for the computation of At*y product
        logger.subinfo('A\' operator', indent_char='*', with_progress=True )
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
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
                    logger.error( 'Too many threads for the EC compartments to evaluate; try decreasing the number', prefix='\n' )
            else :
                self.THREADS['ECt'] = None

            if self.DICTIONARY['nV'] > 0 :
                self.THREADS['ISOt'] = np.zeros( n+1, dtype=np.uint32 )
                N = np.floor( self.DICTIONARY['ISO']['nV']/n )
                
                for i in xrange(1,n) :
                    self.THREADS['ISOt'][i] = self.THREADS['ISOt'][i-1] + N
                self.THREADS['ISOt'][n] = self.DICTIONARY['ISO']['nV']
                
                # check if some threads are not assigned any segment
                if np.count_nonzero( np.diff( self.THREADS['ISOt'].astype(np.int32) ) <= 0 ) :
                    self.THREADS = None
                    logger.error( 'Too many threads for the ISO compartments to evaluate; try decreasing the number', prefix='\n' )
            else :
                self.THREADS['ISOt'] = None

        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )


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
            logger.error( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.KERNELS is None :
            logger.error( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first' )
        if self.THREADS is None :
            logger.error( 'Threads not set; call "set_threads()" first' )

        if self.DICTIONARY['IC']['nF'] <= 0 :
            logger.error( 'No streamline found in the dictionary; check your data' )
        if self.DICTIONARY['EC']['nE'] <= 0 and self.KERNELS['wmh'].shape[0] > 0 :
            logger.error( 'The selected model has EC compartments, but no peaks have been provided; check your data' )

        tic = time.time()
        logger.subinfo('')
        logger.info( 'Building linear operator A:' )

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
                    logger.error( '\nbuild_dir is not empty, unsafe build option.' )
                elif config.nTHREADS is not None:
                    logger.error( '\nThe parameter build_dir has changed, unsafe build option.' )
                else:
                    logger.warning( '\nUsing build_dir, always quit your python console between COMMIT Evaluation.' )

            config.nTHREADS   = self.THREADS['n']
            config.model      = self.model.id
            config.nIC        = self.KERNELS['wmr'].shape[0]
            config.nEC        = self.KERNELS['wmh'].shape[0]
            config.nISO       = self.KERNELS['iso'].shape[0]
            config.build_dir  = build_dir

            sys.dont_write_bytecode = True
            pyximport.install( reload_support=True, language_level=3, build_dir=build_dir, build_in_temp=True, inplace=False )

        if 'commit.operator.operator' in sys.modules :
            del sys.modules['commit.operator.operator']
        import commit.operator.operator

        self.A = commit.operator.operator.LinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS )

        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )


    def get_y( self ):
        """
        Returns a numpy array that corresponds to the 'y' vector of the optimisation problem.
        NB: this can be run only after having loaded the dictionary and the data.
        """
        if self.DICTIONARY is None :
            logger.error( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.niiDWI is None :
            logger.error( 'Data not loaded; call "load_data()" first' )

        y = self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float64)
        # y[y < 0] = 0
        return y


    def set_regularisation(self, regularisers=(None, None, None), lambdas=(None, None, None), is_nonnegative=(True, True, True), params=(None, None, None)):
        """
        Set the regularisation parameters for the optimisation problem.

        Parameters
        ----------
        regularisers - tuple :
            sets the penalty term to be used for each compartment:
                regularisers[0] corresponds to the Intracellular compartment
                regularisers[1] corresponds to the Extracellular compartment
                regularisers[2] corresponds to the Isotropic compartment
            Each regularisers[k] must be one of: {None, 'lasso', 'group_lasso', 'sparse_group_lasso'}:
                'lasso' penalises with the 1-norm of the coefficients
                    corresponding to the compartment.
                'group_lasso' penalises according to the following formulation (see [1]):
                    $\Omega(x) = \lambda \sum_{g\in G} w_g |x_g|
                    Considers both the non-overlapping and the hierarchical formulations.
                    NB: this option is allowed only in the IC compartment.
                'sparse_group_lasso' combines lasso and group lasso penalisations
                    NB: this option is allowed only in the IC compartment.
            Default = (None, None, None).

        lambdas - tuple :
            percentage of the maximum regularisation parameter for each compartment:
                lambdas[0] corresponds to the Intracellular compartment
                lambdas[1] corresponds to the Extracellular compartment
                lambdas[2] corresponds to the Isotropic compartment
            The lambdas correspond to the ones described in the mathematical formulation of the regularisation term
            $\Omega(x) = lambdas[0]*regularisers[0](x) + lambdas[1]*regularisers[1](x) + lambdas[2]*regularisers[2](x)$
            The maximum regularisation parameter is the value of lambda above which it is guaranteed that the optimal 
                solution is zero (computed as in [3] for lasso and as in [4] for group lasso).
            NB: if regularisers[k] is None, then lambdas[k] is ignored.
            NB: lambdas[k] must be a float greater than 0.
            NB: if regularisers[0] is 'sparse_group_lasso', then lambdas[k] must be a tuple of two elements,
                the first corresponding to the l1 penalty and the second to the l2 penalty.
            Default = (None, None, None).

        is_nonnegative - tuple :
            impose a non negativity constraint for each compartment:
                is_nonnegative[0] corresponds to the Intracellular compartment
                is_nonnegative[1] corresponds to the Extracellular compartment
                is_nonnegative[2] corresponds to the Isotropic compartment
            Default = (True, True, True).

        params - tuple :
            dictionary of additional parameters for the regularisation term for each compartment:
                params[0] corresponds to the Intracellular compartment
                params[1] corresponds to the Extracellular compartment
                params[2] corresponds to the Isotropic compartment
            Default = (None, None, None).
            Available options:
                'group_idx' - np.array(np.int32) :
                    group indices for the IC compartment (not implemented for EC and ISO).
                    This field is necessary only if regularisers[0] is 'group_lasso' or 'sparse_group_lasso'.
                    Example:
                        structureIC = np.array([[0,2,5],[1,3,4],[0,1,2,3,4,5],[6]], dtype=np.object_)
                        that is equivalent to
                                    [0,1,2,3,4,5]        [6]
                                    /       \
                                [0,2,5]       [1,3,4]
                        which has two non-overlapping groups, one of which is the union
                        of two other non-overlapping groups.
                'group_weights_type' - string :
                    type of group weights to be used for the IC compartment if regularisers[0] is 'group_lasso' 
                        or 'sparse_group_lasso'
                    Available options are: {'standard', 'adaptive'}:
                        'standard' - each group has as weight the square root of the group size.
                        'adaptive' - the weights are computed as in [2].
                    Default = 'adaptive'.
                'group_weights_prior' - np.array(np.float64) :
                    weights associated to each group of the IC compartment, based on prior knowledge.
                    If None, then the weights are computed as in [2], assuming that the fit has already been 
                        performed without regularisation. Otherwise, the provided weights are used as they are.
                    This field can be specified only if regularisers[0] is 'group_lasso' or 'sparse_group_lasso' 
                        and group_weights_type is 'adaptive'.
                    Default = None.
                'coeff_weights' - np.array(np.float64) :
                    weights associated to each individual element of the compartment (implemented for all compartments).
                    This field can be specified only if the chosen regulariser is 'lasso' or 'sparse_group_lasso'.
                    NB: the weights must have the same size as the number of elements in the compartment.

        References:
            [1] Jenatton et al. - 'Proximal Methods for Hierarchical Sparse Coding'
            [2] Schiavi et al. - 'A new method for accurate in vivo mapping of human brain connections using 
                microstructural and anatomical information'
            [3] Kim et al. - 'An interior-point method for large-scale l1-regularized logistic regression'
            [4] Yuan, Lin - 'Model selection and estimation in regression with grouped variables'
        """

        # functions to compute the maximum value of the regularisation parameter (lambda)

        def compute_lambda_max_lasso(start, size, w_coeff): 
            # Ref. Kim et al. - 'An interior-point method for large-scale l1-regularized logistic regression'
            At = self.A.T
            y  = self.get_y()
            Aty = np.asarray(At.dot(y))
            return np.max(np.abs(Aty[start:start+size]) / w_coeff)

        def compute_lambda_max_group(w_group, idx_group): 
            # Ref. Yuan, Lin - 'Model selection and estimation in regression with grouped variables'
            At = self.A.T
            y  = self.get_y()
            Aty = np.asarray(At.dot(y))
            norm_group = np.zeros( w_group.shape, dtype=np.float64 )
            for g in range(w_group.size):
                norm_group[g] = np.sqrt(np.sum(Aty[idx_group[g]]**2)) / w_group[g]
            return np.max(norm_group)
            

        regularisation = {}

        regularisation['startIC']  = 0
        regularisation['sizeIC']   = int( self.DICTIONARY['IC']['nF'] * self.KERNELS['wmr'].shape[0] )
        regularisation['startEC']  = int( regularisation['sizeIC'] )
        regularisation['sizeEC']   = int( self.DICTIONARY['EC']['nE'] * self.KERNELS['wmh'].shape[0] )
        regularisation['startISO'] = int( regularisation['sizeIC'] + regularisation['sizeEC'] )
        regularisation['sizeISO']  = int( self.DICTIONARY['nV'] * self.KERNELS['iso'].shape[0] )

        regularisation['regIC']  = regularisers[0]
        regularisation['regEC']  = regularisers[1]
        regularisation['regISO'] = regularisers[2]

        regularisation['nnIC']  = is_nonnegative[0]
        regularisation['nnEC']  = is_nonnegative[1]
        regularisation['nnISO'] = is_nonnegative[2]

        dictIC_params, dictEC_params, dictISO_params = params

        tr = time.time()
        logger.subinfo('')
        logger.info( 'Setting regularisation:' )
        
        # set regularisation for the intracellular compartment
        if regularisation['regIC'] == 'lasso':
            if lambdas[0] is None:
                logger.error('Missing regularisation parameter for the IC compartment')
            elif lambdas[0] < 0:
                logger.error('Regularisation parameter for the IC compartment must be greater than 0')
            elif lambdas[0] == 0:
                logger.warning('Regularisation parameter for the IC compartment is 0, the solution will be the same as the one without regularisation')
                regularisation['lambdaIC_perc'] = lambdas[0]
            else:
                regularisation['lambdaIC_perc'] = lambdas[0]
        elif regularisation['regIC'] == 'smoothness':
            logger.error('Not yet implemented')
        elif regularisation['regIC'] == 'group_lasso':
            if lambdas[0] is None:
                logger.error('Missing regularisation parameter for the IC compartment')
            elif lambdas[0] < 0:
                logger.error('Regularisation parameter for the IC compartment must be greater than 0')
            elif lambdas[0] == 0:
                logger.warning('Regularisation parameter for the IC compartment is 0, the solution will be the same as the one without regularisation')
                regularisation['lambdaIC_perc'] = lambdas[0]
            else:
                regularisation['lambdaIC_perc'] = lambdas[0]
            if dictIC_params is None:
                logger.error('Dictionary of additional parameters for the IC compartment not provided')
            if dictIC_params['group_idx'] is None:
                logger.error('Group structure for the IC compartment not provided')
            if dictIC_params['group_weights_type'] not in ['standard', 'adaptive']:
                logger.error('Type of group weights not among the available options, i.e. {standard, adaptive}')
        elif regularisation['regIC'] == 'sparse_group_lasso':
            if len(lambdas[0]) != 2:
                logger.error('Regularisation parameters for the IC compartment ara not exactly two')
            elif lambdas[0][0] < 0 or lambdas[0][1] < 0:
                logger.error('Regularisation parameters for the IC compartment must be greater than 0')
            elif lambdas[0][0] == 0 or lambdas[0][1] == 0:
                logger.warning('Regularisation parameters for the IC compartment are both 0, the solution will be the same as the one without regularisation')
                regularisation['lambdaIC_perc'] = lambdas[0]
            else:
                regularisation['lambdaIC_perc'] = lambdas[0]
            if dictIC_params is None:
                logger.error('Dictionary of additional parameters for the IC compartment not provided')
            if dictIC_params['group_idx'] is None:
                logger.error('Group structure for the IC compartment not provided')
            if dictIC_params['group_weights_type'] not in ['standard', 'adaptive']:
                logger.error('Type of group weights not among the available options, i.e. {standard, adaptive}')

        # check if group_weights_prior is consistent with group_weights_type and group_idx
        if regularisation['regIC'] == 'group_lasso' or regularisation['regIC'] == 'sparse_group_lasso':
            if dictIC_params['group_weights_type'] == 'standard' and 'group_weights_prior' in dictIC_params:
                logger.error('Group weights prior knowledge is not allowed for standard group weights')
            if dictIC_params['group_weights_type'] == 'adaptive' and 'group_weights_prior' in dictIC_params:
                if dictIC_params['group_weights_prior'].size != dictIC_params['group_idx'].size:
                    logger.error('Group weights and group indices must have the same size')

        # check if coeff_weights is consistent with the compartment size
        if regularisation['regIC'] == 'lasso' or regularisation['regIC'] == 'sparse_group_lasso':
            if dictIC_params is not None and 'coeff_weights' in dictIC_params:
                if dictIC_params['coeff_weights'].size != len(self.DICTIONARY['TRK']['kept']):
                    logger.error(f'"coeff_weights" must have the same size as the number of elements in the IC compartment (got {dictIC_params["coeff_weights"].size} but {len(self.DICTIONARY["TRK"]["kept"])} expected)')
                dictIC_params['coeff_weights'] = dictIC_params['coeff_weights'][self.DICTIONARY['TRK']['kept']==1]

        # set regularisation for the extracellular compartment
        if regularisation['regEC'] == 'lasso':
            if regularisation['sizeEC'] == 0:
                logger.error('No extracellular compartment found in the dictionary. Unable to set regularisation for the EC compartment.')
            if lambdas[1] is None:
                logger.error('Missing regularisation parameter for the EC compartment')
            elif lambdas[1] < 0:
                logger.error('Regularisation parameter for the EC compartment must be greater than 0')
            elif lambdas[1] == 0:
                logger.warning('Regularisation parameter for the EC compartment is 0, the solution will be the same as the one without regularisation')
                regularisation['lambdaEC_perc'] = lambdas[1]
            else:
                regularisation['lambdaEC_perc'] = lambdas[1]
            if dictEC_params is not None and 'coeff_weights' in dictEC_params:
                if dictEC_params['coeff_weights'].size != regularisation['sizeEC']:
                    logger.error(f'"coeff_weights" must have the same size as the number of elements in the EC compartment (got {dictEC_params["coeff_weights"].size} but {regularisation["sizeEC"]} expected)')
        elif regularisation['regEC'] == 'smoothness':
            logger.error('Not yet implemented')
        elif regularisation['regEC'] == 'group_lasso':
            logger.error('Not yet implemented')
        elif regularisation['regEC'] == 'sparse_group_lasso':
            logger.error('Not yet implemented')

        # set regularisation for the isotropic compartment
        if regularisation['regISO'] == 'lasso':
            if regularisation['sizeISO'] == 0:
                logger.error('No isotropic compartment found in the dictionary. Unable to set regularisation for the ISO compartment.')
            if lambdas[2] is None:
                logger.error('Missing regularisation parameter for the ISO compartment')
            elif lambdas[2] < 0:
                logger.error('Regularisation parameter for the ISO compartment must be greater than 0')
            elif lambdas[2] == 0:
                logger.warning('Regularisation parameter for the ISO compartment is 0, the solution will be the same as the one without regularisation')
                regularisation['lambdaISO_perc'] = lambdas[2]
            else:
                regularisation['lambdaISO_perc'] = lambdas[2]
            if dictISO_params is not None and 'coeff_weights' in dictISO_params:
                if dictISO_params['coeff_weights'].size != regularisation['sizeISO']:
                    logger.error(f'"coeff_weights" must have the same size as the number of elements in the ISO compartment (got {dictISO_params["coeff_weights"].size} but {regularisation["sizeISO"]} expected)')
        elif regularisation['regISO'] == 'smoothness':
            logger.error('Not yet implemented')
        elif regularisation['regISO'] == 'group_lasso':
            logger.error('Not yet implemented')
        elif regularisation['regISO'] == 'sparse_group_lasso':
            logger.error('Not yet implemented')


        # Check if group indices need to be updated in case of 'group_lasso' or 'sparse_group_lasso'
        if regularisation['regIC'] == 'group_lasso' or regularisation['regIC'] == 'sparse_group_lasso' and (0 in self.DICTIONARY['TRK']['kept']) :
            # update the group indices considering only the kept elements
            dictionary_TRK_kept = self.DICTIONARY['TRK']['kept']
            if 'group_weights_prior' in dictIC_params:
                weightsIC_group = dictIC_params['group_weights_prior']
            else:
                weightsIC_group = None
            idx_in_kept = np.zeros(dictionary_TRK_kept.size, dtype=np.int32) - 1  # -1 is used to flag indices for removal
            idx_in_kept[dictionary_TRK_kept==1] = list(range(self.DICTIONARY['IC']['nF']))

            newICgroup_idx = []
            newweightsIC_group = []
            ICgroup_idx = dictIC_params['group_idx']
            for count, group in enumerate(ICgroup_idx):
                group = idx_in_kept[group]
                idx_to_delete = np.where(group==-1)[0]
                if idx_to_delete.size>0:
                    group = np.delete(group,idx_to_delete)
                    if(group.size>0):
                        newICgroup_idx.append(group)
                else:
                    newICgroup_idx.append(group)
                if weightsIC_group is not None and group.size>0:
                    newweightsIC_group.append(weightsIC_group[count])

            dictIC_params['group_idx'] = np.array(newICgroup_idx, dtype=np.object_)
        else:
            if regularisation['regIC'] == 'group_lasso' or regularisation['regIC'] == 'sparse_group_lasso' and 'group_weights_prior' in dictIC_params:
                newweightsIC_group = dictIC_params['group_weights_prior']

        # check if group weights need to be updated in case of 'group_lasso' or 'sparse_group_lasso'
        if regularisation['regIC'] == 'group_lasso' or regularisation['regIC'] == 'sparse_group_lasso':
            # set the group weights
            if dictIC_params['group_weights_type'] is 'standard':
                group_size = np.array([g.size for g in dictIC_params['group_idx']], dtype=np.int32)
                dictIC_params['group_weights'] = np.sqrt(group_size)
            elif dictIC_params['group_weights_type'] is 'adaptive':
                if 'group_weights_prior' not in dictIC_params: # default weights (both cardinality and x_nnls, like in wiki)
                    # check if fit has been performed
                    if self.x is None:
                        logger.error('Group weights cannot be computed if the fit (without regularisation) has not been performed yet')
                    group_size = np.array([g.size for g in dictIC_params['group_idx']], dtype=np.int32)
                    x_nnls, _, _ = self.get_coeffs(get_normalized=False)
                    group_x_norm = np.array([np.linalg.norm(x_nnls[g])+1e-12 for g in dictIC_params['group_idx']], dtype=np.float64)
                    dictIC_params['group_weights'] = np.sqrt(group_size) / group_x_norm
                else:
                    dictIC_params['group_weights'] = np.array(newweightsIC_group)
            else:
                logger.error('Type of group weights not among the available options, i.e. {standard, adaptive}')

        regularisation['dictIC_params']  = dictIC_params
        regularisation['dictEC_params']  = dictEC_params
        regularisation['dictISO_params'] = dictISO_params

        # update lambdas using lambda_max
        if regularisation['regIC'] == 'lasso':
            if dictIC_params is not None and 'coeff_weights' in dictIC_params:
                regularisation['lambdaIC_max'] = compute_lambda_max_lasso(regularisation['startIC'], regularisation['sizeIC'], dictIC_params['coeff_weights'])
            else:
                regularisation['lambdaIC_max'] = compute_lambda_max_lasso(regularisation['startIC'], regularisation['sizeIC'], np.ones(regularisation['sizeIC'], dtype=np.float64))
            # regularisation['lambdaIC_max'] = compute_lambda_max_lasso(regularisation['startIC'], regularisation['sizeIC'])
            regularisation['lambdaIC'] = regularisation['lambdaIC_perc'] * regularisation['lambdaIC_max']
        if regularisation['regIC'] == 'group_lasso':
            regularisation['lambdaIC_max'] = compute_lambda_max_group(dictIC_params['group_weights'], dictIC_params['group_idx'])
            regularisation['lambdaIC'] = regularisation['lambdaIC_perc'] * regularisation['lambdaIC_max']
        if regularisation['regIC'] == 'sparse_group_lasso':
            regularisation['lambdaIC_max'] = ( compute_lambda_max_lasso(regularisation['startIC'], regularisation['sizeIC']), compute_lambda_max_group(dictIC_params['group_weights'], dictIC_params['group_idx']) )
            regularisation['lambdaIC'] = ( regularisation['lambdaIC_perc'][0] * regularisation['lambdaIC_max'][0], regularisation['lambdaIC_perc'][1] * regularisation['lambdaIC_max'][1] )
        if regularisation['regEC'] == 'lasso':
            if dictEC_params is not None and 'coeff_weights' in dictEC_params:
                regularisation['lambdaEC_max'] = compute_lambda_max_lasso(regularisation['startEC'], regularisation['sizeEC'], dictEC_params['coeff_weights'])
            else:
                regularisation['lambdaEC_max'] = compute_lambda_max_lasso(regularisation['startEC'], regularisation['sizeEC'], np.ones(regularisation['sizeEC'], dtype=np.float64))
            # regularisation['lambdaEC_max'] = compute_lambda_max_lasso(regularisation['startEC'], regularisation['sizeEC'])
            regularisation['lambdaEC'] = regularisation['lambdaEC_perc'] * regularisation['lambdaEC_max']
        if regularisation['regISO'] == 'lasso':
            if dictISO_params is not None and 'coeff_weights' in dictISO_params:
                regularisation['lambdaISO_max'] = compute_lambda_max_lasso(regularisation['startISO'], regularisation['sizeISO'], dictISO_params['coeff_weights'])
            else:
                regularisation['lambdaISO_max'] = compute_lambda_max_lasso(regularisation['startISO'], regularisation['sizeISO'], np.ones(regularisation['sizeISO'], dtype=np.float64))
            # regularisation['lambdaISO_max'] = compute_lambda_max_lasso(regularisation['startISO'], regularisation['sizeISO'])
            regularisation['lambdaISO'] = regularisation['lambdaISO_perc'] * regularisation['lambdaISO_max']

        self.regularisation_params = commit.solvers.init_regularisation(regularisation)

        logger.subinfo( 'IC compartment:', indent_char='*' )
        if (regularisation['regIC'] == 'lasso' or regularisation['regIC'] == 'sparse_group_lasso') and dictIC_params is not None and 'coeff_weights' in dictIC_params:
                logger.subinfo( f'Regularisation type: {regularisation["regIC"]} (weighted version)', indent_lvl=1, indent_char='-' )
        else:
            logger.subinfo( f'Regularisation type: {regularisation["regIC"]}', indent_lvl=1, indent_char='-' )
        logger.subinfo( f'Non-negativity constraint: {regularisation["nnIC"]}', indent_lvl=1, indent_char='-' )
        if regularisation['regIC'] is not None:
            logger.subinfo( f'Lambda max: {regularisation["lambdaIC_max"]}', indent_lvl=1, indent_char='-' )
            logger.subinfo( f'% lambda: {regularisation["lambdaIC_perc"]}', indent_lvl=1, indent_char='-' )
            logger.subinfo( f'Lambda used: {regularisation["lambdaIC"]}', indent_lvl=1, indent_char='-' )
        if regularisation['regIC'] == 'group_lasso' or regularisation['regIC'] == 'sparse_group_lasso':
            logger.subinfo( f'Number of groups: {len(dictIC_params["group_idx"])}', indent_lvl=1, indent_char='-' )
            if 'group_weights_prior' in dictIC_params:
                logger.subinfo( f'Type of group weights: {dictIC_params["group_weights_type"]}, with prior knowledge', indent_lvl=1, indent_char='-' )
            else:
                logger.subinfo( f'Type of group weights: {dictIC_params["group_weights_type"]}', indent_lvl=1, indent_char='-' )

        logger.subinfo( 'EC compartment:', indent_char='*' )
        if regularisation['regEC'] == 'lasso' and dictEC_params is not None and 'coeff_weights' in dictEC_params:
            logger.subinfo( f'Regularisation type: {regularisation["regEC"]} (weighted version)', indent_lvl=1, indent_char='-' )
        else:
            logger.subinfo( f'Regularisation type: {regularisation["regEC"]}', indent_lvl=1, indent_char='-' )
        logger.subinfo( f'Non-negativity constraint: {regularisation["nnEC"]}', indent_lvl=1, indent_char='-' )
        if regularisation['regEC'] is not None:
            logger.subinfo( f'Lambda max: {regularisation["lambdaEC_max"]}', indent_lvl=1, indent_char='-' )
            logger.subinfo( f'% lambda: {regularisation["lambdaEC_perc"]}', indent_lvl=1, indent_char='-' )
            logger.subinfo( f'Lambda used: {regularisation["lambdaEC"]}', indent_lvl=1, indent_char='-' )

        logger.subinfo( 'ISO compartment:', indent_char='*' )
        if regularisation['regISO'] == 'lasso' and dictISO_params is not None and 'coeff_weights' in dictISO_params:
            logger.subinfo( f'Regularisation type: {regularisation["regISO"]} (weighted version)', indent_lvl=1, indent_char='-' )
        else:
            logger.subinfo( f'Regularisation type: {regularisation["regISO"]}', indent_lvl=1, indent_char='-' )
        logger.subinfo( f'Non-negativity constraint: {regularisation["nnISO"]}', indent_lvl=1, indent_char='-' )
        if regularisation['regISO'] is not None: 
            logger.subinfo( f'Lambda max: {regularisation["lambdaISO_max"]}', indent_lvl=1, indent_char='-' )
            logger.subinfo( f'% lambda: {regularisation["lambdaISO_perc"]}', indent_lvl=1, indent_char='-' )
            logger.subinfo( f'Lambda used: {regularisation["lambdaISO"]}', indent_lvl=1, indent_char='-' )

        logger.subinfo( f'[ {time.time() - tr:.1f} seconds ]' )


    def fit( self, tol_fun=1e-3, tol_x=1e-6, max_iter=100, x0=None, confidence_map_filename=None, confidence_map_rescale=False ) :
        """Fit the model to the data.

        Parameters
        ----------
        tol_fun : float
            Tolerance on the objective function (default : 1e-3)
        max_iter : integer
            Maximum number of iterations (default : 100)
        x0 : np.array
            Initial guess for the solution of the problem (default : None)
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
            logger.error( 'Data not loaded; call "load_data()" first' )
        if self.DICTIONARY is None :
            logger.error( 'Dictionary not loaded; call "load_dictionary()" first' )
        if self.KERNELS is None :
            logger.error( 'Response functions not generated; call "generate_kernels()" and "load_kernels()" first' )
        if self.THREADS is None :
            logger.error( 'Threads not set; call "set_threads()" first' )
        if self.A is None :
            logger.error( 'Operator not built; call "build_operator()" first' )

        if self.regularisation_params is None:
            self.set_regularisation()

        logger.subinfo('')
        logger.info( 'Fit model:' )

        # Confidence map
        self.confidence_map_img = None
        self.set_config('confidence_map_filename', None)
        confidence_array = None
        confidence_array_changed = False

        if confidence_map_filename is not None:
            # Loading confidence map
            tic = time.time()
            logger.subinfo( 'Loading confidence map:' )

            if not exists( pjoin( self.get_config('DATA_path'), confidence_map_filename)  ) :
                logger.error( 'Confidence map not found' )

            self.set_config('confidence_map_filename', confidence_map_filename)
            confidence_map  = nibabel.load( pjoin( self.get_config('DATA_path'), confidence_map_filename) )
            self.confidence_map_img = np.asanyarray( confidence_map.dataobj ).astype(np.float32)

            if self.confidence_map_img.ndim not in [3,4]:
                logger.error( 'Confidence map must be 3D or 4D dataset' )

            if self.confidence_map_img.ndim == 3:
                logger.subinfo('Extending the confidence map volume to match the DWI signal volume(s)... ', indent_lvl=1, indent_char='-' )
                self.confidence_map_img = np.repeat(self.confidence_map_img[:, :, :, np.newaxis], self.niiDWI_img.shape[3], axis=3)
            hdr = confidence_map.header if nibabel.__version__ >= '2.0.0' else confidence_map.get_header()
            confidence_map_dim = self.confidence_map_img.shape[0:3]
            confidence_map_pixdim = tuple( hdr.get_zooms()[:3] )
            logger.subinfo('dim    : %d x %d x %d x %d' % self.confidence_map_img.shape, indent_lvl=2, indent_char='-')
            logger.subinfo('pixdim : %.3f x %.3f x %.3f' % confidence_map_pixdim, indent_lvl=2, indent_char='-')

            logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )

            if ( self.get_config('dim') != confidence_map_dim ):
                logger.error( 'Dataset does not have the same geometry (number of voxels) as the DWI signal' )

            if (self.get_config('pixdim') != confidence_map_pixdim ):
                logger.error( 'Dataset does not have the same geometry (voxel size) as the DWI signal' )

            if (self.confidence_map_img.shape != self.niiDWI_img.shape):
                logger.error( 'Dataset does not have the same geometry as the DWI signal' )

            confidence_array = self.confidence_map_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32)

            if(np.isnan(confidence_array).any()):
                confidence_array[np.isnan(confidence_array)] = 0.
                confidence_array_changed = True
                logger.warning('Confidence map contains NaNs. Those values were changed to 0.')

            cMAX = np.max(confidence_array)
            cMIN = np.min(confidence_array)
            if(cMIN == cMAX):
                self.confidence_map_img = None
                self.set_config('confidence_map_filename', None)
                confidence_array = None
                confidence_array_changed = False
                logger.warning('All voxels in the confidence map have the same value. The confidence map will not be used')

            elif(confidence_map_rescale):
                confidence_array = ( confidence_array - cMIN ) / ( cMAX - cMIN )
                logger.subinfo ( '[Confidence map interval was scaled from the original [%.1f, %.1f] to the intended [%.1f, %.1f] linearly]' % ( cMIN, cMAX, np.min(confidence_array), np.max(confidence_array) ), indent_lvl=2 )
                confidence_array_changed = True

            if(confidence_array_changed):
                nV = self.DICTIONARY['nV']
                self.confidence_map_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ] = np.reshape( confidence_array, (nV,-1) ).astype(np.float32)

        if x0 is not None :
            if x0.shape[0] != self.A.shape[1] :
                logger.error( 'x0 dimension does not match the number of columns of the dictionary' )


        self.CONFIG['optimization']                   = {}
        self.CONFIG['optimization']['tol_fun']        = tol_fun
        self.CONFIG['optimization']['tol_x']          = tol_x
        self.CONFIG['optimization']['max_iter']       = max_iter
        self.CONFIG['optimization']['regularisation'] = self.regularisation_params

        # run solver
        t = time.time()
        with ProgressBar(disable=self.verbose!=3, hide_on_exit=True) as pb:
            self.x, opt_details = commit.solvers.solve(self.get_y(), self.A, self.A.T, tol_fun=tol_fun, tol_x=tol_x, max_iter=max_iter, verbose=self.verbose, x0=x0, regularisation=self.regularisation_params, confidence_array=confidence_array)

        self.CONFIG['optimization']['fit_details'] = opt_details
        self.CONFIG['optimization']['fit_time'] = time.time()-t

        logger.subinfo( f'[ {time.strftime("%Hh %Mm %Ss", time.gmtime(self.CONFIG["optimization"]["fit_time"]) )} ]' )


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
            logger.error( 'Model not fitted to the data; call "fit()" first' )

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


    def save_results( self, path_suffix=None, coeffs_format='%.5e', stat_coeffs='sum', save_est_dwi=False, do_reweighting=True ) :
        """Save the output (coefficients, errors, maps etc).

        Parameters
        ----------
        path_suffix : string
            Text to be appended to "Results" to create the output path (default : None)
        stat_coeffs : string
            Stat to be used if more coefficients are estimated for each streamline.
            Options: 'sum', 'mean', 'median', 'min', 'max', 'all' (default : 'sum')
        coeffs_format : string
            Format for saving the coefficients to `streamline_weights.txt` (default: '%.5e')
        save_est_dwi : boolean
            Save the estimated DW-MRI signal (default : False)
        """
        RESULTS_path = 'Results_' + self.model.id
        if path_suffix :
            self.set_config('path_suffix', path_suffix)
            RESULTS_path = RESULTS_path + path_suffix

        logger.subinfo('')
        logger.info( 'Saving results to "%s/*":' % RESULTS_path )
        tic = time.time()

        if self.x is None :
            logger.error( 'Model not fitted to the data; call "fit()" first' )

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
        logger.subinfo('Fitting errors:', indent_lvl=2, indent_char='-' )

        niiMAP_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
        affine = self.niiDWI.affine if nibabel.__version__ >= '2.0.0' else self.niiDWI.get_affine()
        niiMAP     = nibabel.Nifti1Image( niiMAP_img, affine )
        niiMAP_hdr = niiMAP.header if nibabel.__version__ >= '2.0.0' else niiMAP.get_header()
        niiMAP_hdr['descrip'] = 'Created with COMMIT %s'%self.get_config('version')

        y_mea = np.reshape( self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32), (nV,-1) )
        y_est = np.reshape( self.A.dot(self.x), (nV,-1) ).astype(np.float32)

        tmp = np.sqrt( np.mean((y_mea-y_est)**2,axis=1) )
        logger.subinfo(f'RMSE:  {tmp.mean():.3f} +/- {tmp.std():.3f}', indent_lvl=2, indent_char='-')
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP_hdr['cal_min'] = 0
        niiMAP_hdr['cal_max'] = tmp.max()
        nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_RMSE.nii.gz') )

        tmp = np.sum(y_mea**2,axis=1)
        idx = np.where( tmp < 1E-12 )
        tmp[ idx ] = 1
        tmp = np.sqrt( np.sum((y_mea-y_est)**2,axis=1) / tmp )
        tmp[ idx ] = 0
        logger.subinfo(f'NRMSE: {tmp.mean():.3f} +/- {tmp.std():.3f}', indent_lvl=2, indent_char='-')
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP_hdr['cal_min'] = 0
        niiMAP_hdr['cal_max'] = 1
        nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_NRMSE.nii.gz') )

        if self.confidence_map_img is not None:
            confidence_array = np.reshape( self.confidence_map_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32), (nV,-1) )

            tmp = np.sum(confidence_array,axis=1)
            idx = np.where( tmp < 1E-12 )
            tmp[ idx ] = 1
            tmp = np.sqrt( np.sum(confidence_array*(y_mea-y_est)**2,axis=1) / tmp )
            tmp[ idx ] = 0
            logger.subinfo(f'RMSE considering the confidence map:  {tmp.mean():.3f} +/- {tmp.std():.3f}', indent_lvl=2, indent_char='-')
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
            niiMAP_hdr['cal_min'] = 0
            niiMAP_hdr['cal_max'] = tmp.max()
            nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_RMSE_adjusted.nii.gz') )

            tmp = np.sum(confidence_array*y_mea**2,axis=1)
            idx = np.where( tmp < 1E-12 )
            tmp[ idx ] = 1
            tmp = np.sqrt( np.sum(confidence_array*(y_mea-y_est)**2,axis=1) / tmp )
            tmp[ idx ] = 0
            logger.subinfo(f'NRMSE considering the confidence map: {tmp.mean():.3f} +/- {tmp.std():.3f}', indent_lvl=2, indent_char='-')
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
            niiMAP_hdr['cal_min'] = 0
            niiMAP_hdr['cal_max'] = 1
            nibabel.save( niiMAP, pjoin(RESULTS_path,'fit_NRMSE_adjusted.nii.gz') )
            confidence_array = None

        # Map of compartment contributions
        logger.subinfo('Voxelwise contributions:', indent_char='*')

        logger.subinfo('Intra-axonal', indent_lvl=2, indent_char='-', with_progress=True)
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            niiIC_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
            if len(self.KERNELS['wmr']) > 0 :
                offset = nF * self.KERNELS['wmr'].shape[0]
                tmp = ( x[:offset].reshape( (-1,nF) ) * norm_fib.reshape( (-1,nF) ) ).sum( axis=0 )
                xv = np.bincount( self.DICTIONARY['IC']['v'], minlength=nV,
                    weights=tmp[ self.DICTIONARY['IC']['fiber'] ] * self.DICTIONARY['IC']['len']
                ).astype(np.float32)
                niiIC_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv

        logger.subinfo('Extra-axonal', indent_lvl=2, indent_char='-', with_progress=True)
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            niiEC_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
            if len(self.KERNELS['wmh']) > 0 :
                offset = nF * self.KERNELS['wmr'].shape[0]
                tmp = x[offset:offset+nE*len(self.KERNELS['wmh'])].reshape( (-1,nE) ).sum( axis=0 )
                xv = np.bincount( self.DICTIONARY['EC']['v'], weights=tmp, minlength=nV ).astype(np.float32)
                niiEC_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv

        logger.subinfo('Isotropic   ', indent_lvl=2, indent_char='-', with_progress=True)
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            niiISO_img = np.zeros( self.get_config('dim'), dtype=np.float32 )
            if len(self.KERNELS['iso']) > 0 :
                offset = nF * self.KERNELS['wmr'].shape[0] + nE * self.KERNELS['wmh'].shape[0]
                offset_iso = offset + self.DICTIONARY['ISO']['nV']
                tmp = x[offset:offset_iso].reshape( (-1,self.DICTIONARY['ISO']['nV']) ).sum( axis=0 )
                xv = np.bincount( self.DICTIONARY['ISO']['v'], weights=tmp, minlength=nV ).astype(np.float32)
                niiISO_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv

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
        logger.subinfo('Configuration and results:', indent_char='*')
        logger.subinfo('Saving streamline_weights.txt', indent_lvl=2, indent_char='-', with_progress=True)
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
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
                    logger.error( 'Stat not allowed. Possible values: sum, mean, median, min, max, all', prefix='\n' )

            # scale output weights if blur was used
            dictionary_info = load_dictionary_info( pjoin(self.get_config('TRACKING_path'), 'dictionary_info.pickle') )
            if dictionary_info['blur_gauss_extent'] > 0 or dictionary_info['blur_core_extent'] > 0 :
                if stat_coeffs == 'all' :
                    logger.error( 'Not yet implemented. Unable to account for blur in case of multiple streamline constributions.' )
            if "tractogram_centr_idx" in dictionary_info.keys():
                ordered_idx = dictionary_info["tractogram_centr_idx"].astype(np.int64)
                unravel_weights = np.zeros( dictionary_info['n_count'], dtype=np.float64)
                unravel_weights[ordered_idx] = self.DICTIONARY['TRK']['kept'].astype(np.float64)
                temp_weights = unravel_weights[ordered_idx] 
                if dictionary_info['blur_gauss_extent'] > 0 or dictionary_info['blur_core_extent'] > 0:
                    temp_weights[temp_weights>0] = xic[self.DICTIONARY['TRK']['kept']>0] * self.DICTIONARY['TRK']['lenTot'] / self.DICTIONARY['TRK']['len']
                    unravel_weights[ordered_idx] = temp_weights
                    xic = unravel_weights
                else:
                    temp_weights[temp_weights>0] = xic[self.DICTIONARY['TRK']['kept']>0]
                    unravel_weights[ordered_idx] = temp_weights
                    xic = unravel_weights

            else:
                if dictionary_info['blur_gauss_extent'] > 0 or dictionary_info['blur_core_extent'] > 0:
                    xic[ self.DICTIONARY['TRK']['kept']==1 ] *= self.DICTIONARY['TRK']['lenTot'] / self.DICTIONARY['TRK']['len']

            
            self.temp_data['DICTIONARY'] = self.DICTIONARY
            self.temp_data['niiIC_img'] = niiIC_img
            self.temp_data['niiEC_img'] = niiEC_img
            self.temp_data['niiISO_img'] = niiISO_img
            self.temp_data['streamline_weights'] = xic
            self.temp_data['RESULTS_path'] = RESULTS_path

            if hasattr(self.model, '_postprocess') and do_reweighting:
                self.model._postprocess(self.temp_data, verbose=self.verbose)

            np.savetxt( pjoin(RESULTS_path,'streamline_weights.txt'), xic, fmt=coeffs_format )
            self.set_config('stat_coeffs', stat_coeffs)

        # Save to a pickle file the following items:
        #   item 0: dictionary with all the configuration details
        #   item 1: np.array obtained through the optimisation process with the normalised kernels
        #   item 2: np.array renormalisation of coeffs in item 1
        logger.subinfo('results.pickle:              ', indent_char='-', indent_lvl=2, with_progress=True)
        with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:
            xic, xec, xiso = self.get_coeffs()
            x = self.x
            if self.get_config('doNormalizeKernels') :
                x = x * np.hstack( (norm1*norm_fib,norm2,norm3) )
            with open( pjoin(RESULTS_path,'results.pickle'), 'wb+' ) as fid :
                self.CONFIG['optimization']['regularisation'].pop('omega', None)
                self.CONFIG['optimization']['regularisation'].pop('prox', None)
                pickle.dump( [self.CONFIG, x, self.x], fid, protocol=2 )

        if save_est_dwi :
            logger.subinfo('Estimated signal:', indent_char='-', indent_lvl=2, with_progress=True)
            with ProgressBar(disable=self.verbose < 3, hide_on_exit=True, subinfo=True) as pbar:                    
                self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ] = y_est
                nibabel.save( nibabel.Nifti1Image( self.niiDWI_img , affine ), pjoin(RESULTS_path,'fit_signal_estimated.nii.gz') )
                self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ] = y_mea

        logger.subinfo( f'[ {(time.time() - tic):.1f} seconds ]' )
