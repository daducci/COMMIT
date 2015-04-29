import numpy as np
import re
import time
import glob
import sys
import os.path
import nibabel
import cPickle
from dipy.data.fetcher import dipy_home
from dipy.core.geometry import cart2sphere
from dipy.reconst.shm import real_sym_sh_basis
import commit.models
import commit.solvers
import pyximport
pyximport.install( reload_support=True )


class Evaluation :
    """
    Class to hold all the information (data and parameters) when performing an
    evaluation with the COMMIT framework.
    """

    def __init__( self, study_path, subject ) :
        """
        Setup the data structure with default values.

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
        self.CONFIG['study_path'] = study_path
        self.CONFIG['subject']    = subject
        self.CONFIG['DATA_path']  = os.path.join( study_path, subject )

        self.CONFIG['doNormalizeSignal']  = True
        self.CONFIG['doMergeB0']	      = True
        self.CONFIG['doNormalizeKernels'] = True
        self.CONFIG['doDemean']		      = False


    def load_data( self, dwi_filename = 'DWI.nii', scheme_filename = 'DWI.scheme', b0_thr = 0 ) :
        """
        Load the diffusion signal and its corresponding acquisition scheme.

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
        self.CONFIG['dwi_filename']    = dwi_filename
        self.niiDWI  = nibabel.load( os.path.join( self.CONFIG['DATA_path'], dwi_filename) )
        self.niiDWI_img = self.niiDWI.get_data().astype(np.float32)
        self.CONFIG['dim']    = self.niiDWI_img.shape[0:3]
        self.CONFIG['pixdim'] = tuple( self.niiDWI.get_header().get_zooms()[:3] )
        print '\t\t- dim    = %d x %d x %d x %d' % self.niiDWI_img.shape
        print '\t\t- pixdim = %.3f x %.3f x %.3f' % self.CONFIG['pixdim']

        print '\t* Acquisition scheme...'
        self.CONFIG['scheme_filename'] = scheme_filename
        self.CONFIG['b0_thr'] = b0_thr
        self.scheme = Scheme( os.path.join( self.CONFIG['DATA_path'], scheme_filename), b0_thr )
        print '\t\t- %d samples, %d shells' % ( self.scheme.nS, len(self.scheme.shells) )
        print '\t\t- %d @ b=0' % ( self.scheme.b0_count ),
        for i in xrange(len(self.scheme.shells)) :
            print ', %d @ b=%.1f' % ( len(self.scheme.shells[i]['idx']), self.scheme.shells[i]['b'] ),
        print

        if self.scheme.nS != self.niiDWI_img.shape[3] :
            raise ValueError( 'Scheme does not match with DWI data' )

        print '   [ %.1f seconds ]' % ( time.time() - tic )

        # Preprocessing
        tic = time.time()
        print '\n-> Preprocessing:'

        if self.CONFIG['doNormalizeSignal'] :
            print '\t* Normalizing to b0...',
            sys.stdout.flush()
            mean = np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 )
            idx = mean < 1
            mean[ idx ] = 1
            mean = 1 / mean
            mean[ idx ] = 0
            for i in xrange(self.scheme.nS) :
                self.niiDWI_img[:,:,:,i] *= mean
            print '[ min=%.2f,  mean=%.2f, max=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.mean(), self.niiDWI_img.max() )

        if self.CONFIG['doMergeB0'] :
            print '\t* Merging multiple b0 volume(s)...',
            mean = np.expand_dims( np.mean( self.niiDWI_img[:,:,:,self.scheme.b0_idx], axis=3 ), axis=3 )
            self.niiDWI_img = np.concatenate( (self.niiDWI_img[:,:,:,self.scheme.dwi_idx], mean), axis=3 )
        else :
            print '\t* Keeping all b0 volume(s)...',
        print '[ %d x %d x %d x %d ]' % self.niiDWI_img.shape

        if self.CONFIG['doDemean'] :
            print '\t* Demeaning signal...',
            sys.stdout.flush()
            mean = np.repeat( np.expand_dims(np.mean(self.niiDWI_img,axis=3),axis=3), self.niiDWI_img.shape[3], axis=3 )
            self.niiDWI_img = self.niiDWI_img - mean
            print '[ min=%.2f,  mean=%.2f, max=%.2f ]' % ( self.niiDWI_img.min(), self.niiDWI_img.mean(), self.niiDWI_img.max() )

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def set_model( self, model_name ) :
        """
        Set the model to use to describe the signal contributions in each voxel.

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

        self.CONFIG['ATOMS_path'] = os.path.join( self.CONFIG['study_path'], 'kernels', self.model.id )


    def generate_kernels( self, regenerate = False, lmax = 12 ) :
        """
        Generate the high-resolution response functions for each compartment.
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

        print '\n-> Simulating with "%s" model:' % self.model.name

        # check if kernels were already generated
        tmp = glob.glob( os.path.join(self.CONFIG['ATOMS_path'],'A_*.npy') )
        if len(tmp)>0 and not regenerate :
            print '   [ Kernels already computed. Call "generate_kernels( regenerate=True )" to force regeneration. ]'
            return

        # create folder or delete existing files (if any)
        if not os.path.exists( self.CONFIG['ATOMS_path'] ) :
            os.makedirs( self.CONFIG['ATOMS_path'] )
        else :
            for f in glob.glob( os.path.join(self.CONFIG['ATOMS_path'],'*') ) :
                os.remove( f )

        # load precomputed rotation matrices
        self.CONFIG['lmax'] = lmax
        filename = os.path.join( dipy_home, 'COMMIT_aux_matrices_lmax=%d.pickle'%lmax )
        if not os.path.exists( filename ) :
            raise RuntimeError( 'Auxiliary matrices not found; call "precompute_rotation_matrices()" first.' )
        aux = cPickle.load( open(filename,'rb') )

        # calculate indices within structures
        nSH = (lmax+1)*(lmax+2)/2
        idx_IN  = []
        idx_OUT = []
        for s in xrange( len(self.scheme.shells) ) :
            idx_IN.append( range(500*s,500*(s+1)) )
            idx_OUT.append( range(nSH*s,nSH*(s+1)) )

        # Dispatch to the right handler for each model
        tic = time.time()
        self.model.generate( self.CONFIG['ATOMS_path'], self.scheme, aux, idx_IN, idx_OUT )
        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def load_kernels( self ) :
        """
        Load rotated kernels and project to the specific gradient scheme of this subject.
        Dispatch to the proper function, depending on the model.
        """
        if self.model is None :
            raise RuntimeError( 'Model not set; call "set_model()" method first.' )
        if self.scheme is None :
            raise RuntimeError( 'Scheme not loaded; call "load_data()" first.' )

        tic = time.time()
        print '\n-> Resampling kernels for subject "%s":' % self.CONFIG['subject']

        # load precomputed rotation matrices
        lmax = self.CONFIG['lmax']
        filename = os.path.join( dipy_home, 'COMMIT_aux_matrices_lmax=%d.pickle'%lmax )
        if not os.path.exists( filename ) :
            raise RuntimeError( 'Auxiliary matrices not found; call "precompute_rotation_matrices()" first.' )
        aux = cPickle.load( open(filename,'rb') )

        # calculate auxiliary structures
        nSH = (lmax+1)*(lmax+2)/2
        idx_OUT = np.zeros( self.scheme.dwi_count, dtype=np.int32 )
        Ylm_OUT = np.zeros( (self.scheme.dwi_count,nSH*len(self.scheme.shells)), dtype=np.float32 ) # matrix from SH to real space
        idx = 0
        for s in xrange( len(self.scheme.shells) ) :
            nS = len( self.scheme.shells[s]['idx'] )
            idx_OUT[ idx:idx+nS ] = self.scheme.shells[s]['idx']
            _, theta, phi = cart2sphere( self.scheme.shells[s]['grad'][:,0], self.scheme.shells[s]['grad'][:,1], self.scheme.shells[s]['grad'][:,2] )
            tmp, _, _ = real_sym_sh_basis( lmax, theta, phi )
            Ylm_OUT[ idx:idx+nS, nSH*s:nSH*(s+1) ] = tmp
            idx += nS

        # Dispatch to the right handler for each model
        self.KERNELS = self.model.resample( self.CONFIG['ATOMS_path'], idx_OUT, Ylm_OUT )
        nIC  = self.KERNELS['wmr'].shape[0]
        nEC  = self.KERNELS['wmh'].shape[0]
        nISO = self.KERNELS['iso'].shape[0]

        # Remove multiple b0(s)
        if self.CONFIG['doMergeB0'] :
            print '\t* Merging multiple b0 volume(s)...',
            ones = np.expand_dims( np.ones(self.KERNELS['wmr'].shape[0:3],dtype=np.float32), axis=3 )
            self.KERNELS['wmr'] = np.concatenate( (self.KERNELS['wmr'][:,:,:,self.scheme.dwi_idx], ones), axis=3 )
            ones = np.expand_dims( np.ones(self.KERNELS['wmh'].shape[0:3],dtype=np.float32), axis=3 )
            self.KERNELS['wmh'] = np.concatenate( (self.KERNELS['wmh'][:,:,:,self.scheme.dwi_idx], ones), axis=3 )
            ones = np.expand_dims( np.ones(self.KERNELS['iso'].shape[0:1],dtype=np.float32), axis=1 )
            self.KERNELS['iso'] = np.concatenate( (self.KERNELS['iso'][:,self.scheme.dwi_idx], ones), axis=1 )
        else :
            print '\t* Keeping all b0 volume(s)...',

        # ensure contiguous arrays for C part
        self.KERNELS['wmr'] = np.ascontiguousarray( self.KERNELS['wmr'] )
        self.KERNELS['wmh'] = np.ascontiguousarray( self.KERNELS['wmh'] )
        self.KERNELS['iso'] = np.ascontiguousarray( self.KERNELS['iso'] )

        print '[ OK ]'

        # De-mean kernels
        if self.CONFIG['doDemean'] :
            print '\t* Demeaning signal...',
            for j in xrange(181) :
                for k in xrange(181) :
                    for i in xrange(nIC) :
                        self.KERNELS['wmr'][i,j,k,:] -= self.KERNELS['wmr'][i,j,k,:].mean()
                    for i in xrange(nEC) :
                        print nEC
                        self.KERNELS['wmh'][i,j,k,:] -= self.KERNELS['wmh'][i,j,k,:].mean()
            for i in xrange(nISO) :
                self.KERNELS['iso'][i] -= self.KERNELS['iso'][i].mean()
            print '[ OK ]'

        # Normalize atoms
        if self.CONFIG['doNormalizeKernels'] :
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


    def load_dictionary( self, path ) :
        """
        Load the sparse structure previously created with "trk2dictionary" script.

        Parameters
        ----------
        path : string
            Folder containing the output of the trk2dictionary script (relative to subject path)
        """
        if self.niiDWI is None :
            raise RuntimeError( 'Data not loaded; call "load_data()" first.' )

        tic = time.time()
        print '\n-> Loading the dictionary:'
        self.DICTIONARY = {}

        self.CONFIG['TRACKING_path'] = os.path.join(self.CONFIG['DATA_path'],path)
        mask_filename = os.path.join(self.CONFIG['TRACKING_path'],'dictionary_tdi.nii')
        if not os.path.exists( mask_filename ) :
            mask_filename += '.gz'
            if not os.path.exists( mask_filename ) :
                raise RuntimeError( 'Dictionary not found. Execute ''trk2dictionary'' script first.' );

        niiMASK = nibabel.load( mask_filename )
        self.DICTIONARY['MASK'] = (niiMASK.get_data() > 0).astype(np.uint8)

        # segments from the tracts
        # ------------------------
        print '\t* segments from the tracts...',
        sys.stdout.flush()

        self.DICTIONARY['IC'] = {}

        self.DICTIONARY['IC']['trkLen'] = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_trkLen.dict'), dtype=np.float32 )
        self.DICTIONARY['IC']['nF'] = self.DICTIONARY['IC']['trkLen'].size

        self.DICTIONARY['IC']['fiber'] = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_f.dict'), dtype=np.uint32 )
        self.DICTIONARY['IC']['n'] = self.DICTIONARY['IC']['fiber'].size

        vx = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_vx.dict'), dtype=np.uint8 ).astype(np.uint32)
        vy = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_vy.dict'), dtype=np.uint8 ).astype(np.uint32)
        vz = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_vz.dict'), dtype=np.uint8 ).astype(np.uint32)
        self.DICTIONARY['IC']['v'] = vx + self.CONFIG['dim'][0] * ( vy + self.CONFIG['dim'][1] * vz )
        del vx, vy, vz

        ox = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_ox.dict'), dtype=np.uint8 ).astype(np.uint16)
        oy = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_oy.dict'), dtype=np.uint8 ).astype(np.uint16)
        self.DICTIONARY['IC']['o'] = oy + 181*ox
        del ox, oy

        self.DICTIONARY['IC']['len'] = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_IC_len.dict'), dtype=np.float32 )

        if self.CONFIG['doNormalizeKernels'] :
            # divide the length of each segment by the fiber length so that all the columns of the libear operator will have same length
            # NB: it works in conjunction with the normalization of the kernels
            sl = self.DICTIONARY['IC']['len']
            tl = self.DICTIONARY['IC']['trkLen']
            f  = self.DICTIONARY['IC']['fiber']
            for s in xrange(self.DICTIONARY['IC']['n']) :
                sl[s] /= tl[ f[s] ]

        # reorder the segments based on the "v" field
        idx = np.argsort( self.DICTIONARY['IC']['v'], kind='mergesort' )
        self.DICTIONARY['IC']['v']     = self.DICTIONARY['IC']['v'][ idx ]
        self.DICTIONARY['IC']['o']     = self.DICTIONARY['IC']['o'][ idx ]
        self.DICTIONARY['IC']['fiber'] = self.DICTIONARY['IC']['fiber'][ idx ]
        self.DICTIONARY['IC']['len']   = self.DICTIONARY['IC']['len'][ idx ]
        del idx

        print '[ %d fibers and %d segments ]' % ( self.DICTIONARY['IC']['nF'], self.DICTIONARY['IC']['n'] )

        # segments from the peaks
        # -----------------------
        print '\t* segments from the peaks...',
        sys.stdout.flush()

        self.DICTIONARY['EC'] = {}

        vx = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_EC_vx.dict'), dtype=np.uint8 ).astype(np.uint32)
        vy = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_EC_vy.dict'), dtype=np.uint8 ).astype(np.uint32)
        vz = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_EC_vz.dict'), dtype=np.uint8 ).astype(np.uint32)
        self.DICTIONARY['EC']['v'] = vx + self.CONFIG['dim'][0] * ( vy + self.CONFIG['dim'][1] * vz )
        del vx, vy, vz

        self.DICTIONARY['EC']['nE'] = self.DICTIONARY['EC']['v'].size

        ox = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_EC_ox.dict'), dtype=np.uint8 ).astype(np.uint16)
        oy = np.fromfile( os.path.join(self.CONFIG['TRACKING_path'],'dictionary_EC_oy.dict'), dtype=np.uint8 ).astype(np.uint16)
        self.DICTIONARY['EC']['o'] = oy + 181*ox
        del ox, oy

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
        self.DICTIONARY['ISO']['v'] = vx + self.CONFIG['dim'][0] * ( vy + self.CONFIG['dim'][1] * vz )
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

        lut = np.zeros( self.CONFIG['dim'], dtype=np.uint32 ).ravel()
        for i in xrange(idx.size) :
            lut[ idx[i] ] = i
        self.DICTIONARY['IC'][ 'v'] = lut[ self.DICTIONARY['IC'][ 'v'] ]
        self.DICTIONARY['EC'][ 'v'] = lut[ self.DICTIONARY['EC'][ 'v'] ]
        self.DICTIONARY['ISO']['v'] = lut[ self.DICTIONARY['ISO']['v'] ]

        print '         [ OK ]'

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def set_threads( self, n = None ) :
        """
        Set the number of threads to use for the matrix-vector operations with A and A'.

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
                for c in np.bincount( self.DICTIONARY['IC']['v'] ) :
                    tot += c
                    if tot >= N :
                        self.THREADS['IC'][t] = self.THREADS['IC'][t-1] + tot
                        t += 1
                        tot = 0
            self.THREADS['IC'][-1] = self.DICTIONARY['IC']['n']
        else :
            self.THREADS['IC'] = None

        if self.DICTIONARY['EC']['nE'] > 0 :
            self.THREADS['EC'] = np.zeros( n+1, dtype=np.uint32 )
            for i in xrange(n) :
                self.THREADS['EC'][i] = np.searchsorted( self.DICTIONARY['EC']['v'], self.DICTIONARY['IC']['v'][ self.THREADS['IC'][i] ] )
            self.THREADS['EC'][-1] = self.DICTIONARY['EC']['nE']
        else :
            self.THREADS['EC'] = None

        if self.DICTIONARY['nV'] > 0 :
            self.THREADS['ISO'] = np.zeros( n+1, dtype=np.uint32 )
            for i in xrange(n) :
                self.THREADS['ISO'][i] = np.searchsorted( self.DICTIONARY['ISO']['v'], self.DICTIONARY['IC']['v'][ self.THREADS['IC'][i] ] )
            self.THREADS['ISO'][-1] = self.DICTIONARY['nV']
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

                t = 0
                tot = 0
                i1 = 0
                i2 = 0
                N = np.floor(self.DICTIONARY['IC']['n']/n)
                for i in xrange(C.size) :
                    i2 += C[i]
                    tot += C[i]
                    if tot >= N :
                        self.THREADS['ICt'][ i1:i2 ] = t
                        t += 1
                        if t==n-1 :
                            break
                        i1 = i2
                        tot = C[i]
                self.THREADS['ICt'][idx] = self.THREADS['ICt'].copy()

        else :
            self.THREADS['ICt'] = None

        if self.DICTIONARY['EC']['nE'] > 0 :
            self.THREADS['ECt'] = np.zeros( n+1, dtype=np.uint32 )
            N = np.floor( self.DICTIONARY['EC']['nE']/n )
            for i in xrange(1,n) :
                self.THREADS['ECt'][i] = self.THREADS['ECt'][i-1] + N
            self.THREADS['ECt'][-1] = self.DICTIONARY['EC']['nE']
        else :
            self.THREADS['ECt'] = None

        if self.DICTIONARY['nV'] > 0 :
            self.THREADS['ISOt'] = np.zeros( n+1, dtype=np.uint32 )
            N = np.floor( self.DICTIONARY['nV']/n )
            for i in xrange(1,n) :
                self.THREADS['ISOt'][i] = self.THREADS['ISOt'][i-1] + N
            self.THREADS['ISOt'][-1] = self.DICTIONARY['nV']
        else :
            self.THREADS['ISOt'] = None

        print '[ OK ]'

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def build_operator( self ) :
        """
        Compile/build the operator for computing the matrix-vector multiplications by A and A'
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
        config.nIC      = self.KERNELS['wmr'].shape[0]
        config.nEC      = self.KERNELS['wmh'].shape[0]
        config.nISO     = self.KERNELS['iso'].shape[0]
        if not 'commit.operator.operator' in sys.modules :
            import commit.operator.operator
        else :
            reload( sys.modules['commit.operator.operator'] )
        self.A = sys.modules['commit.operator.operator'].LinearOperator( self.DICTIONARY, self.KERNELS, self.THREADS )

        print '   [ %.1f seconds ]' % ( time.time() - tic )


    def fit( self, tol_fun = 1e-3, max_iter = 100, verbose = 1 ) :
        """
        Fit the model to the data.

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

        self.CONFIG['optimization'] = {}
        self.CONFIG['optimization']['tol_fun']  = tol_fun
        self.CONFIG['optimization']['max_iter'] = max_iter
        self.CONFIG['optimization']['verbose']  = verbose

        # run solver
        t = time.time()
        print '\n-> Fit model using "nnls":'
        Y = self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float64)
        self.x = commit.solvers.nnls( Y, self.A, tol_fun=tol_fun, max_iter=max_iter, verbose=verbose )
        self.CONFIG['optimization']['fit_time'] = time.time()-t
        print '   [ %s ]' % ( time.strftime("%Hh %Mm %Ss", time.gmtime(self.CONFIG['optimization']['fit_time']) ) )


    def save_results( self, path_suffix = None ) :
        """
        Save the output (coefficients, errors, maps etc).

        Parameters
        ----------
        path_suffix : string
            Text to be appended to "Results" to create the output path (default : None)
        """
        if self.x is None :
            raise RuntimeError( 'Model not fitted to the data; call "fit()" first.' )

        RESULTS_path = 'Results_' + self.model.id
        if path_suffix :
            self.CONFIG['path_suffix'] = path_suffix
            RESULTS_path = RESULTS_path +'_'+ path_suffix

        print '\n-> Saving results to "%s/*":' % RESULTS_path
        tic = time.time()

        # create folder or delete existing files (if any)
        RESULTS_path = os.path.join( self.CONFIG['TRACKING_path'], RESULTS_path )
        if not os.path.exists( RESULTS_path ) :
            os.makedirs( RESULTS_path )
        else :
            for f in glob.glob( os.path.join(RESULTS_path,'*') ) :
                os.remove( f )
        self.CONFIG['RESULTS_path'] = RESULTS_path

        # Configuration and results
        print '\t* configuration and results...',
        sys.stdout.flush()
        nF = self.DICTIONARY['IC']['nF']
        nE = self.DICTIONARY['EC']['nE']
        nV = self.DICTIONARY['nV']
        if self.CONFIG['doNormalizeKernels'] :
            # renormalize the coefficients
            norm1 = np.tile(self.KERNELS['wmr_norm'],nF)
            norm2 = np.tile(self.KERNELS['wmh_norm'],nE)
            norm3 = np.tile(self.KERNELS['iso_norm'],nV)
            x = self.x / np.hstack( (norm1,norm2,norm3) )
        else :
            x = self.x
        with open( os.path.join(RESULTS_path,'results.pickle'), 'wb+' ) as fid :
            cPickle.dump( [self.CONFIG, x], fid, protocol=2 )
        print '[ OK ]'

        # Map of wovelwise errors
        print '\t* fitting errors:'

        niiMAP_img = np.zeros( self.CONFIG['dim'], dtype=np.float32 )
        niiMAP = nibabel.Nifti1Image( niiMAP_img, self.niiDWI.affine )

        y_mea = np.reshape( self.niiDWI_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'], : ].flatten().astype(np.float32), (nV,-1) )
        y_est = np.reshape( self.A.dot(self.x), (nV,-1) ).astype(np.float32)

        print '\t\t- RMSE...',
        sys.stdout.flush()
        tmp = np.sqrt( np.mean((y_mea-y_est)**2,axis=1) )
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP.get_header()['cal_min'] = 0
        niiMAP.get_header()['cal_max'] = tmp.max()
        nibabel.save( niiMAP, os.path.join(RESULTS_path,'fit_RMSE.nii.gz') )
        print ' [ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() )

        print '\t\t- NRMSE...',
        sys.stdout.flush()
        tmp = np.sum(y_mea**2,axis=1)
        idx = np.where( tmp < 1E-12 )
        tmp[ idx ] = 1
        tmp = np.sqrt( np.sum((y_mea-y_est)**2,axis=1) / tmp )
        tmp[ idx ] = 0
        niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = tmp
        niiMAP.get_header()['cal_min'] = 0
        niiMAP.get_header()['cal_max'] = 1
        nibabel.save( niiMAP, os.path.join(RESULTS_path,'fit_NRMSE.nii.gz') )
        print '[ %.3f +/- %.3f ]' % ( tmp.mean(), tmp.std() )

        # Map of compartment contributions
        print '\t* voxelwise contributions:'

        print '\t\t- intra-axonal',
        sys.stdout.flush()
        niiMAP_img[:] = 0
        if len(self.KERNELS['wmr']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0]
            tmp = x[:offset].reshape( (-1,nF) ).sum( axis=0 )
            xv = np.bincount( self.DICTIONARY['IC']['v'], minlength=nV,
                weights=tmp[ self.DICTIONARY['IC']['fiber'] ] * self.DICTIONARY['IC']['len']
            ).astype(np.float32)
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        nibabel.save( niiMAP, os.path.join(RESULTS_path,'compartment_IC.nii.gz') )
        print '[ OK ]'

        print '\t\t- extra-axonal',
        sys.stdout.flush()
        niiMAP_img[:] = 0
        if len(self.KERNELS['wmh']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0]
            tmp = x[offset:offset+nE*len(self.KERNELS['wmh'])].reshape( (-1,nE) ).sum( axis=0 )
            xv = np.bincount( self.DICTIONARY['EC']['v'], weights=tmp, minlength=nV ).astype(np.float32)
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        nibabel.save( niiMAP, os.path.join(RESULTS_path,'compartment_EC.nii.gz') )
        print '[ OK ]'

        print '\t\t- isotropic',
        sys.stdout.flush()
        niiMAP_img[:] = 0
        if len(self.KERNELS['iso']) > 0 :
            offset = nF * self.KERNELS['wmr'].shape[0] + nE * self.KERNELS['wmh'].shape[0]
            xv = x[offset:].reshape( (-1,nV) ).sum( axis=0 )
            niiMAP_img[ self.DICTIONARY['MASK_ix'], self.DICTIONARY['MASK_iy'], self.DICTIONARY['MASK_iz'] ] = xv
        nibabel.save( niiMAP, os.path.join(RESULTS_path,'compartment_ISO.nii.gz') )
        print '   [ OK ]'

        print '   [ %.1f seconds ]' % ( time.time() - tic )


def precompute_rotation_matrices( lmax = 12 ) :
    """
    Precompute the rotation matrices to rotate the high-resolution kernels (500 directions per shell)

    Parameters
    ----------
    lmax : int
        Maximum SH order to use for the rotation phase (default : 12)
    """
    filename = os.path.join( dipy_home, 'COMMIT_aux_matrices_lmax=%d.pickle'%lmax )
    if os.path.isfile( filename ) :
        return

    print '\n-> Precomputing rotation matrices for l_max=%d:' % lmax
    AUX = {}
    AUX['lmax'] = lmax

    # load file with 500 directions
    grad = np.loadtxt( os.path.join(os.path.dirname(commit.__file__), '500_dirs.txt') )
    for i in xrange(grad.shape[0]) :
        grad[i,:] /= np.linalg.norm( grad[i,:] )
        if grad[i,1] < 0 :
            grad[i,:] = -grad[i,:] # to ensure they are in the spherical range [0,180]x[0,180]

    # matrix to fit the SH coefficients
    _, theta, phi = cart2sphere( grad[:,0], grad[:,1], grad[:,2] )
    tmp, _, _ = real_sym_sh_basis( lmax, theta, phi )
    AUX['fit'] = np.dot( np.linalg.pinv( np.dot(tmp.T,tmp) ), tmp.T )

    # matrices to rotate the functions in SH space
    AUX['Ylm_rot'] = np.zeros( (181,181), dtype=np.object )
    for ox in xrange(181) :
        for oy in xrange(181) :
            tmp, _, _ = real_sym_sh_basis( lmax, ox/180.0*np.pi, oy/180.0*np.pi )
            AUX['Ylm_rot'][ox,oy] = tmp

    # auxiliary data to perform rotations
    AUX['const'] = np.zeros( AUX['fit'].shape[0], dtype=np.float64 )
    AUX['idx_m0'] = np.zeros( AUX['fit'].shape[0], dtype=np.int32 )
    i = 0
    for l in xrange(0,AUX['lmax']+1,2) :
        const  = np.sqrt(4.0*np.pi/(2.0*l+1.0))
        idx_m0 = (l*l + l + 2.0)/2.0 - 1
        for m in xrange(-l,l+1) :
            AUX['const'][i]  = const
            AUX['idx_m0'][i] = idx_m0
            i += 1


    with open( filename, 'wb+' ) as fid :
        cPickle.dump( AUX, fid, protocol=2 )

    print '   [ DONE ]'


def rotate_kernel( K, AUX, idx_IN, idx_OUT, is_isotropic ) :
    """
    Rotate a spherical function.

    Parameters
    ----------
    K : numpy.ndarray
        Spherical function (in signal space) to rotate
    AUX : dictionary
        Auxiliary data structures needed to rotate functions in SH space
    idx_IN : list of list
        Index of samples in input kernel (K) belonging to each shell
    idx_OUT : list of list
        Index of samples in output kernel (K) belonging to each shell
    is_isotropic : boolean
        Indentifies whether K is an isotropic function or not

    Returns
    -------
    KRlm = numpy.array
        Spherical function (in SH space) rotated to 181x181 directions distributed
        on a hemisphere
    """

    # project kernel K to SH space
    Klm = []
    for s in xrange(len(idx_IN)) :
        Klm.append( np.dot( AUX['fit'], K[ idx_IN[s] ] ) )

    n = len(idx_IN)*AUX['fit'].shape[0]

    if is_isotropic == False :
        # fit SH and rotate kernel to 181*181 directions
        KRlm = np.zeros( (181,181,n), dtype=np.float32 )
        for ox in xrange(181) :
            for oy in xrange(181) :
                Ylm_rot = AUX['Ylm_rot'][ox,oy]
                for s in xrange(len(idx_IN)) :
                    KRlm[ox,oy,idx_OUT[s]] = AUX['const'] * Klm[s][AUX['idx_m0']] * Ylm_rot
    else :
        # simply fit SH
        KRlm = np.zeros( n, dtype=np.float32 )
        for s in xrange(len(idx_IN)) :
            KRlm[idx_OUT[s]] = Klm[s].astype(np.float32)

    return KRlm


def resample_kernel( KRlm, nS, idx_out, Ylm_out, is_isotropic ) :
    """
    Resample a spherical function

    Parameters
    ----------
    KRlm : numpy.array
        Rotated spherical functions (in SH space) to project
    nS : integer
        Number of samples in the subject's acquisition scheme
    idx_out : list of list
        Index of samples in output kernel
    Ylm_out : numpy.array
        Matrix to project back all shells from SH space to signal space (of the subject)
    is_isotropic : boolean
        Indentifies whether Klm is an isotropic function or not

    Returns
    -------
    KR = numpy.array
        Rotated spherical functions projected to signal space of the subject
    """
    if is_isotropic == False :
        KR = np.ones( (181,181,nS), dtype=np.float32 )
        for ox in xrange(181) :
            for oy in xrange(181) :
                KR[ox,oy,idx_out] = np.dot( Ylm_out, KRlm[ox,oy,:] ).astype(np.float32)
    else :
        KR = np.ones( nS, dtype=np.float32 )
        KR[idx_out] = np.dot( Ylm_out, KRlm ).astype(np.float32)
    return KR



class Scheme :
    """
    A class to hold information about an acquisition scheme.
    The scheme can be specified in two formats:
    - as a Nx4 matrix: the first three columns are the gradient directions and
      the fourth is the b-value (s/mm^2).
    - as a Nx7 matrix: the first three columns are the gradient directions, and
      the remaining four are: the gradient strength (T/m), big delta (s),
      small delta (s) and echo time (s), respectively.
    The "Camino header line" (eg. VERSION: BVECTOR) is optional.
    """

    def __init__( self, data, b0_thr = 0 ) :
        """
        Initialize the acquisition scheme.

        Parameters
        ----------
        data : string or numpy.ndarray
            The filename of the scheme or a matrix containing the actual values
        b0_thr : float
            The threshold on the b-values to identify the b0 images (default: 0)
        """
        if type(data) is str :
            # try loading from file
            try :
                n = 0 # headers lines to skip to get to the numeric data
                with open(data) as fid :
                    for line in fid :
                        if re.match( r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?', line.strip() ) :
                            break
                        n += 1
                tmp = np.loadtxt( data, skiprows=n )
            except :
                raise IOError( 'Unable to open scheme file' )
            self.load_from_table( tmp, b0_thr )
        else :
            # try loading from matrix
            self.load_from_table( data, b0_thr )


    def load_from_table( self, data, b0_thr = 0 ) :
        """
        Build the structure from an input matrix.

        The first three columns represent the gradient directions.
        Then, we accept two formats to describe each gradient:
            - if the shape of data is Nx4, the 4^ column is the b-value;
            - if the shape of data is Nx7, the last 4 columns are, respectively, the gradient strength, big delta, small delta and TE.

        Parameters
        ----------
        data : numpy.ndarray
            Matrix containing tall the values.
        b0_thr : float
            The threshold on the b-values to identify the b0 images (default: 0)
        """
        self.raw = data

        # number of samples
        self.nS = self.raw.shape[0]

        # set/calculate the b-values
        if self.raw.shape[1] == 4 :
            self.version = 0
            self.b = self.raw[:,3]
        elif self.raw.shape[1] == 7 :
            self.version = 1
            self.b = ( 267.513e6 * self.raw[:,3] * self.raw[:,5] )**2 * (self.raw[:,4] - self.raw[:,5]/3.0) * 1e-6 # in mm^2/s
        else :
            raise ValueError( 'Unrecognized scheme format' )

        # store information about the volumes
        self.b0_thr    = b0_thr
        self.b0_idx    = np.where( self.b <= b0_thr )[0]
        self.b0_count  = len( self.b0_idx )
        self.dwi_idx   = np.where( self.b > b0_thr )[0]
        self.dwi_count = len( self.dwi_idx )

        # ensure the directions are in the spherical range [0,180]x[0,180]
        idx = np.where( self.raw[:,1] < 0 )[0]
        self.raw[idx,0:3] = -self.raw[idx,0:3]

        # store information about each shell in a dictionary
        self.shells = []
        tmp = np.ascontiguousarray( self.raw[:,3:] )
        schemeUnique, schemeUniqueInd = np.unique( tmp.view([('', tmp.dtype)]*tmp.shape[1]), return_index=True )
        schemeUnique = schemeUnique.view(tmp.dtype).reshape((schemeUnique.shape[0], tmp.shape[1]))
        schemeUnique = schemeUnique[np.argsort(schemeUniqueInd)]
        bUnique = self.b[ schemeUniqueInd ]
        for i in xrange(schemeUnique.shape[0]) :
            if bUnique[i] <= b0_thr :
                continue
            shell = {}
            shell['b'] = bUnique[i]
            if self.version == 0 :
                shell['G']     = None
                shell['Delta'] = None
                shell['delta'] = None
                shell['TE']    = None
            else :
                shell['G']     = schemeUnique[i,0]
                shell['Delta'] = schemeUnique[i,1]
                shell['delta'] = schemeUnique[i,2]
                shell['TE']    = schemeUnique[i,3]

            shell['idx']  = np.where( self.b==bUnique[i] )[0]
            shell['grad'] = self.raw[shell['idx'],0:3]
            self.shells.append( shell )


    @property
    def nS( self ) :
        return self.b0_count + self.dwi_count


    def create_high_resolution( self, b_scale = 1 ) :
        """
        Create an high-resolution version (500 directions per shell).
        All other parameters remain the same.

        Parameters
        ----------
        b_scale : float
            If needed, apply a scaling to the b-values (default : 1)
        """
        # load HR directions and ensure they are in the spherical range [0,180]x[0,180]
        grad = np.loadtxt( os.path.join(os.path.dirname(commit.__file__), '500_dirs.txt') )
        for i in xrange(grad.shape[0]) :
            grad[i,:] /= np.linalg.norm( grad[i,:] )
            if grad[i,1] < 0 :
                grad[i,:] = -grad[i,:]

        n = len( self.shells )
        raw = np.zeros( (500*n, 4 if self.version==0 else 7) )
        row = 0
        for i in xrange(n) :
            raw[row:row+500,0:3] = grad
            if self.version == 0 :
                raw[row:row+500,3] = self.shells[i]['b'] * b_scale
            else :
                raw[row:row+500,3] = self.shells[i]['G']
                raw[row:row+500,4] = self.shells[i]['Delta']
                raw[row:row+500,5] = self.shells[i]['delta']
                raw[row:row+500,6] = self.shells[i]['TE']
            row += 500

        return Scheme( raw )
