import numpy as np
from os.path import exists, join as pjoin
from os import remove
import sys
import subprocess
import tempfile
import commit.core
from dipy.core.gradients import gradient_table
from dipy.sims.voxel import single_tensor


class StickZeppelinBall :
    """
    Simulate the response functions according to the Stick-Zeppelin-Ball model.

    The contributions of the tracts are modeled as "sticks", i.e. tensors with
    a given axial diffusivity (d_par) but null perpendicular diffusivity.
    Extra-cellular contributions are modeled as tensors with the same axial diffusivity
    as the sticks (d_par) and whose perpendicular diffusivities are calculated with a tortuosity
    model as a function of the intra-cellular volume fractions (ICVFs).
    Isotropic contributions are modeled as tensors with isotropic diffusivities (d_ISOs).

    References
    ----------
    .. [1] Panagiotaki et al. (2012) Compartment models of the diffusion MR signal
           in brain white matter: A taxonomy and comparison. NeuroImage, 59: 2241-54
    """

    def __init__( self ) :
        self.id     = 'StickZeppelinBall'
        self.name   = 'Stick-Zeppelin-Ball'
        self.d_par  = 1.7E-3              # Parallel diffusivity [mm^2/s]
        self.ICVFs  = [ 0.7 ]             # Intra-cellular volume fraction(s) [0..1]
        self.d_ISOs = [ 1.7E-3, 3.0E-3 ]  # Isotropic diffusivitie(s) [mm^2/s]


    def set( self, d_par, ICVFs, d_ISOs ) :
        self.d_par  = d_par
        self.ICVFs  = ICVFs
        self.d_ISOs = d_ISOs


    def generate( self, out_path, scheme, aux, idx_in, idx_out ) :
        print '\t* 1 stick, %d extra-cellular and %d isotropic' % ( len(self.ICVFs), len(self.d_ISOs) )

        # store the number of samples in the original scheme
        self.nS = scheme.nS

        # create an high-resolution version of the scheme
        scheme_high = scheme.create_high_resolution( b_scale=1 )
        gtab = gradient_table( scheme_high.b, scheme_high.raw[:,0:3] )

        # Stick
        print '\t* A_001...',
        sys.stdout.flush()
        signal = single_tensor( gtab, evals=[0, 0, self.d_par] )
        lm = commit.core.rotate_kernel( signal, aux, idx_in, idx_out, False )
        np.save( pjoin( out_path, 'A_001.npy' ), lm )
        print '[ OK ]'

        # Zeppelin(s)
        idx = 2
        for d in [ self.d_par*(1.0-ICVF) for ICVF in self.ICVFs] :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            signal = single_tensor( gtab, evals=[d, d, self.d_par] )
            lm = commit.core.rotate_kernel( signal, aux, idx_in, idx_out, False )
            np.save( pjoin( out_path, 'A_%03d.npy'%idx ), lm )
            idx += 1
            print '[ OK ]'

        # Ball(s)
        for d in self.d_ISOs :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            signal = single_tensor( gtab, evals=[d, d, d] )
            lm = commit.core.rotate_kernel( signal, aux, idx_in, idx_out, True )
            np.save( pjoin( out_path, 'A_%03d.npy'%idx ), lm )
            idx += 1
            print '[ OK ]'


    def resample( self, in_path, idx_out, Ylm_out ) :
        KERNELS = {}
        KERNELS['model'] = self.id
        KERNELS['wmr']   = np.zeros( (1,181,181,self.nS), dtype=np.float32 )
        KERNELS['wmh']   = np.zeros( (len(self.ICVFs),181,181,self.nS), dtype=np.float32 )
        KERNELS['iso']   = np.zeros( (len(self.d_ISOs),self.nS), dtype=np.float32 )

        # Stick
        print '\t* A_001...',
        sys.stdout.flush()
        lm = np.load( pjoin( in_path, 'A_001.npy' ) )
        KERNELS['wmr'][0,...] = commit.core.resample_kernel( lm, self.nS, idx_out, Ylm_out, False )
        print '[ OK ]'

        # Zeppelin(s)
        idx = 2
        for i in xrange(len(self.ICVFs)) :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            lm = np.load( pjoin( in_path, 'A_%03d.npy'%idx ) )
            KERNELS['wmh'][i,...] = commit.core.resample_kernel( lm, self.nS, idx_out, Ylm_out, False )
            idx += 1
            print '[ OK ]'

        # Ball(s)
        for i in xrange(len(self.d_ISOs)) :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            lm = np.load( pjoin( in_path, 'A_%03d.npy'%idx ) )
            KERNELS['iso'][i,...] = commit.core.resample_kernel( lm, self.nS, idx_out, Ylm_out, True )
            idx += 1
            print '[ OK ]'

        return KERNELS



class CylinderZeppelinBall :
    """
    Simulate the response functions according to the Cylinder-Zeppelin-Ball model.

    The contributions of the tracts are modeled as "cylinders" with specific radii (Rs)
    and a given axial diffusivity (d_par).
    Extra-cellular contributions are modeled as tensors with the same axial diffusivity
    as the sticks (d_par) and whose perpendicular diffusivities are calculated with a tortuosity
    model as a function of the intra-cellular volume fractions (ICVFs).
    Isotropic contributions are modeled as tensors with isotropic diffusivities (d_ISOs).

    NB: this models works only with schemes containing the full specification of
        the diffusion gradients (eg gradient strength, small delta etc).

    NB: this model requires Camino to be installed and properly configured
        in the system; in particular, the script "datasynth" must be placed
        in your system path.

    References
    ----------
    .. [1] Panagiotaki et al. (2012) Compartment models of the diffusion MR signal
           in brain white matter: A taxonomy and comparison. NeuroImage, 59: 2241-54
    """

    def __init__( self ) :
        self.id     = 'CylinderZeppelinBall'
        self.name   = 'Cylinder-Zeppelin-Ball'
        self.d_par  = 1.7E-3              # Parallel diffusivity [mm^2/s]
        self.Rs     = [ 4E-6, 16E-6 ]     # Radii of the axons [micrometers]
        self.ICVFs  = [ 0.7 ]             # Intra-cellular volume fraction(s) [0..1]
        self.d_ISOs = [ 1.7E-3, 3.0E-3 ]  # Isotropic diffusivitie(s) [mm^2/s]


    def set( self, d_par, Rs, ICVFs, d_ISOs ) :
        self.d_par  = d_par
        self.Rs  = Rs
        self.ICVFs  = ICVFs
        self.d_ISOs = d_ISOs


    def generate( self, out_path, scheme, aux, idx_in, idx_out ) :
        if scheme.version != 1 :
            raise RuntimeError( 'This model requires a "VERSION: STEJSKALTANNER" scheme.' )

        print '\t* %d restricted, %d hindered and %d isotropic' % ( len(self.Rs), len(self.ICVFs), len(self.d_ISOs) )

        # store the number of samples in the original scheme
        self.nS = scheme.nS

        # create a high-resolution scheme to pass to 'datasynth'
        scheme_high = scheme.create_high_resolution( b_scale=1E6 )
        filename_scheme = pjoin( out_path, 'scheme.txt' )
        np.savetxt( filename_scheme, scheme_high.raw, fmt='%15.8e', delimiter=' ', header='VERSION: STEJSKALTANNER', comments='' )

        # temporary file where to store "datasynth" output
        filename_signal = pjoin( tempfile._get_default_tempdir(), next(tempfile._get_candidate_names())+'.Bfloat' )

        # Cylinder(s)
        idx = 1
        for R in self.Rs :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()

            CMD = 'datasynth -synthmodel compartment 1 CYLINDERGPD %E 0 0 %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null' % ( self.d_par*1E-6, R, filename_scheme, filename_signal )
            subprocess.call( CMD, shell=True )
            if not exists( filename_signal ) :
                raise RuntimeError( 'Problems generating the signal with "datasynth"' )
            signal  = np.fromfile( filename_signal, dtype='>f4' )
            if exists( filename_signal ) :
                remove( filename_signal )

            lm = commit.core.rotate_kernel( signal, aux, idx_in, idx_out, False )
            np.save( pjoin( out_path, 'A_%03d.npy'%idx ), lm )
            idx += 1
            print '[ OK ]'

        # Zeppelin(s)
        for d in [ self.d_par*(1.0-ICVF) for ICVF in self.ICVFs] :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()

            CMD = 'datasynth -synthmodel compartment 1 ZEPPELIN %E 0 0 %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null' % ( self.d_par*1E-6, d*1e-6, filename_scheme, filename_signal )
            subprocess.call( CMD, shell=True )
            if not exists( filename_signal ) :
                raise RuntimeError( 'Problems generating the signal with "datasynth"' )
            signal  = np.fromfile( filename_signal, dtype='>f4' )
            if exists( filename_signal ) :
                remove( filename_signal )

            lm = commit.core.rotate_kernel( signal, aux, idx_in, idx_out, False )
            np.save( pjoin( out_path, 'A_%03d.npy'%idx ), lm )
            idx += 1
            print '[ OK ]'

        # Ball(s)
        for d in self.d_ISOs :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()

            CMD = 'datasynth -synthmodel compartment 1 BALL %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null' % ( d*1e-6, filename_scheme, filename_signal )
            subprocess.call( CMD, shell=True )
            if not exists( filename_signal ) :
                raise RuntimeError( 'Problems generating the signal with "datasynth"' )
            signal  = np.fromfile( filename_signal, dtype='>f4' )
            if exists( filename_signal ) :
                remove( filename_signal )

            lm = commit.core.rotate_kernel( signal, aux, idx_in, idx_out, True )
            np.save( pjoin( out_path, 'A_%03d.npy'%idx ), lm )
            idx += 1
            print '[ OK ]'


    def resample( self, in_path, idx_out, Ylm_out ) :
        KERNELS = {}
        KERNELS['model'] = self.id
        KERNELS['wmr']   = np.zeros( (len(self.Rs),181,181,self.nS), dtype=np.float32 )
        KERNELS['wmh']   = np.zeros( (len(self.ICVFs),181,181,self.nS), dtype=np.float32 )
        KERNELS['iso']   = np.zeros( (len(self.d_ISOs),self.nS), dtype=np.float32 )

        # Cylinder(s)
        idx = 1
        for i in xrange(len(self.Rs)) :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            lm = np.load( pjoin( in_path, 'A_%03d.npy'%idx ) )
            KERNELS['wmr'][i,...] = commit.core.resample_kernel( lm, self.nS, idx_out, Ylm_out, False )
            idx += 1
            print '[ OK ]'

        # Zeppelin(s)
        for i in xrange(len(self.ICVFs)) :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            lm = np.load( pjoin( in_path, 'A_%03d.npy'%idx ) )
            KERNELS['wmh'][i,...] = commit.core.resample_kernel( lm, self.nS, idx_out, Ylm_out, False )
            idx += 1
            print '[ OK ]'

        # Ball
        for i in xrange(len(self.d_ISOs)) :
            print '\t* A_%03d...' % idx,
            sys.stdout.flush()
            lm = np.load( pjoin( in_path, 'A_%03d.npy'%idx ) )
            KERNELS['iso'][i,...] = commit.core.resample_kernel( lm, self.nS, idx_out, Ylm_out, True )
            idx += 1
            print '[ OK ]'

        return KERNELS