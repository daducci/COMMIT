import os, sys
import pickle
import amico

from commit import trk2dictionary
import commit
commit.core.setup()


# get path to the current directory
# local_path = os.path.dirname( os.path.realpath( __file__ ) )

def run_commit_StickZeppelinBall():
    local_path = "/media/full/DATA/Software/COMMIT_versions/demo_data"
    trk2dictionary.run(
        filename_tractogram = os.path.join( local_path, 'demo01_fibers.tck' ),
        filename_peaks      = os.path.join( local_path, 'peaks.nii.gz' ),
        filename_mask       = os.path.join( local_path, 'WM.nii.gz' ), #'WM.nii.gz',
        fiber_shift         = 0.5,
        peaks_use_affine    = True
    )

    amico.util.fsl2scheme( os.path.join( local_path, 'bvals.txt'), os.path.join( local_path, 'bvecs.txt'), os.path.join( local_path, 'DWI.scheme') )

    mit = commit.Evaluation( local_path, '.' )
    mit.load_data( os.path.join( local_path, 'DWI.nii.gz'), os.path.join( local_path, 'DWI.scheme') )

    mit.set_model( 'StickZeppelinBall' )
    d_par   = 1.7E-3             # Parallel diffusivity [mm^2/s]
    d_perps = [ 0.51E-3 ]        # Perpendicular diffusivity(s) [mm^2/s]
    # d_perps = []        # Perpendicular diffusivity(s) [mm^2/s]
    d_isos  = [ 1.7E-3, 3.0E-3 ] # Isotropic diffusivity(s) [mm^2/s]
    mit.model.set( d_par, d_perps, d_isos )
    mit.generate_kernels( regenerate=True )
    mit.load_kernels()

    mit.load_dictionary( os.path.join( local_path, 'COMMIT') )

    mit.set_threads()
    mit.build_operator()

    mit.fit( tol_fun=1e-3, max_iter=1000 )
    mit.save_results()


def run_commit_BallandStick():
    local_path = "/media/full/DATA/Software/COMMIT_versions/demo_data"
    trk2dictionary.run(
        filename_tractogram = os.path.join( local_path, 'demo01_fibers.tck' ),
        filename_peaks      = os.path.join( local_path, 'peaks.nii.gz' ),
        filename_mask       = os.path.join( local_path, 'WM.nii.gz' ), #'WM.nii.gz',
        fiber_shift         = 0.5,
        peaks_use_affine    = True
    )

    amico.util.fsl2scheme( os.path.join( local_path, 'bvals.txt'), os.path.join( local_path, 'bvecs.txt'), os.path.join( local_path, 'DWI.scheme') )

    mit = commit.Evaluation( '.', '.' )
    mit.load_data( os.path.join( local_path, 'DWI.nii.gz'), os.path.join( local_path, 'DWI.scheme') )

    mit.set_model( 'StickZeppelinBall' )
    d_par   = 1.7E-3             # Parallel diffusivity [mm^2/s]
    d_perps = []        # Perpendicular diffusivity(s) [mm^2/s]
    d_isos  = [ 1.7E-3, 3.0E-3 ] # Isotropic diffusivity(s) [mm^2/s]
    mit.model.set( d_par, d_perps, d_isos )
    mit.generate_kernels( regenerate=True )
    mit.load_kernels()

    mit.load_dictionary( os.path.join( local_path, 'COMMIT') )

    mit.set_threads()
    mit.build_operator()

    mit.fit( tol_fun=1e-3, max_iter=1000 )
    mit.save_results()


def run_commit_VolumeFractions():
    local_path = "/media/full/DATA/Software/COMMIT_versions/demo_data"
    trk2dictionary.run(
        filename_tractogram = os.path.join( local_path, 'demo01_fibers.tck' ),
        filename_peaks      = os.path.join( local_path, 'peaks.nii.gz' ),
        filename_mask       = os.path.join( local_path, 'WM.nii.gz' ), #'WM.nii.gz',
        fiber_shift         = 0.5,
        peaks_use_affine    = True
    )

    amico.util.fsl2scheme( os.path.join( local_path, 'bvals.txt'), os.path.join( local_path, 'bvecs.txt'), os.path.join( local_path, 'DWI.scheme') )

    mit = commit.Evaluation( '.', '.' )
    mit.load_data( os.path.join( local_path, 'DWI.nii.gz'), os.path.join( local_path, 'DWI.scheme') )

    mit.set_model( 'StickZeppelinBall' )
    d_par   = 1.7E-3             # Parallel diffusivity [mm^2/s]
    d_perps = []        # Perpendicular diffusivity(s) [mm^2/s]
    d_isos  = [ 1.7E-3, 3.0E-3 ] # Isotropic diffusivity(s) [mm^2/s]
    mit.model.set( d_par, d_perps, d_isos )
    mit.generate_kernels( regenerate=True )
    mit.load_kernels()

    mit.load_dictionary( os.path.join( local_path, 'COMMIT') )

    mit.set_threads()
    mit.build_operator()

    mit.fit( tol_fun=1e-3, max_iter=1000 )
    mit.save_results()


# def check_config(result_config, ref_config):
    

def check_results(pickle_result, ref_pickle):
    with open(pickle_result, 'rb') as f:
        data = pickle.load(f)
    with open(ref_pickle, 'rb') as f:
        ref_data = pickle.load(f)

    result_optimization = data[0]["optimization"]
    ref_optimization = ref_data[0]["optimization"]
    try:
        assert result_optimization["fit_parameters"] == ref_optimization["fit_parameters"]
    except AssertionError:
        sys.exit(1)


def run_tests():
    local_path = os.path.dirname( os.path.realpath( __file__ ) )
    out_path = "/media/full/DATA/Software/COMMIT_versions/demo_data"
    ref_pickle_StickZeppelinBall = os.path.join( local_path, 'demo_data', 'ref_results', 'ref_results_StickZeppelinBall.pickle' )
    ref_pickle_BallandStick = os.path.join( local_path, 'demo_data', 'ref_results', 'ref_results_BallandStick.pickle' )
    # ref_pickle_VolumeFractions = os.path.join( local_path, 'demo_data',  'results_VolumeFractions.pickle' )
    run_commit_StickZeppelinBall()
    results_pickle = os.path.join( out_path, 'COMMIT', 'Results_StickZeppelinBall', 'results.pickle' )
    check_results(results_pickle, ref_pickle_StickZeppelinBall)
    run_commit_BallandStick()
    results_pickle = os.path.join( out_path, 'COMMIT', 'Results_StickZeppelinBall', 'results.pickle' )
    check_results(results_pickle, ref_pickle_BallandStick)
    # run_commit_VolumeFractions()
    # results_pickle = os.path.join( local_path, 'demo_data', 'COMMIT', 'Results_VolumeFractions', 'results.pickle' )
    # check_results(results_pickle, ref_pickle_VolumeFractions)
    sys.exit(0)

if __name__ == "__main__":
    run_tests()
