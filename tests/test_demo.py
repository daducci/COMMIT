import os
import pickle
import sys

import amico

from commit import trk2dictionary
import commit

from dicelib.ui import setup_logger


logger = setup_logger('test_demo')

commit.core.setup()


def run_commit_StickZeppelinBall(local_path):
    trk2dictionary.run(
        filename_tractogram=os.path.join(local_path, 'demo01_fibers.tck'),
        filename_peaks=os.path.join(local_path, 'peaks.nii.gz'),
        filename_mask=os.path.join(local_path, 'WM.nii.gz'), 
        fiber_shift=0.5,
        peaks_use_affine=True
    )

    amico.util.fsl2scheme(os.path.join(local_path, 'bvals.txt'), 
                          os.path.join(local_path, 'bvecs.txt'), 
                          os.path.join(local_path, 'DWI.scheme'))

    mit = commit.Evaluation(local_path, '.')
    mit.load_data(os.path.join(local_path, 'DWI.nii.gz'),
                  os.path.join(local_path, 'DWI.scheme'))

    mit.set_model('StickZeppelinBall')
    d_par   = 1.7E-3            # Parallel diffusivity [mm^2/s]
    d_perps = [0.51E-3]         # Perpendicular diffusivity(s) [mm^2/s]
    d_isos  = [1.7E-3, 3.0E-3]  # Isotropic diffusivity(s) [mm^2/s]
    mit.model.set(d_par, d_perps, d_isos)
    mit.generate_kernels(regenerate=True)
    mit.load_kernels()

    mit.load_dictionary()

    mit.set_threads()
    mit.build_operator()

    mit.fit(tol_fun=1e-3, max_iter=1000)
    mit.save_results()


def run_commit_BallandStick(local_path):
    trk2dictionary.run(
        filename_tractogram=os.path.join(local_path, 'demo01_fibers.tck'),
        filename_peaks=os.path.join(local_path, 'peaks.nii.gz'),
        filename_mask=os.path.join(local_path, 'WM.nii.gz'), 
        fiber_shift=0.5,
        peaks_use_affine=True
    )

    amico.util.fsl2scheme(os.path.join(local_path, 'bvals.txt'), 
                          os.path.join(local_path, 'bvecs.txt'), 
                          os.path.join(local_path, 'DWI.scheme'))

    mit = commit.Evaluation(local_path, '.')
    mit.load_data(os.path.join(local_path, 'DWI.nii.gz'),
                  os.path.join(local_path, 'DWI.scheme'))

    mit.set_model('StickZeppelinBall')
    d_par   = 1.7E-3            # Parallel diffusivity [mm^2/s]
    d_perps = []                # Perpendicular diffusivity(s) [mm^2/s]
    d_isos  = [1.7E-3, 3.0E-3]  # Isotropic diffusivity(s) [mm^2/s]
    mit.model.set(d_par, d_perps, d_isos)
    mit.generate_kernels(regenerate=True)
    mit.load_kernels()

    mit.load_dictionary()

    mit.set_threads()
    mit.build_operator()

    mit.fit(tol_fun=1e-3, max_iter=1000)
    mit.save_results()


def run_commit_VolumeFractions(local_path):
    trk2dictionary.run(
        filename_tractogram=os.path.join(local_path, 'demo01_fibers.tck'),
        filename_peaks=os.path.join(local_path, 'peaks.nii.gz'),
        filename_mask=os.path.join(local_path, 'WM.nii.gz'),
        fiber_shift=0.5,
        peaks_use_affine=True
    )

    amico.util.fsl2scheme(os.path.join(local_path, 'bvals.txt'), 
                          os.path.join(local_path, 'bvecs.txt'), 
                          os.path.join(local_path, 'DWI.scheme'))

    mit = commit.Evaluation(local_path, '.')
    mit.load_data(os.path.join(local_path, 'DWI.nii.gz'),
                  os.path.join(local_path, 'DWI.scheme'))

    mit.set_model('StickZeppelinBall')
    d_par   = 1.7E-3            # Parallel diffusivity [mm^2/s]
    d_perps = []                # Perpendicular diffusivity(s) [mm^2/s]
    d_isos  = [1.7E-3, 3.0E-3]  # Isotropic diffusivity(s) [mm^2/s]
    mit.model.set(d_par, d_perps, d_isos)
    mit.generate_kernels(regenerate=True)
    mit.load_kernels()

    mit.load_dictionary()

    mit.set_threads()
    mit.build_operator()

    mit.fit(tol_fun=1e-3, max_iter=1000)
    mit.save_results()


def check_results(pickle_result, ref_pickle):
    with open(pickle_result, 'rb') as f:
        data = pickle.load(f)
    with open(ref_pickle, 'rb') as f:
        ref_data = pickle.load(f)

    result_optimization = data[0]["optimization"]
    ref_optimization = ref_data[0]["optimization"]

    try:
        assert abs(result_optimization["fit_details"]["residual"] - ref_optimization["fit_details"]["residual"]) < 1e-4
        assert abs(result_optimization["fit_details"]["regterm"] - ref_optimization["fit_details"]["regterm"]) < 1e-4
        assert abs(result_optimization["fit_details"]["cost_function"] - ref_optimization["fit_details"]["cost_function"]) < 1e-4
        assert abs(result_optimization["fit_details"]["abs_cost"] - ref_optimization["fit_details"]["abs_cost"]) < 1e-4
        assert abs(result_optimization["fit_details"]["rel_cost"] - ref_optimization["fit_details"]["rel_cost"]) < 1e-4
        assert abs(result_optimization["fit_details"]["abs_x"] - ref_optimization["fit_details"]["abs_x"]) < 1e-4
        assert abs(result_optimization["fit_details"]["rel_x"] - ref_optimization["fit_details"]["rel_x"]) < 1e-4
        assert result_optimization["fit_details"]["iterations"] == ref_optimization["fit_details"]["iterations"]

    except AssertionError:
        logger.error("Results do not match")
        sys.exit(1)


def run_tests():
    local_path = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), 'demo_data')
    ref_pickle_StickZeppelinBall = os.path.join(
        local_path, 'ref_results', 'ref_results_StickZeppelinBall.pickle')
    ref_pickle_BallandStick = os.path.join(
        local_path, 'ref_results', 'ref_results_BallandStick.pickle')
    # ref_pickle_VolumeFractions = os.path.join( local_path, 'results_VolumeFractions.pickle' )
    run_commit_StickZeppelinBall(local_path)
    results_pickle = os.path.join(
        local_path, 'COMMIT', 'Results_StickZeppelinBall', 'results.pickle')
    check_results(results_pickle, ref_pickle_StickZeppelinBall)
    run_commit_BallandStick(local_path)
    results_pickle = os.path.join(
        local_path, 'COMMIT', 'Results_StickZeppelinBall', 'results.pickle')
    check_results(results_pickle, ref_pickle_BallandStick)
    # run_commit_VolumeFractions()
    # results_pickle = os.path.join( local_path, 'demo_data', 'COMMIT', 'Results_VolumeFractions', 'results.pickle' )
    # check_results(results_pickle, ref_pickle_VolumeFractions)
    sys.exit(0)


if __name__ == "__main__":
    run_tests()
