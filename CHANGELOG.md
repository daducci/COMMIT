# Change Log
### All notable changes to `COMMIT` will be documented in this file.

## `v2.0.1`<br>_2023-10-17_
### 🛠️Changed
- Bump `dmri-dicelib` to `v1.0.2`
- Bump `dmri-amico` to `v2.0.1`

### 🐛Fixed
- Add `dtype=np.object_` to inhomogemeous shape arrays
- Fix `egg_base` option in `setup.cfg` (create the `build` folder before running `setup()`)

---
---

## `v2.0.0`<br>_2023-09-14_
### 🛠️Changed
- Default `ndirs=500` in `core.generate_kernels()` and in `trk2dictionary.run()`
- Expire the deprecated `ndirs` parameter in `amico.core.setup()`
- Expire the deprecated `filename_trk` and `gen_trk` parameters in `trk2dictionary.run()`
- Removed unused parameter `ndirs` from `trk2dictionary_c.cpp()`
- Build output goes into `build`
- Switched to proprietary license (see `LICENSE` file)

### ✨Added
- Added clustering in trk2dictionary.run()
    * geometry-based clustering using `blur_clust_thr` parameter
    * anatomical information streamlines clustering based on their endpoints using `blur_clust_groupby` parameter

- Added trk2dictionary.run() parallel computation using `n_threads` parameter
- Changed internal streamlines representation using `blur_core_extent` and `blur_gauss_extent` parameters
- Added possibility to keep temporal files created by parallel trk2dictionary.run() and clustering using `keep_temp` parameter
- Parallel compilation

### 🐛Fixed
- Bugfixes

---
---

## [1.6.5] - 2023-08-03

### Changed
- Build output goes into 'build'
- Require `dmri-amico>=1.3.2,<2.0.0`

### Fixed
- Solve [#121](https://github.com/daducci/COMMIT/issues/121)

## [1.6.4] - 2023-02-14

### Fixed
- Invalid command 'bdist_wheel' error when install with pip<=20.0.2

## [1.6.3] - 2022-10-21

### Changed
- Ensure non-negative values in the 'y' vector before fit

### Added
- 'replace_bad_voxels' parameter in 'load_data()' to replace NaN and Inf values in the signal

## [1.6.2] - 2022-05-05

### Fixed
- Problem with very anisotropic voxelsizes

## [1.6.1] - 2022-01-31

### Changed
- generate_kernels: forcing 'ndirs' to 1 in case of 'VolumeFraction' model
- core.setup: removed 'ndirs' parameter, now precomputes all directions

### Added
- Shortcut 'commit.core.setup' to 'commit.setup'

## [1.6.0] - 2021-12-03

### Changed
- Improved the blur functionality

## [1.5.3] - 2021-12-03

### Added
- Option 'coeffs_format' to 'save_results()'

### Changed
- Install information are stored (and taken from) commit/info.py

## [1.5.2] - 2021-10-27

### Fixed
- Missing dependencies

## [1.5.1] - 2021-10-11

### Added
- Option 'get_normalized' to 'save_coeffs()'


## [1.5.0] - 2021-06-19

### Added
- Possibility to specify a voxel confidence map

### Fixed
- Set license: BSD-3-Clause License
- Indices of groups in case of streamlines that were pre-filtered during trk2dictionary

### Changed
- Loading of nii data using (np.asanyarray( nii.dataobj )) in core and trk2dictionary
- confidence_map from float64 to float32

## [1.4.6] - 2021-03-25

### Fixed
- Length of short segments to neglect
- Streamline weights, in case of blur, are properly scaled

## [1.4.5] - 2021-02-08

### Fixed
- operator.pyxbld: Changed the condition to create a new operator
- trk2dictionary.pyx: Check that the tractogram exists before trying to
            load it and remove the try section
- trk2dictionary.run(): fixed bug with blur parameters and computing the blur

### Added
- core.pyx: Add to the function build_operator the parameter build_dir

### Changed
- core.pyx: The function build_operator checks if the LUT configuration
            changed before build a new operator
- verbose variables in core.pyx and solvers.py changed to be boolean
- trk2dictionary.run(): removed 'points_to_skip' option

## [1.4.4] - 2020-10-28

### Changed
- Option to set one single direction in the resolution of the LUT

## [1.4.3] - 2020-10-22

### Added
- store model parameters to results.pickle

## [1.4.2] - 2020-10-22

### Fixed
- trk2dictionary.run(): check for invalid parameters passed to the blur

## [1.4.1] - 2020-10-21

### Fixed
- operator.pyxbld: Changed the condition to create a new operator

### Added
- COMMIT version is stored in results.pickle
- COMMIT version is stored in output NIFTI files

## [1.4.0.4] - 2020-09-24

### Fixed
- trk2dictionary.run(): bug in the blurring functionality
- trk2dictionary.run(): 'blur_sigma' defaults to 0

## [1.4.0.3] - 2020-08-07

### Fixed
- COMMIT_debugger: compilation problem
- COMMIT_debugger: wrong visualization in Linux

## [1.4.0.2] - 2020-08-07

### Changed
- Moved the documentation to the Wiki

## [1.4.0.1] - 2020-08-03

### Changed
- Updated the installation guide

## [1.4.0.0] - 2020-07-30

### Changed
- trk2dictionary.run(): removed 'gen_trk' option
- save_results(): removed 'save_coeff' and 'save_opt_details' parameters
- save_results(): now saving only streamline_weights.txt (not anymore xic.txt, xec.txt, xiso.txt)
- load_dictionary(): renamed 'use_mask' to 'use_all_voxels_in_mask'
- Removed unused 'dictionary_ndirs.dict' file
- trk2dictionary.run(): 'min_fiber_len' defaults to 0.0 for backward compatibility

### Added
- added 'get_coeffs()' function to get all estimated coefficients
- save_results(): added 'stat_coeffs' parameter for saving streamline weights
- trk2dictionary.run(): added 'max_fiber_len' parameter to discard long streamlines
- load_data(): added 'b0_min_signal' to discard voxels with very low signal

## [1.3.9] - 2020-06-09

### Changed
- Modify setup.py and fix spams dependencies

## [1.3.8] - 2020-05-12

### Changed
- Improvements to the COMMIT_debugger.

## [1.3.7] - 2020-04-25

### Changed
- Adapt demos to use d_perps instead of ICVFs for setting model parameters.

## [1.3.6] - 2020-04-22

### Fixed
- Bug when the selected model has EC compartments but no peaks are provided (in trk2dictionary).

## [1.3.5] - 2020-04-08

### Added
- Parameter 'min_fiber_len' in trk2dictionary to discard streamlines shorter than a given length in mm.

### Fixed
- Bug when 'points_to_skip' was higher then streamline length.
- Few corrections to docstring of trk2dictionary.

## [1.3.4] - 2020-04-02

### Changed
- Added colorized output. NB: needs AMICO 1.2.0 or above.

## [1.3.3] - 2020-03-31

### Added
- Added possibility to save the predicted DW-MR signal in save_results.

### Fixed
- Minor cleanup.


## [1.3.2] - 2020-03-27

### Added
- Check if dictionary (upon loading) and data have the same geometry.

### Fixed
- Bug while saving coefficients in save_results.


## [1.3.1] - 2020-03-27

### Fixed
- Improved the loading of the streamlines in trk2dictionary


## [1.3] - 2019-10-30

This version of COMMIT *is not compatible* with [AMICO](https://github.com/daducci/AMICO) v1.0.1 of below. If you update COMMIT to this version, please update AMICO to version 1.1.0 or above.

### Added
- Changelog file to keep tracking of the COMMIT versions.

### Changed
- Added compatibility with low resolution LUTs.

### Fixed
- Nothing.
