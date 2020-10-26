
# Change Log
All notable changes to COMMIT will be documented in this file.


## [1.5.0] - 2020-10-23

### Added
- Added the possibility to specify a voxel confidence map

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
