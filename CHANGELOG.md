
# Change Log
All notable changes to COMMIT will be documented in this file.

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
