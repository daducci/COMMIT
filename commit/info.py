# -*- coding: UTF-8 -*-

# Format version as expected by setup.py (string of form "X.Y.Z")
_version_major = 1
_version_minor = 6
_version_micro = 1
_version_extra = '' #'.dev'
__version__    = "%s.%s.%s%s" % (_version_major,_version_minor,_version_micro,_version_extra)

NAME                = 'dmri-commit'
DESCRIPTION         = 'Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)'
LONG_DESCRIPTION    = """
========
 COMMIT
========

The reconstructions recovered with existing tractography algorithms are not really quantitative even though diffusion MRI is a quantitative modality by nature. As a matter of fact, several techniques have been proposed in recent years to estimate, at the voxel level, intrinsic micro-structural features of the tissue, such as axonal density and diameter, by using multi-compartment models.
Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT) implements a novel framework to re-establish the link between tractography and tissue micro-structure.

Starting from an input set of candidate fiber-tracts, which can be estimated using standard fiber-tracking techniques, COMMIT models the diffusion MRI signal in each voxel of the image as a linear combination of the restricted and hindered contributions generated in every location of the brain by these candidate tracts. Then, COMMIT seeks for the effective contribution of each of them such that they globally fit the measured signal at best.

These weights can be easily estimated by solving a convenient global convex optimization problem and using efficient algorithms. Results clearly demonstrated the benefits of the proposed formulation, opening new perspectives for a more quantitative and biologically-plausible assessment of the structural connectivity in the brain.
"""
URL                 = 'https://github.com/daducci/COMMIT'
DOWNLOAD_URL        = "N/A"
LICENSE             = 'BSD license'
AUTHOR              = 'Alessandro Daducci'
AUTHOR_EMAIL        = 'alessandro.daducci@univr.it'
PLATFORMS           = "OS independent"
MAJOR               = _version_major
MINOR               = _version_minor
MICRO               = _version_micro
VERSION             = __version__