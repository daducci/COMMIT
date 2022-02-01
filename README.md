[![PyPI](https://img.shields.io/pypi/v/dmri-commit)](https://pypi.org/project/dmri-commit/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/dmri-commit)](#)
[![LICENSE](https://img.shields.io/github/license/daducci/commit)](https://github.com/daducci/COMMIT/blob/master/LICENSE)
[![GitHub top language](https://img.shields.io/github/languages/top/daducci/commit?color=lightgray)](#)
[![reference](https://img.shields.io/badge/DOI-10.1109/TMI.2014.2352414-red.svg)](https://ieeexplore.ieee.org/document/6884830)

[![GitHub stars](https://img.shields.io/github/stars/daducci/COMMIT?style=social)](#)
[![GitHub forks](https://img.shields.io/github/forks/daducci/COMMIT?style=social)](#)
[![GitHub watchers](https://img.shields.io/github/watchers/daducci/COMMIT?style=social)](#)
[![GitHub followers](https://img.shields.io/github/followers/daducci?style=social)](#)
[![GitHub contributors](https://img.shields.io/github/contributors-anon/daducci/COMMIT?style=social)](#)
[![Twitter Follow](https://img.shields.io/twitter/follow/ADaducci)](https://twitter.com/intent/follow?screen_name=ADaducci)

# COMMIT

The reconstructions recovered with existing tractography algorithms are not really quantitative even though diffusion MRI is a quantitative modality. COMMIT stands for *Convex Optimization for Microstructure Informed Tractography* and is a **powerful framework for enhancing the anatomical accuracy of the reconstructions** by combining tractography with microstructural features of the neuronal tissue.

<img align="right" src="https://github.com/daducci/COMMIT/wiki/images/filtering_methods.png" height="250">

**How?** Starting from an input set of candidate fiber-tracts estimated using standard fiber-tracking techniques, COMMIT models the diffusion MRI signal in each voxel of the image as a *linear combination* of the restricted and hindered contributions generated in every location of the brain by these candidate tracts. Then, COMMIT seeks for the effective contribution of each of them such that they globally fit the measured signal at best.
These weights can be *efficiently estimated by solving a convenient linear system*.

Results clearly demonstrated the benefits of the proposed formulation, opening new perspectives for a more quantitative and biologically-plausible assessment of the structural connectivity of the brain. See the [references](https://github.com/daducci/COMMIT/wiki/References) for more information.

<p align="center">
<img src="https://github.com/daducci/COMMIT/wiki/images/COMMIT_example.png" height="450">
</p>

## Main features

- Very efficient: COMMIT is implemented in Python but the core of the algorithm is implemented in C++ and using **multi-thread programming** for efficient parallel computation.
- Accepts and works with **any input tractogram** (i.e. set of fiber tracts).
- Can easily implement and consider **any multi-compartment model** available in the literature: possibility to account for restricted, hindered as well as isotropic contributions into the signal forward model.
- **Low memory** consumption using optimized sparse data structures, e.g. it can easily run on a standard laptop with 8GB RAM a full-brain tractogram from the HCP data (1M fibers, 3 shells, 1.25 mm^3 resolution).
- **Soon**: **GPU implementation** for even faster model fitting.


## Documentation

More information/documentation, as well as a series of tutorials, can be found in the [wiki pages](https://github.com/daducci/COMMIT/wiki/Home).

### Installation
To install COMMIT, refer to the [installation guide](https://github.com/daducci/COMMIT/wiki/Installation).

### Getting started

To get started with the COMMIT framework, have a look at [this tutorial](https://github.com/daducci/COMMIT/wiki/Getting-started), which will guide you through the main steps of the processing.

