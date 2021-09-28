# COMMIT

The reconstructions recovered with existing tractography algorithms are *not really quantitative* even though diffusion MRI is a quantitative modality by nature. As a matter of fact, several techniques have been proposed in recent years to estimate, at the voxel level, intrinsic micro-structural features of the tissue, such as axonal density and diameter, by using multi-compartment models. COMMIT implements a novel framework to **re-establish the link between tractography and tissue micro-structure**.

Starting from an input set of candidate fiber-tracts, which can be estimated using standard fiber-tracking techniques, COMMIT models the diffusion MRI signal in each voxel of the image as a *linear combination* of the restricted and hindered contributions generated in every location of the brain by these candidate tracts. Then, COMMIT seeks for the effective contribution of each of them such that they globally fit the measured signal at best.

These weights can be easily estimated by solving a convenient **global convex optimization problem** and using efficient algorithms. Results clearly demonstrated the benefits of the proposed formulation, opening new perspectives for a more quantitative and biologically-plausible assessment of the structural connectivity in the brain.


## Main features

- Accepts and works with **any input tractogram** (i.e. set of fiber tracts).
- Can easily implement and consider **any multi-compartment model** available in the literature: possibility to account for restricted, hindered as well as isotropic contributions into the signal forward model.
- Very efficient: the core of the algorithm is implemented in C++ and using **multi-threading programming** for efficient parallel computation.
- **Low memory** consumption using optimized sparse data structures, e.g. it can easily run on a standard laptop with 8GB RAM a full-brain tractogram from the HCP data (1M fibers, 3 shells, 1.25 mm^3 resolution).
- **Soon**: **GPU implementation** for even faster model fitting.


## Documentation

More information/documentation, as well as a series of tutorials, can be found in the [wiki pages](https://github.com/daducci/COMMIT/wiki/Home).

### Installation

To install COMMIT, refer to the [installation guide](https://github.com/daducci/COMMIT/wiki/Installation).

### Getting started

To get started with the COMMIT framework, have a look at [this tutorial](https://github.com/daducci/COMMIT/wiki/Getting-started), which will guide you through the main steps of the processing.

