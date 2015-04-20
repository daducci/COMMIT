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

## How to cite COMMIT

**COMMIT: Convex Optimization Modeling for Microstructure Informed Tractography**  
Alessandro Daducci, Alessandro Dal Palú, Alia Lemkaddem, Jean-Philippe Thiran  
*IEEE Transactions on Medical Imaging* 34(1) 246-257, 2015  
[Link to publisher](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6884830)

**A convex optimization framework for global tractography**  
Alessandro Daducci, Alessandro Dal Palú, Alia Lemkaddem, Jean-Philippe Thiran  
*IEEE 10th International Symposium on Biomedical Imaging (ISBI)* 524-527, 2013  
[Link to publisher](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6556527)

## Installation
To install COMMIT, please refer to the [installation guide](doc/install.md).

More information/documentation can be found in the [`doc`](doc/) folder.

## Getting started

Tutorials are provided in the [`doc/tutorials`](doc/tutorials/) folder to help you get started with the COMMIT framework.

