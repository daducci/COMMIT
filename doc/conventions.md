# DWI data

The **diffusion MRI signal** is stored as a 4D [NIFTI](http://nifti.nimh.nih.gov/) file, where:

- the *first three dimensions* define the spatial locations of voxels, ie *x*, *y* and *z*;
- the *fourth dimension* contains the diffusion signal for each voxel (*x*,*y*,*z*).

A **companion scheme file** defines all the information about the diffusion acquisition protocol. The scheme is a a text file and can be specified in two formats:

- as a *Nx4 matrix*, where the first three columns are the gradient directions and the fourth contains their b-value (*s/mm^2*).
- as a *Nx7 matrix*, where the first three columns are the gradient directions and the remaining four define the gradient strength (*T/m*), big delta (*s*), small delta (*s*) and echo time (*s*), respectively.

The scheme files follow the [Camino conventions](http://cmic.cs.ucl.ac.uk/camino/index.php?n=Docs.SchemeFiles); here, the header line (eg. `VERSION: BVECTOR`) can be omitted.

# Acquisition protocol

At the moment, COMMIT assumes that the data is acquired with a **multi-shell acquisition protocol** (in addition to an arbitrary number of b0 images). The reason for this is that COMMIT creates internal lookup-tables (LUT) for efficiently computing the matrix-vector multiplications with the linear operator **A**. These LUT are created by performing the rotations in spherical harmonics (SH) space, as this procedure is much faster than generating the single response-functions in all possible orientations. As a consequence, the current version of COMMIT works only with data acquired on (multiple) shells.

For other non shell-like acquisition schemes, e.g. DSI, this procedure is not possible. At the moment, the software threats all the diffusion gradients as belonging to distinct shells (very inefficient).

# Folder structure

Usually, all subjects belonging to the same study are acquired with the *same acquisition scheme*. For this reason, the software implicitly assumes a folder structure as follows:

```
    ├── data
        ├── Study_01    --> subjects acquired with protocol "Study_01"
            ├── Subject_01
            ├── Subject_02
            ├── ...
        ├── Study_02    --> subjects acquired with protocol "Study_02"
            ├── Subject_01
            ├── Subject_02
            ├── ...
        ├── ...
```

For data following this convention, the internal LUT are generated only once for the acquisition scheme of the study and then, for each subject, they are simply projected back very efficiently from the SH space to the specific subject space. In all other cases, the LUT must be regenerated each time.

