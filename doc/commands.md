#  C++ commands

- [trk2dictionary](#trk2dictionary)

- [COMMIT_debugger](#COMMIT_debugger)


## trk2dictionary

Giver an input tractogram, and possibly a peaks file representing the major directions of the extra-celullar water pool in each voxel, it constructs the sparse data-structure that is needed to compute the matrix-vector multiplications **Ax** and **A'y**. The fibers and the extra-cellular directions can be estimated/construcuted with any tool of choice.

```bash
trk2dictionary \
    -i <path/filename> \
    -o <path> \
    [-p <path/filename>] \
    [-w <path/filename>] \
    [-x] [-y] [-z] \
    [-t <float>] \
    [-n <number>] \
    [-s <offset>] \
    [-c] \
    [-h]
```

- `-i`: .trk file with the input tractogram (*mandatory*)
- `-o`: folder where to store the sparse data structure containing the information to build the linear operator **A** (*mandatory*)
- `-p`: .nii file containing the peaks to build the extra-cellular compartment in each voxel (*default: no extra-cellular compartment*)
- `-w`: .nii file containing a binary mask to restrict the analysis only in a portion of the brain (*default: use all voxels covered by fibers*)
- `-x`, `-y`, `-z`, flipz the peaks data along the three axes (*default: no flipping*)
- `-t`: threshold to discard small peaks (*default: 0.1*)
- `-n`: skip *n* points at beginning/end of the fiber (*default: 0*)
- `-s`: apply a shift (in voxel units) to fiber coordinates (*default: 0*)
- `-c`: compute intersections of each fiber segment with the voxel lattice (*default: use centroid of the segment*)


## COMMIT_debugger

This tool allows to display all ingredients used in the COMMIT model in a common 3D space: DWI data (NIFTI format), its acquisition scheme (CAMINO format), extra-cellular directions (NIFTI file) and tractogram (TrackVis format).

**NB**: please note that this tools is released only for debugging purposes and it is very rudimental. For this reason, it is *not compiled by default*; to compile it, be sure to install/configure the OpengGL libraries for your platform, open the file `CMakeLists.txt` and uncomment the corresponding line, ie remove the '#' character from the line `#ADD_SUBDIRECTORY( COMMIT_debugger )`.

```bash
COMMIT_debugger \
    <dwi> \
    <scheme> \
    <peaks> \
    <tracts>
```

- `dwi`, `scheme`: DWI data and corresponding scheme file (NIFTI format). Only the signal from the outer shell is plotted.
- `peaks`: directions to be interpreted as the major directions of the *hindered* water pools in each voxel (NIFTI format).
- `tracts`: tractogram to generate the *restricted* contributions of the tracts in each voxel.

### Shortcuts 

- `1`, `2`, `3`: show/hide the x-, y-, z-plane, respectively. All plotting is performed on these planes.
- `o`, `O`: decrease/increase the opacity of the three planes
- `m`, `M`: decrease/increase the maximum intensity in the image
- `r`: reset the view to its inital settings

- `s`: show/hide the DWI signal
- `X`, `Y`, `Z`: flip the signal glyphs along each axis
- `b`, `B`: change the portion of voxels where to plot the glyphs

- `p`: show/hide the peaks
- `x`, `y`, `z`: flip the peaks along each axis
- `t`, 'T': decrease/increase the peaks threshold
- `n`: normalize the length of the peaks
- `w`, `W`: decrease/increase the line-width of the peaks/glyphs

- `f`: show/hide the tracts
- `c`, `C`: decreases/increases the slab region in which plotting the fibers
- `space`: alternate between "crop fibers to slab around the planes" and "crop the fibers in the current voxel (yellow)"
