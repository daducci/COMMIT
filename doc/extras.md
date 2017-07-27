# External C++ tools

This is a collection of additional tools that might be helpful during a tractogram evaluation with the COMMIT framework. At the moment, the following tool(s) are available:

- [COMMIT_debugger](#commit_debugger)

# How to compile the code

## Install dependencies

You will still need to install the following libraries:

- [CMake](http://www.cmake.org/) to allow cross-platform compilation;
- [Niftilib](https://sourceforge.net/projects/niftilib/) for reading/writing NIFTI files;
- [Blitz++](http://sourceforge.net/projects/blitz/) for efficiently manipulating multi-dimensional arrays;
- [OpenGL](https://www.opengl.org/) and [GLUT](https://www.opengl.org/resources/libraries/glut/) for 3D visualization.

Please follow the corresponding documentation to install these libraries on your platform. Our code was successfully tested on Linux (Ubuntu 14.04) and OSX (10.9 and 10.10) systems.

##  Compile, build and install

Open the terminal and type:

```bash
cd extras
mkdir build
cd build
ccmake ..
```

Hit `c` (twice) and then `g`. This will create the required makefiles for the compilation.
Once back to the terminal, type:

```bash
make
```

This will install the binaries into the `build` folder.


#  List of tools

## COMMIT_debugger

This tool allows one to display, in a common 3D world, all ingredients considered in the COMMIT framework:

- DWI data (4D NIFTI file)
- its acquisition scheme (Camino format)
- voxelwise main diffusion orientations (4D NIFTI file)
- tractogram (TrackVis format).

**NB**: please note that this tool is very rudimental and is released only for debugging purposes.

![Application screenshot](https://github.com/daducci/COMMIT/blob/master/doc/COMMIT_debugger.jpg)

### Synopsis

```bash
COMMIT_debugger \
    <dwi> \
    <scheme> \
    <peaks> \
    [-f <tracts>] \
    [-m <map>]
```

- `dwi`, `scheme`: DWI data and corresponding scheme file. Only the signal from the outer shell is plotted.
- `peaks`: directions to be interpreted as the major directions of the *hindered* water pools in each voxel.
- `tracts`: tractogram to generate the *restricted* contributions of the tracts in each voxel.
- `map`: background map; default is the average b0 computed from the data.

### Shortcuts

- `1`, `2`, `3`: show/hide the x-, y-, z-plane, respectively. All plotting is performed on these planes.
- `o`, `O`: decrease/increase the opacity of the three planes
- `m`, `M`: decrease/increase the maximum intensity in the image
- `r`: reset the view to its inital settings

- `s`: show/hide the DWI signal
- `S`: change shell to display
- `X`, `Y`, `Z`: flip the signal glyphs along each axis
- `b`, `B`: change the portion of voxels where to plot the glyphs

- `p`: show/hide the peaks
- `a`: do/don't use the affine transformation
- `x`, `y`, `z`: flip the peaks along each axis
- `t`, `T`: decrease/increase the peaks threshold
- `n`: normalize the length of the peaks
- `w`, `W`: decrease/increase the line-width of the peaks/glyphs

- `k`, `K`: decrease/increase the threshold for coloring the peaks depending on their norm (range is hard-coded to 0-20). If k=0, then the peaks are colored depending on their directions as usual.
- `l`: change the colormap (lut) for the peaks

- `f`: show/hide the tracts
- `c`, `C`: decrease/increase the slab region in which plotting the fibers
- `space`: alternate between "crop fibers to slab around the planes" and "crop fibers in the current voxel (yellow)"
