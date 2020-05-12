# COMMIT_debugger

This tool allows one to display, in a common 3D world, all objects considered in the COMMIT framework:

- DWI data (4D NIFTI file)
- its acquisition scheme (Camino format)
- voxelwise main diffusion orientations (4D NIFTI file)
- tractogram (TrackVis format only).

![Application screenshot](https://github.com/daducci/COMMIT/blob/master/doc/COMMIT_debugger.jpg)

## Synopsis

```bash
COMMIT_debugger \
    <dwi> \
    <scheme> \
    <peaks> \
    [-f <tracts>] \
    [-m <map>]
```

- `dwi`, `scheme`: DWI data and corresponding scheme file.
- `peaks`: major directions of the *hindered* water pools in each voxel.
- `tracts`: tractogram to generate the *restricted* contributions of the tracts in each voxel.
- `map`: background map; default is one b0 computed from the DWI data.

A **dropdown menu** will appear with right-click of the mouse.

**NB**: please note that this tool is very rudimental and is released only for debugging purposes.




## How to compile the code

### Install dependencies

You need to install the following libraries:

- [CMake](http://www.cmake.org/) to allow cross-platform compilation;
- [Niftilib](https://sourceforge.net/projects/niftilib/) for reading/writing NIFTI files;
- [Blitz++](http://sourceforge.net/projects/blitz/) for efficient manipulation of multi-dimensional arrays;
- [OpenGL](https://www.opengl.org/) and [GLUT](https://www.opengl.org/resources/libraries/glut/) for 3D visualization.

Please follow the corresponding documentation to install these libraries on your platform. Our code was successfully tested on Linux (Ubuntu 14.04) and OSX (10.9 up to 10.15) systems.

#### OSX with homebrew

If your're using [homebrew](https://brew.sh), then the following should work:

```bash
brew tap brewsci/science
brew install cmake
brew install blitz
brew install niftilib
```

The `OpenGL` and `GLUT` libraries are already provided by the operating system.

###  Compile, build and install

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
make install
```

This will install the binaries into the `/usr/local/bin` folder. This installation path can be changed by rerunning `ccmake ..` and then modifying the `CMAKE_INSTALL_PREFIX` parameter to suit your custom needs.
