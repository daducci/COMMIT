# External C++ tools

This is a collection of additional tools that might be helpful during a tractogram evaluation with the COMMIT framework. At the moment, the following tool(s) are availablae:

- [COMMIT_debugger](#commit_debugger)

# How to compile the code

## Install dependencies

You will still need to install the following libraries:

- [CMake](http://www.cmake.org/) to allow cross-platform compilation;
- [Niftilib](https://sourceforge.net/projects/niftilib/) for reading/writing NIFTI files;
- [Blitz++](http://sourceforge.net/projects/blitz/) for efficiently manipulating multi-dimensional arrays;
- [OpenGL](https://www.opengl.org/) and [GLUT](https://www.opengl.org/resources/libraries/glut/) for 3D visualization.

Please follow the corresponding documentation to install these libraries on your platform. Our code was successfully tested on Linux (Ubuntu 14.04) and OSX (10.9 up to 10.15) systems.

### OSX with homebrew

If your're using [homebrew](https://brew.sh), then the following should do:

```bash
brew tap brewsci/science
brew install cmake
brew install blitz
brew install niftilib
```

The `OpenGL` and `GLUT` libraries are already provided by the operating system.

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
make install
```

This will install the binaries into the `/usr/local/bin` folder. This installation path can be changed by rerunning `ccmake ..` and then modifying the `CMAKE_INSTALL_PREFIX` parameter to suits your custom needs.


#  List of tools

## COMMIT_debugger

This tool allows one to display, in a common 3D world, all ingredients considered in the COMMIT framework:

- DWI data (4D NIFTI file)
- its acquisition scheme (Camino format)
- voxelwise main diffusion orientations (4D NIFTI file)
- tractogram (TrackVis format).

**NB**: please note that this tool is very rudimental and is released only for debugging purposes.

![Application screenshot](https://github.com/daducci/COMMIT/
ob/master/doc/COMMIT_debugger.jpg)

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

A **dropdown menu** will appear with right-click of the mouse.
