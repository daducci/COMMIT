# How to "debug" your data

This tool allows one to display in a common 3D space all the objects used by COMMIT (DWI data, streamlines etc...) in order to **spot possible incosistencies between the conventions** of COMMIT and the software that generated the data, e.g. flip in some axes in the DWI data or in the peaks, spatial shift of the streamlines, whether the affine transformation was already applied to the data etc.

**NB**: please note that this tool is very rudimental and is released only for debugging purposes.

![Application screenshot](https://github.com/daducci/COMMIT/blob/master/docs/COMMIT_debugger.jpg)

## Synopsis

```bash
COMMIT_debugger \
    <dwi> \
    <scheme> \
    [-p <peaks> ] \
    [-f <tracts>] \
    [-m <map>]
```

- `dwi`: DWI data (4D NIFTI file);
- `scheme`: corresponding acquisition scheme (Camino format);
- `peaks`: major directions of the *hindered* water pools in each voxel (4D NIFTI file);
- `tracts`: tractogram to generate the *restricted* contributions of the tracts in each voxel (.TRK or .TCK file);
- `map`: background map; default is a b0 computed from the DWI data.

A **dropdown menu** will appear with right-click of the mouse.

## Install dependencies

You need to install the following libraries:

- [CMake](http://www.cmake.org/) to allow cross-platform compilation;
- [Niftilib](https://sourceforge.net/projects/niftilib/) for reading/writing NIFTI files;
- [Blitz++](http://sourceforge.net/projects/blitz/) for efficient manipulation of multi-dimensional arrays;
- [OpenGL](https://www.opengl.org/) and [GLUT](https://www.opengl.org/resources/libraries/glut/) for 3D visualization.

Please follow the corresponding documentation to install these libraries on your platform. Our code was successfully tested on Linux (Ubuntu 14.04) and OSX (10.9 up to 10.15) systems.

### OSX with homebrew

If your're using [homebrew](https://brew.sh), then the following should work:

```bash
brew tap brewsci/science
brew install cmake
brew install blitz
brew install niftilib
```

The `OpenGL` and `GLUT` libraries are already provided by the operating system.

### LINUX with apt-get

If your're using LINUX with apt-get, then the following should work:

```bash
sudo apt-get install nifti-bin
sudo apt-get install libnifti2
sudo apt-get install libnifti-dev
sudo apt-get install freeglut3-dev
sudo apt-get install libblitz0-dev
```

The `OpenGL` libraries are already provided by the operating system.

If you have problems installing `Blitz++`, follow the instructions at https://github.com/blitzpp/blitz

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

This will install the binaries into the `/usr/local/bin` folder. This installation path can be changed by rerunning `ccmake ..` and then modifying the `CMAKE_INSTALL_PREFIX` parameter to suit your custom needs.
