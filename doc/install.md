# Compile C++ code

##  Reading a .trk file ([TrackVis](http://www.trackvis.org/docs/?subsect=fileformat) format)

Open the terminal and type:
```bash
cd c++
mkdir build
cd build
ccmake ..
```
Hit `c` (twice) and then `g`. This will create the required makefiles for the compilation.
Once back to the terminal, type:

```bash
make
sudo make install
```
This will install the binaries (e.g. `trk2dictionary`) in your filesystem. Please specify a folder in your `$PATH`, e.g. `/usr/bin`.

## Creating the linear operator **A**

To compile the mex files that implement the linear operator **A**, you need to setup the MEX-compiler in MATLAB as follows:
```matlab
mex -setup
```
In case of issues, please refer to the [MATLAB documentation](http://www.mathworks.ch/ch/help/matlab/matlab_external/what-you-need-to-build-mex-files.html).

# Install the Camino toolkit

`COMMIT` uses the [Camino](http://camino.org.uk) toolkit to generate the kernels.
To install it, please follow the corresponding [documentation](http://cmic.cs.ucl.ac.uk/camino//index.php?n=Main.Installation).

# Setup paths/variables in MATLAB

Before using COMMIT, you need to copy the file `COMMIT_Setup.txt` and rename it to `COMMIT_Setup.m`.
Modify its content to set the paths to your specific needs, e.g. set `CAMINO_path` to the folder where you installed *Camino* in the previous step.

`COMMIT` is ready to use!
