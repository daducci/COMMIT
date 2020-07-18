# Installation


## Install dependencies

### Python and DIPY

COMMIT is written in [Python](https://www.python.org/) and it internally makes use of the [DIPY](http://dipy.org) library. Also, some parts of the code require to be compiled and this is done via the [Cython](http://cython.org/) module.
Please install and configure all these packages by following the guidelines on the corresponding websites.

> COMMIT was **succesfully tested** on:  
  - OSX 10.10, [Anaconda](http://docs.continuum.io/anaconda/) Python distribution and DIPY 0.9.0dev.

### AMICO

COMMIT shares the code for the generation/rotation of the response-function lookup tables with [AMICO](https://github.com/daducci/AMICO). Please install AMICO following the instructions [here](https://github.com/daducci/AMICO).

> NB: in order to use COMMIT, it is only necessary to install the Python code; no additional modules (e.g. SPAMS and NODDI) are required.

### Camino toolkit

Depending on the forward-model employed, COMMIT can require the [Camino](http://camino.org.uk) toolkit to generate the response functions, e.g. in case of the `Cylinder-Zeppelin-Ball` model.

Please follow the corresponding [documentation](http://cmic.cs.ucl.ac.uk/camino//index.php?n=Main.Installation) to install Camino and make sure to include the folder containing the script `datasynth` in your system path.

### Cuda toolkit (optional)

COMMIT has GPU acceleration support for fast model fitting. In order to use COMMIT with GPU acceleration, it is necessary to install the [CUDA toolkit](https://docs.nvidia.com/cuda/cuda-quick-start-guide/index.html#introduction). After the installation of CUDA, make sure to add CUDA and the CUDA libraries to the PATH. For example, if the CUDA Toolkit 11.0 was installed in the default folder (`/usr/local/cuda-11.0/`), run in a system shell

```bash
export PATH=/usr/local/cuda-11.0/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-11.0/lib64\
                         ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

COMMIT uses the `CUDAHOME` variable to compile some parts of the code. In a system shell, under the previous context, run the command:

```bash
export CUDAHOME=/usr/local/cuda-11.0/
```

**NOTE:** Only NVIDIA GPUs with compute capability >= 5.0 are supported.

**NOTE:** It is recommended to have the latest [NVIDIA drivers](https://www.nvidia.com/Download/index.aspx?lang=en-us) installed.

## Install COMMIT

Open the system shell, go to the folder where you downloaded this repository and run:

```bash
pip install .
```

COMMIT is now available in your Python interpreter and can be imported as usual:

```python
import commit
```

### Uninstall COMMIT

Open the system shell and run:

```bash
pip uninstall commit
```
