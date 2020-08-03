# Installation

COMMIT is written in [Python](https://www.python.org/) and requires Python 3; version 2 is no longer supported. Also, COMMIT internally makes use of several libraries, e.g.:
- [AMICO](https://github.com/daducci/AMICO) for the generation/rotation of the response functions;
- [DIPY](http://dipy.org) to handle the diffusion data;
- [NiBabel](https://nipy.org/nibabel/) to read/write neuroimaging file formats;
- [Cython](http://cython.org/) to compile parts of the code.

The following installation procedures will install all these required libraries. In case of problems with a package, please refer to the corresponding website or write us for assistance.


## Install from PYPI

Open the system shell and run:

```bash
pip install dmri-commit
```

This will download and install COMMIT from the [Python Package Index](https://pypi.org).

### Import from sources

To be sure of having the lastest updates, download the source code from this repository, open the system shell and go to the folder where you downloaded this repository, and run:

```bash
pip install .
```

## Import COMMIT

COMMIT is now available in your Python interpreter and can be imported as usual:

```python
import commit
```

## Uninstall COMMIT

Open the system shell and run:

```bash
pip uninstall dmri-commit
```

## Notes

1. In order to use COMMIT, it is only necessary to install the Python code; no additional modules (e.g. SPAMS and NODDI) are required.

2. Depending on the forward-model employed, COMMIT may require the [Camino](http://camino.org.uk) toolkit to generate the response functions, e.g., in case of the `Cylinder-Zeppelin-Ball` model. Please follow the corresponding [documentation](http://cmic.cs.ucl.ac.uk/camino//index.php?n=Main.Installation) to install Camino and make sure to include the folder containing the script `datasynth` in your system path.

