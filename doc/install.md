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

## Install COMMIT

Open the system shell, go to the folder where you downloaded this repository and run:

```bash
python setup.py build --force
python setup.py install
```

COMMIT is now available in your Python interpreter and can be imported as usual:

```python
import commit
```
