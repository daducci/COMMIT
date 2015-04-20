# Forward models

COMMIT is *not* a model, but a *framework*: it allows the combination of a tractogram with any generic multi-compartment model, accounting for possible signal contributions arising from **restricted**, **hindered** and **isotropic** water pools.

Two classical models are already included in the current version of the software, i.e. `Stick-Zeppelin-Ball` and `Cylinder-Zeppelin-Ball` defined in [(Panagiotaki et al., Neuroimage, 2012)](http://www.sciencedirect.com/science/article/pii/S1053811911011566). Each compartment can be selectively enabled/disabled; this means, for example, that many other models are already implicitly included, such as the `Ball&Stick` that can be obtained from the `Stick-Zeppelin-Ball` model by disabling the `Zeppelin` contributions.

Additional multi-compartment models can be easily added. A model is defined as a *class* in the file `models.py` and must expose (at least) the following methods for:

1) Setting the specific **parameters of the model**; the method must be named `set`, but the actual signature is model-dependent:

```python
def set( self, ... ) :
```

2) **Generating high-resolution response functions** and rotate them (in SH space), with the following signature:

```python
def generate( self, out_path, scheme, aux, idx_in, idx_out ) :

Parameters
----------
out_path : string
    The path where to store the rotated kernels.

scheme : Scheme class
    The original acquisition scheme.

aux : dictionary
    Auxiliary data structures needed to rotate functions in SH space.

idx_in : list of list
    Index of samples in input kernel belonging to each shell.

idx_out : list of list
    Index of samples in output kernel belonging to each shell.
```

3) **Projecting the response functions** from SH space to signal space of the subject, with the following signature:

```python
def resample( self, in_path, idx_out, Ylm_out ) :

Parameters
----------
in_path : string
    The path where the rotated kernels in SH space are stored

idx_out : list of list
    Index of samples in output kernel belonging to each shell

Ylm_out : numpy.array
    Matrix to project back all shells from SH space to signal space (of the subject)

Returns
-------
KERNELS : dict
    Contains all the response functions projected to the signal space of the subject
```