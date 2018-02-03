
You can find the ipython notebook version of this tutorial [at this link](tutorial_solvers.ipynb).

# Advanced solvers

This tutorial shows how to exploit the advanced features of the COMMIT framework from the side of the **optimisation problem**. The general formulation is the following:
\begin{equation}
x^* = \arg\min_{x\in R^n_+} \frac12 \|Ax-y\|_2^2 + \lambda_{IC}\Omega_{IC}(x) + \lambda_{EC}\Omega_{EC}(x) + \lambda_{ISO}\Omega_{ISO}(x),
\end{equation}
where $A$ is the COMMIT dictionary, $n$ is defined in such a way that the product $Ax$ makes sense and $y$ is the datum that we want to fit. The three regularisation terms allow us to exploit ***distinct penalties for each compartment***.

*Note*: before exploring this tutorial, you should follow the [Getting Started](https://github.com/daducci/COMMIT/tree/master/doc/tutorials/GettingStarted) tutorial.


### Download and unpack the data

Download and extract the **example dataset** from the following [ZIP archive](http://hardi.epfl.ch/static/data/COMMIT_demos/LausanneTwoShell.zip), which contains the following files:

- `DWI.nii`: a diffusion MRI dataset with 100 measurements distributed on 2 shells, respectively at b=700 s/mm^2 and b=2000 s/mm^2;
- `DWI.scheme`: its corresponding acquisition scheme;
- `peaks.nii.gz`: main diffusion orientations estimated with CSD;
- `fibers.trk`: tractogram with about 280K fibers estimated using a streamline-based algorithm;
- `WM.nii.gz`: white-matter mask extracted from an anatomical T1w image.


<span style="color:crimson">**Make sure that your working directory is the folder where you unzipped the downloaded archive.**</span>


```python
path_to_the_directory_with_the_unzipped_archive = '.' # edit this
cd path_to_the_directory_with_the_unzipped_archive
```

### Load the usual COMMIT structure


```python
from commit import trk2dictionary

trk2dictionary.run(
    filename_trk   = 'LausanneTwoShell/fibers.trk',
    path_out       = 'LausanneTwoShell/CommitOutput',
    filename_peaks = 'LausanneTwoShell/peaks.nii.gz',
    filename_mask  = 'LausanneTwoShell/WM.nii.gz',
    fiber_shift    = 0.5,
    peaks_use_affine = True
)

import commit
mit = commit.Evaluation( '.', 'LausanneTwoShell' )
mit.load_data( 'DWI.nii', 'DWI.scheme' )

mit.set_model( 'StickZeppelinBall' )

d_par = 1.7E-3              # Parallel diffusivity [mm^2/s]
ICVFs = [ 0.7 ]             # Intra-cellular volume fraction(s) [0..1]
d_ISOs = [ 1.7E-3, 3.0E-3 ] # Isotropic diffusivitie(s) [mm^2/s]

mit.model.set( d_par, ICVFs, d_ISOs )
mit.generate_kernels( regenerate=True )
mit.load_kernels()

mit.load_dictionary( 'CommitOutput' )
mit.set_threads()
mit.build_operator()
```

### Perform clustering of the streamlines

You will need `dipy`, which is among the requirements of COMMIT, hence there should be no problem.

The `threshold` parameter has to be tuned for each brain. Do not consider our choice as a standard one.


```python
from nibabel import trackvis as tv
fname='LausanneTwoShell/fibers.trk'
streams, hdr = tv.read(fname)
streamlines = [i[0] for i in streams]

from dipy.segment.clustering import QuickBundles
threshold = 15.0
qb = QuickBundles(threshold=threshold)
clusters = qb.cluster(streamlines)

import numpy as np
structureIC = np.array([c.indices for c in clusters])
weightsIC   = np.array([1.0/np.sqrt(len(c)) for c in structureIC])
```

Notice that we defined `structure_IC` as a `numpy.array` that contains a list of lists containing the indices associated to each group. We know it sounds a little bit bizarre but it computationally convenient.

### Define the regularisation term
Each compartment must be regularised separately. The user can choose among the following penalties:

- $\sum_{g\in G}w_g\|x_g\|_k$ : `commit.solvers.group_sparsity` with $k\in \{2, \infty\}$ (only for IC compartment)

- $\|x\|_1$ : `commit.solvers.norm1`

- $\|x\|_2$ : `commit.solvers.norm2`

- $\iota_{\ge 0}(x)$ : `commit.solvers.non_negative` (Default for all compartments)

If the chosen regularisation for the IC compartment is $\sum_{g\in G}\|x_g\|_k$, we can define $k$ via the `group_norm` field, which must be one between

- $\|x\|_2$ : `commit.solvers.norm2` (Default)

- $\|x\|_\infty$ : `commit.solvers.norminf`

In this example we consider the following penalties:

- Intracellular: group sparsity with 2-norm of each group

- Extracellular: 2-norm

- Isotropic: 1-norm


```python
regnorms = [commit.solvers.group_sparsity, commit.solvers.norm2, commit.solvers.norm1]

group_norm = 2 # each group is penalised with its 2-norm
```

The regularisation parameters are specified within the lambdas field. Again, do not consider our choice as a standard one.


```python
lambdas = [10.,10.,10.]
```

### Call the constructor of the data structure


```python
regterm = commit.solvers.init_regularisation(mit,
                                             regnorms    = regnorms,
                                             structureIC = structureIC,
                                             weightsIC   = weightsIC,
                                             group_norm  = group_norm,
                                             lambdas     = lambdas)
```

### Call the fit function to perform the optimisation


```python
mit.fit(regularisation=regterm, max_iter=1000)
```

### Save the results


```python
suffix = 'IC'+str(regterm[0])+'EC'+str(regterm[1])+'ISO'+str(regterm[2])
mit.save_results(path_suffix=suffix)
```
