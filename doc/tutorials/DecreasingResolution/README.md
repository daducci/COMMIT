# Decreasing Resolution

This tutorial illustrates how to **drecrease the number of directions per shell** in COMMIT.

## Starting Script

Check the **getting started tutorial** from the following [link](https://github.com/daducci/COMMIT/tree/master/doc/tutorials/GettingStarted), which contains an example dataset and the following code:

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
commit.core.setup()

mit = commit.Evaluation( '.', 'LausanneTwoShell' )
mit.load_data( 'DWI.nii', 'DWI.scheme' )

mit.set_model( 'StickZeppelinBall' )
d_par = 1.7E-3                          # Parallel diffusivity [mm^2/s]
ICVFs = [ 0.7 ]                         # Intra-cellular volume fraction(s) [0..1]
d_ISOs = [ 1.7E-3, 3.0E-3 ]             # Isotropic diffusivitie(s) [mm^2/s]
mit.model.set( d_par, ICVFs, d_ISOs )

mit.generate_kernels( regenerate=True )
mit.load_kernels()
mit.load_dictionary( 'CommitOutput' )

mit.set_threads()
mit.build_operator()

mit.fit(tol_fun=1e-3, max_iter=200)
mit.save_results()
```

## Change the resolution

In order to change the number of used directions per shell in COMMIT, specify the desired number of directions by changing the value of the parameter `ndirs` in the functions:

- trk2dictionary.run()
- commit.core.setup()
- mit.generate_kernels()

If the value of `ndirs` is not specified, COMMIT will use as default `ndirs = 181*181 = 32,761` directions per shell. 

In this example, we modify the above code to use 500 directions. The modified code should be like this:

```python
from commit import trk2dictionary
trk2dictionary.run(
    filename_trk   = 'LausanneTwoShell/fibers.trk',
    path_out       = 'LausanneTwoShell/CommitOutput',
    filename_peaks = 'LausanneTwoShell/peaks.nii.gz',
    filename_mask  = 'LausanneTwoShell/WM.nii.gz',
    fiber_shift    = 0.5,
    peaks_use_affine = True,
    ndirs = 500
)

import commit
commit.core.setup( ndirs=500 )

mit = commit.Evaluation( '.', 'LausanneTwoShell' )
mit.load_data( 'DWI.nii', 'DWI.scheme' )

mit.set_model( 'StickZeppelinBall' )
d_par = 1.7E-3
ICVFs = [ 0.7 ]
d_ISOs = [ 1.7E-3, 3.0E-3 ]
mit.model.set( d_par, ICVFs, d_ISOs )

mit.generate_kernels( regenerate=True, ndirs=500 )
mit.load_kernels()
mit.load_dictionary( 'CommitOutput' )

mit.set_threads()
mit.build_operator()

mit.fit(tol_fun=1e-3, max_iter=200)
mit.save_results()
```

Please note that the value of `ndirs` must match in all the functions. COMMIT does not support an arbitrary number of directions. The value of `ndirs` has to be one the values in the set: {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 32761 (default)}.