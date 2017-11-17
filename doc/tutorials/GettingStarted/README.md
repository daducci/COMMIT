# Getting started

This tutorial illustrates the basics for using the COMMIT framework to **evaluate the evidence of a tractogram**.

## Download data

Download and extract the **example dataset** from the following [ZIP archive](http://hardi.epfl.ch/static/data/COMMIT_demos/LausanneTwoShell.zip), which contains the following files:

- `DWI.nii`: a diffusion MRI dataset with 100 measurements distributed on 2 shells, respectively at b=700 s/mm^2 and b=2000 s/mm^2;
- `DWI.scheme`: its corresponding acquisition scheme;
- `peaks.nii.gz`: main diffusion orientations estimated with CSD;
- `fibers.trk`: tractogram with about 280K fibers estimated using a streamline-based algorithm;
- `WM.nii.gz`: white-matter mask extracted from an anatomical T1w image.


## Convert the tractogram

Open the *Python interpreter* and go to the folder where you downloaded/unzipped the archive. Then run the following commands:

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
```

The output should be something like this:

```
-> Creating the dictionary from tractogram:
	* Segment position = COMPUTE INTERSECTIONS
	* Fiber shift X    = 0.500 (voxel-size units)
	* Fiber shift Y    = 0.500 (voxel-size units)
	* Fiber shift Z    = 0.500 (voxel-size units)
	* Points to skip   = 0
	* Loading data:
		* tractogram
			- 106 x 106 x 60
			- 2.0000 x 2.0000 x 2.0000
			- 283522 fibers
		* filtering mask
			- 106 x 106 x 60
			- 2.0000 x 2.0000 x 2.0000
		* EC orientations
			- 106 x 106 x 60 x 9
			- 2.0000 x 2.0000 x 2.0000
			- ignoring peaks < 0.10 * MaxPeak
			- flipping axes : [ x=True, y=True, z=False ]
		* output written to "LausanneTwoShell/CommitOutput"
	* Exporting IC compartments:
          [ 283522 fibers, 24388967 segments ]
	* Exporting EC compartments:
          [ 53021 voxels, 145472 segments ]
   [ 44.6 seconds ]
```

Please note that, in this particular example, in order to have all the data in the same reference system we had to:

- apply a translation of half voxel to the fibers.

![Flipping in the data](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/GettingStarted/debugger_screenshot2.jpg)

## Load the diffusion data

Precompute the rotation matrices used internally by COMMIT to create the lookup-tables for the response functions:

```python
import commit
commit.core.setup()
```

Now, load the data:

```python
mit = commit.Evaluation( '.', 'LausanneTwoShell' )
mit.load_data( 'DWI.nii', 'DWI.scheme' )
```

The output should be something like:

```
-> Loading data:
	* DWI signal...
		- dim    = 106 x 106 x 60 x 100
		- pixdim = 2.000 x 2.000 x 2.000
	* Acquisition scheme...
		- 100 samples, 2 shells
		- 10 @ b=0 , 30 @ b=700.0 , 60 @ b=2000.0
   [ 0.2 seconds ]

-> Preprocessing:
	* Normalizing to b0... [ min=0.00,  mean=0.64, max=36.15 ]
	* Merging multiple b0 volume(s)... [ 106 x 106 x 60 x 91 ]
   [ 0.5 seconds ]
```

## Set the forward-model

For this example we made use of the `Stick-Zeppelin-Ball` model described in [(Panagiotaki et al., NeuroImage, 2012)](http://www.sciencedirect.com/science/article/pii/S1053811911011566):

- the contributions of the tracts are modeled as "sticks", i.e. tensors with a given axial diffusivity (`1.7*10^-3 mm^2/s`) but null perpendicular diffusivity;
- extra-cellular contributions are modeled as tensors with the same axial diffusivity as the sticks (1.7*10^-3 mm^2/s) and whose perpendicular diffusivities are calculated with a tortuosity model as a function of the intra-cellular volume fractions (`0.7`);
- isotropic contributions are modeled as tensors with isotropic diffusivities (`1.7*10^-3 mm^2/s` and `3.0*10^-3 mm^2/s`).

Setup the parameters of the model and **generate the lookup-tables**:

```python
mit.set_model( 'StickZeppelinBall' )

d_par = 1.7E-3                          # Parallel diffusivity [mm^2/s]
ICVFs = [ 0.7 ]                         # Intra-cellular volume fraction(s) [0..1]
d_ISOs = [ 1.7E-3, 3.0E-3 ]             # Isotropic diffusivitie(s) [mm^2/s]

mit.model.set( d_par, ICVFs, d_ISOs )
mit.generate_kernels( regenerate=True )
mit.load_kernels()
```

and the output should look like:

```
-> Simulating with "Stick-Zeppelin-Ball" model:
	* 1 stick, 1 extra-cellular and 2 isotropic
	* A_001... [ OK ]
	* A_002... [ OK ]
	* A_003... [ OK ]
	* A_004... [ OK ]
   [ 1.5 seconds ]

-> Resampling kernels for subject "LausanneTwoShell":
	* A_001... [ OK ]
	* A_002... [ OK ]
	* A_003... [ OK ]
	* A_004... [ OK ]
	* Merging multiple b0 volume(s)... [ OK ]
	* Normalizing... [ OK ]
   [ 1.0 seconds ]
```

## Load the sparse data-structure

Load in memory the sparse data-structure previously created with `trk2dicitonary.run()`:

```python
mit.load_dictionary( 'CommitOutput' )
```

The output should show that around 280K fibers have been loaded, in addition to 145K segments for the extra-cellular contributions in the 53K voxels of the white matter:

```
-> Loading the dictionary:
	* segments from the tracts... [ 283522 fibers and 24388967 segments ]
	* segments from the peaks...  [ 145472 segments ]
	* isotropic contributions...  [ 53021 voxels ]
	* post-processing...          [ OK ]
   [ 14.8 seconds ]
```

## Build the linear operator A

Now it's time to build the linear operator **A** to compute the matrix-vector multiplications for solving the linear system. This operator uses information from the segments loaded in the previous step and the lookup-tables for the response functions; it also needs to know the workload to be assigned to each thread durint the multiplications. To this aim, run the following commands:

```python
mit.set_threads()
mit.build_operator()
```

The output should be something similar to this:

```
-> Distributing workload to different threads:
	* number of threads : 4
	* A operator...  [ OK ]
	* A' operator... [ OK ]
   [ 3.5 seconds ]

-> Building linear operator A:
   [ 2.1 seconds ]
```

NB: the *number of threads* is automatically set to the maximum number of cores in the system (4 in this example), but this setting can be manually set.

## Fit the model to the data

To fit the model (`Stick-Zeppelin-Ball` in this case) to the data, simply run:

```python
mit.fit( tol_fun = 1e-3, max_iter = 200 )
```

The optimization progress is displayed by default:

```
-> Fit model using "nnls":
|     ||Ax-y||     |  Cost function    Abs error      Rel error    |     Abs x          Rel x
------|------------------|-----------------------------------------------|------------------------------
1  |   7.5552614e+02  |  2.8540987e+05  4.0602923e+05  1.4226180e+00  |  5.4262515e+01  1.0000000e+00
2  |   6.7997468e+02  |  2.3118278e+05  5.4227093e+04  2.3456372e-01  |  1.6229691e+01  2.6520302e-01
3  |   6.2490484e+02  |  1.9525303e+05  3.5929749e+04  1.8401635e-01  |  1.4457099e+01  2.0335528e-01
...
...
...
137  |   1.4197542e+02  |  1.0078510e+04  1.0588051e+01  1.0505571e-03  |  1.5019784e+00  4.0796383e-03
138  |   1.4190279e+02  |  1.0068201e+04  1.0309090e+01  1.0239257e-03  |  1.4936457e+00  4.0495040e-03
139  |   1.4183213e+02  |  1.0058177e+04  1.0024696e+01  9.9667126e-04  |  1.4848343e+00  4.0182480e-03
< Stopping criterion: REL_OBJ >
[ 00h 07m 04s ]
```

where the columns report, respectively, the *iteration number*, the *cost function* and its *relative change*.

## Storing the results

The results and the output maps can be stored to files as follows:

```python
mit.save_results()
```

As shown in the output, the results are saved in the folder `Results_StickZeppelinBall`:

```
-> Saving results to "Results_StickZeppelinBall/*":
	* configuration and results... [ OK ]
	* fitting errors:
		- RMSE...  [ 0.059 +/- 0.018 ]
		- NRMSE... [ 0.117 +/- 0.037 ]
	* voxelwise contributions:
		- intra-axonal [ OK ]
		- extra-axonal [ OK ]
		- isotropic    [ OK ]
   [ 2.0 seconds ]
```

The following figure shows the **density of the tracts** [(Calamante et al., NeuroImage, 2010)](http://www.sciencedirect.com/science/article/pii/S1053811910009766) of the original tractogram (left) and of its optimized version (right):

![Track-density](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/GettingStarted/density.png)

It is also possible to visualize voxelwise maps of the corresponding contributions of the **extra-cellular space** (left) and other **isotropic contaminations** (right):

![Compartments](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/GettingStarted/compartments.png)

Finally, the **fitting error** in each voxel can also be inspected:

![fitting error](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/GettingStarted/NRMSE.png)
