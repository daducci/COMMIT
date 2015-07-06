# Comparison between COMMIT and LIFE (on STN96 data from the LiFE original publication)

In this example, we show the **importance of using adequate multi-compartment models** to be able to effectively evaluate the evidence of a tractogram, i.e. set of fiber tracts. For more information, please refer to the following abstract (#3148):

> **On evaluating the accuracy and biological plausibility of diffusion MRI tractograms**  
> *David Romascano, Alessandro Dal Palú, Jean-Philippe Thiran, and Alessandro Daducci*

that recently has been **specially selected for a power pitch presentation** (less than 3% of submitted papers) at the annual *International Society for Magnetic Resonance in Medicine* (ISMRM) meeting in Toronto (30/05-05/06 2015)!

To this aim, we evaluate the performance of the the *LiFE* model that was recently described in [(Pestilli et al, Nat Methods, Sep 2014)](http://www.nature.com/nmeth/journal/v11/n10/abs/nmeth.3098.html). To model the diffusion MR signal in each voxel, **LiFE considers only contributions arising from the tracts** crossing a particular voxel (i.e. restricted diffusion). Notably, *LiFE* does not consider the extra-cellular space around the fibers (i.e. hindered diffusion) and all the partial volume that can occur with gray matter and CSF. On the other hand, *COMMIT* can account for all possible compartments that contribute to the signal in a voxel.

As a matter of fact, *LiFE* can be considered a **special case** of our *COMMIT* framework
 [(Daducci et al, IEEE Trans Med Imaging, Aug 2014)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6884830); in fact, it corresponds to the preliminary formulation we had proposed in [(Daducci et al, ISBI, Apr 2013)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6556527). Hence, the *COMMIT* framework can be used to evaluate both approaches.
 

## Download the data

1. Create the folder `STN96/scan1` in your data directory.

2. Download the original DWI data from [here](https://stacks.stanford.edu/file/druid:cs392kv3054/life_demo_data.tar.gz).

3. Extract the file `life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz` from the archive, unzip it and move it to the `scan1` folder with the name `DWI.nii`, i.e.

 ```bash
 gunzip life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz
 mv life_demo_scan1_subject1_b2000_150dirs_stanford.nii STN96/scan1/DWI.nii
 ```

4. Download precomputed reconstructions from [here](http://hardi.epfl.ch/static/data/COMMIT_demos/STN96_scan1.zip). This archive contains a CSD reconstruction + probabilistic tracking performed according to the experimental setting used in the corresponding publication (e.g. CSD implemented in *MrTrix* and probabilistic tracking with 500000 tracts). 

5. Unzip the file content into the `STN96/scan1` folder in your data directory.

## Convert the tracts to the internal data structure

Use the module `commit.trk2dictionary` to convert the tracts contained in the file `STN96/scan1/Tracking/PROB/fibers.trk` to the internal sparse data structure used by the COMMIT framework. In the Python shell run:

```python
from commit import trk2dictionary
trk2dictionary.run(
    filename_trk   = 'STN96/scan1/Tracking/PROB/fibers.trk',
    path_out       = 'STN96/scan1/Tracking/PROB',
    filename_peaks = 'STN96/scan1/CSD/CSD_FODsh_peaks.nii',
    filename_mask  = 'STN96/scan1/WM.nii'
)
```
This will create the necessary data structure (`STN96/scan1/Tracking/PROB/dictionary_*`) containing all the details of the tracts.

NB: the output tractogram from *MrTrix* has been already converted to the format accepted by *COMMIT*, i.e. [TrackVis format](http://www.trackvis.org/docs/?subsect=fileformat) with fiber coordinates (in mm) in image space, where the coordinate (0,0,0) corresponds to the corner of first voxel.

## Process data with COMMIT

Setup a *Microstructure Informed Tractography (mit)* experiment and load the data:

```python
import commit
commit.core.setup()

mit = commit.Evaluation( 'STN96', 'scan1' )
mit.CONFIG['doNormalizeSignal'] = False
mit.CONFIG['doDemean'] = False
mit.load_data( 'DWI.nii', 'DWI.scheme' )
```

Calculate the **kernels** corresponding to the different compartments. In this example, we use 1 kernel for intra-axonal compartment (i.e. Stick), 1 for extra-axonal space (i.e. Zeppelin) and 2 to model partial volume with gray matter and CSF:

```python
mit.set_model( 'StickZeppelinBall' )
mit.model.set( 1.7E-3, [ 0.7 ], [ 1.7E-3, 3.0E-3 ] )
mit.generate_kernels( regenerate=True )
mit.load_kernels()
```

Load the *sparse data structure* that represents the linear operators **A** and **A'**:

```python
mit.load_dictionary( 'Tracking/PROB' )
```

**Solve** the inverse problem according to the *COMMIT* model:

```python
mit.set_threads()
mit.build_operator()
mit.fit()
mit.save_results( 'COMMIT' )
```

Result will be stored in `STN96/scan1/Tracking/PROB/Results_StickZeppelinBall_COMMIT/`.


## Process data with LiFE

Setup and load the data; this time, however, we will apply the *demeaning procedure* used in *LiFE* to both data and kernels:

```python
mit = commit.Evaluation( 'STN96', 'scan1' )
mit.CONFIG['doNormalizeSignal'] = False
mit.CONFIG['doDemean'] = True
mit.load_data( 'DWI.nii', 'DWI.scheme' )
```

Calculate the **kernel** corresponding to the intra-cellular compartment (the only one considered in *LiFE*); in this example, thus, we use only 1 kernel for intra-axonal compartment (i.e. Stick):

```python
mit.set_model( 'StickZeppelinBall' )
mit.model.set( 1.7E-3, [], [] )
mit.generate_kernels( regenerate=True )
mit.load_kernels()
```

Load the *sparse data structure* that represents the linear operators **A** and **At**:

```python
mit.load_dictionary( 'Tracking/PROB' )
```

**Solve** the inverse problem according to the *LiFE* model:

```python
mit.set_threads()
mit.build_operator()
mit.fit()
mit.save_results( 'LIFE' )
```

Result will be stored in `STN96/scan1/Tracking/PROB/Results_StickZeppelinBall_LIFE/`.


## Compare the two models

Let's first analyze the performance of the two approaches in the **native space in which the two models perform the fitting**. In fact, *LiFE* does not fit the model to the acquired diffusion MR signal, but rather to the signal after removing the mean value in each voxel, i.e. demeaned signal.

It is important to note that as the two models actually work in different spaces (different values), a normalization of the error metrics is required in order to compare their accuracy in explaining the measured diffusion MR data. To this aim, we use the *Normalized RMSE (NRMSE)* as quality measure. Please note that the normalization constant used in each voxel quantifies the magnitude of the data in that voxel, hence the values are expressed as *percentage error* with respect to the actual measurements considered in the voxel, i.e. measured diffusion MR signal for *COMMIT* and demeaned signal for *LiFE*.

We then load the *NRMSE* fit error of the two models, as follows:

```python
import nibabel, numpy, pylab

# load error maps
niiERR_L = nibabel.load( 'STN96/scan1/Tracking/PROB/Results_StickZeppelinBall_LIFE/fit_NRMSE.nii.gz' );
niiERR_L_img = 100.0 * niiERR_L.get_data()
niiERR_C = nibabel.load( 'STN96/scan1/Tracking/PROB/Results_StickZeppelinBall_COMMIT/fit_NRMSE.nii.gz' );
niiERR_C_img = 100.0 * niiERR_C.get_data()

# load mask
niiMASK = nibabel.load( 'STN96/scan1/Tracking/PROB/dictionary_tdi.nii.gz' );
niiMASK_img = niiMASK.get_data()
```

Then we plot the fitting error with *LiFE* in a representative slice of the brain where two important fiber bundles cross (CST and CC):

```python
# plot the NRMSE with LiFE
pylab.figure(1,facecolor='white')
h = pylab.imshow( numpy.rot90(niiERR_L_img[:,69,:].squeeze()), interpolation='nearest', cmap='hot' )
h.set_clim(0.0,100.0)
pylab.colorbar()
pylab.axis('off')
h.axes.get_xaxis().set_visible(False)
h.axes.get_yaxis().set_visible(False)
yL = niiERR_L_img[ niiMASK_img>0 ]
pylab.title( 'LiFE : %.1f%% +/- %.1f%%' % ( yL.mean(), yL.std() ) )
```

![NRMSE for LiFE](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig1.png)

The average fitting error is, in this case, pretty high, i.e. **69.7% ± 16.5%**. Also, we see that *LiFE* shows the highest errors in regions with crossing fibers and close to gray matter, as expected (see [this abstract](ISMRM_3148.pdf)).

We plot now the fitting error with *COMMIT*:

```python
# plot the NRMSE with COMMIT
pylab.figure(2,facecolor='white')
h = pylab.imshow( numpy.rot90(niiERR_C_img[:,69,:].squeeze()), interpolation='nearest', cmap='hot' )
h.set_clim(0.0,100.0)
pylab.colorbar()
pylab.axis('off')
h.axes.get_xaxis().set_visible(False)
h.axes.get_yaxis().set_visible(False)
yC = niiERR_C_img[ niiMASK_img>0 ]
pylab.title( 'COMMIT : %.1f%% +/- %.1f%%' % ( yC.mean(), yC.std() ) )
```

![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig2.png)

The average fitting error is drastically reduced with *COMMIT*, i.e. (**19.3% ± 4.7%**). Also, a more homogeneous distribution of the errors can be observed, notably in crossing regions and in proximity to gray matter.
 
Now we can directly compare the *fitting error distributions* of the two models:

```python
# direct comparison of the NRMSE of LiFE and COMMIT
pylab.figure(3,facecolor='white')
pylab.clf()
pylab.hold(True)
N = numpy.count_nonzero(niiMASK_img>0)
hL, _ = numpy.histogram( yL, bins=60, range=(0,100), density=True )
hC, _ = numpy.histogram( yC, bins=60, range=(0,100), density=True )
pylab.plot( hL, '-', color=[.8,0,0], linewidth=3, label='LiFE' )
pylab.plot( hC, '-', color=[0,.8,0], linewidth=3, label='COMMIT' )
pylab.xticks( numpy.linspace(0,60,11,dtype=numpy.uint8), numpy.linspace(0,100,11,dtype=numpy.uint8) )
pylab.grid(True)
pylab.xlabel( 'NRMSE [%]' )
pylab.ylabel( 'percentage of voxels' )
pylab.legend()
pylab.title( 'Error distributions' )
```

![Histograms comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig3.png)

Also, we can directly compare their fitting errors *voxel-by-voxel* with the following scatter-plot:

```python
# voxelwise comparison of the NRMSE of LiFE and COMMIT
pylab.figure(4,facecolor='white')
pylab.clf()
pylab.hold(True)
pylab.plot( yL, yC, 'bx' )
pylab.plot( [0,100], [0,100], 'k--', linewidth=2 )
pylab.grid(True)
pylab.axis([0,100,0,100])
pylab.xlabel( 'NRMSE [%] with LiFE' )
pylab.ylabel( 'NRMSE [%] with COMMIT' )
pylab.title( 'Error scatterplot' )
```

![Scatterplot comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig4.png)

As we can see, in all voxels the *COMMIT* model **always explains the data much better** than the *LiFE* model.


## Compare the two models (continued)

One might also want to **evaluate how well both models explain the measured diffusion MRI signal** acquired with the scanner.
To this end, we need to *add back the mean* to the data used by the *LiFE* model and utilize the previously estimated fiber weights. By doing this we can directly compare the two models with respect to the same common quantity, i.e. the acquired diffusion MRI signal.
No normalization is needed in this case and we can then use the *RMSE* (expressed in raw signal units) to compare **the accuracy of the fit** of the two approaches.

To this aim, it is simply necessary to perform the following operations after processing the data with *LiFE*:

```python
# reload the DWI data and KERNELS (LUT) and DO NOT remove the mean
mit.CONFIG['doDemean'] = False
mit.load_data( 'DWI.nii', 'DWI.scheme' )
mit.load_kernels()
mit.build_operator()

# recompute the error metrics
mit.save_results( 'LIFE_2' )
```

By doing this, both the measurements **y** and the signal **Ax** predicted by the *LiFE* model will be compared using the *NMSE* error metric to evaluate how well the *LiFE* model actually explains the measured diffusion MRI signal.
We then load the *RMSE* errors and compare the accuracy of the two models, as follows:

```python
# load error maps
niiERR_L = nibabel.load( 'STN96/scan1/Tracking/PROB/Results_StickZeppelinBall_LIFE_2/fit_RMSE.nii.gz' );
niiERR_L_img = niiERR_L.get_data()
niiERR_C = nibabel.load( 'STN96/scan1/Tracking/PROB/Results_StickZeppelinBall_COMMIT/fit_RMSE.nii.gz' );
niiERR_C_img = niiERR_C.get_data()

# plot the RMSE with LiFE
pylab.figure(5,facecolor='white')
h = pylab.imshow( numpy.rot90(niiERR_L_img[:,69,:].squeeze()), interpolation='nearest', cmap='hot' )
h.set_clim(0.0,200.0)
pylab.colorbar()
pylab.axis('off')
h.axes.get_xaxis().set_visible(False)
h.axes.get_yaxis().set_visible(False)
yL = niiERR_L_img[ niiMASK_img>0 ]
pylab.title( 'LiFE : %.1f +/- %.1f' % ( yL.mean(), yL.std() ) )

# plot the RMSE with COMMIT
pylab.figure(6,facecolor='white')
h = pylab.imshow( numpy.rot90(niiERR_C_img[:,69,:].squeeze()), interpolation='nearest', cmap='hot' )
h.set_clim(0.0,200.0)
pylab.colorbar()
pylab.axis('off')
h.axes.get_xaxis().set_visible(False)
h.axes.get_yaxis().set_visible(False)
yC = niiERR_C_img[ niiMASK_img>0 ]
pylab.title( 'COMMIT : %.1f +/- %.1f' % ( yC.mean(), yC.std() ) )

# direct comparison of the RMSE of LiFE and COMMIT
pylab.figure(7,facecolor='white')
pylab.clf()
pylab.hold(True)
N = numpy.count_nonzero(niiMASK_img>0)
hL, _ = numpy.histogram( yL, bins=100, range=(0,300), density=True )
hC, _ = numpy.histogram( yC, bins=100, range=(0,300), density=True )
pylab.plot( hL, '-', color=[.8,0,0], linewidth=3, label='LiFE' )
pylab.plot( hC, '-', color=[0,.8,0], linewidth=3, label='COMMIT' )
pylab.xticks( numpy.linspace(0,100,7,dtype=numpy.uint16), numpy.linspace(0,300,7,dtype=numpy.uint16) )
pylab.grid(True)
pylab.xlabel( 'RMSE [raw signal units]' )
pylab.ylabel( 'percentage of voxels' )
pylab.legend()
pylab.title( 'Error distributions' )

# voxelwise comparison of the NRMSE of LiFE and COMMIT
pylab.figure(8,facecolor='white')
pylab.clf()
pylab.hold(True)
pylab.plot( yL, yC, 'bx' )
pylab.plot( [0,350], [0,350], 'k--', linewidth=2 )
pylab.grid(True)
pylab.axis([0,350,0,350])
pylab.xlabel( 'RMSE [raw signal units] with LiFE' )
pylab.ylabel( 'RMSE [raw signal units] with COMMIT' )
pylab.title( 'Error scatterplot' )
```

![RMSE for LiFE](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig5.png)

![RMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig6.png)

![Histogram comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig7.png)

![Scatterplot comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/tutorials/LiFE_STN96/RESULTS_Fig8.png)

As we can see, the results essentially lead to the the same results, as previously highlighted using the *NRMSE* metric, de facto showing the **superiority of the *COMMIT* model in explaining the measured diffusion MRI signal** with respect to *LiFE*.
