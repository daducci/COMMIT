# Comparison between COMMIT and LIFE (on STN96 data from the LiFE original publication)

The *LiFE* method recently described in [(Pestilli et al, Nat Methods, Sep 2014)](http://www.nature.com/nmeth/journal/v11/n10/abs/nmeth.3098.html) can be considerd a **special case** of our *COMMIT* framework
 [(Daducci et al, IEEE Trans Med Imaging, Aug 2014)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6884830); notably, it corresponds to the preliminary version of it that we had proposed in [(Daducci et al, ISBI, Apr 2013)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6556527) back in 2013.
 
To model the signal in each voxel, *LiFE* considers only the restricted contributions arising from all the tracts crossing a particular voxel (i.e. restricted diffusion). However, *LiFE* does not consider the extra-cellular space around the fibers (i.e. hindered diffusion) and all the partial volume that can occur with gray matter and CSF. On the other hand, *COMMIT* can account for all possible compartments that are present in a voxel and that contribute to the diffusion MR signal.

In this demo we will show the **importance of using adequate multi-compartment models** to be able to effectively evaluate the evidence of a tractogram, i.e. set of fiber tracts.

For more information, please refer to [this abstract](ISMRM_3148.pdf), that recently has been **specially selected for a power pitch presentation** (less than 3% of submitted papers) at the annual *International Society for Magnetic Resonance in Medicine* (ISMRM) meeting in Toronto (30/05-05/06 2015)!
 

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

Use the script `trk2dictionary` in the `c++` folder to convert the tracts contained in the file `STN96/scan1/Tracking/PROB/fibers.trk` to the internal sparse data structure used by *COMMIT*:

```bash
trk2dictionary \
   -i STN96/scan1/Tracking/PROB/fibers.trk \
   -p STN96/scan1/CSD/CSD_FODsh_peaks.nii \
   -w STN96/scan1/WM.nii \
   -o STN96/scan1/Tracking/PROB/ \
   -c
```
This will create the necessary data structure (`STN96/scan1/Tracking/PROB/dictionary_*`) containing all the details if the tracts.

Please be sure to first compile the code in the `c++` folder and to have the resulting binaries (in particular `trk2dictionary`) in your path.
NB: the output tractogram from *MrTrix* has been already converted to the format accepted by *COMMIT*, i.e. [TrackiVis format](http://www.trackvis.org/docs/?subsect=fileformat) with fiber coordinates (in mm) in image space, where the coordinate (0,0,0) corresponds to the corner of first voxel.

## Process data with COMMIT

Start MATLAB.

Setup and load the data:

```matlab
clear, clc
COMMIT_Setup

CONFIG = COMMIT_Config( 'STN96', 'scan1' );
CONFIG.doDemean	= false;
DATA_Load
```

Calculate the **kernels** corresponding to the different compartments. In this example, we use 1 kernel for intra-axonal compartment (i.e. Stick), 1 for extra-axonal space (i.e. Zeppelin) and 2 to model partial volume with gray matter and CSF:

```matlab
CONFIG.kernels.namePostfix  = 'COMMIT';
CONFIG.kernels.d            = 1.7;
CONFIG.kernels.Rs           = [ 0 ];
CONFIG.kernels.ICVFs        = [ 0.7 ];
CONFIG.kernels.dISOs        = [ 3.0 1.7 ];

KERNELS_CreateFolderForHighResolutionKernels( CONFIG );
KERNELS_PrecomputeRotationMatrices();
KERNELS_StickZeppelinBall_Generate( CONFIG );
KERNELS_ActiveAx_RotateAndSave( CONFIG );

KERNELS = KERNELS_Load( CONFIG );
KERNELS_ProcessAtoms
```

Compile the linear operators **A** and **At**:

```matlab
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup
```

**Solve** the inverse problem according to the  *COMMIT* model:

```matlab
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults
```

Result will be stored in `STN96/scan1/Tracking/PROB/Results_COMMIT/`.


## Process data with LiFE

Setup and load the data; this time, however, we will apply the *demeaning procedure* used in *LiFE* to both data and kernels:

```matlab
clear, clc
COMMIT_Setup

CONFIG = COMMIT_Config( 'STN96', 'scan1' );
CONFIG.doDemean	= true;
DATA_Load
```

Calculate the **kernel** corresponding to the intra-cellular compartment (the only one considered in *LiFE*); in this example, thus, we use only 1 kernel for intra-axonal compartment (i.e. Stick):

```matlab
CONFIG.kernels.namePostfix  = 'LIFE';
CONFIG.kernels.d            = 1.7;
CONFIG.kernels.Rs           = [ 0 ];
CONFIG.kernels.ICVFs        = [  ];
CONFIG.kernels.dISOs        = [  ];

KERNELS_CreateFolderForHighResolutionKernels( CONFIG );
KERNELS_PrecomputeRotationMatrices();
KERNELS_StickZeppelinBall_Generate( CONFIG );
KERNELS_ActiveAx_RotateAndSave( CONFIG );

KERNELS = KERNELS_Load( CONFIG );
KERNELS_ProcessAtoms
```

Compile the linear operators **A** and **At**:

```matlab
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup
```

**Solve** the inverse problem according to the *LiFE* model:

```matlab
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults
```

Result will be stored in `STN96/scan1/Tracking/PROB/Results_LIFE/`.


## Compare the two models

Let's first analyze the performances of the two approaches in the **native space in which the two models perform the fit**. In fact, *LiFE* does not fit the model to the acquired diffusion MRI signal, but rather on the signal after removing the mean value in each voxel, i.e. demeaned signal.

It is important to note that as the two models actually work in different spaces (different values as measurements), a normalization of the error metrics is required in order to compare their accuracy in explaining the input data. To this aim, we will use the *Normalized RMSE (NRMSE)* as quality measure. Please note that the normalization constant used in each voxel quantifies the magnitude of the data in that voxel, so the values are expressed as *percentage error* with respect to the actual measurements used in the voxel, i.e. measured diffusion MRI signal for *COMMIT* and demeaned signal for *LiFE*.

We then load the *NRMSE* fit error of the two models, as follows:

```matlab
niiERR_L = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_NRMSE.nii') );
niiERR_C = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_NRMSE.nii') );
```

Then we plot the fitting error with *LiFE* in a representative slice of the brain where two important fiber bundles cross (CST and CC):

```matlab
% plot the NRMSE of the LiFE model
figure(1), clf, imagesc( rot90(squeeze(100*niiERR_L.img(:,70,:))), [0 100] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = 100*niiERR_L.img( niiMASK.img>0 );
title( sprintf('LiFE : %.1f%% +/- %.1f%%', mean(yL), std(yL) ))
```

![NRMSE for LiFE](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig1.png)

The average fitting error is, in this case, pretty high, i.e. **68.9% ± 17.9%**. Also, we see that *LiFE* shows the highest errors in regions with crossing fibers and close to gray matter, as expected (see [this abstract](ISMRM_3148.pdf)).

We plot now the fitting error with *COMMIT*:

```matlab
% plot the NRMSE of the COMMIT model
figure(2), clf, imagesc( rot90(squeeze(100*niiERR_C.img(:,70,:))), [0 100] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = 100*niiERR_C.img( niiMASK.img>0 );
title( sprintf('COMMIT : %.1f%% +/- %.1f%%', mean(yL), std(yL) ))
```

![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig2.png)

The average fitting error is drastically reduced with *COMMIT*, i.e. (**18.3% ± 4.7%**). Also, a more homogeneous distribution of the errors can be observed, notably in crossing regions and in proximity with gray matter.
 
We now directly compare the *fitting error distributions* of the two models:

```matlab
% direct comparison of the NRMSE of LiFE and COMMIT
figure(3), clf, hold on
x = linspace(0,100,60);
yL = hist( 100*niiERR_L.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
yC = hist( 100*niiERR_C.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on, axis tight
xlabel( 'NRMSE [%]' ), ylabel( 'percentage of voxels' )
pbaspect([1 1 1])
legend('LiFE','COMMIT')
title('Error distributions')
```

![Histograms comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig3.png)

Also, we can directly compare *voxel-by-voxel* the fitting errors of the two models:

```matlab
% voxelwise comparison of the NRMSE of LiFE and COMMIT
figure(4), clf, hold on
yL = 100*niiERR_L.img( niiMASK.img>0 );
yC = 100*niiERR_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 100], [0 100], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 100 0 100])
xlabel( 'NRMSE [%] with LiFE' ), ylabel( 'NRMSE [%] with COMMIT' )
title('Error scatterplot')
```

![Scatterplot comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig4.png)

As we can see, in all voxels the *COMMIT* model **always explains the data much better** than the *LiFE* model.


## Compare the two models (2/2)

However, one might also want to **evaluate how well both models explain the measured diffusion MRI signal** acquired with the scanner.
Thus, we need to *add back the mean* to the data used by the *LiFE* model and utilize the previously estimated fiber weights. this way, we can compare the two models with respect to the same common quantity, i.e. the acquired diffusion MRI signal.
In this case, no normalization is needed and to compare **the accuracy of the fit** we can use the *RMSE* (expressed in raw signal units). 

To perform this comparison, it is simply necessary to perform the following operations:

```matlab
% reload the DWI data and DO NOT remove the mean
CONFIG.doDemean	= false;
DATA_Load

% reload kernel LUT and, again, DO NOT remove the mean
KERNELS_ActiveAx_RotateAndSave( CONFIG );
KERNELS = KERNELS_Load( CONFIG );
KERNELS_ProcessAtoms

% recompute the error metrics
OPTIMIZATION_SaveResults
```

At this point, we can load the *RMSE* fit errors and compare the results, as follows:

```matlab
niiERR_L  = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_RMSE.nii') );
niiERR_C  = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_RMSE.nii') );

% plot the RMSE of the LiFE model
figure(5), clf, imagesc( rot90(squeeze(niiERR_L.img(:,70,:))), [0 200] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = niiERR_L.img( niiMASK.img>0 );
title( sprintf('LiFE : %.1f +/- %.1f', mean(yL), std(yL) ))
saveas(gcf,'RESULTS_Fig5.png')

% plot the RMSE of the COMMIT model
figure(6), clf, imagesc( rot90(squeeze(niiERR_C.img(:,70,:))), [0 200] )
axis ij image off, cm = hot(256); cm(1,:) = 0; colormap(cm); colorbar
yL = niiERR_C.img( niiMASK.img>0 );
title( sprintf('COMMIT : %.1f +/- %.1f', mean(yL), std(yL) ))
saveas(gcf,'RESULTS_Fig6.png')

% direct comparison of the RMSE of LiFE and COMMIT
figure(7), clf, hold on
x = linspace(0,300,100);
yL = hist( niiERR_L.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
yC = hist( niiERR_C.img(niiMASK.img>0), x ) / nnz(niiMASK.img>0);
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on, axis tight
xlabel( 'RMSE [raw signal units]' ), ylabel( 'percentage of voxels' )
pbaspect([1 1 1])
legend('LiFE','COMMIT')
title('Error distributions')
saveas(gcf,'RESULTS_Fig7.png')

% voxelwise comparison of the RMSE of LiFE and COMMIT
figure(8), clf, hold on
yL = niiERR_L.img( niiMASK.img>0 );
yC = niiERR_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 260], [0 260], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 260 0 260])
xlabel( 'RMSE [raw signal units] with LiFE' ), ylabel( 'RMSE [raw signal units] with COMMIT' )
title('Error scatterplot')
saveas(gcf,'RESULTS_Fig8.png')
```

![RMSE for LiFE](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig5.png)

![RMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig6.png)

![Histogram comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig7.png)

![Scatterplot comparison LiFE vs COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig8.png)

As we can see, the results essentially show the same results as highlightes using the *NRMSE* metric.
