# Comparison between COMMIT and LIFE (on STN96 data from the LiFE original publication)

The *LiFE* method recently described in [(Pestilli et al, Nat Methods, 2014)](http://www.nature.com/nmeth/journal/v11/n10/abs/nmeth.3098.html) can be considerd a **special case** of our *COMMIT* framework
 [(Daducci et al, IEEE Trans Med Imaging, 2014)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6884830); notably, it corresponds to the preliminary version of it that we had proposed in [(Daducci et al, ISBI, 2013)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6556527) back in 2013.
 
To model the signal in each voxel, *LiFE* considers only the restricted contributions arising from all the tracts crossing a paarticular voxel (i.e. restricted diffusion). However, *LiFE* does not consider the extra-cellular space around the fibers (i.e. hindered diffusion) and all the partial volume that can occur with gray matter and CSF. On the other hand, *COMMIT* can account for all possible compartments that are present in a voxel and that contribute to the diffusion MR signal.

In this demo we will show the **importance of using adequate compartment models** in order to be able to effectively evaluate the evidence of a tractogram, i.e. set of fiber tracts.
 

## Download data

1. Download precomputed data (CSD reconstruction + probabilistic tracking) from [here](http://hardi.epfl.ch/static/data/COMMIT_demos/STN96_scan1.zip).
2. Unzip the file content into the `STN96/scan1` folder in your data directory.
3. Download the original DWI data from [here](https://stacks.stanford.edu/file/druid:cs392kv3054/life_demo_data.tar.gz).
4. Extract the file `life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz` from the archive, unzip it and move it to the `scan1` folder with the name `DWI.nii`, i.e.

```bash
gunzip life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz
mv life_demo_scan1_subject1_b2000_150dirs_stanford.nii STN96/scan1/DWI.nii
```


## Process data with COMMIT

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
CONFIG.kernels.doNormalize = false;
KERNELS_ProcessAtoms
```

Compile the linear operators **A** and **At**:

```matlab
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup
```

Actually **solve** the inverse problem according to the  *COMMIT* model:

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
CONFIG.kernels.doNormalize = false;
KERNELS_ProcessAtoms
```

Compile the linear operators **A** and **At**:

```matlab
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup
```

Actually **solve** the inverse problem according to the  *LiFE* model:

```matlab
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults
```
Result will be stored in `STN96/scan1/Tracking/PROB/Results_LIFE/`.


## Compare the two models

Load the computed *NRMSE* metric of both models:

```matlab
clear, clc
niiNMSE_L = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_NRMSE.nii') );
niiNMSE_C = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_NRMSE.nii') );
niiMASK   = load_untouch_nii( fullfile('scan1','Tracking','PROB','dictionary_mask.nii') );
```

Plot the fitting error with *LiFE*:

```matlab
figure(1), imagesc( rot90(squeeze(100 * niiNMSE_L.img(:,70,:))), [0 100] )
axis ij image off, colorbar
yL = 100 * niiNMSE_L.img( niiMASK.img>0 );
title( sprintf('LiFE : %.1f%% +/- %.1f%%', mean(yL), std(yL) ))
```

![NRMSE for LiFE](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig1.png)

The average fitting error is, in this case, pretty high, i.e. **68.6% ± 18.0%**. Also, we see that *LiFE* shows the highest errors in regions with crossing fibers and close to gray matter, as expected.

We plot now the fitting error with *COMMIT*:

```matlab
figure(2), imagesc( rot90(squeeze(100 * niiNMSE_C.img(:,70,:))), [0 100] )
axis ij image off, colorbar
yC = 100 * niiNMSE_C.img( niiMASK.img>0 );
title( sprintf('COMMIT : %.1f%% +/- %.1f%%', mean(yC), std(yC) ))
```

![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig2.png)

The average fitting error is drastically reduced with *COMMIT*, i.e. (**19.3% ± 4.6%**). Also, a more homogeneous distribution of the errors can be observed, notably in crossing regions and in proximity with gray matter.
 
Compare the *fitting error histograms* of the two models:

```matlab
figure(3), set(gcf,'Color',[1 1 1]);
cla, set(gca,'Color',[.97 .97 .97]); hold on
x = linspace(0,100,100);
yL = hist( 100 * niiNMSE_L.img(niiMASK.img>0), x );
yC = hist( 100 * niiNMSE_C.img(niiMASK.img>0), x );
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on
legend( 'LiFE', 'COMMIT' )
xlabel( 'NRMSE [%]' ), ylabel('frequency')
title('LiFE vs COMMIT')
```

![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig3.png)

Compare *voxel-by-voxel* the fitting error:

```matlab
figure(4), set(gcf,'Color',[1 1 1]);
cla, set(gca,'Color',[.97 .97 .97]); hold on
yL = 100 * niiNMSE_L.img( niiMASK.img>0 );
yC = 100 * niiNMSE_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 100], [0 100], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 100 0 100])
xlabel( 'NMSE [%] with LiFE' ), ylabel( 'NMSE [%] with COMMIT' )
title('LiFE vs COMMIT')
```

![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig4.png)
