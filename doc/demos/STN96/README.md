# Comparison between COMMIT and LIFE on STN96 data from the LiFE original publication

## Download data

1. Download precomputed data (CSD reconstruction + probabilistic tracking) from:
   (http://hardi.epfl.ch/static/data/COMMIT_demos/STN96_scan1.zip)
2. Unzip the file into the `STN96/scan1` folder in your data directory
3. Download the original DWI data from:
   (https://stacks.stanford.edu/file/druid:cs392kv3054/life_demo_data.tar.gz)
4. Extract the file `life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz` from the archive, unzip it and move it into the `scan1` folder with the name `DWI.nii`, i.e.
```bash 
gunzip life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz
mv life_demo_scan1_subject1_b2000_150dirs_stanford.nii data/STN96/scan1/DWI.nii
```

## Process data with LiFE

```matlab
clearvars, clearvars -global, clc
COMMIT_Setup

%% SETUP & DATA
%  ============
CONFIG = COMMIT_Config( 'STN96', 'scan1' );
CONFIG.doDemean	= true;
DATA_Load

%% KERNELS
%  =======
CONFIG.kernels.namePostfix  = 'LIFE';
CONFIG.kernels.d            = 1.7;
CONFIG.kernels.Rs           = [ 0 ];
CONFIG.kernels.ICVFs        = [  ];
CONFIG.kernels.dISOs        = [  ];

% Calculate the kernels
KERNELS_CreateFolderForHighResolutionKernels( CONFIG );
KERNELS_PrecomputeRotationMatrices();
KERNELS_StickZeppelinBall_Generate( CONFIG );
KERNELS_ActiveAx_RotateAndSave( CONFIG );

% Processing
KERNELS = KERNELS_Load( CONFIG );
CONFIG.kernels.doNormalize = false;
KERNELS_ProcessAtoms

%% LINEAR OPERATORS A and A_t
%  ==========================
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup

%% SOLVE inverse problem
%  =====================
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults
```

## Process data with COMMIT
```matlab
clearvars, clearvars -global, clc
COMMIT_Setup

%% SETUP & DATA
%  ============
CONFIG = COMMIT_Config( 'STN96', 'scan1' );
CONFIG.doDemean	= false;
DATA_Load

%% KERNELS
%  =======
CONFIG.kernels.namePostfix  = 'COMMIT';
CONFIG.kernels.d            = 1.7;
CONFIG.kernels.Rs           = [ 0 ];
CONFIG.kernels.ICVFs        = [ 0.7 ];
CONFIG.kernels.dISOs        = [ 3.0 1.7 ];

% Calculate the kernels
KERNELS_CreateFolderForHighResolutionKernels( CONFIG );
KERNELS_PrecomputeRotationMatrices();
KERNELS_StickZeppelinBall_Generate( CONFIG );
KERNELS_ActiveAx_RotateAndSave( CONFIG );

% Processing
KERNELS = KERNELS_Load( CONFIG );
CONFIG.kernels.doNormalize = false;
KERNELS_ProcessAtoms

%% LINEAR OPERATORS A and A_t
%  ==========================
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup

%% SOLVE inverse problem
%  =====================
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults
```


## Compare the two models

Load the computed NRMSE metric:
```matlab
clearvars, clearvars -global, clc

niiNMSE_L = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_LIFE','fit_NRMSE.nii') );
niiNMSE_C = load_untouch_nii( fullfile('scan1','Tracking','PROB','Results_COMMIT','fit_NRMSE.nii') );
niiMASK   = load_untouch_nii( fullfile('scan1','Tracking','PROB','dictionary_mask.nii') );
```

Plot the **NRMSE** for *LiFE*:
```matlab
figure(1), imagesc( rot90(squeeze(niiNMSE_L.img(:,70,:))), [0 1] )
axis ij image off, colorbar
yL = niiNMSE_L.img( niiMASK.img>0 );
title( sprintf('LiFE : %.2f +/- %.2f', mean(yL), std(yL) ))
```
![NRMSE for LiFE](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig1.png)


Plot the **NRMSE** for *COMMIT*:
```matlab
figure(2), imagesc( rot90(squeeze(niiNMSE_C.img(:,70,:))), [0 1] )
axis ij image off, colorbar
yC = niiNMSE_C.img( niiMASK.img>0 );
title( sprintf('COMMIT : %.2f +/- %.2f', mean(yC), std(yC) ))
```
![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig2.png)

Compare the **NRMSE** histograms of the two models:
```matlab
figure(3), set(gcf,'Color',[1 1 1]);
cla, set(gca,'Color',[.97 .97 .97]); hold on
x = linspace(0,1,100);
yL = hist( niiNMSE_L.img(niiMASK.img>0), x );
yC = hist( niiNMSE_C.img(niiMASK.img>0), x );
plot( x, yL, '- ', 'LineWidth', 3, 'Color',[.8 0 0] )
plot( x, yC, '- ', 'LineWidth', 3, 'Color',[0 .8 0] )
grid on, box on
legend( 'LiFE', 'COMMIT' )
xlabel( 'NRMSE' ), ylabel('frequency')
title('LiFE vs COMMIT')
```
![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig3.png)

Compare the **NRMSE** voxel-by-voxel between the two models:
```matlab
figure(4), set(gcf,'Color',[1 1 1]);
cla, set(gca,'Color',[.97 .97 .97]); hold on
yL = niiNMSE_L.img( niiMASK.img>0 );
yC = niiNMSE_C.img( niiMASK.img>0 );
plot( yL, yC, 'bx' )
plot( [0 1], [0 1], 'k--', 'LineWidth', 2 )
grid on, box on
axis([0 1 0 1])
xlabel( 'LiFE' ), ylabel( 'COMMIT' )
```
![NRMSE for COMMIT](https://github.com/daducci/COMMIT/blob/master/doc/demos/STN96/RESULTS_Fig4.png)
