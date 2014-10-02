clearvars, clearvars -global, clc
COMMIT_Setup


%% SETUP & DATA
CONFIG = COMMIT_Config( 'STN96', 'scan1' );
CONFIG.doDemean	= true;
DATA_Load


%% KERNELS
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
CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','PROB');
DICTIONARY_LoadSegments

CONFIG.OPTIMIZATION.nTHREADS = 4;
OPTIMIZATION_Setup


%% SOLVE inverse problem
OPTIMIZATION_Solve
OPTIMIZATION_SaveResults

