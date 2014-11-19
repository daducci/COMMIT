%
% Set the default parameters for a specific experiment (protocol/subject)
%
% Parameters
% ----------
% protocol : string
%   Folder for a specific acquisition protocol, possible containing different subjects
% subject : string
%   Folder of the specific subject
%
% Returns
% -------
% CONFIG : struct
% 	Struct containing all the default parameters and setup of the experiment
%
function [ CONFIG ] = COMMIT_Config( protocol, subject )

global COMMIT_path DATA_path

CONFIG = [];

CONFIG.protocol			= protocol;
CONFIG.subject			= subject;

CONFIG.DATA_path		= fullfile(DATA_path,CONFIG.protocol,CONFIG.subject);
CONFIG.dwiFilename		= fullfile(CONFIG.DATA_path,'DWI.nii');
CONFIG.normalizeSignal	= true;
CONFIG.doDemean			= false;
CONFIG.useReference		= false;

CONFIG.schemeFilename	    = fullfile(CONFIG.DATA_path,'DWI.scheme');
CONFIG.scheme			    = [];

CONFIG.TRACKING_path		= fullfile(CONFIG.DATA_path,'Tracking','GIBBS');

CONFIG.kernels = [];
CONFIG.kernels.namePostfix	= 'COMMIT';
CONFIG.kernels.d			= 1.7;				% float: parallel diffusivity to use (units of 1e-10, eg specify 17.0 for 17.0 mm^2/s)
CONFIG.kernels.Rs			= [ 0.5 10 ];		% float array: radii of the cylinders for the intra-cellular atoms (units of 1e-3, eg specify 2.0 for 2.0 micrometers)
CONFIG.kernels.ICVFs		= [ 0.7 ];			% float array: the intra-cellular volume fractions (used in the tortuosity model) to simulate the extra-cellular atoms
CONFIG.kernels.dISOs		= [ 3.0 1.7 ];		% float array: free diffusivities to simulate the isotropic atoms (units of 1e-10, eg specify 3.0 for 3.0 mm^2/s)
CONFIG.kernels.doNormalize  = true;

CONFIG.OPTIMIZATION = [];
CONFIG.OPTIMIZATION.nTHREADS = 1;
CONFIG.OPTIMIZATION.verbosity       = 3;
CONFIG.OPTIMIZATION.optTol			= 1e-5;
CONFIG.OPTIMIZATION.maxIter			= 100;
CONFIG.OPTIMIZATION.solver			= 'NNLS';

CONFIG.OPTIMIZATION.SNR_estimated	= NaN;
CONFIG.OPTIMIZATION.RW_maxIter		= 0;
CONFIG.OPTIMIZATION.RW_tauN			= 1;
CONFIG.OPTIMIZATION.RW_tauD			= 1;

