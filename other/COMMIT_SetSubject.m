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
function COMMIT_SetSubject( protocol, subject )

    global COMMIT_data_path
    global CONFIG

    CONFIG = [];

    CONFIG.protocol			= protocol;
    CONFIG.subject			= subject;

    CONFIG.DATA_path		= fullfile( COMMIT_data_path, CONFIG.protocol, CONFIG.subject );
    CONFIG.dwiFilename		= fullfile( CONFIG.DATA_path, 'DWI.nii' );
	CONFIG.dim              = [];
	CONFIG.pixdim           = [];
    CONFIG.normalizeSignal	= true;
    CONFIG.doDemean			= false;
    CONFIG.useReference		= false;
    CONFIG.normalizeKernels = true;

    CONFIG.schemeFilename	= fullfile( CONFIG.DATA_path, 'DWI.scheme' );
    CONFIG.scheme			= [];
    CONFIG.b0_thr           = 1;
