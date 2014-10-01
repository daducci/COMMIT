%
% Precompute rotated versions (180x180 directions) of each atom simulated with the Camino toolkit and
% store them matching the acquisition protocol used for the subject (specific scheme file)
% 
% Parameters
% ----------
% CONFIG : struct
% 	Struct containing all the parameters and setup of the experiment
% lmax : unsigned int
%   Maximum spherical harmonics order to use for the rotation phase
%
function KERNELS_ActiveAx_RotateAndSave( CONFIG, lmax )
	if nargin < 2, lmax = 12; end
	global COMMIT_path

	TIME = tic();	
	fprintf( '\n-> Rotating kernels to 180x180 directions for subject="%s" and d=%4.2f:\n', CONFIG.subject, CONFIG.kernels.d );

	% check if original scheme exists
	if ~exist( CONFIG.schemeFilename, 'file' )
		error( '[KERNELS_ActiveAx_RotateAndSave] File "%s" not found', CONFIG.schemeFilename )
	else
		scheme = KERNELS_LoadScheme( CONFIG.schemeFilename );
	end
	
	% check if auxiliary matrices have been precomputed
	auxFilename = fullfile(COMMIT_path,'code','kernels',sprintf('AUX_matrices__lmax=%d.mat',lmax) );
	if ~exist( auxFilename, 'file' )
		error( '[KERNELS_ActiveAx_RotateAndSave] Auxiliary matrices "%s" not found', auxFilename )
	else
		AUX = load( auxFilename );
	end

	[tmp,~,~] = fileparts( CONFIG.schemeFilename );
	OUTPUT_path = fullfile(tmp,'..','kernels',sprintf('d=%4.2f',CONFIG.kernels.d));
	
	% precompute aux data structures
	% ------------------------------
	idx_IN  = [];
	idx_OUT = [];
	Ylm_OUT = [];
	row = 1;
	for i = 1:numel(scheme.shells)
		idx_IN{end+1}  = row : row+500-1;
		idx_OUT{end+1} = scheme.shells{i}.idx;
		[colatitude, longitude] = cart2sphere( scheme.shells{i}.grad(:,1), scheme.shells{i}.grad(:,2), scheme.shells{i}.grad(:,3) );
		Ylm_OUT{end+1} = createYlm( AUX.lmax, colatitude, longitude ); % matrix from SH to real space
		row = row+500;
	end

	% rotate kernel in all direction and sample according to subject's acquisition scheme
	% -----------------------------------------------------------------------------------
	KERNELS = {};
	KERNELS.d_par     = CONFIG.kernels.d;
	KERNELS.nS        = scheme.nS;
	KERNELS.wmr       = [];
	KERNELS.wmr_radii = [];
	KERNELS.wmr_norm  = [];
	KERNELS.wmh       = [];
	KERNELS.wmh_icvf  = [];
	KERNELS.wmh_norm  = [];
	KERNELS.iso       = [];
	KERNELS.iso_d     = [];
	KERNELS.iso_norm  = [];

	% intra-cellular
	fprintf( '\t- %d intra-cellular\n', numel(CONFIG.kernels.Rs) );
	for R = CONFIG.kernels.Rs
		fprintf( '\t\t* R = %5.2f micrometers...', R );
		TIME2 = tic();
		KERNELS.wmr{end+1} = rotate_kernel( fullfile(OUTPUT_path,sprintf('Er_R=%05.2f.Bfloat',R)), scheme, AUX, idx_IN, idx_OUT, Ylm_OUT );
		KERNELS.wmr_radii(end+1) = R;
		KERNELS.wmr_norm(end+1) = norm( KERNELS.wmr{end}(:,1,1) );
		fprintf( ' [%.1f seconds]\n', toc(TIME2) );
	end

	% extra-cellular
	fprintf( '\t- %d extra-cellular\n', numel(CONFIG.kernels.ICVFs) );
	for ICVF = CONFIG.kernels.ICVFs
		fprintf( '\t\t* ICVF = %.2f...', ICVF );
		TIME2 = tic();
		KERNELS.wmh{end+1} = rotate_kernel( fullfile(OUTPUT_path,sprintf('Eh_icvf=%05.2f.Bfloat',ICVF)), scheme, AUX, idx_IN, idx_OUT, Ylm_OUT );
		KERNELS.wmh_icvf(end+1) = ICVF;
		KERNELS.wmh_norm(end+1) = norm( KERNELS.wmh{end}(:,1,1) );
		fprintf( ' [%.1f seconds]\n', toc(TIME2) );
	end
	
	% isotropic
	fprintf( '\t- %d isotropic\n', numel(CONFIG.kernels.dISOs) );
	for dISO = CONFIG.kernels.dISOs
		fprintf( '\t\t* dISO = %5.2f ...', dISO );
		TIME2 = tic();
		KERNELS.iso{end+1} = resample_iso_kernel( fullfile(OUTPUT_path,sprintf('Ei_dIso=%05.2f.Bfloat',dISO)), scheme, AUX, idx_IN, idx_OUT, Ylm_OUT );
		KERNELS.iso_d(end+1) = dISO;
		KERNELS.iso_norm(end+1) = norm( KERNELS.iso{end}(:,1,1) );
		fprintf( ' [%.1f seconds]\n', toc(TIME2) );
	end

	% save
	fprintf( '\t- saving to file...' );
	TIME2 = tic();
	kernelsFilename = fullfile( CONFIG.DATA_path, sprintf('kernels_d=%4.2f',CONFIG.kernels.d) );
	if isempty( CONFIG.kernels.namePostfix )
		save( [ kernelsFilename '.mat'], 'KERNELS')
	else
		save( [ kernelsFilename '_' CONFIG.kernels.namePostfix '.mat'], 'KERNELS')
	end
	fprintf( ' [%.1f seconds]\n', toc(TIME2) );

	fprintf( '  [ %.1f seconds ]\n', toc(TIME) );
end
