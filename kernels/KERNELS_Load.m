% Load precomputed kernels
% 
% Parameters
% ----------
% CONFIG : struct
% 	Struct containing all the parameters and setup of the experiment
% 
% Returns
% -------
% KERNELS : struct
% 	Atoms for each compartment rotated along 181x181 directions
function [ KERNELS ] = KERNELS_Load( CONFIG )

	dictFilename = fullfile( CONFIG.DATA_path, sprintf('kernels_d=%4.2f',CONFIG.kernels.d) );
	if isempty( CONFIG.kernels.namePostfix )
		dictFilename = [ dictFilename '.mat'];
	else
		dictFilename = [ dictFilename '_' CONFIG.kernels.namePostfix '.mat'];
	end

	fprintf( '\n-> Loading precomputed kernels "%s":\n', dictFilename );
		
	% check if kernels file exists
	if ~exist( dictFilename, 'file' )
		error( '\t[KERNELS_Load] folder "%s" not found', dictFilename )
	end

	load( dictFilename );
	
	fprintf( '   [ %d intra-cellular, %d extra-cellular, %d isotropic ]\n', numel(KERNELS.wmr), numel(KERNELS.wmh), numel(KERNELS.iso) );
end
