%
% Simulate data with the Camino toolkit (its path must be set in following CAMINO_path variable)
% 
% Parameters
% ----------
% CONFIG : struct
% 	Struct containing all the parameters and setup of the experiment
%
function KERNELS_ActiveAx_Generate( CONFIG )
	TIME = tic();	
	global CAMINO_path

	fprintf( '\n-> Simulating high-resolution "ActiveAx" kernels with Camino for protocol "%s":\n', CONFIG.protocol );
	
	% check if original scheme exists
	if ~exist( CONFIG.schemeFilename, 'file' )
		error( '\t[SimulateDataWithCamino] File "%s" not found', CONFIG.schemeFilename )
	end

	% check if high-resolution scheme has been created
	[tmp,~,~] = fileparts( CONFIG.schemeFilename );
	schemeHrFilename = fullfile(tmp,'..','kernels','protocol_HR.scheme');
	if ~exist( schemeHrFilename, 'file' )
		error( '[SimulateDataWithCamino] File "protocol_HR.scheme" not found in folder "%s"', tmp )
	end
	
	% create folder where to store simulated kernels for this diffusivity value (d)
	OUTPUT_path = fullfile(tmp,'..','kernels',sprintf('d=%4.2f',CONFIG.kernels.d));
	[~,~,~] = mkdir( OUTPUT_path );
	
	% Simulate intra-cellular compartments
	fprintf( '\t- intra-cellular (varying the radius)\n' );
	for R = CONFIG.kernels.Rs
		fprintf( '\t\t- R = %5.2f micrometers\n', R );
		filename = fullfile(OUTPUT_path,sprintf('Er_R=%05.2f.Bfloat',R));
		if exist( filename, 'file' ), delete( filename ); end
		CMD = sprintf( '%s/datasynth -synthmodel compartment 1 CYLINDERGPD %E 0 0 %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null', CAMINO_path, CONFIG.kernels.d*1e-9, R*1e-6, schemeHrFilename, filename );
		[status result] = system( CMD );
		if status>0
			disp(result)
			error( '\t\t- ERROR: problems generating the signal\n' );
		end
	end
	fprintf( '\t  [ %.1f seconds ]\n', toc(TIME) );

	% Simulate extra-cellular compartments
	fprintf( '\t- extra-cellular (varying the intra-cellular volume fraction)\n' );
	TIME = tic();
	for ICVF = CONFIG.kernels.ICVFs
		d_perp = CONFIG.kernels.d * ( 1.0 - ICVF );
		fprintf( '\t\t- ICVF = %3.2f (d_perp = %5.2f)\n', ICVF, d_perp );
		filename = fullfile(OUTPUT_path,sprintf('Eh_icvf=%05.2f.Bfloat',ICVF));
		if exist( filename, 'file' ), delete( filename ); end
		CMD = sprintf( '%s/datasynth -synthmodel compartment 1 ZEPPELIN %E 0 0 %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null', CAMINO_path, CONFIG.kernels.d*1e-9, d_perp*1e-9, schemeHrFilename, filename );
		[status result] = system( CMD );
		if status>0
			disp(result)
			error( '\t\t- ERROR: problems generating the signal\n' );
		end
	end
	fprintf( '\t  [ %.1f seconds ]\n', toc(TIME) );
	
	% Simulate isotropic compartments
	fprintf( '\t- isotropic (varying the diffusivity)\n' );
	TIME = tic();
	for dISO = CONFIG.kernels.dISOs
		fprintf( '\t\t- d_iso = %5.2f\n', dISO );
		filename = fullfile(OUTPUT_path,sprintf('Ei_dIso=%05.2f.Bfloat',dISO));
		if exist( filename, 'file' ), delete( filename ); end
		CMD = sprintf( '%s/datasynth -synthmodel compartment 1 BALL %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null', CAMINO_path, dISO*1e-9, schemeHrFilename, filename );
		[status result] = system( CMD );
		if status>0
			disp(result)
			error( '\t\t- ERROR: problems generating the signal\n' );
		end
	end
	fprintf( '\t  [ %.1f seconds ]\n', toc(TIME) );
	
	fprintf( '  [ DONE ]\n' );
end
