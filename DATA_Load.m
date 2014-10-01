%
% Load data and perform some preprocessing
%
global niiSIGNAL

fprintf( '\n-> Loading and setup:\n' );

fprintf( '\t* Loading data...\n' );
niiSIGNAL = load_untouch_nii( CONFIG.dwiFilename );
niiSIGNAL.img = single(niiSIGNAL.img);
CONFIG.scheme = KERNELS_LoadScheme( CONFIG.schemeFilename );
fprintf( '\t\t* dim    = %d x %d x %d x %d\n' , niiSIGNAL.hdr.dime.dim(2:5) );
fprintf( '\t\t* pixdim = %.3f x %.3f x %.3f\n', niiSIGNAL.hdr.dime.pixdim(2:4) );
if CONFIG.scheme.nS == niiSIGNAL.hdr.dime.dim(5)
	fprintf( '\t\t* %d measurements divided in %d shells (%d b=0)\n', CONFIG.scheme.nS, numel(CONFIG.scheme.shells), CONFIG.scheme.b0_count );
	fprintf( '\t  [ OK ]\n' );
else
	error( '[main] Data and scheme do not match\n' );
end

if ( CONFIG.normalizeSignal )
	fprintf( '\t* Normalizing to b0...' );
	meanB0 = mean( niiSIGNAL.img(:,:,:,CONFIG.scheme.b0_idx), 4 );
	idx = meanB0 < 1;
	meanB0 = 1 ./ meanB0;
	meanB0( idx ) = 0;
	niiSIGNAL.img = niiSIGNAL.img .* repmat( meanB0, [1 1 1 niiSIGNAL.hdr.dime.dim(5)] );
	clear meanB0 idx
	fprintf( ' [ min=%.3f,  max=%.3f ]\n', min(niiSIGNAL.img(:)), max(niiSIGNAL.img(:)) );
end

fprintf( '\t* Removing b0 volume(s)...' );
niiSIGNAL.img = niiSIGNAL.img(:,:,:,CONFIG.scheme.dwi_idx);
niiSIGNAL.hdr.dime.dim(5) = CONFIG.scheme.dwi_count;
fprintf( ' [ %d x %d x %d x %d ]\n' , niiSIGNAL.hdr.dime.dim(2:5) );

if ( CONFIG.useReference )
	fprintf( '\t* Adding reference...' );
	niiSIGNAL.img = niiSIGNAL.img(:,:,:,[1 1:niiSIGNAL.hdr.dime.dim(5)]);
	niiSIGNAL.img(:,:,:,1) = 1;
	niiSIGNAL.hdr.dime.dim(5) = niiSIGNAL.hdr.dime.dim(5) + 1;
	fprintf( ' [ %d x %d x %d x %d ]\n' , niiSIGNAL.hdr.dime.dim(2:5) );
end

if ( CONFIG.doDemean )
	fprintf( '\t* Demean signal as LiFE...' );
	niiSIGNAL.img = niiSIGNAL.img - repmat( mean(niiSIGNAL.img,4), [1 1 1 niiSIGNAL.hdr.dime.dim(5)] );
	fprintf( ' [ min=%.3f,  max=%.3f ]\n', min(niiSIGNAL.img(:)), max(niiSIGNAL.img(:)) );
end


% prepare data for the optimization
% ---------------------------------
fprintf( '\t* Permute axes to improve performance...' );
niiSIGNAL.img = permute( niiSIGNAL.img, [4 1 2 3] );
niiSIGNAL.hdr.dime.dim(2:5) = [ size(niiSIGNAL.img,1) size(niiSIGNAL.img,2) size(niiSIGNAL.img,3) size(niiSIGNAL.img,4) ];
niiSIGNAL.hdr.dime.pixdim(2:5) = niiSIGNAL.hdr.dime.pixdim([5 2:4]);
fprintf( ' [ %d x %d x %d x %d ]\n' , niiSIGNAL.hdr.dime.dim(2:5) );

fprintf( '   [ DONE ]\n' );
