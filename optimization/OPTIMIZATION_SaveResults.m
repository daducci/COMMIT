ticID = tic;

% add a postfix description to the folder name
if isfield( CONFIG.kernels, 'namePostfix' ) 
	RESULTS_path = [ 'Results_' CONFIG.kernels.namePostfix ];
else
	RESULTS_path = 'Results';
end

[~,tracking,~] = fileparts(CONFIG.TRACKING_path);
fprintf( '\n-> saving results to "%s/%s/Tracking/%s/%s/":\n', CONFIG.protocol, CONFIG.subject, tracking, RESULTS_path );
clear tracking

RESULTS_path = fullfile( CONFIG.TRACKING_path, RESULTS_path );
[~,~,~] = mkdir( RESULTS_path ); 


% configuration and results
% -------------------------
fprintf( '\t- configuration\n');
save( fullfile(RESULTS_path,'CONFIG.mat'), 'CONFIG' )
print( gcf, fullfile( RESULTS_path, sprintf('summary.png') ), '-dpng' );


% map of wovelwise errors
% -----------------------
fprintf( '\t- voxelwise map of fitting errors\n');
y_mea = niiSIGNAL.img;
y_est = zeros( niiSIGNAL.hdr.dime.dim(2:5) );
y_est( DICTIONARY.MASKidx ) = A * CONFIG.OPTIMIZATION.x;
if ( CONFIG.useReference )
	y_mea = y_mea(2:end,:,:,:);
	y_est = y_est(2:end,:,:,:);
end

% Reconstruction NRMSE
niiERR = niiSIGNAL;
niiERR.img = squeeze( sqrt(sum((y_mea-y_est).^2,1) ./ sum(y_mea.^2,1)) );
niiERR.img( DICTIONARY.MASK==0 | isnan(niiERR.img) | isinf(niiERR.img) ) = 0;
niiERR.hdr.dime.datatype = 16;
niiERR.hdr.dime.bitpix   = 32;
niiERR.hdr.dime.dim(1:5) = [ 3 niiERR.hdr.dime.dim(3:5) 1 ];
niiERR.hdr.dime.pixdim(2:5) = niiERR.hdr.dime.pixdim(3:6);
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = 1;
niiERR.hdr.dime.cal_min = 0;
niiERR.hdr.dime.cal_max = 1;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_NRMSE.nii') )

figure(99), clf, hist( niiERR.img(niiERR.img>0), 200 ), axis tight, grid on, title('NRMSE'), drawnow
print( gcf, fullfile( RESULTS_path, sprintf('fit_NRMSE.png') ), '-dpng' );
close( figure(99) )

% Reconstruction RMSE
niiERR.img = squeeze( sqrt(mean((y_mea-y_est).^2,1)) );
niiERR.img( DICTIONARY.MASK==0 ) = 0;
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = max( niiERR.img(:) );
niiERR.hdr.dime.cal_min = 0;
niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_RMSE.nii') )

figure(99), clf, hist( niiERR.img(niiERR.img>0), 200 ), axis tight, grid on, title('RMSE'), drawnow
print( gcf, fullfile( RESULTS_path, sprintf('fit_RMSE.png') ), '-dpng' );
close( figure(99) )

% Reconstruction SNR
niiERR.img = squeeze( 10*log10( mean(y_mea.^2,1) ./ mean((y_mea-y_est).^2,1) ) );
niiERR.img( DICTIONARY.MASK==0 | isnan(niiERR.img) | isinf(niiERR.img) ) = 0;
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = max( niiERR.img(:) );
niiERR.hdr.dime.cal_min = 0;
niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_SNR.nii') )

figure(99), clf, hist( niiERR.img(niiERR.img>0), 200 ), axis tight, grid on, title('SNR [dB]'), drawnow
print( gcf, fullfile( RESULTS_path, sprintf('fit_SNR.png') ), '-dpng' );
close( figure(99) )

% Volume fractions
% ----------------
fprintf( '\t- volume fractions\n' );

if ( CONFIG.kernels.doNormalize )
	xW = CONFIG.OPTIMIZATION.x ./ [ repmat( KERNELS.wmr_norm, 1, A.nF ) repmat( KERNELS.wmh_norm, 1, A.nE ) repmat( KERNELS.iso_norm, 1, A.nV ) ]';
else
	xW = CONFIG.OPTIMIZATION.x;
end

niiERR.hdr.dime.dim(1) = 4;
niiERR.hdr.dime.dim(5) = 3;
niiERR.hdr.dime.glmin = 0;
niiERR.hdr.dime.glmax = 1;
niiERR.hdr.dime.cal_min = niiERR.hdr.dime.glmin;
niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
niiERR.img = zeros(niiERR.hdr.dime.dim(2:5),'single');

% isotropic
fprintf( '\t\t- isotropic\n' );
IMG  = zeros( DICTIONARY.dim );
base = A.nF*A.nR + A.nE*A.nT;
for s = 1:numel(DICTIONARY.ISO.v)
	v = DICTIONARY.ISO.v(s) + 1;
	xx = xW(base+s : A.nV : base+A.nV*A.nI);
	for k = 1:A.nI
		IMG(v) = IMG(v) + xx(k);
	end
end
niiERR.img( :,:,:,3 ) = IMG;

% extra-axonal
if numel(KERNELS.wmh) > 0
	fprintf( '\t\t- extra-axonal\n' );
	IMG  = zeros( DICTIONARY.dim );
	base = A.nF*A.nR;
	for s = 1:numel(DICTIONARY.EC.v)
		v = DICTIONARY.EC.v(s) + 1;
		xx = xW( base+s: A.nE : base+A.nE*A.nT);
		for k = 1:A.nT
			IMG(v) = IMG(v) + xx(k);
		end
	end
	niiERR.img( :,:,:,2 ) = IMG;
end

% intra-axonal
fprintf( '\t\t- intra-axonal\n' );
IMG   = zeros( DICTIONARY.dim );
IMG2  = zeros( DICTIONARY.dim );
for s = 1:numel(DICTIONARY.IC.v)
	f = DICTIONARY.IC.fiber(s) + 1;
	v = DICTIONARY.IC.v(s) + 1;
	for k = 1:A.nR
		w = xW(f + (k-1)*A.nF) * DICTIONARY.IC.len(s);
		IMG(v) = IMG(v) + w;
	end
	IMG2(v) = IMG2(v) + DICTIONARY.IC.len(s);
end

niiERR.img( :,:,:,1 ) = IMG;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'vf.nii') )

niiERR.img( :,:,:,1 ) = IMG2;
save_untouch_nii( niiERR, fullfile(RESULTS_path,'vf__raw.nii') )


fprintf( '   [ %.2f seconds ]\n', toc(ticID) );

clear ticID RESULTS_path niiERR y_est y_mea IMG IMG2 w xW xx f v s k base

