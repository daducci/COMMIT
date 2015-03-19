%
% Save to file the output maps and the configuration
%
% Parameters
% ----------
% namePostfix : string
%   Append a given string/tag at the end of the folder where output data is saved (default: '')
%
function COMMIT_SaveResults( namePostfix )
    global CONFIG DICTIONARY KERNELS niiSIGNAL A
    ticID = tic;
    
    if ( nargin < 1 ), namePostfix = ''; end
    CONFIG.namePostfix = namePostfix;

    RESULTS_path = [ 'Results_' CONFIG.model.id ];
    if ( ~isempty(CONFIG.namePostfix) ), RESULTS_path = [ RESULTS_path '_' CONFIG.namePostfix]; end
    fprintf( '\n-> Saving results to "%s/*":\n', RESULTS_path );
    RESULTS_path = fullfile( CONFIG.TRACKING_path, RESULTS_path );
    [~,~,~] = mkdir( RESULTS_path );
    delete( fullfile(RESULTS_path,'*') );


    % Configuration and results
    % =========================
    fprintf( '\t- configuration... ');
    save( fullfile(RESULTS_path,'CONFIG.mat'), 'CONFIG', '-v7' )
    print( gcf, fullfile( RESULTS_path, sprintf('summary.png') ), '-dpng' );
    fprintf( '[OK]\n');


    % Map of wovelwise errors
    % =======================
    fprintf( '\t- fitting errors:\n');

    y_mea = niiSIGNAL.img;
    y_est = zeros( niiSIGNAL.hdr.dime.dim(2:5), 'single' );
    y_est( DICTIONARY.MASKidx ) = A * CONFIG.OPTIMIZATION.x;
    if ( CONFIG.useReference )
        y_mea = y_mea(2:end,:,:,:);
        y_est = y_est(2:end,:,:,:);
    end

    niiERR = [];
    niiERR.hdr = niiSIGNAL.hdr;
    niiERR.untouch = 1;
    niiERR.hdr.dime.dim(1:5)    = [ 3 CONFIG.dim 1 ];
    niiERR.hdr.dime.pixdim(2:5) = [ CONFIG .pixdim 1 ];
    niiERR.hdr.dime.datatype = 16;
    niiERR.hdr.dime.bitpix   = 32;
    niiERR.hdr.dime.glmin    = 0;
    niiERR.hdr.dime.cal_min  = 0;

    % Reconstruction RMSE
    fprintf( '\t\t- RMSE ' );
    niiERR.img = squeeze( sqrt(mean((y_mea-y_est).^2,1)) );
    niiERR.img( DICTIONARY.MASK==0 ) = 0;
    niiERR.hdr.dime.glmax = median( niiERR.img(:) );
    niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
    save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_RMSE.nii') )
    tmp = niiERR.img( DICTIONARY.MASK>0 );
    fprintf( ' [ %.3f +/- %.3f ]\n', mean(tmp), std(tmp) );

    % Reconstruction NRMSE
    fprintf( '\t\t- NRMSE ' );
    niiERR.img = squeeze( sqrt(sum((y_mea-y_est).^2,1) ./ sum(y_mea.^2,1)) );
    niiERR.img( DICTIONARY.MASK==0 | isnan(niiERR.img) | isinf(niiERR.img) ) = 0;
    niiERR.hdr.dime.glmax = 1;
    niiERR.hdr.dime.cal_max = 1;
    save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_NRMSE.nii') )
    tmp = niiERR.img( DICTIONARY.MASK>0 );
    fprintf( '[ %.3f +/- %.3f ]\n', mean(tmp), std(tmp) );

    % Reconstruction SNR
    fprintf( '\t\t- SNR  ' );
    niiERR.img = squeeze( 10*log10( mean(y_mea.^2,1) ./ mean((y_mea-y_est).^2,1) ) );
    niiERR.img( DICTIONARY.MASK==0 | isnan(niiERR.img) | isinf(niiERR.img) ) = 0;
    niiERR.hdr.dime.glmax = max( niiERR.img(:) );
    niiERR.hdr.dime.cal_max = niiERR.hdr.dime.glmax;
    save_untouch_nii( niiERR, fullfile(RESULTS_path,'fit_SNR.nii') )
    tmp = niiERR.img( DICTIONARY.MASK>0 );
    fprintf( ' [ %.3f +/- %.3f ]\n', mean(tmp), std(tmp) );

    clear y_est y_mea tmp
    fprintf( '\t  [OK]\n');


    % Volume fractions
    % ================
    fprintf( '\t- compartment maps:\n' );

    niiERR.img = [];
    niiERR.hdr.dime.glmax   = 1;
    niiERR.hdr.dime.cal_max = 1;
    niiERR.img = zeros( CONFIG.dim, 'single' );

    % renormalize the coefficients
    if ( CONFIG.normalizeKernels )
        xW = CONFIG.OPTIMIZATION.x ./ [ reshape(repmat(KERNELS.wmr_norm,A.nF,1),1,[]) reshape(repmat(KERNELS.wmh_norm,A.nE,1),1,[]) reshape(repmat(KERNELS.iso_norm,A.nV,1),1,[]) ]';
    else
        xW = CONFIG.OPTIMIZATION.x;
    end

    % isotropic
    fprintf( '\t\t- isotropic\n' );
    if numel(KERNELS.iso) > 0
        tmp = zeros( DICTIONARY.nV, 1, 'single' );
        base = A.nF*A.nR + A.nE*A.nT;
        for s = 1:numel(DICTIONARY.ISO.v)
            v = DICTIONARY.ISO.v(s)/KERNELS.nS + 1;
            xx = xW(base+s : A.nV : base+A.nV*A.nI);
            tmp(v) = tmp(v) + sum(xx);
        end
        niiERR.img( DICTIONARY.MASK>0 ) = tmp;
    else
        niiERR.img( DICTIONARY.MASK>0 ) = 0;
    end
    save_untouch_nii( niiERR, fullfile(RESULTS_path,'compartment_ISO.nii') )

    % extra-axonal
    fprintf( '\t\t- extra-axonal\n' );
    if numel(KERNELS.wmh) > 0
        tmp  = zeros( DICTIONARY.nV, 1, 'single' );
        base = A.nF*A.nR;
        for s = 1:numel(DICTIONARY.EC.v)
            v = DICTIONARY.EC.v(s)/KERNELS.nS + 1;
            xx = xW( base+s : A.nE : base+A.nE*A.nT );
            tmp(v) = tmp(v) + sum(xx);
        end
        niiERR.img( DICTIONARY.MASK>0 ) = tmp;
    else
        niiERR.img( DICTIONARY.MASK>0 ) = 0;
    end
    save_untouch_nii( niiERR, fullfile(RESULTS_path,'compartment_EC.nii') )

    % intra-axonal
    fprintf( '\t\t- intra-axonal\n' );
    xx = sum( reshape( xW( 1:A.nF*A.nR ), A.nF, A.nR ), 2 );
    tmp = zeros( DICTIONARY.nV, 1, 'single' );
    for s = 1:DICTIONARY.IC.n
        f = double(DICTIONARY.IC.fiber(s) + 1);
        v = DICTIONARY.IC.v(s)/KERNELS.nS + 1;
        tmp(v) = tmp(v) + xx(f) * DICTIONARY.IC.len(s);
    end
    niiERR.img( DICTIONARY.MASK>0 ) = tmp;
    save_untouch_nii( niiERR, fullfile(RESULTS_path,'compartment_IC.nii') )

    fprintf( '\t  [OK]\n');


    fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
end
