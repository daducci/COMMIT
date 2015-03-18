%
% Resample rotated kernels to the acquisition protocol used for the specific subject
%
% Parameters
% ----------
% lmax : unsigned int
%   Maximum spherical harmonics order to use for the rotation phase
%
function COMMIT_ResampleKernels( lmax )
    if nargin < 1, lmax = 12; end
    global CONFIG COMMIT_data_path KERNELS

    TIME = tic();
    fprintf( '\n-> Resampling rotated kernels for subject "%s":\n', CONFIG.subject );

    % check if original scheme exists
    if ~exist( CONFIG.schemeFilename, 'file' ) || isempty(CONFIG.scheme)
        error( '[COMMIT_ResampleKernels] File "%s" not found', CONFIG.schemeFilename )
    end

    % check if auxiliary matrices have been precomputed
    auxFilename = fullfile(COMMIT_data_path,sprintf('AUX_matrices__lmax=%d.mat',lmax) );
    if ~exist( auxFilename, 'file' )
        error( '[COMMIT_ResampleKernels] Auxiliary matrices "%s" not found', auxFilename )
    else
        AUX = load( auxFilename );
    end


    % Precompute aux data structures
    % ==============================
    nSH = (AUX.lmax+1)*(AUX.lmax+2)/2;
    idx_OUT = zeros(CONFIG.scheme.dwi_count,1,'single');
    Ylm_OUT = zeros(CONFIG.scheme.dwi_count,nSH*numel(CONFIG.scheme.shells),'single');
    idx = 1;
    for s = 1:numel(CONFIG.scheme.shells)
        nS = numel(CONFIG.scheme.shells{s}.idx);
        idx_OUT(idx:idx+nS-1) = CONFIG.scheme.shells{s}.idx;
        [colatitude, longitude] = COMMIT_Cart2sphere( CONFIG.scheme.shells{s}.grad(:,1), CONFIG.scheme.shells{s}.grad(:,2), CONFIG.scheme.shells{s}.grad(:,3) );
        Ylm_OUT(idx:idx+nS-1, [1:nSH]+(s-1)*nSH) = COMMIT_CreateYlm( AUX.lmax, colatitude, longitude ); % matrix from SH to real space
        idx = idx + nS;
    end


    % Dispatch to the right handler for each model
    % ============================================
    if ~isempty(CONFIG.model)
        CONFIG.model.ResampleKernels( fullfile(COMMIT_data_path,CONFIG.protocol,'kernels',CONFIG.model.id), idx_OUT, Ylm_OUT );
    else
        error( '[COMMIT_ResampleKernels] Model not set' )
    end


    % Ensure single data type for mex code
    % ------------------------------------
    fprintf('\t- convert to single data type...');
    for i = 1:numel(KERNELS.wmr)
        KERNELS.wmr{i} = single( KERNELS.wmr{i}(CONFIG.scheme.dwi_idx,:,:) );
    end
    for i = 1:numel(KERNELS.wmh)
        KERNELS.wmh{i} = single( KERNELS.wmh{i}(CONFIG.scheme.dwi_idx,:,:) );
    end
    for i = 1:numel(KERNELS.iso)
        KERNELS.iso{i} = single( KERNELS.iso{i}(CONFIG.scheme.dwi_idx,:,:) );
    end
    KERNELS.nS = CONFIG.scheme.dwi_count;
    fprintf( ' [ OK ]\n' );


    % De-mean kernels as LiFE
    % -----------------------
    if ( CONFIG.doDemean )
        fprintf('\t - de-meaning atoms...');
        for j = 1:181
        for k = 1:181
            for i = 1:numel(KERNELS.wmr)
                KERNELS.wmr{i}(:,j,k) = KERNELS.wmr{i}(:,j,k) - mean(KERNELS.wmr{i}(:,j,k),1);
            end
            for i = 1:numel(KERNELS.wmh)
                KERNELS.wmh{i}(:,j,k) = KERNELS.wmh{i}(:,j,k) - mean(KERNELS.wmh{i}(:,j,k),1);
            end
            for i = 1:numel(KERNELS.iso)
                KERNELS.iso{i}(:,j,k) = KERNELS.iso{i}(:,j,k) - mean(KERNELS.iso{i}(:,j,k),1);
            end
        end
        end
        fprintf( ' [ OK ]\n' );
    end


    % Add b0 reference (promote unitary sum)
    % --------------------------------------
    if ( CONFIG.useReference )
        fprintf('\t- adding reference...');
        for i = 1:numel(KERNELS.wmr)
            KERNELS.wmr{i} = KERNELS.wmr{i}([1 1:KERNELS.nS],:,:);
            KERNELS.wmr{i}(1,:,:) = 1;
        end
        for i = 1:numel(KERNELS.wmh)
            KERNELS.wmh{i} = KERNELS.wmh{i}([1 1:KERNELS.nS],:,:);
            KERNELS.wmh{i}(1,:,:) = 1;
        end
        for i = 1:numel(KERNELS.iso)
            KERNELS.iso{i} = KERNELS.iso{i}([1 1:KERNELS.nS],:,:);
            KERNELS.iso{i}(1,:,:) = 1;
        end
        KERNELS.nS = KERNELS.nS + 1;
        fprintf( ' [ OK ]\n' );
    end


    % Normalize atoms
    % ---------------
    if ( CONFIG.normalizeKernels )
        fprintf('\t- normalizing to have unitary norm...');
        KERNELS.wmr_norm = [];
        for i = 1:numel(KERNELS.wmr)
            KERNELS.wmr_norm(i) = norm( KERNELS.wmr{i}(:,1,1) );
            for j = 1:181
            for k = 1:181
                KERNELS.wmr{i}(:,j,k) = KERNELS.wmr{i}(:,j,k) / KERNELS.wmr_norm(i);
            end
            end
        end
        KERNELS.wmh_norm = [];
        for i = 1:numel(KERNELS.wmh)
            KERNELS.wmh_norm(i) = norm( KERNELS.wmh{i}(:,1,1) );
            for j = 1:181
            for k = 1:181
                KERNELS.wmh{i}(:,j,k) = KERNELS.wmh{i}(:,j,k) / KERNELS.wmh_norm(i);
            end
            end
        end
        KERNELS.iso_norm = [];
        for i = 1:numel(KERNELS.iso)
            KERNELS.iso_norm(i) = norm( KERNELS.iso{i}(:,1,1) );
            KERNELS.iso{i} = KERNELS.iso{i} / KERNELS.iso_norm(i);
        end
        fprintf( ' [ OK ]\n' );
    end

    fprintf( '   [ %.1f seconds ]\n', toc(TIME) );
end
