classdef COMMIT_ACTIVEAX < handle

properties
    id, name   % id and name of the model
    dPar       % parallel diffusivity [units of mm^2/s]
    IC_Rs      % radii of the axons [units of 1E-6 (micrometers)]
    IC_VFs     % volume fractions of the axons
    dISOs      % isotropic diffusivities [units of mm^2/s]
    bscale     % scaling factor for b-values
end


methods

    % ===========================
    % Initialization of the model
    % ===========================
    function obj = COMMIT_ACTIVEAX()
        global CONFIG

        if ( CONFIG.scheme.version ~= 1 )
            error( '[COMMIT_ACTIVEAX] The scheme must be "STEJSKALTANNER" to use this model' )
        end

        obj.id        = 'ACTIVEAX';
        obj.name      = 'ActiveAx';
        obj.dPar      = 1.7 * 1E-3;
        obj.dISOs     = [2.0 3.0] * 1E-3;
        obj.IC_Rs     = [4 16] * 1e-6;
        obj.IC_VFs    = [0.7];
        obj.bscale    = 1;
    end


    % ===============================
    % Set the parameters of the model
    % ===============================
    function obj = Set( obj, dPar, IC_Rs, IC_VFs, dISOs )
        obj.dPar      = dPar;
        obj.dISOs     = dISOs;
        obj.IC_Rs     = IC_Rs;
        obj.IC_VFs    = IC_VFs;
    end


    % ==================================================================
    % Generate high-resolution kernels and rotate them in harmonic space
    % ==================================================================
    function GenerateKernels( obj, ATOMS_path, schemeHR, AUX, idx_IN, idx_OUT )
        global CONFIG COMMIT_data_path CAMINO_path

        % check if high-resolution scheme has been created
        schemeHrFilename = fullfile(ATOMS_path,'protocol_HR.scheme');
        if ~exist( schemeHrFilename, 'file' )
            error( '[COMMIT_ACTIVEAX.GenerateKernels] File "protocol_HR.scheme" not found in folder "%s"', ATOMS_path )
        end

        filenameHr = [tempname '.Bfloat'];

        % Restricted
        % ==========
        idx = 1;
        for R = obj.IC_Rs
            TIME = tic();
            fprintf( '\t\t- A_%03d... ', idx );

            % generate
            if exist( filenameHr, 'file' ), delete( filenameHr ); end
            CMD = sprintf( '%s/datasynth -synthmodel compartment 1 CYLINDERGPD %E 0 0 %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null', CAMINO_path, obj.dPar*1e-6, R, schemeHrFilename, filenameHr );
            [status result] = system( CMD );
            if status>0
                disp(result)
                error( '[COMMIT_ACTIVEAX.GenerateKernels] Problems generating the signal with datasynth' );
            end

            % rotate and save
            fid = fopen( filenameHr, 'r', 'b' );
            signal = fread(fid,'float');
            fclose(fid);
            delete( filenameHr );
            lm = COMMIT_RotateKernel( signal, AUX, idx_IN, idx_OUT, false );
            save( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), '-v6', 'lm' )

            idx = idx + 1;
            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end


        % Hindered
        % ========
        for ICVF = obj.IC_VFs
            TIME = tic();
            fprintf( '\t\t- A_%03d... ', idx );

            % generate
            d_perp = obj.dPar * ( 1.0 - ICVF );
            if exist( filenameHr, 'file' ), delete( filenameHr ); end
            CMD = sprintf( '%s/datasynth -synthmodel compartment 1 ZEPPELIN %E 0 0 %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null', CAMINO_path, obj.dPar*1e-6, d_perp*1e-6, schemeHrFilename, filenameHr );
            [status result] = system( CMD );
            if status>0
                disp(result)
                error( '[COMMIT_ACTIVEAX.GenerateKernels] problems generating the signal' );
            end

            % rotate and save
            fid = fopen( filenameHr, 'r', 'b' );
            signal = fread(fid,'float');
            fclose(fid);
            delete( filenameHr );
            lm = COMMIT_RotateKernel( signal, AUX, idx_IN, idx_OUT, false );
            save( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), '-v6', 'lm' )

            idx = idx + 1;
            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end


        % Isotropic
        % =========
        for dISO = obj.dISOs
            TIME = tic();
            fprintf( '\t\t- A_%03d... ', idx );

            % generate
            if exist( filenameHr, 'file' ), delete( filenameHr ); end
            CMD = sprintf( '%s/datasynth -synthmodel compartment 1 BALL %E -schemefile %s -voxels 1 -outputfile %s 2> /dev/null', CAMINO_path, dISO*1e-6, schemeHrFilename, filenameHr );
            [status result] = system( CMD );
            if status>0
                disp(result)
                error( '[COMMIT_ACTIVEAX.GenerateKernels] problems generating the signal' );
            end

            % resample and save
            fid = fopen( filenameHr, 'r', 'b' );
            signal = fread(fid,'float');
            fclose(fid);
            delete( filenameHr );
            lm = COMMIT_RotateKernel( signal, AUX, idx_IN, idx_OUT, true );
            save( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), '-v6', 'lm' )

            idx = idx + 1;
            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end
    end


    % ==============================================
    % Project kernels from harmonic to subject space
    % ==============================================
    function ResampleKernels( obj, ATOMS_path, idx_OUT, Ylm_OUT )
        global CONFIG COMMIT_data_path KERNELS

        % Setup the KERNELS structure
        % ===========================
        KERNELS = {};
        KERNELS.model = obj.id;
        KERNELS.nS    = CONFIG.scheme.nS;
        KERNELS.dPar  = obj.dPar;
        KERNELS.wmr   = [];
        KERNELS.wmh   = [];
        KERNELS.iso   = [];

        % Restricted
        % ==========
        idx = 1;
        KERNELS.wmr_radii = obj.IC_Rs;
        for R = obj.IC_Rs
            TIME = tic();
            fprintf( '\t- A_%03d...  ', idx );

            load( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), 'lm' );
            KERNELS.wmr{end+1}       = COMMIT_ResampleKernel( lm, idx_OUT, Ylm_OUT, false );
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end


        % Hindered
        % ========
        KERNELS.wmh_icvf = obj.IC_VFs;
        for ICVF = obj.IC_VFs
            TIME = tic();
            fprintf( '\t- A_%03d...  ', idx );

            load( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), 'lm' );
            KERNELS.wmh{end+1}      = COMMIT_ResampleKernel( lm, idx_OUT, Ylm_OUT, false );
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end


        % Isotropic
        % =========
        KERNELS.iso_d = obj.dISOs;
        for dISO = obj.dISOs
            TIME = tic();
            fprintf( '\t- A_%03d...  ', idx );

            load( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), 'lm' );
            KERNELS.iso{end+1}      = COMMIT_ResampleKernel( lm, idx_OUT, Ylm_OUT, true );
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end

    end

end

end
