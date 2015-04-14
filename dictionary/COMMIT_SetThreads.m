%
% Set the number of threads to use for the computation of A*x and A'y products
%
% Parameters
% ----------
% nTHREADS : integer
%   Number of threads to use (default: 1)
%
function COMMIT_SetThreads( nTHREADS )
    global COMMIT_code_path CONFIG niiSIGNAL
    global DICTIONARY KERNELS THREADS

    if ( nargin < 1 || nTHREADS < 1 ), nTHREADS = 1; end
    THREADS = [];
    THREADS.n = nTHREADS;

    ticID = tic;
    fprintf( '\n-> Computing the load for each thread\n   ==================================\n' );
    fprintf( '\t- number of threads : %d\n', THREADS.n );

    % Distribute VOXELS to DIFFERENT THREADS to avoid conflicts in A*x product
    % ========================================================================
    fprintf( '\t- A*x product...  ' );

    % intra-cellular segments
    if DICTIONARY.IC.n > 0
        THREADS.IC = zeros( nTHREADS+1, 1, 'uint32' );
        if nTHREADS>1
            N = floor(DICTIONARY.IC.n/nTHREADS);
            C = histc(DICTIONARY.IC.v,[0:max(DICTIONARY.IC.v)]);
            segPerThread = zeros( nTHREADS+1, 1, 'uint32' );
            tID = 1;
            tot = 0;
            for i = 1:numel(C)
                tot = tot + C(i);
                if tot >= N
                    segPerThread( tID+1 ) = tot;
                    tID   = tID+1;
                    tot   = 0;
                end
            end
            segPerThread( end ) = tot;
            THREADS.IC = cumsum(segPerThread);
        end
        THREADS.IC(end) = DICTIONARY.IC.n;
    else
        THREADS.IC = [];
    end

    % extra-cellular segments
    if DICTIONARY.EC.nE > 0
        THREADS.EC = zeros( nTHREADS+1, 1, 'uint32' );
        for i = 1:nTHREADS
            THREADS.EC(i) = find( DICTIONARY.EC.v == DICTIONARY.IC.v(THREADS.IC(i)+1), 1, 'first' ) - 1;
        end
        THREADS.EC(end) = DICTIONARY.EC.nE;
    else
        THREADS.EC = [];
    end

    % isotropic segments
    if DICTIONARY.nV > 0
        THREADS.ISO = zeros( nTHREADS+1, 1, 'uint32' );
        for i = 1:nTHREADS
            THREADS.ISO(i) = find( DICTIONARY.ISO.v == DICTIONARY.IC.v(THREADS.IC(i)+1), 1, 'first' ) - 1;
        end
        THREADS.ISO(end) = DICTIONARY.nV;
    else
        THREADS.ISO = [];
    end

    fprintf( '[OK]\n' );


    % Distribute COMPARTMENTS to DIFFERENT THREADS to avoid conflicts in At*y product
    % ===============================================================================

    fprintf( '\t- A''*y product... ' );

    % intra-cellular segments
    % -----------------------
%     THREADS.ICt = (nTHREADS-1)*ones( DICTIONARY.IC.n , 1, 'uint8' );
%     if nTHREADS>1
%         [~,~,Uidx] = unique(DICTIONARY.IC.fiber,'stable');
%         C = accumarray(Uidx,1);
%         tID = 0;
%         tot = 0;
%         iPrev = 1;
%         for i = 1:numel(C)
%             if tot >= floor(DICTIONARY.IC.n/nTHREADS)
%                 THREADS.ICt( Uidx>=iPrev & Uidx<i ) = tID;
%                 tID   = tID+1;
%                 if tID==nTHREADS-1, break, end
%                 iPrev = i;
%                 tot   = C(i);
%             else
%                 tot = tot + C(i);
%             end
%         end
%     end

    THREADS.ICt = (nTHREADS-1)*ones( DICTIONARY.IC.n , 1, 'uint8' );
    if nTHREADS>1
        [~,idx] = sort(DICTIONARY.IC.fiber);
        C = histc(DICTIONARY.IC.fiber,[0:max(DICTIONARY.IC.fiber)]);
        tID = 0;
        tot = 0;
        i1 = 1;
        i2 = 1;
        N = floor(DICTIONARY.IC.n/nTHREADS);
        for i = 1:numel(C)
            i2 = i2 + C(i);
            tot = tot + C(i);
            if tot >= N
                THREADS.ICt(i1:i2) = tID;
                tID = tID+1;
                if tID==nTHREADS-1, break, end
                i1 = i2+1;
                tot = C(i);
            end
        end
    end
    THREADS.ICt(idx) = THREADS.ICt;


    % extra-cellular segments
    % -----------------------
    if DICTIONARY.EC.nE > 0
        THREADS.ECt = zeros( nTHREADS+1, 1, 'uint32' );
        N = floor(DICTIONARY.EC.nE/nTHREADS);
        for i = 1:nTHREADS-1
            THREADS.ECt(i+1) = THREADS.ECt(i) + N;
        end
        THREADS.ECt(end) = DICTIONARY.EC.nE;
    else
        THREADS.ECt = [];
    end

    % isotropic segments
    % ------------------
    if DICTIONARY.nV > 0
        THREADS.ISOt = zeros( nTHREADS+1, 1, 'uint32' );
        N = floor(DICTIONARY.nV/nTHREADS);
        for i = 1:nTHREADS-1
            THREADS.ISOt(i+1) = THREADS.ISOt(i) + N;
        end
        THREADS.ISOt(end) = DICTIONARY.nV;
    else
        THREADS.ISO = [];
    end
    fprintf( '[OK]\n' );

    fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
end
