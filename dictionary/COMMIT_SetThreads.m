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
    CONFIG.OPTIMIZATION = [];
    CONFIG.OPTIMIZATION.nTHREADS = nTHREADS;
    THREADS = [];

    ticID = tic;
    fprintf( '\n-> Computing the load for each thread\n   ==================================\n' );
    fprintf( '\t- number of threads : %d\n', CONFIG.OPTIMIZATION.nTHREADS );

    % Distribute VOXELS to DIFFERENT THREADS to avoid conflicts in A*x product
    % ========================================================================
    fprintf( '\t- A*x product...  ' );

    % intra-cellular segments
    % -----------------------
    tmp = (CONFIG.OPTIMIZATION.nTHREADS-1)*ones( DICTIONARY.IC.n , 1, 'uint8' );
    if CONFIG.OPTIMIZATION.nTHREADS>1
        [~,~,Uidx] = unique(DICTIONARY.IC.v,'stable');
        C = accumarray(Uidx,1);
        tID = 0;
        tot = 0;
        iPrev = 1;
        for i = 1:numel(C)
            if tot >= floor(DICTIONARY.IC.n/CONFIG.OPTIMIZATION.nTHREADS)
                tmp( Uidx>=iPrev & Uidx<i ) = tID;
                tID   = tID+1;
                if tID==CONFIG.OPTIMIZATION.nTHREADS-1, break, end
                iPrev = i;
                tot   = C(i);
            else
                tot = tot + C(i);
            end
        end
    end

    % store only the pointers of the start of each thread's block
    THREADS.IC = zeros( CONFIG.OPTIMIZATION.nTHREADS+1, 1, 'uint32' );
    for i = 0:CONFIG.OPTIMIZATION.nTHREADS-1
        THREADS.IC(i+1) = find( tmp==i, 1, 'first' ) - 1;
    end
    THREADS.IC(end) = DICTIONARY.IC.n;
    clear Uidx tID tot iPrev C

    % extra-cellular segments
    % -----------------------
    if DICTIONARY.EC.nE > 0
        tmp = (CONFIG.OPTIMIZATION.nTHREADS-1)*ones( DICTIONARY.EC.nE, 1, 'uint8' );
        for i = 1:CONFIG.OPTIMIZATION.nTHREADS
            tmp( DICTIONARY.EC.v >= DICTIONARY.IC.v(THREADS.IC(i)+1) & DICTIONARY.EC.v <= DICTIONARY.IC.v(THREADS.IC(i+1)) ) = i-1;
        end

        % store only the pointers of the start of each thread's block
        THREADS.EC = zeros( CONFIG.OPTIMIZATION.nTHREADS+1, 1, 'uint32' );
        for i = 0:CONFIG.OPTIMIZATION.nTHREADS-1
            THREADS.EC(i+1) = find( tmp==i, 1, 'first' ) - 1;
        end
        THREADS.EC(end) = DICTIONARY.EC.nE;
    else
        THREADS.EC = [];
    end

    % isotropic segments
    % ------------------
    tmp = (CONFIG.OPTIMIZATION.nTHREADS-1)*ones( DICTIONARY.nV, 1, 'uint8' );
    for i = 1:CONFIG.OPTIMIZATION.nTHREADS
        tmp( DICTIONARY.ISO.v >= DICTIONARY.IC.v(THREADS.IC(i)+1) & DICTIONARY.ISO.v <= DICTIONARY.IC.v(THREADS.IC(i+1)) ) = i-1;
    end

    % store only the pointers of the start of each thread's block
    THREADS.ISO = zeros( CONFIG.OPTIMIZATION.nTHREADS+1, 1, 'uint32' );
    for i = 0:CONFIG.OPTIMIZATION.nTHREADS-1
        THREADS.ISO(i+1) = find( tmp==i, 1, 'first' ) - 1;
    end
    THREADS.ISO(end) = DICTIONARY.nV;
    clear tmp i
    fprintf( '[OK]\n' );


    % Distribute COMPARTMENTS to DIFFERENT THREADS to avoid conflicts in At*y product
    % ===============================================================================

    fprintf( '\t- A''*y product... ' );

    % intra-cellular segments
    % -----------------------
    tmp = (CONFIG.OPTIMIZATION.nTHREADS-1)*ones( DICTIONARY.IC.n , 1, 'uint8' );
    if CONFIG.OPTIMIZATION.nTHREADS>1
        [~,~,Uidx] = unique(DICTIONARY.IC.fiber,'stable');
        C = accumarray(Uidx,1);
        tID = 0;
        tot = 0;
        iPrev = 1;
        for i = 1:numel(C)
            if tot >= floor(DICTIONARY.IC.n/CONFIG.OPTIMIZATION.nTHREADS)
                tmp( Uidx>=iPrev & Uidx<i ) = tID;
                tID   = tID+1;
                if tID==CONFIG.OPTIMIZATION.nTHREADS-1, break, end
                iPrev = i;
                tot   = C(i);
            else
                tot = tot + C(i);
            end
        end
    end

    % store only the pointers of the start of each thread's block
    THREADS.ICt = tmp;
    clear Uidx tID tot iPrev C i tmp

    % extra-cellular segments
    % -----------------------
    if DICTIONARY.EC.nE > 0
        THREADS.ECt = zeros( CONFIG.OPTIMIZATION.nTHREADS+1, 1, 'uint32' );
        n =    floor(DICTIONARY.EC.nE/CONFIG.OPTIMIZATION.nTHREADS);
        for i = 1:CONFIG.OPTIMIZATION.nTHREADS-1
            THREADS.ECt(i+1) = THREADS.ECt(i) + n;
        end
        THREADS.ECt(end) = DICTIONARY.EC.nE;
    else
        THREADS.ECt = [];
    end

    % isotropic segments
    % ------------------
    THREADS.ISOt = zeros( CONFIG.OPTIMIZATION.nTHREADS+1, 1, 'uint32' );
    n = floor(DICTIONARY.nV/CONFIG.OPTIMIZATION.nTHREADS);
    for i = 1:CONFIG.OPTIMIZATION.nTHREADS-1
        THREADS.ISOt(i+1) = THREADS.ISOt(i) + n;
    end
    THREADS.ISOt(end) = DICTIONARY.nV;
    fprintf( '[OK]\n' );

    fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
end
