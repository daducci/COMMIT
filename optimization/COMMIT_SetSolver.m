%
%
%
% Parameters
% ----------
%
function COMMIT_SetSolver( solver, optTol, maxIter, verbosity )
    global COMMIT_code_path CONFIG niiSIGNAL
    global DICTIONARY KERNELS Y A At

    ticID = tic;
    fprintf( '\n-> Setting the problem/solver\n   ==========================\n' );

    if ( nargin < 4), verbosity = 1; end
    if ( nargin < 3), maxIter = 100; end
    if ( nargin < 2), optTol = 1e-3; end
    if ( nargin < 1), solver = 'NNLS'; end
    
    CONFIG.OPTIMIZATION = [];
    CONFIG.OPTIMIZATION.verbosity = verbosity;
    CONFIG.OPTIMIZATION.optTol    = optTol;
    CONFIG.OPTIMIZATION.maxIter   = maxIter;
    CONFIG.OPTIMIZATION.solver    = solver;

    fprintf( '\t- Solver    : "%s"\n', CONFIG.OPTIMIZATION.solver );
    fprintf( '\t- Max iter  : %d\n',   CONFIG.OPTIMIZATION.maxIter );
    fprintf( '\t- Tol       : %.2e\n', CONFIG.OPTIMIZATION.optTol );


    % Compiling the mex-files encoding the linear operator A
    % ======================================================
    fprintf( '\t- compiling mex files... ' );
    mex( ...
        sprintf('-DnIC=%d',numel(KERNELS.wmr)), sprintf('-DnEC=%d',numel(KERNELS.wmh)), sprintf('-DnISO=%d',numel(KERNELS.iso)), ...
        '-O', '-silent', '-outdir', fullfile(COMMIT_code_path,'dictionary'), fullfile(COMMIT_code_path,'dictionary','COMMIT_MakeOperator_A.cpp') ...
    );

    mex( ...
        sprintf('-DnIC=%d',numel(KERNELS.wmr)), sprintf('-DnEC=%d',numel(KERNELS.wmh)), sprintf('-DnISO=%d',numel(KERNELS.iso)), ...
        '-O', '-silent', '-outdir', fullfile(COMMIT_code_path,'dictionary'), fullfile(COMMIT_code_path,'dictionary','COMMIT_MakeOperator_At.cpp') ...
    );

    A  = [];
    At = [];
    A  = COMMIT_MakeOperator();
    At = A';
    fprintf( '[OK]\n' );


    % Preparing array Y with the signal in each voxel
    % ===============================================
    fprintf( '\t- preparing Y array... ' );
    Y = double( niiSIGNAL.img( DICTIONARY.MASKidx ) );
    fprintf( '[OK]\n' );

    fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
end
