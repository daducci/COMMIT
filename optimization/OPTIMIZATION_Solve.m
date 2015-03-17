ticID = tic;
fprintf( '\n-> Solving inverse problem\n   =======================\n' );
fprintf( '\t - Solver    : "%s"\n', CONFIG.OPTIMIZATION.solver );
fprintf( '\t - Max iter  : %d\n',   CONFIG.OPTIMIZATION.maxIter );
fprintf( '\t - Tol       : %.2e\n', CONFIG.OPTIMIZATION.optTol );
fprintf( '\t - # threads : %d\n\n', CONFIG.OPTIMIZATION.nTHREADS );

switch ( CONFIG.OPTIMIZATION.solver )

    case 'NNLS'
		param = [];
		param.verbose  = CONFIG.OPTIMIZATION.verbosity;
		param.max_iter = CONFIG.OPTIMIZATION.maxIter;
		param.rel_obj  = CONFIG.OPTIMIZATION.optTol;
		x = solve_nnls( Y, A, param );
		CONFIG.OPTIMIZATION.x = x;
		clear param x
        
    otherwise
        error( '[OPTIMIZATION_Solve] Solver not recognized' );

end

CONFIG.OPTIMIZATION.time = toc(ticID);
fprintf( '   [ %.0fh %.0fm %.0fs ]\n', floor(CONFIG.OPTIMIZATION.time/3600), floor(mod(CONFIG.OPTIMIZATION.time/60,60)), mod(CONFIG.OPTIMIZATION.time,60) )
clear ticID

