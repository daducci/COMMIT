%
% Solve the inverse model
%
function COMMIT_Fit( )
    global CONFIG Y A At

    ticID = tic;
    fprintf( '\n-> Solving inverse problem\n   =======================\n' );

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
            error( '[COMMIT_Fit] Solver not recognized' );

    end

    CONFIG.OPTIMIZATION.time = toc(ticID);
    fprintf( '%dh %02dm %02ds\n', floor(mod(CONFIG.OPTIMIZATION.time, [0, 3600, 60]) ./ [3600, 60, 1]) )
end
