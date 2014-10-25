ticID = tic;
fprintf( '\n-> Solving inverse problem\n   =======================\n' );
fprintf( '\t - Solver    : "%s"\n', CONFIG.OPTIMIZATION.solver );
fprintf( '\t - Max iter  : %d\n', CONFIG.OPTIMIZATION.maxIter );
fprintf( '\t - Tol       : %.2e\n', CONFIG.OPTIMIZATION.optTol );
fprintf( '\t - # threads : %d\n\n', CONFIG.OPTIMIZATION.nTHREADS );

switch ( CONFIG.OPTIMIZATION.solver )

	case { 'BPDN', 'LASSO' }
		spgParam = spgSetParms( 'verbosity',CONFIG.OPTIMIZATION.verbosity, 'optTol',CONFIG.OPTIMIZATION.optTol, 'lsTol', 1e-4, 'bpTol', 1e-4, 'iterations',CONFIG.OPTIMIZATION.maxIter );
		spgParam.project     = @(x,weights,tau) NormL1NN_project(x,weights,tau);
		spgParam.primal_norm = @(x,weights    ) NormL1NN_primal(x,weights);
		spgParam.dual_norm   = @(x,weights    ) NormL1NN_dual(x,weights);
		spgParam.weights = 1;
		CONFIG.OPTIMIZATION.spgParam = spgParam;
		clear spgParam

		if ( strcmp(CONFIG.OPTIMIZATION.solver,'BPDN') )
% 			sigma_noise = sqrt(2/(4-pi)) / CONFIG.OPTIMIZATION.SNR_estimated;
			sigma_noise = 1 / CONFIG.OPTIMIZATION.SNR_estimated;
			CONFIG.OPTIMIZATION.epsilon = sqrt( chi2inv(0.99, 2*DICTIONARY.nV*KERNELS.nS) ) * sigma_noise;

			fprintf( '   - epsilon = %.2f\n', CONFIG.OPTIMIZATION.epsilon );
			[x, r, g, info] = spgl1(A, Y, [], CONFIG.OPTIMIZATION.epsilon, [], CONFIG.OPTIMIZATION.spgParam );
			for iter = 1:CONFIG.OPTIMIZATION.RW_maxIter
				x_prev = x;
				CONFIG.OPTIMIZATION.spgParam.weights = CONFIG.OPTIMIZATION.RW_tauN ./ ( abs(x_prev) + CONFIG.OPTIMIZATION.RW_tauD);
				[x, r, g, info] = spgl1(A, Y, [], CONFIG.OPTIMIZATION.epsilon, [], CONFIG.OPTIMIZATION.spgParam );
				if  abs( norm(x_prev,1)-norm(x,1) ) / norm(x_prev,1) < 1e-2, break, end
			end
		else
			[x, r, g, info] = spgl1(A, Y, 100*size(A,2), [], At*Y, CONFIG.OPTIMIZATION.spgParam );
		end
		CONFIG.OPTIMIZATION.x = x;
		CONFIG.OPTIMIZATION.output = info;
		clear sigma_noise x r g info x_prev iter

	case 'NNLS'
		opt = solopt;
		opt.use_tolo = true;
		opt.tolo = CONFIG.OPTIMIZATION.optTol;
		opt.verbose = CONFIG.OPTIMIZATION.verbosity;
		opt.compute_obj = 1;
		opt.maxit = CONFIG.OPTIMIZATION.maxIter;
 		x = At*Y;
		x(x<0) = 0;
		out = sbb( A, Y, x, opt );
		CONFIG.OPTIMIZATION.x = out.x;
		CONFIG.OPTIMIZATION.output = out;
		clear opt x out
end

CONFIG.OPTIMIZATION.time = toc(ticID);
fprintf( '   [ %.2f seconds ]\n', CONFIG.OPTIMIZATION.time );
clear ticID

