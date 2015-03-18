function sol = solve_nnls(y, A, param)
% solve_nnls - Solve non negative least squares problem.
%
% sol = solve_nnls(y, A, At, PARAM) solves:
%
%   min  0.5*||y-A x||_2^2 s.t. x >= 0
%
%
% Y contains the measurements. A is the forward measurement operator and
% At the associated adjoint operator.
% PARAM a Matlab structure containing the following fields:
%
%   General parameters:
%
%   - verbose: 0 no log, 1 print main steps, 2 print all steps.
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - rel_obj: minimum relative change of the objective value (default:
%   1e-4)
%       The algorithm stops if
%           | f(x(t)) - f(x(t-1)) | / f(x(t)) < rel_obj,
%       where x(t) is the estimate of the solution at iteration t.
%
%
%
% The problem is solved using the forward-backward algorithm with
% FISTA-like accelaration
%
%
%
% Author: Rafael Carrillo
% E-mail: carrillo.rafael@epfl.ch
% Date: Oct. 26, 2014
%

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end


% Initialization
if isfield(param,'initsol')
    xhat = param.initsol;
    res = A*xhat - y;
    grad = A'*res;
    prev_obj = 0.5*norm(res)^2;
else
    xhat = zeros(size(A,2),1);
	grad = -(A'*y);
    prev_obj = 0.5*norm(y)^2;
end

iter = 1;
prev_sol = xhat;
told = 1;
beta = 0.9;
qfval = prev_obj;

maxNorm  = norm(y);
objTimes = nan*ones(param.max_iter,1);
xNorm0   = nan*ones(param.max_iter,1);


%Step size computation
L = norm(A*grad)^2 / norm(grad)^2;
mu = 1.9 / L;

% Main loop
while 1

    %Log
    if param.verbose >= 1
        fprintf('%4i   ', iter);
	end

    %Gradient descend step
    sol = xhat - mu*grad;

    %Projection onto the positive orthant
    sol = real(sol);
    sol(sol<0) = 0;

    %Stepsize check
	tmp = sol-xhat;
    q = qfval + real(tmp'*grad)...
        + 0.5/mu*norm(tmp)^2;
    res = A*sol - y;
    curr_obj = 0.5*norm(res)^2;

    while (curr_obj > q)
        %Backtracking rule
        mu = beta*mu;
        %Gradient descend step
        sol = xhat - mu*grad;

        %Projection onto the positive orthant
        sol = real(sol);
        sol(sol<0) = 0;

        %New stepsize check
		tmp = sol-xhat;
        q = qfval + real(tmp'*grad)...
            + 0.5/mu*norm(tmp)^2;
        res = A*sol - y;
        curr_obj = 0.5*norm(res)^2;
    end


    % Global stopping criterion
    rel_obj  = abs(curr_obj - prev_obj)/curr_obj;
    if param.verbose >= 1
        fprintf('%13.7e   %13.7e\n', sqrt(2*curr_obj), rel_obj);
		if param.verbose >= 3
			objTimes(iter) = sqrt(2*curr_obj);
			xNorm0(iter) = nnz( sum(reshape( sol(1:A.nR*A.nF), A.nF, A.nR ),2) );
			plotResult(A, sol, objTimes, xNorm0, iter, maxNorm )
		end
	end

    if (rel_obj < param.rel_obj)
        crit_BPDN = 'REL_OBJ';
        break;
    elseif iter >= param.max_iter
        crit_BPDN = 'MAX_IT';
        break;
    end

    % FISTA update
    t = (1+sqrt(1+4*told^2))/2;
    xhat = sol + (told-1)/t * (sol - prev_sol);

    % Gradient computation
    res = A*xhat - y;
    grad = A'*res;

    % Update variables
    iter = iter + 1;
    prev_obj = curr_obj;
    prev_sol = sol;
    told = t;
    qfval = 0.5*norm(res)^2;

end

% Log
if param.verbose>=1
    % L1 norm
    fprintf('\n Solution found:\n');
    fprintf(' Final objective value: %e\n', curr_obj);

    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit_BPDN);
end

end


%================ PLOT ==================
function plotResult(A, x, objTimes, xNorm0, iter, maxNorm )
	sfigure(1); clf
	subplot(3,1,1), hold on

	% print IC coefficients
	if A.nR>=1
		idx = 1 : A.nF;
		plot( idx, x(idx), '.', 'Color',[1 0 0] )
	end
	if A.nR>=2
		idx = idx(end)+1 : idx(end)+A.nF;
		plot( idx, x(idx), '.', 'Color',[.5 0 0] )
	end
	if A.nR>=3
		idx = idx(end)+1 : idx(end)+A.nF;
		plot( idx, x(idx), '.', 'Color',[1 0 0] )
	end
	if A.nR>=4
		idx = idx(end)+1 : idx(end)+A.nF;
		plot( idx, x(idx), '.', 'Color',[.5 0 0] )
	end

	% print EC coefficients
	if A.nE > 0
		if A.nT>=1
			idx = A.nR*A.nF+1 : A.nR*A.nF+A.nE;
			plot( idx, x(idx), '.', 'Color',[0 1 0] )
		end
		if A.nT>=2
			idx = idx(end)+1 : idx(end)+A.nE;
			plot( idx, x(idx), '.', 'Color',[0 .5 0] )
		end
		if A.nT>=3
			idx = idx(end)+1 : idx(end)+A.nE;
			plot( idx, x(idx), '.', 'Color',[0 1 0] )
		end
		if A.nT>=4
			idx = idx(end)+1 : idx(end)+A.nE;
			plot( idx, x(idx), '.', 'Color',[0 .5 0] )
		end
	end

	% print CSF coefficients
	if A.nI>=1
		idx = idx(end)+1 : idx(end)+A.nV;
		plot( idx, x(idx), '.', 'Color',[0 0 1] )
	end
	if A.nI>=2
		idx = idx(end)+1 : idx(end)+A.nV;
		plot( idx, x(idx), '.', 'Color',[0 0 .5] )
	end
	if A.nI>=3
		idx = idx(end)+1 : idx(end)+A.nV;
		plot( idx, x(idx), '.', 'Color',[0 0 1] )
	end
	if A.nI>=4
		idx = idx(end)+1 : idx(end)+A.nV;
		plot( idx, x(idx), '.', 'Color',[0 0 .5] )
	end

	grid on, box on, axis tight
	title('x')

	subplot(3,1,2), plot( objTimes )
	axis([1 numel(objTimes) 0 maxNorm*1.02])
	grid on, title( sprintf('|| Ax - y ||_2 = %.2f',objTimes(iter)) )

	if A.nF > 0
		subplot(3,1,3), plot( xNorm0 )
		axis([1 numel(objTimes) 0 A.nF])
		grid on, title( sprintf('|| x_{fibers} ||_0 = %d (%.1f%%)', xNorm0(iter),xNorm0(iter)/A.nF*100) )
	end

	drawnow
end
