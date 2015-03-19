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

if param.verbose>=1
    maxNorm  = norm(y);
    objTimes = nan*ones(param.max_iter,1);
    xNorm0   = nan*ones(param.max_iter,1);
end

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
    objTimes(iter) = sqrt(2*curr_obj);
    xNorm0(iter) = nnz( sum(reshape( sol(1:A.nR*A.nF), A.nF, A.nR ),2) );
    if param.verbose >= 1
        fprintf('%13.7e   %13.7e\n', objTimes(iter), rel_obj);
		if param.verbose >= 2
			COMMIT_PlotSolverStatus(A, sol, objTimes, xNorm0, iter, maxNorm );
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
if param.verbose >= 1
    fprintf('\n Solution found:\n');
    fprintf(' Final objective value: %e\n', curr_obj);
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit_BPDN);
    
    if param.verbose >= 2
        COMMIT_PlotSolverStatus(A, sol, objTimes, xNorm0, iter, maxNorm );
    end
end


end
