function out = sbb(A, b, x0, opt)

% function out = sbb(A, b, x0, opt)
% Solve a nonnegative least squares problem, min ||Ax - b||_2, s.t. x \ge 0
% Barzilai-Borwein Gradient Projection.
%
%
% x0 -- Starting vector (useful for warm-starts).
%
% OPTIONS -- This structure contains important opt that control how the
% optimization procedure runs. To obtain a default structure the user can
% use 'opt = solopt'. Use 'help solopt' to get a description of
% what the individual opt mean.
%
% OUT contains the solution and other information about the optimization or
% info indicating whether the method succeeded or failed.
%
% See also: solopt
%
% Version 1.2 (c) 2009  Dongmin Kim  and Suvrit Sra
% 
    global gA gb startTime
	
	maxNorm = norm(b);

    % initialize parameters
    gA = A;
    gb = b;
    startTime = tic;
    
    out = getout();
    out.x = x0;
    out.refx = x0;
    out.obj = func(out.x);
    out.refobj = out.obj;
    out.grad = gradf(out.x);
    out.refgrad = out.grad;
    out.oldgrad = out.grad;
    out.iter = 0;
    out.time = 0;
    out.algo = 'SBB';
    out.startTime = startTime;
    
    out.status = 'Failure';
    out.step = 1;
    out.fcnt = 0;
    out.gcnt = 1;
	
	out.objTimes  = nan*ones(opt.maxit,1);
	out.xNorm0    = nan*ones(opt.maxit,1);

    h = 0;  % don't show status for now,
    beta = 1; % no scaling at the beginning,
    gtrh = 1e-9;  % threshold for gradient scaling,

    % threshold for identifying active vars,
    % opt.tolk * 1e-2 makes sense in practice and results in better
    % identification of final active variables. However, this also
    % entails unnecessary--in my opinion--theoretical complication 
    % hence to keep the brevity of our paper, initialize it to 1e-9.
    % We should recover the original heuristic for the official release.
    ptrh = 1e-9;  

    gp = find(out.x == 0 & out.grad > 0); % initial free vars.
    
    % Main loop
    while 1
        out.iter = out.iter + 1;
		
		% check descent, compute obj, count it!!
		out.obj = func(out.x);
		out.fcnt = out.fcnt + 1;
		out.objTimes(out.iter) = sqrt(2*out.obj);
		out.xNorm0(out.iter) = nnz( sum(reshape( out.x(1:A.nR*A.nF), A.nF, A.nR ),2) );
		fprintf( '%4d\t%.3E', out.iter, out.objTimes(out.iter) );
		if out.iter>1
			fprintf( '\t%.3E', abs(1-out.oldobj/out.obj));
		end
		fprintf('\n');

			
		if ( opt.verbose >= 3 )
			plotResult(A, out, maxNorm )
		end


            if ((out.refobj - out.obj) < opt.sigma ...
                    * (out.refgrad' * (out.refx - out.x)))
                beta = beta * opt.eta;
                fprintf('diminishing beta!\n');
            else
                out.refx = out.x;
                out.refobj = out.obj;
                out.refgrad = out.grad;
            end
%         end
        
        gp = find(out.x < ptrh & out.grad > ptrh);
        out.x(gp) = 0;
        out.grad(gp) = 0;
        out.oldgrad(gp) = 0;

        Ag = A * out.oldgrad;
        
        % termination
        if norm(Ag) < gtrh
            term_reason = 11;
            break;
        else
            term_reason = check_termination(opt, out);
            if (term_reason > 0)
                break;
            end
        end
        
        % update x & gradient
        if (mod(out.iter, 2) == 0)
            step = (out.oldgrad' * out.oldgrad) / (Ag' * Ag);
        else
            numer = Ag' * Ag;
            Ag = A' * Ag;
            Ag(gp) = 0;
            step = numer / (Ag' * Ag);
        end

        out.x = out.x - beta * step * out.grad;
        out.x(out.x < 0) = 0;                      % projection
        out.oldgrad = out.grad;
        out.grad = gradf(out.x);
        out.gcnt = out.gcnt + 1;
        
        if (opt.compute_obj)
            out.oldobj = out.obj;
            out.obj = func(out.x);
            out.fcnt = out.fcnt + 1;
        end
        
    end % of while

    out.obj = func(out.x);                % objective value

	plotResult(A, out, maxNorm )
	
	
    %% ------------------------------------------------
    %  Final statistics and wrap up
    %  ------------------------------------------------
    out.time = toc(out.startTime);
    out.status = 'Success';
    out.term_reason = set_term_reason(term_reason);
    idx = find(out.x > 0 | (out.x == 0 & out.grad < 0));
    out.gradnorm = norm(abs(out.grad(idx)), 'inf');
    
    if (opt.verbose)
        fprintf('DONE: "%s"\n', out.term_reason);
    end
end
 
%--------------------------------------------------------------------------
% Objective functions and gradients
%
function fc = func(x)
    global gA gb;

    Ax = gA * x;
    fc = 0.5 * norm(Ax - gb)^2;
end

function gc = gradf(x)
    global gA gb;

    Ax = gA * x;
    gc = (gA' * (Ax - gb));
end


%================ PLOT ==================
function plotResult(A, out, maxNorm )
	sfigure(1); clf
	subplot(3,1,1), hold on

	x = out.x;
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

	subplot(3,1,2), plot( out.objTimes )
	axis tight;	if max(out.objTimes)>0, ylim([0 maxNorm*1.02]),	end
	grid on, title( sprintf('|| Ax - y ||_2 = %.2f',out.objTimes(out.iter)) )

	if A.nF > 0
		subplot(3,1,3), plot( out.xNorm0 )
		axis tight, ylim([0 A.nF])
		grid on, title( sprintf('|| x_{fibers} ||_0 = %d (%.1f%%)', out.xNorm0(out.iter),out.xNorm0(out.iter)/A.nF*100) )
	end

	drawnow
end
