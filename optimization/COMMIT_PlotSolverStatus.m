function COMMIT_PlotSolverStatus(A, x, objTimes, xNorm0, iter, maxNorm )
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

	grid on, box on, axis tight, axis([0 numel(x)+1 0 max(x)])
	title('x')

	subplot(3,1,2), plot( objTimes )
    axis([1-1e-3 numel(objTimes)+1e-3 0 maxNorm*1.02])
	grid on, title( sprintf('|| Ax - y ||_2 = %.2f',objTimes(iter)) )

	if A.nF > 0
		subplot(3,1,3), plot( xNorm0 )
		axis([1-1e-3 numel(objTimes)+1e-3 0 A.nF])
		grid on, title( sprintf('|| x_{fibers} ||_0 = %d (%.1f%%)', xNorm0(iter),xNorm0(iter)/A.nF*100) )
	end

	drawnow
end
