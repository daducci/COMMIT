function v = check_termination(options, out)
%% ------------------------------------------------------
%  various termination criteria
%  ------------------------------------------------------

% First check if we are doing termination based on running time
if (options.time_limit)
%out.time = etime(clock, out.start_time);
out.time = toc(out.startTime);
  if (out.time >= options.maxtime)
    v = 1;
    return;
  end
end

% Now check if we are doing break by tolx
if (options.use_tolx)
  if (norm(out.x-out.oldx)/norm(out.oldx) < options.tolx)
    v = 2;
    return;
  end
end

% Are we doing break by tolo (tol obj val)
if (options.use_tolo && out.iter > 2)
  delta = abs(1-out.oldobj/out.obj);
  if (delta < options.tolo)
    v = 3;
    return;
  end
end

% Finally the plain old check if max iter has been achieved
if (out.iter >= options.maxit)
  v = 4;
  return;
end

% KKT violation
if (options.use_kkt)
    idx = find(out.x ~= 0 | (out.x == 0 & out.grad < 0));

    %if abs(out.x' * out.grad) <= options.tolk
    if norm(abs(out.grad(idx)), 'inf') <= options.tolk
        v = 7;
        return;
    end
end

% Gradient check
if (options.use_tolg)
  nr = norm(out.grad, 'inf');
  %  fprintf('%g\n', nr);
  if (nr < options.tolg)
    v = 8;
    return;
  end
end

% fixed-point threshold
if (options.use_tolp)
    tmp = out.x - out.grad;
    tmp(tmp < 0) = 0;
    if norm(tmp - out.x, 'inf') <= options.tolp
        v = 10;
        return;
    end
end


% No condition met, so return false
v = 0;
