function show_status(iter, options, h)

if (options.verbose)
  if (options.asgui)
    waitbar(iter/options.maxit, h);
  else
    fprintf('.');
    if (mod(iter, 30) == 0)
      fprintf('\n');
    end
  end
end
