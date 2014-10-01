function options = solopt(varargin)
% SOLOPT  --  Creates a default options structure
%
% OPTIONS = SOLOPT
%
% The fields are described as follows
%
% ===========================================================================
%                         NLS
% ===========================================================================
% 
% OPTIONS.MAXIT      -- maximum number of iterations. Default = 13000
% OPTIONS.MAXTIME    -- maximum time (in seconds) optimization is allowed to
%                       run. In case you do not want this value to impact the 
%                       code, set the value of options.time_limit = 0
% OPTIONS.MAXMEM     -- maximum number of vectors for P-L-BFGS. Default = 7.
% OPTIONS.TIME_LIMIT -- has a default value of 1
% OPTIONS.USE_TOLX and OPTIONS.TOLX -- the first variable determines if
%                       relative change in the solution variable 'x' should 
%                       be used to determine a stopping criterion. The value 
%                       of this tolerance is 1e-6 by default and can be 
%                       changed by setting OPTIONS.TOLX
% OPTIONS.USE_TOLO and OPTIONS.TOLO -- stopping criterion based on relative 
%                       change in objective function value
% OPTIONS.USE_TOLP and OPTIONS.TOLP -- stopping criterion based on fixed-point 
%                       iteration of x, i.e., ||Proj(x - g) - x||_inf
% OPTIONS.USE_KKT and OPTIONS.TOLK -- stopping criterion based on KKT
%                       violation. Currently only support for nnls problems.
% OPTIONS.VERBOSE    -- whether to print info on screen while optimizing 
%                       (default = true)
% OPTIONS.TAU = 1e-7 -- This is a scaled version of the 'tau' parameter 
%                       usually smaller values lead to better answers.
% OPTIONS.COMPUTE_OBJ -- whether to compute obj. value at each iteration or 
%                       not (default is 1, i.e., yes) set to 0 to save some cpu
% OPTIONS.ASGUI      -- if to report progress of algo by a waitbar (default = 0)
% OPTIONS.PBB_GRADIENT_NORM -- Default is 1e-9
% OPTIONS.MAX_FUNC_EVALS -- maximum number of iterations for
%                       line-search. (default = 100).
% OPTIONS.BETA
% and OPTIONS.SIGMA  -- constants for line-search. Defaults are
%                       0.0498 and 0.298 respectively.
% OPTIONS.ALGO       -- select the underlying algorithm. Available algorithms 
%                       are  
%                      'PQN'  : projected quasi-Newton,
%                      'PLB'  : projected limited momory BFGS,
%                      'BBGP' : projected Barzilai-Borwein. Default is 'BBGP'.
% OPTIONS.OBJFN      -- name of objective function. Values honored by pqn,
%                       plb, etc. are:
%                       'lsq': for ||Ax-b||^2
%                       'kl' : for KL(b, Ax)
%                       ''
%
% ===========================================================================
%                         NMA
% ===========================================================================
%
% OPTIONS.PAD        -- small number to pad the denominator in Lee&Seung algos.
%                      (default = 1e-9).
% OPTIONS.TOLB       -- initial tolerance to optimize the factor B.
%                      (default = 1e-2).
% OPTIONS.TOLC       -- initial tolerance to optimize the factor C.
%                      (default = 1e-2).
% OPTIONS.STEP_B     -- step to decrease the tolerance for the factor B.
%                      (default = 1e-2).
% OPTIONS.STEP_C     -- step to decrease the tolerance for the factor C.
%                      (default = 1e-2).
% OPTIONS.TOLMIN     -- minimum tolerance for the factors B & C.
%                      (default = 1e-10).
% OPTIONS.FN_MAXIT  -- FNMAe requires substantially small number of
%                         iterations.
%                      (default = 30).
% OPTIONS.FN_ALGO   -- select the underlying algorithm for FNMAe. 
%                       Available algorithms are  
%                      'PQN'  : projected quasi-Newton,
%                      'PLB'  : projected limited momory BFGS,
%                      'BBGP' : projected Barzilai-Borwein. 
%                      NOTE that the default is 'PLB' for FNMAe.
% OPTIONS.MAX_FUNC_EVALS -- maximum number of iterations when 'BBGP' is selected
%                    for the underlying NLS solver of FNMAe. NMA_ALS & FNMAi 
%                    also use this value as the maximum number of iterations.
%                   (default = 100).  
%
% ===========================================================================
%                         NTA
% ===========================================================================
%
% OPTIONS.DECOMP     -- select the tensor decomposition.
%                       Available decompositions are  
%                      'PARAFAC' & 'TUCKER',  Default is 'PARAFAC'.
%
% OPTIONS.FN_MAXIT, 
% OPTIONS.FN_ALGO, 
% OPTIONS.TOLB, 
% OPTIONS.STEP_B,
% OPTIONS.TOLMIN.    -- FNTA shares these options with FNMAe.
%
% 


options.maxit = 50000;
options.maxtime = 20000;
options.maxmem = 7;
options.time_limit = 1;
options.use_tolx = 0;                   
options.use_tolo = 0;
options.use_tolg = 0;
options.use_tolp = 0;
options.use_kkt = 1;
options.tolx = 1e-8;                    % who knows things stop before 10sec!
options.tolo = 1e-6;
options.tolp = 1e-4;
options.tolk = 1e-6;
%options.tolg = 1e-4;
options.tolg = .5;
options.verbose = 0;                    % initially
options.tau = 1e-7;             
options.compute_obj = 0;
options.asgui = 0;
options.max_func_evals = 100;
options.pbb_gradient_norm = 1e-9;
options.beta = 0.0498;
options.sigma = 0.298;
options.eta = 0.99;
options.nulliter = 100;
options.objfn = 'lsq';

if nargin == 1
  options.algo = varargin{1};
else   % Default
  options.algo = 'BBGP';
end

% ===========================================================================
options.pad = 1e-8;
options.tolb = 1e-2;
options.tolc = 1e-2;
options.step_b = 1e-2;
options.step_c = 1e-2;
options.tolmin = 1e-8;
options.fn_maxit = 30;
options.fn_algo = 'PLB';

% ===========================================================================
options.decomp = 'PARAFAC';

