import numpy as np
from math import sqrt
from warnings import warn
from commit.proximals import fista, omega_hierarchical, prox_hierarchical
import sys
eps = np.finfo(float).eps

def issolver(solver=None):
    """
    Check if this version of COMMIT is endowed with the wanted solver.

    Parameters
    ----------
    solver: string (default=None)
        name of the solver

    Returns
    ----------
    Boolean if solver is a non-empty string, otherwise it returns the list of
    available solvers

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    """
    list_of_solvers = ['nnls','nnlsl1','gnnls','hnnls']
    if solver==None:
        return list_of_solvers
    elif solver in list_of_solvers:
        return True
    else:
        return False

def init_regularisation(thenorm = 2):
    """Initialise the structure that describes the regularisation term
    $ \Omega(x) = \lambda \sum_{g\inG} w_g \| x_{|g} \|_n $

    Parameters
    ----------
    thenorm : integer
        n-norm to be used in the sum (default: 2).

    Returns
    ----------
    Python dictionary with the following fields
        startIC : integer
            index from which the IntraCellular (IC) compartment starts.
        sizeIC : np.array of integers
            number of coefficients related to the EC compartment OR list of sizes of each group of the IC compartment (check the documentation of commit.solvers.gnnls).
        structureIC : np.array of objects
            list of indices of each group of the structured sparsity penalty term (check the documentation of commit.solvers.{gnnls,hnnls}).
        weightsIC: np.array of float64
            weights associated to each group encoded in sizeIC or structureIC
        lambdaIC : float64
            regularisation parameter for the IC compartment
        normIC : integer
            n-norm for the regularisation of the IC compartment.

        startEC : integer
            index from which the ExtraCellular (EC) compartment starts.
        sizeEC : integer or np.array
            number of coefficients related to the EC compartment.
        lambdaEC : float64
            regularisation parameter for the EC compartment.
        normEC : integer
            n-norm for the regularisation of the EC compartment.

        startISO : integer
            index from which the isotropic (ISO) compartment starts.
        sizeISO : integer
            number of coefficients related to the ISO compartment.
        lambdaISO : float64
            regularisation parameter for the ISO compartment.
        normISO : integer
            n-norm for the regularisation of the ISO compartment.

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    """
    regularisation = {}
    regularisation['startIC'] = 0
    regularisation['sizeIC'] = None # TODO: get the size from the mit object
    regularisation['structureIC'] = None
    regularisation['weightsIC'] = None

    regularisation['startEC'] = None
    regularisation['sizeEC'] = None

    regularisation['startISO'] = None
    regularisation['sizeISO'] = None

    regularisation['lambdaIC'] = .0
    regularisation['lambdaEC'] = .0
    regularisation['lambdaISO'] = .0

    regularisation['normIC'] = int(thenorm)
    regularisation['normEC'] = int(thenorm)
    regularisation['normISO'] = int(thenorm)

    # Solver-specific fields
    regularisation['group_is_ordered'] = False

    return regularisation

## Interface for NNLS solver
def nnls(y, A, At, tol_fun = 1e-4, tol_x = 1e-6, max_iter = 1000, verbose = 1):
    """Solve Non-Negative Least Squares (NNLS)

       min 0.5 * || y - Ax ||_2^2,          s.t. x >= 0

    Parameters
    ----------
    y : 1-d array of doubles.
        Signal to be fit.

    A : matrix or object endowed with the .dot() method.
        Dictionary describing the forward model.

    At : matrix or class exposing the .dot() method.
        Adjoint operator of A.

    tol_fun : double, optional (default: 1e-4).
        Minimum relative change of the objective value. The algorithm stops if:
               | f(x(t)) - f(x(t-1)) | / f(x(t)) < tol_fun,
        where x(t) is the estimate of the solution at iteration t.

    tol_x : double, optional (default: 1e-6).
        Minimum relative change of the solution x. The algorithm stops if:
               || x(t) - x(t-1) || / || x(t) || < tol_x,
        where x(t) is the estimate of the solution at iteration t.

    max_iter : integer, optional (default: 1000).
        Maximum number of iterations.

    verbose : integer, optional (default: 1).
        0 no log, 1 print each iteration results.

    Returns
    -------
    x : 1-d array of doubles.
        Minimizer of the NNLS problem.

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    """
    x0 = np.zeros( A.shape[1], dtype=np.float64 )

    def prox(x):
        x[x<0.0] = 0.0
        return x
    omega = lambda x: 0.0

    return fista(y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, prox)

## Interface for GNNLS solver
def gnnls(y, A, At, tol_fun = 1e-4, tol_x = 1e-6, max_iter = 1000, verbose = 1, regularisation = None) :
    """Solve Group-structured Non-Negative Least Squares (GNNLS)

       min 0.5 * || Ax - y ||_2^2 + \Omega(x),      s.t. x >= 0
        where \Omega(x) is a regularisation functional of the form
            \Omega(x) = \lambda \sum_{g\in G} \|x_{g}\|.
        The groups defined in G must not overlap.
        If lambda=0 the solver is equivalent to NNLS.

    Parameters
    ----------
    y : 1-d array of doubles.
        Signal to be fit.

    A : matrix or object endowed with the .dot() method.
        Dictionary describing the forward model.

    At : matrix or class exposing the .dot() method.
        Adjoint operator of A.

    tol_fun : double, optional (default: 1e-4).
        Minimum relative change of the objective value. The algorithm stops if:
               | f(x(t)) - f(x(t-1)) | / f(x(t)) < tol_fun,
        where x(t) is the estimate of the solution at iteration t.

    tol_x : double, optional (default: 1e-6).
        Minimum relative change of the solution x. The algorithm stops if:
               || x(t) - x(t-1) || / || x(t) || < tol_x,
        where x(t) is the estimate of the solution at iteration t.

    max_iter : integer, optional (default: 1000).
        Maximum number of iterations.

    verbose : integer, optional (default: 1).
        0 = no log, 1 = print each iteration results.

    regularisation : python dictionary
        Structure initialised by init_regularisation(). The necessary fields are:
            *lambdaIC
            *normIC
            *sizeIC or structureIC
            *weightsIC
            *startIC
            *group_is_ordered

        If group_is_ordered==True ,field sizeIC must be a list that encodes the length of each group, e.g.
            sizeIC = np.array([1,23,4,13,48])
        means that starting from startIC=0 we have five groups represented by the
        streamlines associated to indices (python-like)
            0, 1:23, 23:27, 27:40, 40:88
        If group_is_ordered==False, the

    Returns
    -------
    x : 1-d array of doubles.
        Minimizer of the GNNLS problem.

    References
    -----
    [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    """
    x0 = np.ones( A.shape[1], dtype=np.float64 )

    lambdaIC = regularisation.get('lambdaIC')
    if lambdaIC == 0.0:
        return nnls(y, A, At, tol_fun, tol_x, max_iter, verbose)

    normIC      = regularisation.get('normIC')
    sizeIC      = regularisation.get('sizeIC')
    weightsIC   = regularisation.get('weightsIC')
    startIC     = regularisation.get('startIC')
    structureIC = regularisation.get('structureIC')

    if normIC != np.inf and normIC!= 2:
        raise ValueError('Only 2-norm and infty-norm are allowed for the regularisation term.')

    # make the list-like structure suitable for the solver
    if regularisation['group_is_ordered']==True:
        warn('The ordered group structure will be deprecated. Consider using the structureIC field for defining the group structure.', DeprecationWarning)
        if len(sizeIC) != len(weightsIC):
            raise ValueError( 'weightsIC and sizeIC have different lengths.' )

        bundles = np.cumsum(np.insert(sizeIC,0,0))
        structureIC = np.array([range(bundles[k],bundles[k+1]) for k in range(0,len(bundles)-1)])
        regularisation['structureIC'] = structureIC
        del bundles
    elif len(structureIC) != len(weightsIC):
        raise ValueError('weightsIC and structureIC have different lengths.')

    omega = lambda x: omega_hierarchical( x, structureIC, weightsIC, lambdaIC, normIC )
    prox  = lambda x:  prox_hierarchical( x, structureIC, weightsIC, lambdaIC, normIC )

    return  fista(y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, prox)

## Interface for HNNLS solver
def hnnls(y, A, At, tol_fun = 1e-4, tol_x = 1e-6, max_iter = 1000, verbose = 1, regularisation = None) :
    """Solve Hierarchical Non-Negative Least Squares (HNNLS)

            min 0.5 * || y - Ax ||_2^2 + \Omega(x),      s.t. x >= 0
        where \Omega(x) is a regularisation functional of the form
            \Omega(x) = \lambda \sum_{g\in G} \|x_{g}\|
        and the groups in G define an hierarchical structure.

        E.g.: the tree structure
                        []
                      /    \
                    []      []
                   /  \    /  \
                 [a]  [b][c]  [d]
        corresponds to the structure
        G = { ([a],[b],[a,b],[c],[d],[c,d],[a,b,c,d]) , << },
        where << is a large order relation. See [1] for an extended description
        of the hierarchical structure required by the solver.

        The tree has not necessarily to be binary and contains the particular
        case of group structure (see gnnls solver).

        If lambda=0 the solver is equivalent to NNLS.


    Parameters
    ----------
    y : 1-d array of doubles.
        Signal to be fit.

    A : matrix or object endowed with the .dot() method.
        Dictionary describing the forward model.

    At : matrix or class exposing the .dot() method.
        Adjoint operator of A.

    tol_fun : double, optional (default: 1e-4).
        Minimum relative change of the objective value. The algorithm stops if:
               | f(x(t)) - f(x(t-1)) | / f(x(t)) < tol_fun,
        where x(t) is the estimate of the solution at iteration t.

    tol_x : double, optional (default: 1e-6).
        Minimum relative change of the solution x. The algorithm stops if:
               || x(t) - x(t-1) || / || x(t) || < tol_x,
        where x(t) is the estimate of the solution at iteration t.

    max_iter : integer, optional (default: 1000).
        Maximum number of iterations.

    verbose : integer, optional (default: 1).
        0 no log, 1 print each iteration results.

    regularisation : python dictionary
        Structure initialised by init_regularisation(). The necessary fields are:
            *lambdaIC
            *normIC
            *structureIC
            *weightsIC
            *startIC
        Field structureIC is supposed to be given as np.array of lists, which
        will result in a dtype=object np.array.

    TODO
    -------
    Regularisation of EC and ISO compartments

    Returns
    -------
    x : 1-d array of doubles.
        Minimizer of the HNNLS problem.

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    Acknowledgment: Rafael Carrillo - lts5 @ EPFL
    References:
        [1] Jenatton et al. - `Proximal Methods for Hierarchical Sparse Coding`
    """
    x0 = np.ones( A.shape[1], dtype=np.float64 )

    if regularisation is None:
        raise ValueError('The given tree structure is empty. Check the documentation.')

    lambdaIC = regularisation.get('lambdaIC')
    if lambdaIC == 0.0:
        return nnls(y, A, At, tol_fun, tol_x, max_iter, verbose)

    normIC = regularisation.get('normIC')
    structureIC = regularisation.get('structureIC')
    weightsIC = regularisation.get('weightsIC')

    if normIC != 1 and normIC!= 2:
        raise ValueError('Only 1-norm and 2-norm are allowed for the regularisation term.')

    if len(structureIC) != len(weightsIC):
        raise ValueError( 'weightsIC and structureIC have different lengths.' )

    omega = lambda x: omega_hierarchical( x, structureIC, weightsIC, lambdaIC, normIC )
    prox  = lambda x:  prox_hierarchical( x, structureIC, weightsIC, lambdaIC, normIC )

    return fista(y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, prox)

## Interface for NNLSL1 solver
def nnlsl1(y, A, At, tol_fun = 1e-4, tol_x = 1e-6, max_iter = 1000, verbose = 1, regularisation = None) :
    """Solve Non-Negative Least Squares with L1 regularization (NNLSL1)

        min 0.5 * || y - Ax ||_2^2 + lambda_IC ||x_IC||_1
                                   + lambda_EC ||x_EC||_1
                                   + lambda_ISO ||x_ISO||_1,
            s.t. x >= 0.

    If lambda*=0 the solver is equivalent to NNLS.

    The regularisation term must contain all the informations regarding start
    and size of each compartment. The variables related EC and ISO can be set
    to None and this results in an L1 regularisation of only the IC compartment.

    Parameters
    ----------
    y : 1-d array of doubles.
        Signal to be fit.

    A : matrix or object endowed with the .dot() method.
        Dictionary describing the forward model.

    At : matrix or class exposing the .dot() method.
        Adjoint operator of A.

    tol_fun : double, optional (default: 1e-4).
        Minimum relative change of the objective value. The algorithm stops if:
               | f(x(t)) - f(x(t-1)) | / f(x(t)) < tol_fun,
        where x(t) is the estimate of the solution at iteration t.

    tol_x : double, optional (default: 1e-6).
        Minimum relative change of the solution x. The algorithm stops if:
               || x(t) - x(t-1) || / || x(t) || < tol_x,
        where x(t) is the estimate of the solution at iteration t.

    max_iter : integer, optional (default: 1000).
        Maximum number of iterations.

    verbose : integer, optional (default: 1).
        0 no log, 1 print each iteration results.

    regularisation : python dictionary
        Structure initialised by init_regularisation(). The only necessary
        fields are *lambdaIC and *sizeIC

    Returns
    -------
    x : 1-d array of doubles.
        Minimizer of the NNLSL1 problem.

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    Acknowledgment: Rafael Carrillo - lts5 @ EPFL
    References:
        [1] Francis Bach et al. "Convex optimization with sparsity-inducing norms."
    """
    x0 = np.ones( A.shape[1], dtype=np.float64 )

    if regularisation is None:
        raise ValueError('The given tree structure is empty. Check the documentation.')


    # Regularise IC compartment
    lambdaIC = regularisation.get('lambdaIC')
    if lambdaIC == 0.0:
        return nnls(y, A, At, tol_fun, tol_x, max_iter, verbose)
    startIC = regularisation.get('startIC')
    if startIC == None:
        startIC = 0
    sizeIC = regularisation.get('sizeIC')
    if sizeIC == None:
        sizeIC = len(x0)
    omegaIC = lambda x: lambdaIC * sum(x[startIC:sizeIC])
    proxIC  = lambda x: prox_nnl1( x, lambdaIC, startIC, sizeIC )
    # if (not lambdaIC == None) and (not startIC == None) and (not sizeIC == None):
    #     omegaIC = lambda x: lambdaIC * sum(x[startIC:sizeIC])
    #     proxIC  = lambda x: prox_nnl1( x, lambdaIC, startIC, sizeIC )
    # else:
    #     omegaIC = lambda x: 0.0
    #     proxIC  = lambda x: np.zeros_like(x)

    # Regularise EC compartment
    lambdaEC = regularisation.get('lambdaEC')
    startEC = regularisation.get('startEC')
    sizeEC = regularisation.get('sizeEC')
    if (not lambdaEC == None) and (not startEC == None) and (not sizeEC == None):
        omegaEC = lambda x: lambdaEC * sum(x[startEC:sizeEC])
        proxEC  = lambda x: prox_nnl1( x, lambdaEC, startEC, sizeEC )
    else:
        omegaEC = lambda x: 0.0
        proxEC  = lambda x: np.zeros_like(x)

    # Regularise ISO compartment
    lambdaISO = regularisation.get('lambdaISO')
    startISO = regularisation.get('startISO')
    sizeISO = regularisation.get('sizeISO')
    if (not lambdaISO == None) and (not startISO == None) and (not sizeISO == None):
        omegaISO = lambda x: lambdaISO * sum(x[startISO:sizeISO])
        proxISO  = lambda x:   prox_nnl1( x, lambdaISO, startISO, sizeISO )
    else:
        omegaISO = lambda x: 0.0
        proxISO  = lambda x: np.zeros_like(x)

    omega = lambda x: omegaIC(x) + omegaEC(x) + omegaISO(x)
    prox  = lambda x: proxIC(x)  + proxEC(x)  + proxISO(x)


    return  fista(y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, prox)


## # Splitting methods
def  fista( y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, proximal) :
    """
    Solve the regularised least squares problem

        argmin_x 0.5*||Ax-y||_2^2 + Omega(x)

    with the FISTA algorithm described in [1].

    Notes
    -----
    Author: Matteo Frigo - athena @ INRIA
    Acknowledgment: Rafael Carrillo - lts5 @ EPFL
    References:
        [1] Beck & Teboulle - `A Fast Iterative Shrinkage Thresholding
            Algorithm for Linear Inverse Problems`
    """


    # Initialization
    res = -y.copy()
    xhat = x0.copy()
    x = np.zeros_like(xhat)
    res += A.dot(xhat)
    xhat = proximal( xhat )
    reg_term = omega( xhat )
    prev_obj = 0.5 * np.linalg.norm(res)**2 + reg_term

    told = 1
    beta = 0.9
    prev_x = xhat.copy()
    grad = np.asarray(At.dot(res))
    qfval = prev_obj

    # Step size computation
    L = ( np.linalg.norm( A.dot(grad) ) / np.linalg.norm(grad) )**2
    mu = 1.9 / L

    # Main loop
    if verbose >= 1 :
        print
        print "      |     ||Ax-y||     |  Cost function    Abs error      Rel error    |     Abs x          Rel x"
        print "------|------------------|-----------------------------------------------|------------------------------"
    iter = 1
    while True :
        if verbose >= 1 :
            print "%4d  |" % iter,
            sys.stdout.flush()

        # Smooth step
        x = xhat - mu*grad

        # Non-smooth step
        x = proximal( x )
        reg_term_x = omega( x )

        # Check stepsize
        tmp = x-xhat
        q = qfval + np.real( np.dot(tmp,grad) ) + 0.5/mu * np.linalg.norm(tmp)**2 + reg_term_x
        res = A.dot(x) - y
        res_norm = np.linalg.norm(res)
        curr_obj = 0.5 * res_norm**2 + reg_term_x

        # Backtracking
        while curr_obj > q :
            # Smooth step
            mu = beta*mu
            x = xhat - mu*grad

            # Non-smooth step
            x = proximal( x )
            reg_term_x = omega( x )

            # Check stepsize
            tmp = x-xhat
            q = qfval + np.real( np.dot(tmp,grad) ) + 0.5/mu * np.linalg.norm(tmp)**2 + reg_term_x
            res = A.dot(x) - y
            res_norm = np.linalg.norm(res)
            curr_obj = 0.5 * res_norm**2 + reg_term_x

        # Global stopping criterion
        abs_obj = abs(curr_obj - prev_obj)
        rel_obj = abs_obj / curr_obj
        abs_x   = np.linalg.norm(x - prev_x)
        rel_x   = abs_x / ( np.linalg.norm(x) + eps )
        if verbose >= 1 :
            print "  %13.7e  |  %13.7e  %13.7e  %13.7e  |  %13.7e  %13.7e" % ( res_norm, curr_obj, abs_obj, rel_obj, abs_x, rel_x )

        if abs_obj < eps :
            criterion = "Absolute tolerance on the objective"
            break
        elif rel_obj < tol_fun :
            criterion = "Relative tolerance on the objective"
            break
        elif abs_x < eps :
            criterion = "Absolute tolerance on the unknown"
            break
        elif rel_x < tol_x :
            criterion = "Relative tolerance on the unknown"
            break
        elif iter >= max_iter :
            criterion = "Maximum number of iterations"
            break

        # FISTA update
        t = 0.5 * ( 1 + sqrt(1+4*told**2) )
        xhat = x + (told-1)/t * (x - prev_x)

        # Gradient computation
        res = A.dot(xhat) - y
        xarr = np.asarray(x)

        grad = np.asarray(<np.ndarray[np.float64_t]> At.dot(res))

        # Update variables
        iter += 1
        prev_obj = curr_obj
        prev_x = x.copy()
        told = t
        qfval = 0.5 * np.linalg.norm(res)**2


    if verbose >= 1 :
        print "< Stopping criterion: %s >" % criterion

    opt_details = {}
    opt_details['residual'] = res_norm
    opt_details['cost function'] = curr_obj
    opt_details['abs cost'] = abs_obj
    opt_details['rel cost'] = rel_obj
    opt_details['abs x'] = abs_x
    opt_details['rel x'] = rel_x
    opt_details['iteration'] = iter

    return x, opt_details
