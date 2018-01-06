import numpy as np
import sys

eps = np.finfo(float).eps


def nnls( y, A, At = None, tol_fun = 1e-4, tol_x = 1e-9, max_iter = 100, verbose = 1, x0 = None ) :
    """Solve non negative least squares problem of the following form:

       min 0.5*||y-A x||_2^2 s.t. x >= 0

    The problem is solved using the forward-backward algorithm with FISTA-like acceleration.

    Parameters
    ----------
    y : 1-d array of doubles.
        Contains the measurements.

    A : matrix or class exposing the .dot() method.
        It is the forward measurement operator.

    At : matrix or class exposing the .dot() method.
        It is the corresponding adjoint operator (default: computed as A.T).

    tol_fun : double, optional (default: 1e-4).
        Minimum relative change of the objective value. The algorithm stops if:
               | f(x(t)) - f(x(t-1)) | / f(x(t)) < tol_fun,
        where x(t) is the estimate of the solution at iteration t.

    tol_x : double, optional (default: 1e-9).
        Minimum relative change of the solution x. The algorithm stops if:
               || x(t) - x(t-1) || / || x(t) || < tol_x,
        where x(t) is the estimate of the solution at iteration t.

    max_iter : integer, optional (default: 100).
        Maximum number of iterations.

    verbose : integer, optional (default: 1).
        0 no log, 1 print each iteration results.

    x0 : 1-d array of double, optional (default: automatically computed).
        Initial solution.

    Returns
    -------
    x : 1-d array of doubles.
        Best solution in the least-squares sense.

    Notes
    -----
    Author: Rafael Carrillo
    E-mail: rafael.carrillo@epfl.ch
    """
    # Initialization
    CONFIG = {}
    if At is None :
        At = A.T

    if x0 is not None :
        xhat = x0
        res = A.dot(xhat) - y
    else :
        xhat = np.zeros( A.shape[1], dtype=np.float64 )
        res = -y
    grad = At.dot(res)
    prev_obj = 0.5 * np.linalg.norm(res)**2
    iter = 1
    told = 1
    prev_x = xhat
    beta = 0.9
    qfval = prev_obj

    # Step size computation
    L = np.linalg.norm( A.dot(grad) )**2 / np.linalg.norm(grad)**2
    mu = 1.9 / L

    # Main loop
    if verbose >= 1 :
        print "      |     ||Ax-y||     |  Cost function    Abs error      Rel error    |     Abs x          Rel x"
        print "------|------------------|-----------------------------------------------|------------------------------"

    while True :
        if verbose >= 1 :
            print "%4d  |" % iter,
            sys.stdout.flush()

        # Gradient descend step
        x = xhat - mu*grad

        # Projection onto the positive orthant
        x = np.real( x )
        x[ x<0 ] = 0

        # Stepsize check
        tmp = x-xhat
        q = qfval + np.real( np.dot(tmp,grad) ) + 0.5/mu * np.linalg.norm(tmp)**2
        res = A.dot(x) - y
        curr_obj = 0.5 * np.linalg.norm(res)**2

        # Backtracking
        while curr_obj > q :
            # Gradient descend step
            mu = beta*mu
            x = xhat - mu*grad

            # Projection onto the positive orthant
            x = np.real( x )
            x[ x<0 ] = 0

            # New stepsize check
            tmp = x-xhat
            q = qfval + np.real( np.dot(tmp,grad) ) + 0.5/mu * np.linalg.norm(tmp)**2
            res = A.dot(x) - y
            curr_obj = 0.5 * np.linalg.norm(res)**2

        # Global stopping criterion
        abs_obj = np.abs(curr_obj - prev_obj)
        rel_obj = abs_obj / curr_obj
        abs_x   = np.linalg.norm(x - prev_x)
        rel_x   = abs_x / ( np.linalg.norm(x) + eps )
        if verbose >= 1 :
            print "  %13.7e  |  %13.7e  %13.7e  %13.7e  |  %13.7e  %13.7e" % ( np.sqrt(2.0*curr_obj), curr_obj, abs_obj, rel_obj, abs_x, rel_x )

        if abs_obj < eps :
            criterion = "ABS_OBJ"
            break
        elif rel_obj < tol_fun :
            criterion = "REL_OBJ"
            break
        elif abs_x < eps :
            criterion = "ABS_X"
            break
        elif rel_x < tol_x :
            criterion = "REL_X"
            break
        elif iter >= max_iter :
            criterion = "MAX_IT"
            break

        # FISTA update
        t = 0.5 * ( 1 + np.sqrt(1+4*told**2) )
        xhat = x + (told-1)/t * (x - prev_x)

        # Gradient computation
        res = A.dot(xhat) - y
        grad = At.dot(res)

        # Update variables
        iter += 1
        prev_obj = curr_obj
        prev_x = x
        told = t
        qfval = 0.5 * np.linalg.norm(res)**2

    if verbose >= 1 :
        print "< Stopping criterion: %s >" % criterion

    CONFIG['||Ax-y||'] = round(np.sqrt(2.0*curr_obj), 5)
    CONFIG['Cost function'] = round(curr_obj, 5)
    CONFIG['Abs error'] = round(abs_obj, 5)
    CONFIG['Rel error'] = round(rel_obj, 5)
    CONFIG['Abs x'] = round(abs_x, 5)
    CONFIG['Rel x'] = round(rel_x, 5)
    CONFIG['iteration'] = iter

    return x
