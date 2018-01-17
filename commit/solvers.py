import numpy as np
from math import sqrt
import sys
eps = np.finfo(float).eps

from commit.proximals import (non_negativity,
                             omega_group_sparsity,
                             prox_group_sparsity,
                             soft_thresholding,
                             projection_onto_l2_ball)
group_sparsity = -1
non_negative = 0
norm1 = 1
norm2 = 2
list_regnorms = [group_sparsity, non_negative, norm1, norm2]
list_group_sparsity_norms = [2, np.inf]

def init_regularisation(commit_evaluation,
                        regnorms = (non_negative, non_negative, non_negative),
                        structureIC = None, weightsIC = None, group_norm = 2,
                        group_is_ordered = False,
                        lambdas = (.0,.0,.0) ):

    regularisation = {}

    regularisation['startIC']  = 0
    regularisation['sizeIC']   = int( commit_evaluation.DICTIONARY['IC']['nF'])#*len(commit_evaluation.model.ICVFs) )
    regularisation['startEC']  = int( regularisation['sizeIC'] )
    regularisation['sizeEC']   = int( commit_evaluation.DICTIONARY['EC']['nE'] )
    regularisation['startISO'] = int( regularisation['sizeIC'] + regularisation['sizeEC'] )
    regularisation['sizeISO']  = int( commit_evaluation.DICTIONARY['nV'])#*len(commit_evaluation.model.d_ISOs) )

    regularisation['normIC']  = regnorms[0]
    regularisation['normEC']  = regnorms[1]
    regularisation['normISO'] = regnorms[2]

    regularisation['lambdaIC']  = lambdas[0]
    regularisation['lambdaEC']  = lambdas[1]
    regularisation['lambdaISO'] = lambdas[2]

    # Solver-specific fields
    regularisation['group_is_ordered'] = group_is_ordered  # This option will be deprecated in future release
    regularisation['structureIC']      = structureIC
    regularisation['weightsIC']        = weightsIC
    regularisation['group_norm']       = group_norm

    return regularisation


def regularisation2omegaprox(regularisation):
    lambdaIC  = float(regularisation.get('lambdaIC'))
    lambdaEC  = float(regularisation.get('lambdaEC'))
    lambdaISO = float(regularisation.get('lambdaISO'))
    if lambdaIC < 0.0 or lambdaEC < 0.0 or lambdaISO < 0.0:
        raise ValueError('Negative regularisation parameters are not allowed')

    normIC  = regularisation.get('normIC')
    normEC  = regularisation.get('normEC')
    normISO = regularisation.get('normISO')
    if not normIC in list_regnorms:
        raise ValueError('normIC must be one of commit.solvers.{group_sparsity,non_negative,norm1,norm2}')
    if not normEC in list_regnorms:
        raise ValueError('normEC must be one of commit.solvers.{group_sparsity,non_negative,norm1,norm2}')
    if not normISO in list_regnorms:
        raise ValueError('normISO must be one of commit.solvers.{group_sparsity,non_negative,norm1,norm2}')

    ## NNLS case
    if (lambdaIC == 0.0 and lambdaEC == 0.0 and lambdaISO == 0.0) or (normIC == non_negative and normEC == non_negative and normISO == non_negative):
        omega = lambda x: 0.0
        prox  = lambda x: non_negativity(x, 0, len(x))
        return omega, prox

    ## All other cases
    #
    # Intracellular Compartment
    startIC = regularisation.get('startIC')
    sizeIC  = regularisation.get('sizeIC')
    if lambdaIC == 0.0:
        omegaIC = lambda x: 0.0
        proxIC  = lambda x: non_negativity(x, startIC, sizeIC)
    elif normIC == norm2:
        omegaIC = lambda x: lambdaIC * np.linalg.norm(x[startIC:sizeIC])
        proxIC  = lambda x: projection_onto_l2_ball(x, lambdaIC, startIC, sizeIC)
    elif normIC == norm1:
        omegaIC = lambda x: lambdaIC * sum( x[startIC:sizeIC] )
        proxIC  = lambda x: soft_thresholding(x, lambdaIC, startIC, sizeIC)
    elif normIC == non_negative:
        omegaIC = lambda x: 0.0
        proxIC  = lambda x: non_negativity(x, startIC, sizeIC)
    elif normIC == group_sparsity:
        weightsIC   = regularisation.get('weightsIC')
        structureIC = regularisation.get('structureIC')
        if regularisation.get('group_is_ordered'): # This option will be deprecated in future release
            raise DeprecationWarning('The ordered group structure will be deprecated. Consider using the structureIC field for defining the group structure.')
            bundles = np.cumsum(np.insert(sizeIC,0,0))
            structureIC = np.array([range(bundles[k],bundles[k+1]) for k in range(0,len(bundles)-1)])
            regularisation['structureIC'] = structureIC
            del bundles
        if not len(structureIC) == len(weightsIC):
            raise ValueError('Number of groups and weights do not coincide.')
        group_norm = regularisation.get('group_norm')
        if not group_norm in list_group_sparsity_norms:
            raise ValueError('Wrong norm in the structured sparsity term. Choose between %s.' % str(list_group_sparsity_norms))

        omegaIC = lambda x: omega_group_sparsity( x, structureIC, weightsIC, lambdaIC, normIC )
        proxIC  = lambda x:  prox_group_sparsity( x, structureIC, weightsIC, lambdaIC, normIC )
    else:
        raise ValueError('Type of regularisation for IC compartment not recognized.')


    # Extracellular Compartment
    startEC = regularisation.get('startEC')
    sizeEC  = regularisation.get('sizeEC')
    if lambdaEC == 0.0:
        omegaEC = lambda x: 0.0
        proxEC  = lambda x: non_negativity(x, startEC, sizeEC)
    elif normEC == norm2:
        omegaEC = lambda x: lambdaEC * np.linalg.norm(x[startEC:sizeEC])
        proxEC  = lambda x: projection_onto_l2_ball(x, lambdaEC, startEC, sizeEC)
    elif normEC == norm1:
        omegaEC = lambda x: lambdaEC * sum( x[startEC:sizeEC] )
        proxEC  = lambda x: soft_thresholding(x, lambdaEC, startEC, sizeEC)
    elif normEC == non_negative:
        omegaEC = lambda x: 0.0
        proxEC  = lambda x: non_negativity(x, startEC, sizeEC)
    else:
        raise ValueError('Type of regularisation for EC compartment not recognized.')

    # Isotropic Compartment
    startISO = regularisation.get('startISO')
    sizeISO  = regularisation.get('sizeISO')
    if lambdaISO == 0.0:
        omegaISO = lambda x: 0.0
        proxISO  = lambda x: non_negativity(x, startISO, sizeISO)
    elif normISO == norm2:
        omegaISO = lambda x: lambdaISO * np.linalg.norm(x[startISO:sizeISO])
        proxISO  = lambda x: projection_onto_l2_ball(x, lambdaISO, startISO, sizeISO)
    elif normISO == norm1:
        omegaISO = lambda x: lambdaISO * sum( x[startISO:sizeISO] )
        proxISO  = lambda x: soft_thresholding(x, lambdaISO, startISO, sizeISO)
    elif normISO == non_negative:
        omegaISO = lambda x: 0.0
        proxISO  = lambda x: non_negativity(x, startISO, sizeISO)
    else:
        raise ValueError('Type of regularisation for ISO compartment not recognized.')

    omega = lambda x: omegaIC(x) + omegaEC(x) + omegaISO(x)
    prox = lambda x: proxIC(proxEC(proxISO(x)))

    return omega, prox

def solver(y, A, At, tol_fun = 1e-4, tol_x = 1e-6, max_iter = 1000, verbose = 1, x0 = None, regularisation = None):
    if regularisation is None:
        omega = lambda x: 0.0
        prox  = lambda x: non_negativity(x, 0, len(x))
    else:
        omega, prox = regularisation2omegaprox(regularisation)

    if x0 is None:
        x0 = np.ones(A.shape[1])

    return fista( y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, prox)

def fista( y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, proximal) :
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

        grad = np.asarray(At.dot(res))

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
    opt_details['cost_function'] = curr_obj
    opt_details['abs_cost'] = abs_obj
    opt_details['rel_cost'] = rel_obj
    opt_details['abs_x'] = abs_x
    opt_details['rel _x'] = rel_x
    opt_details['iterations'] = iter
    opt_details['stopping_criterion'] = criterion

    return x, opt_details
