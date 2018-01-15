import numpy as np
from math import sqrt
import sys
eps = np.finfo(float).eps

   ##   #    # #    # # #      #   ##   #####  #   #
  #  #  #    #  #  #  # #      #  #  #  #    #  # #
 #    # #    #   ##   # #      # #    # #    #   #
 ###### #    #   ##   # #      # ###### #####    #
 #    # #    #  #  #  # #      # #    # #   #    #
 #    #  ####  #    # # ###### # #    # #    #   #
from commit.proximals import non_negativity,
                             omega_hierarchical,
                             prox_hierarchical,
                             soft_thresholding,
                             projection_onto_l2_ball
hierarchical = -1
non_negative = 0
norm1 = 1
norm2 = 2
list_regnorms = [hierarchical, non_negative, norm1, norm2]

def init_regularisation(commit_evaluation = None,
                        regnorms = (non_negative, non_negative, non_negative),
                        structureIC = None, weightsIC = None,
                        group_is_ordered = False
                        lambdas = (.0,.0,.0) ):
    if commit_evaluation is None:
        raise ValueError('commit.Evaluation object expected but not found')

    regularisation = {}

    regularisation['startIC']  = 0
    regularisation['sizeIC']   = int( commit_evaluation.DICTIONARY['IC']['nF']*len(commit_evaluation.model.ICVFs) )
    regularisation['startEC']  = int( regularisation['sizeIC'] )
    regularisation['sizeEC']   = int( commit_evaluation.DICTIONARY['EC']['nE'] )
    regularisation['startISO'] = int( regularisation['sizeIC'] + regularisation['sizeEC'] )
    regularisation['sizeISO']  = int( commit_evaluation.DICTIONARY['nV']*len(commit_evaluation.model.d_ISOs) )

    regularisation['normIC']  = regIC
    regularisation['normEC']  = regEC
    regularisation['normISO'] = regISO

    regularisation['lambdaIC']  = lambdas[0]
    regularisation['lambdaEC']  = lambdas[1]
    regularisation['lambdaISO'] = lambdas[2]

    # Solver-specific fields
    regularisation['group_is_ordered'] = group_is_ordered
    regularisation['structureIC']      = structureIC
    regularisation['weightsIC']        = weightsIC

    return regularisation


def regularisation2omegaprox(regularisation):
    lambdaIC  = regularisation.get('lambdaIC')
    lambdaEC  = regularisation.get('lambdaEC')
    lambdaISO = regularisation.get('lambdaISO')

    normIC  = regularisation.get('normIC')
    normEC  = regularisation.get('normEC')
    normISO = regularisation.get('normISO')

    ## NNLS case
    if (lambdaIC == 0.0 and lambdaEC == 0.0 and lambdaISO = 0.0) or (normIC == 0 and normEC == 0 and normISO == 0):
        omega = lambda x: 0.0
        prox  = lambda x: non_negativity(x)
        return omega, prox

    ## All other cases
    #
    # Intracellular Compartment
    startIC = regularisation.get('startIC')
    sizeIC  = regularisation.get('sizeIC')
    if lambdaIC > 0.0:
        if normIC == norm2:
            omegaIC = lambda x: lambdaIC * np.linalg.norm(x[startIC:sizeIC])
            proxIC  = lambda x: projection_onto_l2_ball(x, lambdaIC, startIC, sizeIC)
        elif normIC == norm1:
            omegaIC = lambda x: lambdaIC * sum( x[startIC:sizeIC] )
            proxIC  = lambda x: soft_thresholding(x, lambdaIC, startIC, sizeIC)
        elif normIC == non_negative:
            omegaIC = lambda x: 0.0
            proxIC  = lambda x: non_negativity(x, startIC, sizeIC)
        elif normIC == hierarchical:
            weightsIC   = regularisation.get('weightsIC')
            structureIC = regularisation.get('structureIC')
            if regularisation.get('group_is_ordered'): # make the list-like structure suitable for the solver
                raise DeprecationWarning('The ordered group structure will be deprecated. Consider using the structureIC field for defining the group structure.')
                bundles = np.cumsum(np.insert(sizeIC,0,0))
                structureIC = np.array([range(bundles[k],bundles[k+1]) for k in range(0,len(bundles)-1)])
                regularisation['structureIC'] = structureIC
                del bundles

            if len(structureIC) != len(weightsIC):
                raise ValueError('Number of groups and weights do not coincide.')

            omegaIC = lambda x: omega_hierarchical( x, structureIC, weightsIC, lambdaIC, normIC )
            proxIC  = lambda x: prox_hierarchical( x, structureIC, weightsIC, lambdaIC, normIC )
        else:
            raise ValueError('Type of regularisation for IC compartment not recognized.')
    else:
        omegaIC = lambda x: 0.0
        proxIC  = lambda x: non_negativity(x, startIC, sizeIC)

    # Extracellular Compartment
    if lambdaEC > 0.0:
        if normEC == norm2:
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
    else:
        omegaEC = lambda x: 0.0
        proxEC  = lambda x: non_negativity(x, startEC, sizeEC)

    # Isotropic Compartment
    if lambdaISO == 0.0:
        normISO = non_negative
        omegaISO = lambda x: 0.0
        proxISO  = lambda x: non_negativity(x, startISO, sizeISO)

    if normISO not in list_regnorms:
        raise ValueError('Type of regularisation for ISO compartment not recognized.')
    elif lambdaISO > 0.0:
        if normISO == norm2:
            omegaISO = lambda x: lambdaISO * np.linalg.norm(x[startISO:sizeISO])
            proxISO  = lambda x: projection_onto_l2_ball(x, lambdaISO, startISO, sizeISO)
        elif normISO == norm1:
            omegaISO = lambda x: lambdaISO * sum( x[startISO:sizeISO] )
            proxISO  = lambda x: soft_thresholding(x, lambdaISO, startISO, sizeISO)
        elif normISO == non_negative:
            omegaISO = lambda x: 0.0
            proxISO  = lambda x: non_negativity(x, startISO, sizeISO)

    omega = lambda x: omegaIC(x) + omegaEC(x) + omegaISO(x)
    prox  = lambda x: proxIC(x)  + proxEC(x)  + proxISO(x)

    return omega, prox

## # Solver wrapper
  #####  ####### #       #     # ####### ######
 #     # #     # #       #     # #       #     #
 #       #     # #       #     # #       #     #
  #####  #     # #       #     # #####   ######
       # #     # #        #   #  #       #   #
 #     # #     # #         # #   #       #    #
  #####  ####### #######    #    ####### #     #
def solver(y, A, At, tol_fun = 1e-4, tol_x = 1e-6, max_iter = 1000, verbose = 1, x0 = None, regularisation = None):
    if regularisation==None:
        omega = lambda x: 0.0
        prox  = lambda x: non_negativity(x)
    else:
        omega, prox = regularisation2omegaprox(regularisation)

    if x0 = None:
        x0 = np.zeros(A.shape[1])

    return fista( y, A, At, tol_fun, tol_x, max_iter, verbose, x0, omega, proximal)

## # Splitting methods
 ###### #  ####  #####   ##
 #      # #        #    #  #
 #####  #  ####    #   #    #
 #      #      #   #   ######
 #      # #    #   #   #    #
 #      #  ####    #   #    #
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
    opt_details['cost function'] = curr_obj
    opt_details['abs cost'] = abs_obj
    opt_details['rel cost'] = rel_obj
    opt_details['abs x'] = abs_x
    opt_details['rel x'] = rel_x
    opt_details['iterations'] = iter
    opt_details['stopping criterion'] = criterion

    return x, opt_details
