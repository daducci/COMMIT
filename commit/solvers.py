import sys

from math import sqrt

import numpy as np

from dicelib.ui import setup_logger

logger = setup_logger('solvers')
eps = np.finfo(float).eps
list_regularizers = [None, 'lasso', 'group_lasso', 'sparse_group_lasso']
from commit.proximals import non_negativity, omega_group_lasso, prox_group_lasso, soft_thresholding, w_soft_thresholding, omega_sparse_group_lasso, omega_w_sparse_group_lasso
# removed, for now, projection_onto_l2_ball


def init_regularisation(regularisation_params):
    """
    Initialize the regularisation parameters.

    References:
        [1] Jenatton et al. - 'Proximal Methods for Hierarchical Sparse Coding'
    """

    # check if regularisations are in the list
    if regularisation_params['regIC'] not in list_regularizers or regularisation_params['regEC'] not in list_regularizers or regularisation_params['regISO'] not in list_regularizers:
        logger.error('Regularisation not in the list')
    
    startIC  = regularisation_params.get('startIC')
    sizeIC   = regularisation_params.get('sizeIC')
    startEC  = regularisation_params.get('startEC')
    sizeEC   = regularisation_params.get('sizeEC')
    startISO = regularisation_params.get('startISO')
    sizeISO  = regularisation_params.get('sizeISO')

    # create array with al coeff_weights for weighted version of 'lasso'
    all_coeff_weights = np.ones(sizeIC+sizeEC+sizeISO, dtype=np.float64)
    if regularisation_params.get('dictIC_params') is not None and "coeff_weights_kept" in regularisation_params['dictIC_params'].keys():
        all_coeff_weights[startIC:(startIC+sizeIC)] = regularisation_params['dictIC_params']["coeff_weights_kept"]
    # if regularisation_params.get('dictEC_params') is not None and "coeff_weights" in regularisation_params['dictEC_params'].keys():
    #     all_coeff_weights[startEC:(startEC+sizeEC)] = regularisation_params['dictEC_params']["coeff_weights"]
    # if regularisation_params.get('dictISO_params') is not None and "coeff_weights" in regularisation_params['dictISO_params'].keys():
    #     all_coeff_weights[startISO:(startISO+sizeISO)] = regularisation_params['dictISO_params']["coeff_weights"]

    ############################
    # INTRACELLULAR COMPARTMENT#
    ############################

    dictIC_params = regularisation_params.get('dictIC_params')

    if regularisation_params['regIC'] is None:
        omegaIC = lambda x: 0.0
        if regularisation_params.get('nnIC')==True:
            proxIC = lambda x, _: non_negativity(x,startIC,sizeIC)
        else:
            proxIC = lambda x, _: x

    elif regularisation_params['regIC'] == 'lasso':
        lambdaIC = regularisation_params['lambdaIC']
        # check if weights are provided
        if dictIC_params is not None and "coeff_weights_kept" in dictIC_params.keys():
            # w = dictIC_params["coeff_weights"]
            omegaIC = lambda x: lambdaIC * np.linalg.norm(all_coeff_weights[startIC:sizeIC]*x[startIC:sizeIC],1)
            if regularisation_params.get('nnIC'):
                proxIC = lambda x, scaling: non_negativity(w_soft_thresholding(x,all_coeff_weights,scaling*lambdaIC,startIC,sizeIC),startIC,sizeIC)
            else:
                proxIC = lambda x, scaling: w_soft_thresholding(x,all_coeff_weights,scaling*lambdaIC,startIC,sizeIC)
        else:
            omegaIC = lambda x: lambdaIC * np.linalg.norm(x[startIC:sizeIC],1)
            if regularisation_params.get('nnIC'):
                proxIC = lambda x, scaling: non_negativity(soft_thresholding(x,scaling*lambdaIC,startIC,sizeIC),startIC,sizeIC)
            else:
                proxIC = lambda x, scaling: soft_thresholding(x,scaling*lambdaIC,startIC,sizeIC)

    # elif regularisation_params['regIC'] == 'smoothness':
    #     lambdaIC = regularisation_params.get('lambdaIC')
    #     omegaIC = lambda x: lambdaIC * np.linalg.norm(x[startIC:sizeIC])
    #     proxIC  = lambda x: projection_onto_l2_ball(x, lambdaIC, startIC, sizeIC)

    elif regularisation_params['regIC'] == 'group_lasso':
        if not len(dictIC_params['group_idx_kept']) == len(dictIC_params['group_weights']):
            logger.error('Number of groups and weights do not match')

        lambda_group_IC = regularisation_params['lambdaIC']

        # convert to new data structure (needed for faster access)
        N = np.sum([g.size for g in dictIC_params['group_idx_kept']])
        groupIdxIC = np.zeros( (N,), dtype=np.int32 )
        groupSizeIC = np.zeros( (dictIC_params['group_idx_kept'].size,), dtype=np.int32 )
        pos = 0
        for i, g in enumerate(dictIC_params['group_idx_kept']) :
            groupSizeIC[i] = g.size
            groupIdxIC[pos:(pos+g.size)] = g[:]
            pos += g.size

        omegaIC = lambda x: omega_group_lasso( x, groupIdxIC, groupSizeIC, dictIC_params['group_weights'], lambda_group_IC )
        if regularisation_params.get('nnIC'):
            proxIC = lambda x, scaling: non_negativity(prox_group_lasso(x,groupIdxIC,groupSizeIC,dictIC_params['group_weights'],scaling*lambda_group_IC),startIC,sizeIC)
        else:
            proxIC = lambda x, scaling: prox_group_lasso(x,groupIdxIC,groupSizeIC,dictIC_params['group_weights'],scaling*lambda_group_IC)
  
    elif regularisation_params['regIC'] == 'sparse_group_lasso':
        if not len(dictIC_params['group_idx_kept']) == len(dictIC_params['group_weights']):
            logger.error('Number of groups and weights do not match')

        lambdaIC = regularisation_params['lambdaIC'][0]
        lambda_group_IC = regularisation_params['lambdaIC'][1]

        # convert to new data structure (needed for faster access)
        N = np.sum([g.size for g in dictIC_params['group_idx_kept']])
        groupIdxIC  = np.zeros( (N,), dtype=np.int32 )
        groupSizeIC = np.zeros( (dictIC_params['group_idx_kept'].size,), dtype=np.int32 )
        pos = 0
        for i, g in enumerate(dictIC_params['group_idx_kept']) :
            groupSizeIC[i] = g.size
            groupIdxIC[pos:(pos+g.size)] = g[:]
            pos += g.size

        if "coeff_weights_kept" in dictIC_params.keys():
            # w = dictIC_params["coeff_weights"]
            omegaIC = lambda x: omega_w_sparse_group_lasso( x, all_coeff_weights, groupIdxIC, groupSizeIC, dictIC_params['group_weights'], lambdaIC, lambda_group_IC)

            if regularisation_params.get('nnIC'):
                proxIC = lambda x, scaling: non_negativity(prox_group_lasso(w_soft_thresholding(x,all_coeff_weights,scaling*lambdaIC,startIC,sizeIC),groupIdxIC,groupSizeIC,dictIC_params['group_weights'],scaling*lambda_group_IC), startIC,sizeIC)
            else:
                proxIC = lambda x, scaling: prox_group_lasso(w_soft_thresholding(x,all_coeff_weights,scaling*lambdaIC,startIC,sizeIC),groupIdxIC,groupSizeIC,dictIC_params['group_weights'],scaling*lambda_group_IC)
        else:
            omegaIC = lambda x: omega_sparse_group_lasso( x, groupIdxIC, groupSizeIC, dictIC_params['group_weights'], lambdaIC, lambda_group_IC)

            if regularisation_params.get('nnIC'):
                proxIC = lambda x, scaling: non_negativity(prox_group_lasso(soft_thresholding(x,scaling*lambdaIC,startIC,sizeIC),groupIdxIC,groupSizeIC,dictIC_params['group_weights'],scaling*lambda_group_IC), startIC,sizeIC)
            else:
                proxIC = lambda x, scaling: prox_group_lasso(soft_thresholding(x,scaling*lambdaIC,startIC,sizeIC),groupIdxIC,groupSizeIC,dictIC_params['group_weights'],scaling*lambda_group_IC)


    ###########################
    # EXTRCELLULAR COMPARTMENT#
    ###########################

    dictEC_params = regularisation_params.get('dictEC_params')

    if regularisation_params['regEC'] is None:
        omegaEC = lambda x: 0.0
        if regularisation_params.get('nnEC')==True:
            proxEC = lambda x, _: non_negativity(x,startEC,sizeEC)
        else:
            proxEC = lambda x, _: x

    elif regularisation_params['regEC'] == 'lasso':
        lambdaEC = regularisation_params.get('lambdaEC')
        # check if weights are provided
        # if dictEC_params is not None and "coeff_weights" in dictEC_params.keys():
        #     omegaEC = lambda x: lambdaEC * np.linalg.norm(all_coeff_weights[startEC:sizeEC]*x[startEC:sizeEC],1)
        #     if regularisation_params.get('nnEC'):
        #         proxEC = lambda x, scaling: non_negativity(w_soft_thresholding(x,all_coeff_weights,scaling*lambdaEC,startEC,sizeEC),startEC,sizeEC)
        #     else:
        #         proxEC = lambda x, scaling: w_soft_thresholding(x,all_coeff_weights,scaling*lambdaEC,startEC,sizeEC)
        # else:
        omegaEC = lambda x: lambdaEC * np.linalg.norm(x[startEC:(startEC+sizeEC)],1)
        if regularisation_params.get('nnEC'):
            proxEC = lambda x, scaling: non_negativity(soft_thresholding(x,scaling*lambdaEC,startEC,sizeEC),startEC,sizeEC)
        else:
            proxEC = lambda x, scaling: soft_thresholding(x,scaling*lambdaEC,startEC,sizeEC)

    # elif regularisation_params['regIC'] == 'smoothness':
    #     lambdaEC = regularisation_params.get('lambdaEC')
    #     omegaEC = lambda x: lambdaEC * np.linalg.norm(x[startEC:(startEC+sizeEC)])
    #     proxEC  = lambda x: projection_onto_l2_ball(x, lambdaEC, startEC, sizeEC)


    ########################
    # ISOTROPIC COMPARTMENT#
    ########################

    dictISO_params = regularisation_params.get('dictISO_params')

    if regularisation_params['regISO'] is None:
        omegaISO = lambda x: 0.0
        if regularisation_params.get('nnISO')==True:
            proxISO = lambda x, _: non_negativity(x,startISO,sizeISO)
        else:
            proxISO  = lambda x, _: x

    elif regularisation_params['regISO'] == 'lasso':
        lambdaISO = regularisation_params.get('lambdaISO')
        # check if weights are provided
        # if dictISO_params is not None and "coeff_weights" in dictISO_params.keys():
        #     omegaISO = lambda x: lambdaISO * np.linalg.norm(all_coeff_weights[startISO:sizeISO]*x[startISO:sizeISO],1)
        #     if regularisation_params.get('nnISO'):
        #         proxISO = lambda x, scaling: non_negativity(w_soft_thresholding(x,all_coeff_weights,scaling*lambdaISO,startISO,sizeISO),startISO,sizeISO)
        #     else:
        #         proxISO = lambda x, scaling: w_soft_thresholding(x,all_coeff_weights,scaling*lambdaISO,startISO,sizeISO)
        # else:
        omegaISO = lambda x: lambdaISO * np.linalg.norm(x[startISO:(startISO+sizeISO)],1)
        if regularisation_params.get('nnISO'):
            proxISO  = lambda x, scaling: non_negativity(soft_thresholding(x,scaling*lambdaISO,startISO,sizeISO),startISO,sizeISO)
        else:
            proxISO  = lambda x, scaling: soft_thresholding(x,scaling*lambdaISO,startISO,sizeISO)

    # elif regularisation_params['regISO'] == 'group_lasso':
    #     lambdaISO = regularisation_params.get('lambdaISO')
    #     omegaISO = lambda x: lambdaISO * np.linalg.norm(x[startISO:(startISO+sizeISO)])
    #     proxISO  = lambda x: projection_onto_l2_ball(x, lambdaISO, startISO, sizeISO)

    omega = lambda x: omegaIC(x) + omegaEC(x) + omegaISO(x)
    prox = lambda x, scaling: proxIC(proxEC(proxISO(x, scaling), scaling), scaling)

    regularisation_params["omega"] = omega
    regularisation_params["prox"] = prox

    return regularisation_params


def evaluate_model(y, A, x, regularisation=None):
    if regularisation is None:
        omega = lambda x: 0.0
    else:
        omega, _ = init_regularisation(regularisation)

    return 0.5*np.linalg.norm(A.dot(x)-y)**2 + omega(x)


def solve(y, A, At, tol_fun=1e-4, tol_x=1e-6, max_iter=1000, verbose=True, x0  =None, regularisation=None, confidence_array=None):
    """
    Solve the regularised least squares problem

        argmin_x 0.5*|| sqrt(W) ( Ax-y ) ||_2^2 + Omega(x)

    with the Omega described by 'regularisation' and W is the confidence_array

    Check the documentation of commit.solvers.init_regularisation to see how to
    solve a specific problem.
    """
    if regularisation is None:
        omega = lambda x: 0.0
        prox  = lambda x, _: x
    else:
        omega, prox = regularisation['omega'], regularisation['prox']

    if x0 is None:
        x0 = np.zeros(A.shape[1])

    if confidence_array is not None:
        confidence_array = np.sqrt(confidence_array)

    return fista(y, A, At, omega, prox, confidence_array, tol_fun, tol_x, max_iter, verbose, x0)


def fista(y, A, At, omega, prox, sqrt_W=None, tol_fun=1e-4, tol_x=1e-6, max_iter=1000, verbose=False, x0=None):
    """
    Solve the regularised least squares problem

        argmin_x 0.5*|| sqrt(W) ( Ax-y ) ||_2^2 + Omega(x)

    with the FISTA algorithm described in [1].

    The penalty term and its proximal operator must be defined in such a way
    that they already contain the regularisation parameter.

    References:
        [1] Beck & Teboulle - `A Fast Iterative Shrinkage Thresholding
            Algorithm for Linear Inverse Problems`
    """

    # Initialization
    if x0 is None:
        x0 = np.zeros(A.shape[1])
    xhat = x0.copy()
    x = np.zeros_like(xhat)

    if sqrt_W is not None:
        res = sqrt_W * (A.dot(xhat) - y)
        grad = np.asarray(At.dot(sqrt_W * res))
    else:
        res = A.dot(xhat) - y
        grad = np.asarray(At.dot(res))
    
    print("grad", grad[grad!=0][:50])
    print("res", res[res!=0][:50])

    prox( xhat, 1.0 )
    reg_term = omega( xhat )
    prev_obj = 0.5 * np.linalg.norm(res)**2 + reg_term

    told = 1
    beta = 0.9
    prev_x = xhat.copy()
    qfval = prev_obj


    # Step size computation
    if sqrt_W is not None:
        L = ( np.linalg.norm( sqrt_W * A.dot(grad) ) / np.linalg.norm(grad) )**2
    else:
        L = ( np.linalg.norm( A.dot(grad) ) / np.linalg.norm(grad) )**2
    step_size = 1.9 / L
    # Main loop
    if verbose==4 :
        print()
        print( "      |  1/2||Ax-y||^2      Omega      |  Cost function    Abs error      Rel error    |      Abs x          Rel x    " )
        print( "------|--------------------------------|-----------------------------------------------|------------------------------" )
    iter = 1
    while True :
        if verbose==4 :
            print( "%4d  |" % iter, end="" )
            sys.stdout.flush()

        # Smooth step
        x = xhat - step_size*grad

        # Non-smooth step
        prox( x, step_size )
        reg_term_x = omega( x )
        # Check stepsize
        tmp = x-xhat
        q = qfval + np.real( np.dot(tmp,grad) ) + 0.5/step_size * np.linalg.norm(tmp)**2 + reg_term_x
        if sqrt_W is not None:
            res = sqrt_W * ( A.dot(x) - y )
        else:
            res = A.dot(x) - y
        res_norm = np.linalg.norm(res)
        curr_obj = 0.5 * res_norm**2 + reg_term_x

        # Backtracking
        while curr_obj > q :
            # Smooth step
            step_size = beta*step_size
            x = xhat - step_size*grad

            # Non-smooth step
            prox( x, step_size )
            reg_term_x = omega( x )

            # Check stepsize
            tmp = x-xhat
            q = qfval + np.real( np.dot(tmp,grad) ) + 0.5/step_size * np.linalg.norm(tmp)**2 + reg_term_x
            if sqrt_W is not None:
                res = sqrt_W * ( A.dot(x) - y )
            else:
                res = A.dot(x) - y
            res_norm = np.linalg.norm(res)
            curr_obj = 0.5 * res_norm**2 + reg_term_x

        # Global stopping criterion
        abs_obj = abs(curr_obj - prev_obj)
        rel_obj = abs_obj / curr_obj
        abs_x   = np.linalg.norm(x - prev_x)
        rel_x   = abs_x / ( np.linalg.norm(x) + eps )
        if verbose==4 :
            print( "  %13.7e  %13.7e  |  %13.7e  %13.7e  %13.7e  |  %13.7e  %13.7e" % ( 0.5 * res_norm**2, reg_term_x, curr_obj, abs_obj, rel_obj, abs_x, rel_x ), flush=True )

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
        if sqrt_W is not None:
            res = sqrt_W * ( A.dot(xhat) - y )
            grad = np.asarray(At.dot(sqrt_W * res))
        else:
            res = A.dot(xhat) - y
            grad = np.asarray(At.dot(res))

        # Update variables
        iter += 1
        prev_obj = curr_obj
        prev_x = x.copy()
        told = t
        qfval = 0.5 * np.linalg.norm(res)**2

    if verbose==4 :
        print( "< Stopping criterion: %s >" % criterion, flush=True )

    opt_details = {}
    opt_details['residual'] = 0.5*res_norm**2
    opt_details['regterm'] = reg_term_x
    opt_details['cost_function'] = curr_obj
    opt_details['abs_cost'] = abs_obj
    opt_details['rel_cost'] = rel_obj
    opt_details['abs_x'] = abs_x
    opt_details['rel_x'] = rel_x
    opt_details['iterations'] = iter
    opt_details['stopping_criterion'] = criterion

    return x, opt_details
