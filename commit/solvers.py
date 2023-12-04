"""
Author: Matteo Frigo - lts5 @ EPFL and Dep. of CS @ Univ. of Verona

This structure is based on the previous work of Rafael Carrillo and was
supported by the LTS5 laboratory at EPFL, Lausanne.
"""
from __future__ import print_function
import numpy as np
from math import sqrt
import sys
import warnings
eps = np.finfo(float).eps

 # removed, for now, projection_onto_l2_ball
from commit.proximals import non_negativity, omega_group_sparsity, prox_group_sparsity, soft_thresholding
list_regnorms = [None, 'l1', 'group_sparsity'] # removed 'l2' as never tested
list_group_sparsity_norms = ['l2'] # removed 'np.inf' because of issue #54


def init_regularisation(
    commit_evaluation,
    regnorms = (None, None, None),
    is_nonnegative = (True, True, True),
    structureIC = None, weightsIC = None, group_norm = 2,
    lambdas = (0.0, 0.0, 0.0)
):
    """
    Initialise the data structure that defines Omega in:
        argmin_x 0.5*||Ax-y||_2^2 + Omega(x)

    Input
    -----
    commit_evaluation - commit.Evaluation object :
        'dictionary' and 'model' have to be loaded beforehand.

    regnorms - tuple :
        sets the penalty term to be used for each compartment:
            regnorms[0] corresponds to the Intracellular compartment
            regnorms[1] corresponds to the Extracellular compartment
            regnorms[2] corresponds to the Isotropic compartment
        Each regnorms[k] must be one of: {None, 'l1', 'group_sparsity'}:
            'l1' penalises with the 1-norm of the coefficients
                corresponding to the compartment.
            'group_sparsity' penalises according to the following formulation (see [1]):
                $\Omega(x) = \lambda \sum_{g\in G} w_g |x_g|
                Considers both the non-overlapping and the hierarchical formulations.
                NB: this option is allowed only in the IC compartment.
        Default = (None, None, None).

    lambdas - tuple :
        regularisation parameter for each compartment.
        The lambdas correspond to the ones described in the mathematical
        formulation of the regularisation term
        $\Omega(x) = lambdas[0]*regnorm[0](x) + lambdas[1]*regnorm[1](x) + lambdas[2]*regnorm[2](x)$
        Default = (0.0, 0.0, 0.0).

    is_nonnegative - tuple :
        impose a non negativity constraint for each compartment:
            is_nonnegative[0] corresponds to the Intracellular compartment
            is_nonnegative[1] corresponds to the Extracellular compartment
            is_nonnegative[2] corresponds to the Isotropic compartment
        Default = (True, True, True).

    structureIC - np.array(list(list), dtype=np.object_) :
        group structure for the IC compartment.
            This field is necessary only if regterm[0]='group_sparsity'.
            Example:
                structureIC = np.array([[0,2,5],[1,3,4],[0,1,2,3,4,5],[6]], dtype=np.object_)

                that is equivalent to
                            [0,1,2,3,4,5]        [6]
                              /       \
                        [0,2,5]       [1,3,4]
                which has two non overlapping groups, one of which is the union
                of two other non-overlapping groups.

    weightsIC - np.array(np.float64) :
        this defines the weights associated to each group of structure IC.

    group_norm - number :
        norm type for the 'group_sparsity' penalisation of the IC compartment.
            To be chosen among { 'l2' }.
            Default: group_norm = 'l2'.

    References:
        [1] Jenatton et al. - 'Proximal Methods for Hierarchical Sparse Coding'
    """
    regularisation = {}

    regularisation['startIC']  = 0
    regularisation['sizeIC']   = int( commit_evaluation.DICTIONARY['IC']['nF'] * commit_evaluation.KERNELS['wmr'].shape[0])
    regularisation['startEC']  = int( regularisation['sizeIC'] )
    regularisation['sizeEC']   = int( commit_evaluation.DICTIONARY['EC']['nE'] * commit_evaluation.KERNELS['wmh'].shape[0])
    regularisation['startISO'] = int( regularisation['sizeIC'] + regularisation['sizeEC'] )
    regularisation['sizeISO']  = int( commit_evaluation.DICTIONARY['nV'] * commit_evaluation.KERNELS['iso'].shape[0])

    regularisation['normIC']  = regnorms[0]
    regularisation['normEC']  = regnorms[1]
    regularisation['normISO'] = regnorms[2]

    regularisation['lambdaIC']  = float( lambdas[0] )
    regularisation['lambdaEC']  = float( lambdas[1] )
    regularisation['lambdaISO'] = float( lambdas[2] )

    regularisation['nnIC']  = is_nonnegative[0]
    regularisation['nnEC']  = is_nonnegative[1]
    regularisation['nnISO'] = is_nonnegative[2]

    # Check if group indices need to be updated in case of 'group_sparsity'
    if (structureIC is not None) and (0 in commit_evaluation.DICTIONARY['TRK']['kept']) :
        dictionary_TRK_kept = commit_evaluation.DICTIONARY['TRK']['kept']

        idx_in_kept = np.zeros(dictionary_TRK_kept.size, dtype=np.int32) - 1  # -1 is used to flag indices for removal
        idx_in_kept[dictionary_TRK_kept==1] = list(range(commit_evaluation.DICTIONARY['IC']['nF']))

        newStructureIC = []
        newWeightsIC = []
        for count, group in enumerate(structureIC):
            group = idx_in_kept[group]
            idx_to_delete = np.where(group==-1)[0]
            if idx_to_delete.size>0:
                group = np.delete(group,idx_to_delete)
                if(group.size>0):
                    newStructureIC.append(group)
                    newWeightsIC.append(weightsIC[count])
            else:
                newStructureIC.append(group)
                newWeightsIC.append(weightsIC[count])

        structureIC = np.array(newStructureIC, dtype=np.object_)
        weightsIC = np.array(newWeightsIC)

    regularisation['structureIC'] = structureIC
    regularisation['weightsIC']   = weightsIC
    regularisation['group_norm']  = group_norm

    return regularisation


def regularisation2omegaprox(regularisation):
    lambdaIC  = float(regularisation.get('lambdaIC'))
    lambdaEC  = float(regularisation.get('lambdaEC'))
    lambdaISO = float(regularisation.get('lambdaISO'))
    if lambdaIC<0.0 or lambdaEC<0.0 or lambdaISO<0.0:
        raise ValueError('Negative regularisation parameters are not allowed')

    normIC  = regularisation.get('normIC')
    normEC  = regularisation.get('normEC')
    normISO = regularisation.get('normISO')
    if not normIC in list_regnorms:
        raise ValueError('normIC not implemented')
    if not normEC in list_regnorms:
        raise ValueError('normEC not implemented')
    if not normISO in list_regnorms:
        raise ValueError('normISO not implemented')

    # Intracellular Compartment
    startIC = regularisation.get('startIC')
    sizeIC  = regularisation.get('sizeIC')
    if normIC is None:
        omegaIC = lambda x: 0.0
        if regularisation.get('nnIC')==True:
            proxIC = lambda x, scaling: non_negativity(x,startIC,sizeIC)
        else:
            proxIC = lambda x, scaling: x
    elif normIC == 'l1':
        omegaIC = lambda x: lambdaIC * np.linalg.norm(x[startIC:sizeIC],1)
        if regularisation.get('nnIC'):
            proxIC = lambda x, scaling: non_negativity(soft_thresholding(x,scaling*lambdaIC,startIC,sizeIC),startIC,sizeIC)
        else:
            proxIC = lambda x, scaling: non_negativity(x,startIC,sizeIC)
    # elif normIC == 'l2':
    #     omegaIC = lambda x: lambdaIC * np.linalg.norm(x[startIC:sizeIC])
    #     proxIC  = lambda x: projection_onto_l2_ball(x, lambdaIC, startIC, sizeIC)
    elif normIC == 'group_sparsity':
        structureIC = regularisation.get('structureIC')
        groupWeightIC = regularisation.get('weightsIC')
        if not len(structureIC) == len(groupWeightIC):
            raise ValueError('Number of groups and weights do not coincide.')
        group_norm = regularisation.get('group_norm')
        if not group_norm in list_group_sparsity_norms:
            raise ValueError('Wrong norm in the structured sparsity term. Choose between %s.' % str(list_group_sparsity_norms))

        # convert to new data structure (needed for faster access)
        N = np.sum([g.size for g in structureIC])
        groupIdxIC  = np.zeros( (N,), dtype=np.int32 )
        groupSizeIC = np.zeros( (structureIC.size,), dtype=np.int32 )
        pos = 0
        for i, g in enumerate(structureIC) :
            groupSizeIC[i] = g.size
            groupIdxIC[pos:(pos+g.size)] = g[:]
            pos += g.size

        omegaIC = lambda x: omega_group_sparsity( x, groupIdxIC, groupSizeIC, groupWeightIC, lambdaIC, group_norm )
        #TODO: verify if COMMIT2 results are better than before
        if regularisation.get('nnIC'):
            proxIC = lambda x, scaling: non_negativity(prox_group_sparsity(x,groupIdxIC,groupSizeIC,groupWeightIC,scaling*lambdaIC,group_norm),startIC,sizeIC)
        else:
            proxIC = lambda x, scaling: prox_group_sparsity(x,groupIdxIC,groupSizeIC,groupWeightIC,scaling*lambdaIC,group_norm)
    else:
        raise ValueError('Type of regularisation for IC compartment not recognized.')

    # Extracellular Compartment
    startEC = regularisation.get('startEC')
    sizeEC  = regularisation.get('sizeEC')
    if normEC is None:
        omegaEC = lambda x: 0.0
        if regularisation.get('nnEC')==True:
            proxEC = lambda x, scaling: non_negativity(x,startEC,sizeEC)
        else:
            proxEC = lambda x, scaling: x
    elif normEC == 'l1':
        omegaEC = lambda x: lambdaEC * np.linalg.norm(x[startEC:(startEC+sizeEC)],1)
        if regularisation.get('nnEC'):
            proxEC = lambda x, scaling: non_negativity(soft_thresholding(x,scaling*lambdaEC,startEC,sizeEC),startEC,sizeEC)
        else:
            proxEC = lambda x, scaling: soft_thresholding(x,scaling*lambdaEC,startEC,sizeEC)
    # elif normEC == 'l2':
    #     omegaEC = lambda x: lambdaEC * np.linalg.norm(x[startEC:(startEC+sizeEC)])
    #     proxEC  = lambda x: projection_onto_l2_ball(x, lambdaEC, startEC, sizeEC)
    else:
        raise ValueError('Type of regularisation for EC compartment not recognized.')

    # Isotropic Compartment
    startISO = regularisation.get('startISO')
    sizeISO  = regularisation.get('sizeISO')
    if normISO is None:
        omegaISO = lambda x: 0.0
        if regularisation.get('nnISO')==True:
            proxISO = lambda x, scaling: non_negativity(x,startISO,sizeISO)
        else:
            proxISO  = lambda x, scaling: x
    elif normISO == 'l1':
        omegaISO = lambda x: lambdaISO * np.linalg.norm(x[startISO:(startISO+sizeISO)],1)
        if regularisation.get('nnISO'):
            proxISO  = lambda x, scaling: non_negativity(soft_thresholding(x,scaling*lambdaISO,startISO,sizeISO),startISO,sizeISO)
        else:
            proxISO  = lambda x, scaling: soft_thresholding(x,scaling*lambdaISO,startISO,sizeISO)
    # elif normISO == 'l2':
    #     omegaISO = lambda x: lambdaISO * np.linalg.norm(x[startISO:(startISO+sizeISO)])
    #     proxISO  = lambda x: projection_onto_l2_ball(x, lambdaISO, startISO, sizeISO)
    else:
        raise ValueError('Type of regularisation for ISO compartment not recognized.')

    omega = lambda x: omegaIC(x) + omegaEC(x) + omegaISO(x)
    prox = lambda x, scaling: proxIC(proxEC(proxISO(x,scaling),scaling),scaling)
    return omega, prox


def evaluate_model(y, A, x, regularisation = None):
    if regularisation is None:
        omega = lambda x: 0.0
    else:
        omega, _ = regularisation2omegaprox(regularisation)

    return 0.5*np.linalg.norm(A.dot(x)-y)**2 + omega(x)


def solve(y, A, At, tol_fun=1e-4, tol_x=1e-6, max_iter=1000, verbose=True, x0=None, regularisation=None, confidence_array=None):
    """
    Solve the regularised least squares problem

        argmin_x 0.5*|| sqrt(W) ( Ax-y ) ||_2^2 + Omega(x)

    with the Omega described by 'regularisation' and W is the confidence_array

    Check the documentation of commit.solvers.init_regularisation to see how to
    solve a specific problem.
    """
    if regularisation is None:
        omega = lambda x: 0.0
        prox  = lambda x, scaling: x
    else:
        omega, prox = regularisation2omegaprox(regularisation)

    if x0 is None:
        x0 = np.zeros(A.shape[1])

    if confidence_array is not None:
        confidence_array = np.sqrt(confidence_array)

    return fista( y, A, At, omega, prox, confidence_array, tol_fun, tol_x, max_iter, verbose, x0)


def fista( y, A, At, omega, prox, sqrt_W=None, tol_fun=1e-4, tol_x=1e-6, max_iter=1000, verbose=False, x0=None) :
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
    if verbose :
        print()
        print( "      |  1/2||Ax-y||^2      Omega      |  Cost function    Abs error      Rel error    |      Abs x          Rel x    " )
        print( "------|--------------------------------|-----------------------------------------------|------------------------------" )
    iter = 1
    while True :
        if verbose :
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
        if verbose :
            print( "  %13.7e  %13.7e  |  %13.7e  %13.7e  %13.7e  |  %13.7e  %13.7e" % ( 0.5 * res_norm**2, reg_term_x, curr_obj, abs_obj, rel_obj, abs_x, rel_x ) )

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

    if verbose :
        print( "< Stopping criterion: %s >" % criterion )

    opt_details = {}
    opt_details['residual'] = 0.5*res_norm**2
    opt_details['regterm'] = reg_term_x
    opt_details['cost_function'] = curr_obj
    opt_details['abs_cost'] = abs_obj
    opt_details['rel_cost'] = rel_obj
    opt_details['abs_x'] = abs_x
    opt_details['rel _x'] = rel_x
    opt_details['iterations'] = iter
    opt_details['stopping_criterion'] = criterion

    return x, opt_details
