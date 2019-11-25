import numpy as np
np.random.seed(90)
import pandas as pd
import time

from parser import read_data
from meta import meta_data


def calc_opt_params(beta,exp,exp_sig,sig):
    ratio = (sig**2.0)/(exp_sig**2.0)
    opt_params = (ratio*(exp-beta))/(1.0+ratio)
    return opt_params


def normal_loglike(x, mu, sig):
    exp_val = -((x - mu)** 2.0)/(2.0 *(sig ** 2.0))
    pre_exp = 1.0/(np.sqrt(2.0*np.pi*(sig ** 2.0)))
    logp = np.log(pre_exp*np.exp(exp_val))
    return logp


def calc_score(beta,exp,exp_sig,sig,opt_params):
    f_q = normal_loglike(opt_params,0,sig)
    err = exp - opt_params - beta
    f_err = normal_loglike(err,0,exp_sig)
    f = f_q + f_err
    f_comps = [f_q, f_err]
    return f, f_comps


def saxs_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp_saxs = exp_data['saxs'].data['value'].values  # shape: (37,)
    exp_sigma = exp_data['saxs'].data['error'].values  # shape: (37,)

    if indices is None:
        bc_saxs = old_vals - (bc_data['saxs'].data.values[popped_structure, :] - bc_data['saxs'].data.values[new_index, :] )/100.
    else:
        bc_ensemble = bc_data['saxs'].data.values[indices, :]  # shape: (100, 37)
        bc_saxs = np.mean(bc_ensemble, axis=0)  # shape: (37,)

    # optimization
    opt_params = calc_opt_params(bc_saxs, exp_saxs, exp_sigma, bc_data['saxs'].sigma)
    f, f_comps = calc_score(bc_saxs, exp_saxs, exp_sigma, bc_data['saxs'].sigma, opt_params)

    sse_saxs = np.sum((exp_saxs - bc_saxs) ** 2.0)
    total_score_saxs = np.sum(f)

    return sse_saxs, total_score_saxs, bc_saxs#{i+1:bc_saxs[i] for i in range(len(bc_saxs))}


def cs_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp_cs = exp_data['cs'].data['value'].values  # shape: (262,)
    exp_sigma = exp_data['cs'].data['error'].values  # shape: (262,)
    atom_types = exp_data['cs'].data['atomname'].values # shape: (262,)

    if indices is None:
        bc_cs = old_vals - (bc_data['cs'].data.values[popped_structure, :] - bc_data['cs'].data.values[new_index, :] )/100.
    else:
        bc_ensemble = bc_data['cs'].data.values[indices, :]  # shape: (100, 262)
        bc_cs = np.mean(bc_ensemble, axis=0)  # shape: (262,)

    bc_sigma = np.array([bc_data['cs'].sigma[atom_type] for atom_type in atom_types])  # shape: (262,)

    # optimization
    opt_params = calc_opt_params(bc_cs, exp_cs, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc_cs, exp_cs, exp_sigma, bc_sigma, opt_params)

    sse = np.sum((exp_cs - bc_cs) ** 2.0)
    total_score = np.sum(f)

    return sse, total_score, bc_cs#{i+1:bc_cs[i] for i in range(len(bc_cs))}


def fret_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp = exp_data['fret'].data  # scalar
    exp_sigma = exp_data['fret'].sigma  # scalar

    if indices is None:
        bc = old_vals - (
                    bc_data['fret'].data.values[popped_structure, :] - bc_data['fret'].data.values[new_index, :]) / 100.
    else:
        bc_ensemble = bc_data['fret'].data.values[indices, :]  # shape: (100, 1)
        bc = np.mean(bc_ensemble, axis=0)  # shape: (1,)

    bc_sigma = bc_data['fret'].sigma # scalar

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_sigma, opt_params)

    sse = np.sum((exp - bc) ** 2.0)
    total_score = np.sum(f)

    return sse, total_score, bc[0]


def determinant3(a):
    """
    Determinant of a 3x3 matrix
    :param a: 3x3 matrix
    :return: determinant of a
    """
    det = a[0, 0] * a[1, 1] * a[2, 2]
    det += a[0, 1] * a[1, 2] * a[2, 0]
    det += a[0, 2] * a[1, 0] * a[2, 1]
    det -= a[0, 2] * a[1, 1] * a[2, 0]
    det -= a[0, 1] * a[1, 0] * a[2, 2]
    det -= a[0, 0] * a[1, 2] * a[2, 1]
    return det


def solve_3_eqs(a, b):
    """
    Solve Ax=b for a 3x3 A and a 3x1 b
    :param a: 3x3 matrix
    :param b: size 3 vector
    :return: x
    """
    ax = np.array([
        [b[0], a[0, 1], a[0, 2]],
        [b[1], a[1, 1], a[1, 2]],
        [b[2], a[2, 1], a[2, 2]]
    ])

    ay = np.array([
        [a[0, 0], b[0], a[0, 2]],
        [a[1, 0], b[1], a[1, 2]],
        [a[2, 0], b[2], a[2, 2]]
    ])

    az = np.array([
        [a[0, 0], a[0, 1], b[0]],
        [a[1, 0], a[1, 1], b[1]],
        [a[2, 0], a[2, 1], b[2]]
    ])

    det = determinant3(a)
    detx = determinant3(ax)
    dety = determinant3(ay)
    detz = determinant3(az)

    x = np.array([detx / det, dety / det, detz / det])
    return x


def calc_opt_params_JC(alpha1, alpha2,exp_j,exp_sig,mus,sigs):
    """
    takes an average phi value and experimental data point and calculates optimal A, B, C
    calculates parameters for a single phi angle
    :return:
    """
    a = np.zeros((3, 3))
    b = np.zeros((3,))

    a[0][0] = 1.0/(sigs[0] ** 2.0) + ((alpha2/exp_sig)**2.0)
    a[1][1] = 1.0/(sigs[1] ** 2.0) + ((alpha1/exp_sig)**2.0)
    a[2][2] = 1.0/(sigs[2] ** 2.0) + 1.0/(exp_sig ** 2.0)

    a[0][1] = alpha1*alpha2/(exp_sig ** 2.0)
    a[1][0] = alpha1*alpha2/(exp_sig ** 2.0)
    a[0][2] = alpha2/(exp_sig ** 2.0)
    a[2][0] = alpha2/(exp_sig ** 2.0)
    a[1][2] = alpha1/(exp_sig ** 2.0)
    a[2][1] = alpha1/(exp_sig ** 2.0)

    b[0] = mus[0]/(sigs[0] ** 2.0) + exp_j*alpha2/(exp_sig ** 2)
    b[1] = mus[1]/(sigs[1] ** 2.0) + exp_j*alpha1/(exp_sig ** 2)
    b[2] = mus[2]/(sigs[2] ** 2.0) + exp_j/(exp_sig ** 2)

    opt_params = solve_3_eqs(a,b)

    return opt_params


def calc_score_JC(alpha1, alpha2, exp_j, exp_sig, opt_params, mus, sigs):
    """
    calculates score for a single phi angle/residue, can be used on a single protein or ensemble
    calculates score given 1. alpha1 and alpha2 2. exptl J and sigma 3. optimized parameters 4. mus and sigs for A,B,C
    :returns: total score f and array of [f_a, f_b, f_c, f_err]
    """
    f_a = normal_loglike(opt_params[0],mus[0],sigs[0])
    f_b = normal_loglike(opt_params[1],mus[1],sigs[1])
    f_c = normal_loglike(opt_params[2],mus[2],sigs[2])
    err = exp_j - opt_params[0]*alpha2 - opt_params[1]*alpha1 - opt_params[2]
    f_err = normal_loglike(err,0,exp_sig)
    f = f_a + f_b + f_c + f_err
    f_comps = [f_a, f_b, f_c, f_err]
    return f, f_comps


def jc_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp = exp_data['jc'].data['value'].values  # shape: (47,)
    exp_sigma = exp_data['jc'].data['error'].values  # shape: (47,)

    if indices is None:
        pop_alpha = bc_data['jc'].data.values[popped_structure, :]
        add_alpha = bc_data['jc'].data.values[new_index, :]

        bc_alpha1 = old_vals[0] - (pop_alpha - add_alpha)/100.
        assert bc_alpha1.shape == (47,)

        bc_alpha2 = old_vals[1] - (np.square(pop_alpha) - np.square(add_alpha))/100.  # shape: (47,)
        assert bc_alpha2.shape == (47,)
    else:
        bc_ensemble_alpha1 = bc_data['jc'].data.values[indices, :]  # shape: (100, 47)
        bc_alpha1 = np.mean(bc_ensemble_alpha1, axis=0)  # shape: (47,)

        bc_ensemble_alpha2 = np.square(bc_data['jc'].data.values[indices, :])  # shape: (100, 47)
        bc_alpha2 = np.mean(bc_ensemble_alpha2, axis=0)  # shape: (47,)

    bc_sigma = [bc_data['jc'].sigma[i] for i in ["A", "B", "C"]]
    bc_mu = [bc_data['jc'].mu[i] for i in ["A", "B", "C"]]

    # optimization
    sse = 0.0
    total_score = 0.0
    bc_alpha_vals = {}
    jcoup_vals = []
    for k in range(len(exp)):
        opt_params = calc_opt_params_JC(bc_alpha1[k], bc_alpha2[k], exp[k], exp_sigma[k], bc_mu, bc_sigma)
        f, f_comps = calc_score_JC(bc_alpha1[k], bc_alpha2[k], exp[k], exp_sigma[k], opt_params, bc_mu, bc_sigma)

        sse += (opt_params[0]*bc_alpha2[k] + opt_params[1]*bc_alpha1[k] + opt_params[2] - exp[k]) ** 2.0
        total_score += f
        bc_alpha_vals[k+1] = [bc_alpha1[k], bc_alpha2[k]]
        jcoup_vals.append(opt_params[0] * bc_alpha2[k] + opt_params[1] * bc_alpha1[k] + opt_params[2])

    return sse, total_score, [bc_alpha1, bc_alpha2], jcoup_vals#bc_alpha_vals


def noe_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    dist_value = exp_data['noe'].data['dist_value'].values  # shape: (93,)
    upper_bound_value = exp_data['noe'].data['upper_bound_value'].values  # shape: (93,)
    lower_bound_value = exp_data['noe'].data['lower_bound_value'].values  # shape: (93,)
    range_val = upper_bound_value - lower_bound_value
    exp_sigma = range_val / 2.0

    exp_distance = []  # shape: (93,)
    for k in range(len(dist_value)):
        if dist_value[k] == upper_bound_value[k] or dist_value[k] == lower_bound_value[k]:
            exp_distance.append((upper_bound_value[k] + lower_bound_value[k]) / 2.0)
        else:
            exp_distance.append(dist_value[k])
    exp_distance = np.array(exp_distance)

    if indices is None:
        popped = np.power(bc_data['noe'].data.values[popped_structure, :], -6.0)
        added = np.power(bc_data['noe'].data.values[new_index, :], -6.0)
        avg_distance = (np.power(old_vals, -6.0)*100. - (popped - added) )/100.
        avg_distance = np.power(avg_distance, (-1./6.))
    else:
        bc_ensemble = np.power(bc_data['noe'].data.values[indices, :], -6.0)  # shape: (100, 93)
        avg_distance = np.power(np.mean(bc_ensemble, axis=0), (-1./6.))  # shape: (93,)

    # optimization
    opt_params = calc_opt_params(avg_distance, exp_distance, exp_sigma, bc_data['noe'].sigma)
    f, f_comps = calc_score(avg_distance, exp_distance, exp_sigma, bc_data['noe'].sigma, opt_params)

    sse = np.sum((exp_distance - avg_distance) ** 2.0, axis=0)
    total_score = np.sum(f)

    return sse, total_score, avg_distance#{i+1:avg_distance[i] for i in range(len(avg_distance))}


def pre_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp_distance = exp_data['pre'].data['value'].values  # shape: (68,)
    upper_bound_value = exp_data['pre'].data['upper'].values  # shape: (68,)
    lower_bound_value = exp_data['pre'].data['lower'].values  # shape: (68,)
    range_val = upper_bound_value + lower_bound_value
    exp_sigma = range_val / 4.0

    if indices is None:
        popped = np.power(bc_data['pre'].data.values[popped_structure, :], -6.0)
        added = np.power(bc_data['pre'].data.values[new_index, :], -6.0)
        avg_distance = (np.power(old_vals, -6.0)*100. - (popped - added) )/100.
        avg_distance = np.power(avg_distance, (-1./6.))
    else:
        bc_ensemble = np.power(bc_data['pre'].data.values[indices, :], -6.0)  # shape: (100, 68)
        avg_distance = np.power(np.mean(bc_ensemble, axis=0), (-1./6.))  # shape: (68,)

    # optimization
    opt_params = calc_opt_params(avg_distance, exp_distance, exp_sigma, bc_data['pre'].sigma)
    f, f_comps = calc_score(avg_distance, exp_distance, exp_sigma, bc_data['pre'].sigma, opt_params)

    sse = np.sum((exp_distance - avg_distance) ** 2.0, axis=0)
    total_score = np.sum(f)

    return sse, total_score, avg_distance#{i+1:avg_distance[i] for i in range(len(avg_distance))}


def rdc_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp = exp_data['rdc'].data['value'].values  # shape: (28,)
    exp_sigma = exp_data['rdc'].data['error'].values  # shape: (28,)

    if indices is None:
        bc = old_vals - (
                    bc_data['rdc'].data.values[popped_structure, :] - bc_data['rdc'].data.values[new_index, :]) / 100.
    else:
        bc_ensemble = bc_data['rdc'].data.values[indices, :]  # shape: (100, 28)
        bc = np.mean(bc_ensemble, axis=0)  # shape: (28,)

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_data['rdc'].sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_data['rdc'].sigma, opt_params)

    sse = np.sum((exp - bc) ** 2.0, axis=0)
    total_score = np.sum(f)

    return sse, total_score, bc#{i+1:bc[i] for i in range(len(bc))}


def rh_optimization_ensemble(exp_data, bc_data, indices, old_vals=None, popped_structure=None, new_index=None):
    # prepare data
    exp = exp_data['rh'].data  # scalar
    exp_sigma = exp_data['rh'].sigma  # scalar

    if indices is None:
        bc = old_vals - (
                    bc_data['rh'].data.values[popped_structure, :] - bc_data['rh'].data.values[new_index, :]) / 100.
    else:
        bc_ensemble = bc_data['rh'].data.values[indices, :]  # shape: (100, 1)
        bc = np.mean(bc_ensemble, axis=0)  # shape: (1,)

    bc_sigma = bc_data['rh'].sigma # scalar

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_sigma, opt_params)

    sse = np.sum((exp - bc) ** 2.0)
    total_score = np.sum(f)

    return sse, total_score, bc[0]


def maximize_score_SAXS(exp_data, bc_data, ens_size, indices, num_structs, old_score_SAXS, old_score_CS, old_score_FRET, old_score_JC,
                        old_score_NOEs, old_score_PREs, old_score_RDCs, old_score_RH, bc_SAXS_vals, bc_shift_vals,
                        bc_FRET_val, bc_alpha_vals, bc_dist_vals_NOE, bc_dist_vals_PRE, bc_RDC_vals, bc_RH_val):
    indices = list(indices)
    old_SAXS_vals = bc_SAXS_vals
    old_shift_vals = bc_shift_vals
    old_FRET_val = bc_FRET_val
    old_alpha_vals = bc_alpha_vals
    old_dist_vals_NOE = bc_dist_vals_NOE
    old_dist_vals_PRE = bc_dist_vals_PRE
    old_RDC_vals = bc_RDC_vals
    old_RH_val = bc_RH_val
    accepted = 0
    for iterations in range(10000):
        pop_index = np.random.randint(0, ens_size, 1)[0]
        popped_structure = indices[pop_index]
        indices.pop(pop_index)
        struct_found = False
        while not struct_found:
            new_index = np.random.randint(0, num_structs, 1)[0]
            if new_index != popped_structure:
                indices.append(new_index)
                struct_found = True

        # SAXS
        restraint_SSE_SAXS, new_score_SAXS, new_SAXS_vals = saxs_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_SAXS_vals, popped_structure,
                                                                                       new_index)

        # CS
        restraint_SSE_CS, new_score_CS, new_shift_vals = cs_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_shift_vals, popped_structure,
                                                                                       new_index)

        # FRET
        restraint_SSE_FRET, new_score_FRET, new_FRET_val = fret_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_FRET_val, popped_structure,
                                                                                       new_index)

        # JC
        restraint_SSE_JC, new_score_JC, new_alpha_vals, jcoup_vals = jc_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_alpha_vals, popped_structure,
                                                                                       new_index)

        # NOE
        restraint_SSE_NOE, new_score_NOEs, new_dist_vals_NOE = noe_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_dist_vals_NOE, popped_structure,
                                                                                       new_index)

        # PRE
        restraint_SSE_PRE, new_score_PREs, new_dist_vals_PRE = pre_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_dist_vals_PRE, popped_structure,
                                                                                       new_index)

        # RDC
        restraint_SSE_RDC, new_score_RDCs, new_RDC_vals = rdc_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_RDC_vals, popped_structure,
                                                                                       new_index)

        # RH
        restraint_SSE_RH, new_score_RH, new_RH_val = rh_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_RH_val, popped_structure,
                                                                                       new_index)

        # if iterations%1000 == 0:
        #     print("\n*** iteration %i\n"%iterations)
        #     print(popped_structure)
        #     print(new_score_SAXS, new_score_CS, new_score_FRET, new_score_JC)
        #     print(new_score_NOEs, new_score_PREs, new_score_RDCs, new_score_RH)

        if old_score_SAXS + old_score_CS + old_score_FRET + old_score_JC + old_score_NOEs + old_score_PREs + \
                old_score_RDCs + old_score_RH > new_score_SAXS + new_score_CS + new_score_FRET + new_score_JC + \
                new_score_NOEs + new_score_PREs + new_score_RDCs + new_score_RH:
            indices.pop(-1)
            indices.append(popped_structure)
            # print repr(i) + ' 0 ' + repr(new_score_SAXS) + ' ' + repr(old_score_SAXS) + ' ' + repr(indices)
        else:
            #print repr(i+1) + ' ' + repr(new_score_SAXS) + ' ' + repr(old_score_SAXS) + ' ' + repr(indices)
            old_score_SAXS = new_score_SAXS
            old_score_CS = new_score_CS
            old_SAXS_vals = new_SAXS_vals
            old_shift_vals = new_shift_vals
            old_score_FRET = new_score_FRET
            old_FRET_val = new_FRET_val
            old_score_JC = new_score_JC
            old_alpha_vals = new_alpha_vals
            old_score_NOEs = new_score_NOEs
            old_dist_vals_NOE = new_dist_vals_NOE
            old_score_PREs = new_score_PREs
            old_dist_vals_PRE = new_dist_vals_PRE
            old_score_RDCs = new_score_RDCs
            old_RDC_vals = new_RDC_vals
            old_score_RH = new_score_RH
            old_RH_val = new_RH_val
            accepted = accepted + 1

            new_SSE_SAXS = restraint_SSE_SAXS
            new_SSE_CS = restraint_SSE_CS
            new_SSE_FRET = restraint_SSE_FRET
            new_SSE_JC = restraint_SSE_JC
            new_SSE_NOE = restraint_SSE_NOE
            new_SSE_PRE = restraint_SSE_PRE
            new_SSE_RDC = restraint_SSE_RDC
            new_SSE_RH = restraint_SSE_RH

    return old_score_SAXS, new_SSE_SAXS, old_score_CS, new_SSE_CS, old_score_FRET, new_SSE_FRET, old_score_JC, new_SSE_JC, old_score_NOEs, new_SSE_NOE, old_score_PREs, new_SSE_PRE, old_score_RDCs, new_SSE_RDC, old_score_RH, new_SSE_RH, accepted, indices, jcoup_vals


def main(exp_data, bc_data, ens_size=100, epochs=250, output_file=None, verbose=False):
    """"""

    #  get size
    pool_size = bc_data['cs'].data.shape[0]
    if verbose: print("\n### pool size: %i"%pool_size)

    # time
    t0 = time.time()

    textfile = open(output_file, 'w')
    for it in range(epochs):
        # random selection with replacement
        indices = np.random.randint(0, pool_size, ens_size)   # shape: (100,)

        # SAXS optimization on ensemble
        sse_saxs, total_score_saxs, bc_saxs = saxs_optimization_ensemble(exp_data, bc_data, indices)

        # CS optimization on ensemble
        sse_cs, total_score_cs, bc_cs = cs_optimization_ensemble(exp_data, bc_data, indices)

        # FRET optimization on ensemble
        sse_fret, total_score_fret, bc_fret  = fret_optimization_ensemble(exp_data, bc_data, indices)

        # JC optimization on ensemble
        sse_jc, total_score_jc, bc_alpha_vals_jc, _  = jc_optimization_ensemble(exp_data, bc_data, indices)

        # NOEs optimization on ensemble
        sse_noe, total_score_noe, bc_dist_vals_noe  = noe_optimization_ensemble(exp_data, bc_data, indices)

        # PREs optimization on ensemble
        sse_pre, total_score_pre, bc_pre  = pre_optimization_ensemble(exp_data, bc_data, indices)

        # RDCs optimization on ensemble
        sse_rdc, total_score_rdc, bc_rdc  = rdc_optimization_ensemble(exp_data, bc_data, indices)

        # RH optimization on ensemble
        sse_rh, total_score_rh, bc_rh  = rh_optimization_ensemble(exp_data, bc_data, indices)

        # Fianl maximization
        max_score_SAXS, opt_SSE_SAXS, max_score_CS, opt_SSE_CS, max_score_FRET, opt_SSE_FRET, max_score_JC, opt_SSE_JC, \
        max_score_NOEs, opt_SSE_NOE, max_score_PREs, opt_SSE_PRE, max_score_RDCs, opt_SSE_RDC, max_score_RH, opt_SSE_RH, \
        accepted, new_indices, best_jcoups = maximize_score_SAXS(exp_data, bc_data, ens_size, indices, pool_size, total_score_saxs,
                                                                 total_score_cs, total_score_fret, total_score_jc,
                                                                 total_score_noe, total_score_pre, total_score_rdc,
                                                                 total_score_rh, bc_saxs, bc_cs,
                                                                 bc_fret,
                                                                 bc_alpha_vals_jc, bc_dist_vals_noe, bc_pre,
                                                                 bc_rdc, bc_rh)
        s = repr(it) + ' ' + repr(accepted) + ' ' + repr(max_score_SAXS) + ' ' + repr(
            max_score_CS) + ' ' + repr(max_score_FRET) + ' ' + repr(max_score_JC) + ' ' + repr(
            max_score_NOEs) + ' ' + repr(max_score_PREs) + ' ' + repr(max_score_RDCs) + ' ' + repr(
            max_score_RH) + ' ' + repr((opt_SSE_SAXS / 37.0) ** 0.5) + ' ' + repr(
            (opt_SSE_CS / 262.0) ** 0.5) + ' ' + repr(opt_SSE_FRET ** 0.5) + ' ' + repr(
            (opt_SSE_JC / 47.0) ** 0.5) + ' ' + repr((opt_SSE_NOE / 93.0) ** 0.5) + ' ' + repr(
            (opt_SSE_PRE / 68.0) ** 0.5) + ' ' + repr((opt_SSE_RDC / 28.0) ** 0.5) + ' ' + repr(
            opt_SSE_RH ** 0.5) + ' ' + repr(new_indices) + ' ' + repr(best_jcoups) + '\n'
        textfile.write(s)
        textfile.flush()

        if verbose: print("\n### iteration: %i  (elapsed time: %f seconds)"%(it+1, time.time()-t0))

    textfile.close()


# read data
filenames = meta_data()
exp_data = read_data(filenames[0], mode='exp')
bc_data = read_data(filenames[1], mode='trades_uf')
main(exp_data, bc_data, epochs=250, output_file='drksh3_unfolded_trades.txt', verbose=True)
#Todo: fix for other modes
