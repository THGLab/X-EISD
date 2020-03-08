"""
This module contains all functions that are required to perform a single property maximum loglike optimization.
"""

import numpy as np


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
    exp_sigma = exp_data['saxs'].data['error'].values * 0.1  # shape: (37,)

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


def vect_calc_opt_params_jc(alpha1, alpha2, exp_j, exp_sig, mus, sigs):
    """
    takes an average phi value and experimental data point and calculates optimal A, B, C
    calculates parameters for a single phi angle
    :return:
    """
    a = np.zeros((alpha1.shape[0], 3, 3))
    b = np.zeros((alpha1.shape[0], 3))

    a[:, 0, 0] = 1.0 / (sigs[0] ** 2.0) + ((alpha2 / exp_sig) ** 2.0)
    a[:, 1, 1] = 1.0 / (sigs[1] ** 2.0) + ((alpha1 / exp_sig) ** 2.0)
    a[:, 2, 2] = 1.0 / (sigs[2] ** 2.0) + 1.0 / (exp_sig ** 2.0)

    a[:, 0, 1] = alpha1 * alpha2 / (exp_sig ** 2.0)
    a[:, 1, 0] = alpha1 * alpha2 / (exp_sig ** 2.0)
    a[:, 0, 2] = alpha2 / (exp_sig ** 2.0)
    a[:, 2, 0] = alpha2 / (exp_sig ** 2.0)
    a[:, 1, 2] = alpha1 / (exp_sig ** 2.0)
    a[:, 2, 1] = alpha1 / (exp_sig ** 2.0)

    b[:, 0] = mus[0] / (sigs[0] ** 2.0) + exp_j * alpha2 / (exp_sig ** 2)
    b[:, 1] = mus[1] / (sigs[1] ** 2.0) + exp_j * alpha1 / (exp_sig ** 2)
    b[:, 2] = mus[2] / (sigs[2] ** 2.0) + exp_j / (exp_sig ** 2)

    opt_params = np.array([np.linalg.solve(a[i], b[i]) for i in range(a.shape[0])])  # shape: (47,3)

    return opt_params


def vect_calc_score_JC(alpha1, alpha2, exp_j, exp_sig, opt_params, mus, sigs):
    """
    calculates score for a single phi angle/residue, can be used on a single protein or ensemble
    calculates score given 1. alpha1 and alpha2 2. exptl J and sigma 3. optimized parameters 4. mus and sigs for A,B,C
    :returns: total score f and array of [f_a, f_b, f_c, f_err]
    """
    f_a = normal_loglike(opt_params[:,0],mus[0],sigs[0]) # shape: (47,)
    f_b = normal_loglike(opt_params[:,1],mus[1],sigs[1]) # shape: (47,)
    f_c = normal_loglike(opt_params[:,2],mus[2],sigs[2]) # shape: (47,)
    err = exp_j - opt_params[:,0]*alpha2 - opt_params[:,1]*alpha1 - opt_params[:,2] # shape: (47,)
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
    opt_params = vect_calc_opt_params_jc(bc_alpha1, bc_alpha2, exp, exp_sigma, bc_mu, bc_sigma) #shape: (47,3)
    f, f_comps = vect_calc_score_JC(bc_alpha1, bc_alpha2, exp, exp_sigma, opt_params, bc_mu, bc_sigma)
    sse = (opt_params[:,0]*bc_alpha2 + opt_params[:,1]*bc_alpha1 + opt_params[:,2] - exp) ** 2.0
    sse = np.sum(sse,axis=0)
    total_score = np.sum(f, axis=0)
    jcoup_vals = list(opt_params[:,0] * bc_alpha2 + opt_params[:,1] * bc_alpha1 + opt_params[:,2])

    return sse, total_score, [bc_alpha1, bc_alpha2], jcoup_vals


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
