import numpy as np
import pandas as pd
import time

from parser import read_data
from meta import meta_data


def calc_opt_params(beta,exp,exp_sig,sig):
    assert beta.shape == exp.shape
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


def saxs_optimization_ensemble(exp_data, bc_data, indices):
    # prepare data
    exp_saxs = exp_data['saxs'].data['value'].values  # shape: (37,)
    exp_sigma = exp_data['saxs'].data['error'].values  # shape: (37,)
    bc_ensemble = bc_data['saxs'].data.values[indices, :]  # shape: (100, 37)
    bc_saxs = np.mean(bc_ensemble, axis=0)  # shape: (37,)

    # optimization
    opt_params = calc_opt_params(bc_saxs, exp_saxs, exp_sigma, bc_data['saxs'].sigma)
    f, f_comps = calc_score(bc_saxs, exp_saxs, exp_sigma, bc_data['saxs'].sigma, opt_params)

    sse_saxs = np.sum((exp_saxs - bc_saxs) ** 2.0)
    total_score_saxs = np.sum(f)

    return sse_saxs, total_score_saxs, bc_saxs


def cs_optimization_ensemble(exp_data, bc_data, indices):
    # prepare data
    exp_cs = exp_data['cs'].data['value'].values  # shape: (262,)
    exp_sigma = exp_data['cs'].data['error'].values  # shape: (262,)
    atom_types = exp_data['cs'].data['atomname'].values # shape: (262,)

    bc_ensemble = bc_data['cs'].data.values[indices, :]  # shape: (100, 262)
    bc_cs = np.mean(bc_ensemble, axis=0)  # shape: (262,)
    bc_sigma = np.array([bc_data['cs'].sigma[atom_type] for atom_type in atom_types])  # shape: (262,)

    # optimization
    opt_params = calc_opt_params(bc_cs, exp_cs, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc_cs, exp_cs, exp_sigma, bc_sigma, opt_params)

    sse = np.sum((exp_cs - bc_cs) ** 2.0)
    total_score = np.sum(f)

    return sse, total_score, bc_cs


def fret_optimization_ensemble(exp_data, bc_data, indices):
    # prepare data
    exp = exp_data['fret'].data  # scalar
    exp_sigma = exp_data['fret'].sigma  # scalar

    bc_ensemble = bc_data['fret'].data.values[indices, :]  # shape: (100, 1)
    bc = np.mean(bc_ensemble, axis=0)  # shape: (1,)
    bc_sigma = bc_data['fret'].sigma # scalar

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_sigma, opt_params)

    sse = np.sum((exp - bc) ** 2.0)
    total_score = np.sum(f)

    return sse, total_score, bc


def main(exp_data, bc_data, ens_size=100, epochs=250, output_file=None, verbose=False):
    """

    Parameters
    ----------
    ens_size
    epochs
    mode: str, optional (default: 'trades_all')
        This parameter must be one of the followings:
            - trades_all
            - trades-unfold
            - ensemble
            - mixed

    output_file
    verbose

    Returns
    -------

    """

    #  get size
    pool_size = bc_data['cs'].shape[0]
    if verbose: print("\n### pool size: %i"%pool_size)

    # time
    t0 = time.time()

    for it in range(epochs):
        if verbose: print("\n### iteration: %i  (elapsed time: %f seconds)"%(it+1, time.time()-t0))

        # random selection with replacement
        indices = np.random.randint(0, pool_size-1, ens_size)   # shape: (100,)

        # SAXS optimization on ensemble
        sse_saxs, total_score_saxs, bc_saxs = saxs_optimization_ensemble(exp_data, bc_data, indices)

        # CS optimization on ensemble
        sse_cs, total_score_cs, bc_cs = cs_optimization_ensemble(exp_data, bc_data, indices)

        # FRET optimization on ensemble
        sse_fret, total_score_fret, bc_fret  = fret_optimization_ensemble(exp_data, bc_data, indices)



# read data
filenames = meta_data()
exp_data = read_data(filenames[0], mode='exp')
bc_data = read_data(filenames[1], mode='bc')
#Todo: fix for other modes
