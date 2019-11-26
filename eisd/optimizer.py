import numpy as np
np.random.seed(90)

import time

from eisd.utils import modes
from eisd.scorers import *


def maximize_score_SAXS(exp_data, bc_data, ens_size, indices, num_structs, flags,
                        old_score_SAXS, old_score_CS, old_score_FRET, old_score_JC,
                        old_score_NOEs, old_score_PREs, old_score_RDCs, old_score_RH,
                        bc_SAXS_vals, bc_shift_vals, bc_FRET_val, bc_alpha_vals,
                        bc_dist_vals_NOE, bc_dist_vals_PRE, bc_RDC_vals, bc_RH_val):
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
        if flags['saxs']:
            restraint_SSE_SAXS, new_score_SAXS, new_SAXS_vals = saxs_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_SAXS_vals, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_SAXS, new_score_SAXS, new_SAXS_vals = [0,0,0]

        # CS
        if flags['cs']:
            restraint_SSE_CS, new_score_CS, new_shift_vals = cs_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_shift_vals, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_CS, new_score_CS, new_shift_vals = [0,0,0]


        # FRET
        if flags['fret']:
            restraint_SSE_FRET, new_score_FRET, new_FRET_val = fret_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_FRET_val, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_FRET, new_score_FRET, new_FRET_val = [0,0,0]

        # JC
        if flags['jc']:
            restraint_SSE_JC, new_score_JC, new_alpha_vals, jcoup_vals = jc_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_alpha_vals, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_JC, new_score_JC, new_alpha_vals, jcoup_vals = [0,0,0,[0]]

        # NOE
        if flags['noe']:
            restraint_SSE_NOE, new_score_NOEs, new_dist_vals_NOE = noe_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_dist_vals_NOE, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_NOE, new_score_NOEs, new_dist_vals_NOE = [0,0,0]


        # PRE
        if flags['pre']:
            restraint_SSE_PRE, new_score_PREs, new_dist_vals_PRE = pre_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_dist_vals_PRE, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_PRE, new_score_PREs, new_dist_vals_PRE = [0,0,0]


        # RDC
        if flags['rdc']:
            restraint_SSE_RDC, new_score_RDCs, new_RDC_vals = rdc_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_RDC_vals, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_RDC, new_score_RDCs, new_RDC_vals = [0,0,0]


        # RH
        if flags['rh']:
            restraint_SSE_RH, new_score_RH, new_RH_val = rh_optimization_ensemble(exp_data, bc_data, None,
                                                                                       old_RH_val, popped_structure,
                                                                                       new_index)
        else:
            restraint_SSE_RH, new_score_RH, new_RH_val = [0,0,0]


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
            old_score_FRET = new_score_FRET
            old_score_JC = new_score_JC
            old_score_NOEs = new_score_NOEs
            old_score_PREs = new_score_PREs
            old_score_RDCs = new_score_RDCs
            old_score_RH = new_score_RH

            old_SAXS_vals = new_SAXS_vals
            old_shift_vals = new_shift_vals
            old_FRET_val = new_FRET_val
            old_alpha_vals = new_alpha_vals
            old_dist_vals_NOE = new_dist_vals_NOE
            old_dist_vals_PRE = new_dist_vals_PRE
            old_RDC_vals = new_RDC_vals
            old_RH_val = new_RH_val

            new_SSE_SAXS = restraint_SSE_SAXS
            new_SSE_CS = restraint_SSE_CS
            new_SSE_FRET = restraint_SSE_FRET
            new_SSE_JC = restraint_SSE_JC
            new_SSE_NOE = restraint_SSE_NOE
            new_SSE_PRE = restraint_SSE_PRE
            new_SSE_RDC = restraint_SSE_RDC
            new_SSE_RH = restraint_SSE_RH

            accepted = accepted + 1

    return old_score_SAXS, new_SSE_SAXS, old_score_CS, new_SSE_CS, old_score_FRET, new_SSE_FRET, old_score_JC, new_SSE_JC, old_score_NOEs, new_SSE_NOE, old_score_PREs, new_SSE_PRE, old_score_RDCs, new_SSE_RDC, old_score_RH, new_SSE_RH, accepted, indices, jcoup_vals


def main(exp_data, bc_data, ens_size=100, epochs=250, mode='all', output_file=None, verbose=False):
    """"""

    #  get size
    pool_size = bc_data['cs'].data.shape[0]
    if verbose: print("\n### pool size: %i"%pool_size)

    # time
    t0 = time.time()

    flags = modes(mode)

    textfile = open(output_file, 'w')
    for it in range(epochs):
        # random selection with replacement
        indices = np.random.randint(0, pool_size, ens_size)   # shape: (100,)

        # SAXS optimization on ensemble
        if flags['saxs']:
            sse_saxs, total_score_saxs, bc_saxs = saxs_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_saxs, total_score_saxs, bc_saxs = [None, 0, None]

        # CS optimization on ensemble
        if flags['cs']:
            sse_cs, total_score_cs, bc_cs = cs_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_cs, total_score_cs, bc_cs = [None, 0, None]

        # FRET optimization on ensemble
        if flags['fret']:
            sse_fret, total_score_fret, bc_fret = fret_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_fret, total_score_fret, bc_fret = [None, 0, None]

        # JC optimization on ensemble
        if flags['jc']:
            sse_jc, total_score_jc, bc_alpha_vals_jc, _  = jc_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_jc, total_score_jc, bc_alpha_vals_jc, _ = [None, 0, None, None]

        # NOEs optimization on ensemble
        if flags['noe']:
            sse_noe, total_score_noe, bc_dist_vals_noe = noe_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_noe, total_score_noe, bc_dist_vals_noe = [None, 0, None]

        # PREs optimization on ensemble
        if flags['pre']:
            sse_pre, total_score_pre, bc_pre  = pre_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_pre, total_score_pre, bc_pre = [None, 0, None]

        # RDCs optimization on ensemble
        if flags['rdc']:
            sse_rdc, total_score_rdc, bc_rdc  = rdc_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_rdc, total_score_rdc, bc_rdc = [None, 0, None]

        # RH optimization on ensemble
        if flags['rh']:
            sse_rh, total_score_rh, bc_rh  = rh_optimization_ensemble(exp_data, bc_data, indices)
        else:
            sse_rh, total_score_rh, bc_rh = [None, 0, None]

        # Fianl maximization
        max_score_SAXS, opt_SSE_SAXS, max_score_CS, opt_SSE_CS, max_score_FRET, opt_SSE_FRET, max_score_JC, opt_SSE_JC, \
        max_score_NOEs, opt_SSE_NOE, max_score_PREs, opt_SSE_PRE, max_score_RDCs, opt_SSE_RDC, max_score_RH, opt_SSE_RH, \
        accepted, new_indices, best_jcoups = maximize_score_SAXS(exp_data, bc_data, ens_size, indices, pool_size, flags,
                                                                 total_score_saxs, total_score_cs, total_score_fret,
                                                                 total_score_jc,
                                                                 total_score_noe, total_score_pre, total_score_rdc,
                                                                 total_score_rh,
                                                                 bc_saxs, bc_cs, bc_fret, bc_alpha_vals_jc,
                                                                 bc_dist_vals_noe, bc_pre, bc_rdc, bc_rh)
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

