"""This is a sample script to run eisd with our local data."""
import numpy as np
np.random.seed(91)
import os
import pandas as pd

from eisd.utils import meta_data
from eisd.parser import read_data
from eisd.utils import make_pairs
from eisd.optimizer import main



if __name__ == '__main__':

    # parameters
    data_path = "../eisd-pkg" # path to experimental data and structure pools (james' compilation)
    structure = 'ensemble'  # ['trades', 'trades_uf', 'mixed', 'ensemble']
    # run_mode: How do you want to run the eisd optimizer?
    #   - all: optimize on all experimental observables together
    #   - single: optimize on individual experimental observable
    #   - dual: optimize on pairs of experimental observables
    # run_mode = 'dual'
    opt_type = 'max' # 'max', 'mc', None
    beta = 0.1
    # in the new_prop_error directory the following properties are affected: CS, FRET, PRE, SAXS

    run_mode = 'pre_indices'
    pre_indices_mode = 'dual'  # 'dual' or 'signle'
    # pre_indices_path = 'local/newrun_2020/new_prop_errors2/positive_mc/%s/single_mode/opt_%s/'%(structure,str(opt_type))
    pre_indices_path = 'local/newrun_2020/new_prop_errors2/positive_mc/%s/dual_mode/opt_%s/'%(structure,str(opt_type))

    # read files
    filenames = meta_data(data_path)
    exp_data = read_data(filenames['exp'], mode='exp')
    bc_data = read_data(filenames[structure], mode=structure)

    # run_mode: all
    if run_mode == 'all':
        if opt_type == 'mc':
            abs_output = "local/newrun_2020/new_prop_errors2/positive_mc/%s/all_mode/opt_%s_b%s"%(structure,str(opt_type), str(beta))
        else:
            abs_output = "local/newrun_2020/new_prop_errors2/positive_mc/%s/all_mode/opt_%s"%(structure,str(opt_type))

        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        # main(exp_data, bc_data, epochs=1000, mode='all', optimization=True, output_dir=abs_output, verbose=False)
        main(exp_data, bc_data, epochs=1000, mode='all', beta=beta, opt_type=opt_type, output_dir=abs_output, verbose=False)

    # run_mode: dual
    elif run_mode == 'dual':
        # pairs = make_pairs()
        # pairs = [['saxs', 'jc'],['cs', 'jc'], ['fret', 'jc'],['jc', 'noe'],['jc', 'pre'],['jc', 'rdc'],['jc', 'rh']]
        pairs = [['fret', 'noe'], ['fret', 'pre'], ['fret', 'rdc'], ['fret', 'rh']]
        for pair in pairs:
            if opt_type == 'mc':
                abs_output = "local/newrun_2020/new_prop_errors2/positive_mc/%s/dual_mode/opt_%s_b%s/%s_%s"%(structure, str(opt_type), str(beta) ,pair[0], pair[1])
            else:
                abs_output = "local/newrun_2020/new_prop_errors2/positive_mc/%s/dual_mode/opt_%s/%s_%s"%(structure, str(opt_type),pair[0], pair[1])

            if not os.path.exists(abs_output):
                os.makedirs(abs_output)
            main(exp_data, bc_data, epochs=1000, mode=pair, beta=beta, opt_type=opt_type, output_dir=abs_output, verbose=False)

    # run_mode: single
    elif run_mode == 'single':
        # single_modes = ['saxs', 'cs', 'fret', 'jc', 'noe', 'pre', 'rdc', 'rh']
        single_modes = ['saxs', 'cs', 'fret']
        for mode in single_modes:
            if opt_type == 'mc':
                abs_output = "local/newrun_2020/new_prop_errors2/positive_mc/%s/single_mode/opt_%s_b%s/%s"%(structure, str(opt_type),str(beta),mode)
            else:
                abs_output = "local/newrun_2020/new_prop_errors2/positive_mc/%s/single_mode/opt_%s/%s"%(structure, str(opt_type), mode)

            if not os.path.exists(abs_output):
                os.makedirs(abs_output)

            main(exp_data, bc_data, epochs=1000, mode=mode, beta=beta, opt_type=opt_type, output_dir=abs_output, verbose=True)

    # run_mode: pre_indices; to calculate scores/RMSDs for the remaining unoptimized properties
    elif run_mode == 'pre_indices':
        if pre_indices_mode == 'single':
            for prop in ['saxs', 'cs', 'jc', 'fret', 'noe', 'pre', 'rh', 'rdc']:
                ind_path = os.path.join(pre_indices_path, prop, 'indices.csv')
                pre_indices = pd.read_csv(ind_path, header=None).values

                abs_output = os.path.join(pre_indices_path, prop, 'all_processed')
                print(abs_output)
                if not os.path.exists(abs_output):
                    os.makedirs(abs_output)
                main(exp_data, bc_data, pre_indices=pre_indices, output_dir=abs_output)

        elif pre_indices_mode == 'dual':
            pairs = make_pairs()
            for pair in pairs:
                ind_path = os.path.join(pre_indices_path, '%s_%s'%(pair[0], pair[1]), 'indices.csv')
                pre_indices = pd.read_csv(ind_path, header=None).values

                abs_output = os.path.join(pre_indices_path, '%s_%s'%(pair[0], pair[1]), 'all_processed')
                print(abs_output)
                if not os.path.exists(abs_output):
                    os.makedirs(abs_output)
                main(exp_data, bc_data, pre_indices=pre_indices, output_dir=abs_output)