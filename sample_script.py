"""This is a sample script to run eisd with our local data."""
import numpy as np
np.random.seed(91)
import os

from eisd.utils import meta_data
from eisd.parser import read_data
from eisd.utils import make_pairs
from eisd.optimizer import main



if __name__ == '__main__':

    # parameters
    data_path = "../eisd-pkg" # path to experimental data and structure pools (james' compilation)
    structure = 'trades_uf'  # ['trades', 'trades_uf', 'mixed', 'ensemble']
    # run_mode: How do you want to run the eisd optimizer?
    #   - all: optimize on all experimental observables together
    #   - single: optimize on individual experimental observable
    #   - dual: optimize on pairs of experimental observables
    run_mode = 'all'

    # read files
    filenames = meta_data(data_path)
    exp_data = read_data(filenames['exp'], mode='exp')
    bc_data = read_data(filenames[structure], mode=structure)

    # run_mode: all
    if run_mode == 'all':
        abs_output = "newrun/unfolded/all_unoptimized"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        # main(exp_data, bc_data, epochs=1000, mode='all', optimization=True, output_dir=abs_output, verbose=False)
        main(exp_data, bc_data, epochs=1000, mode='all', optimization=False, output_dir=abs_output, verbose=False)

    # run_mode: dual
    elif run_mode == 'dual':
        # pairs = make_pairs()
        # pairs = [['saxs', 'jc'],['cs', 'jc'], ['fret', 'jc'],['jc', 'noe'],['jc', 'pre'],['jc', 'rdc'],['jc', 'rh']]
        pairs = [['saxs', 'noe'], ['cs', 'noe'], ['cs', 'pre'], ['cs', 'rdc'],['fret', 'noe'], ['noe', 'pre'], ['noe', 'rdc']]
        for pair in pairs:
            abs_output = "temp/unfolded/dual_mode/%s_%s"%(pair[0], pair[1])
            if not os.path.exists(abs_output):
                os.makedirs(abs_output)
            main(exp_data, bc_data, epochs=1000, mode=pair, output_dir=abs_output, verbose=False)

    # run_mode: single
    elif run_mode == 'single':
        # single_modes = ['saxs', 'cs', 'fret', 'jc', 'noe', 'pre', 'rdc', 'rh']
        single_modes = ['jc']
        for mode in single_modes:
            abs_output = "newrun/unfolded/single_mode/%s"%mode
            if not os.path.exists(abs_output):
                os.makedirs(abs_output)

            main(exp_data, bc_data, epochs=1000, mode=mode, output_dir=abs_output, verbose=True)
