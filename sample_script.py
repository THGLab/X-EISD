"""This is a sample script to run eisd with our local data."""
import numpy as np
np.random.seed(90)
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
    run_mode = 'dual'

    # read files
    filenames = meta_data(data_path)
    exp_data = read_data(filenames['exp'], mode='exp')
    bc_data = read_data(filenames[structure], mode=structure)

    # run_mode: all
    if run_mode == 'all':
        abs_output = "newrun/unfolded/all"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        output_file = os.path.join(abs_output, 'drksh3_folded_trades_all.txt')
        main(exp_data, bc_data, epochs=250, mode='all', output_file=output_file, verbose=True)

    # run_mode: dual
    elif run_mode == 'dual':
        abs_output = "newrun/unfolded/dual_mode"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        # pairs = make_pairs()
        # pairs = [['saxs', 'jc'],['cs', 'jc'], ['fret', 'jc'],['jc', 'noe'],['jc', 'pre'],['jc', 'rdc'],['jc', 'rh']]
        pairs = [['saxs', 'noe'], ['cs', 'noe'], ['cs', 'pre'], ['cs', 'rdc'],['fret', 'noe'], ['noe', 'pre'], ['noe', 'rdc']]
        for pair in pairs[4:]:
            output_file = os.path.join(abs_output, "drksh3_unfolded_trades_%s_%s.txt"%(pair[0], pair[1]))
            main(exp_data, bc_data, epochs=1000, mode=pair, output_file=output_file, verbose=False)

    # run_mode: single
    elif run_mode == 'single':
        abs_output = "newrun/unfolded/single_mode"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        # single_modes = ['saxs', 'cs', 'fret', 'jc', 'noe', 'pre', 'rdc', 'rh']
        single_modes = ['jc']
        for mode in single_modes:
            output_file = os.path.join(abs_output, "drksh3_unfolded_trades_%s.txt"%mode)
            main(exp_data, bc_data, epochs=1000, mode=mode, output_file=output_file, verbose=True)
