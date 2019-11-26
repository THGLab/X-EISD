import os

from eisd.meta import meta_data
from eisd.parser import read_data
from eisd.utils import make_paris
from eisd.optimizer import main


if __name__ == '__main__':

    data_path = "../eisd-pkg" # path to experimental data and structure pools (temporary)
    filenames = meta_data(data_path)
    exp_data = read_data(filenames[0], mode='exp')
    bc_data = read_data(filenames[1], mode='trades_uf')

    # run_mode: How do you want to run the eisd?
    #   - all: optimize on all experimental observables together
    #   - single: optimize on individual experimental observable
    #   - dual: optimize on pairs of experimental observables
    run_mode = 'single'

    # all
    if run_mode == 'all':
        abs_output = "newrun/unfolded/all"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        output_file = os.path.join(abs_output, 'drksh3_folded_trades_all.txt')
        main(exp_data, bc_data, epochs=250, mode='all', output_file=output_file, verbose=True)

    # dual mode
    elif run_mode == 'dual':
        abs_output = "newrun/unfolded/dual_mode"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        pairs = make_paris()
        for pair in pairs[::-1]:
            output_file = os.path.join(abs_output, "drksh3_unfolded_trades_%s_%s.txt"%(pair[0], pair[1]))
            main(exp_data, bc_data, epochs=1000, mode=pair, output_file=output_file, verbose=True)

    # single mode
    elif run_mode == 'single':
        abs_output = "newrun/unfolded/single_mode"
        if not os.path.exists(abs_output):
            os.makedirs(abs_output)

        # single_modes = ['saxs', 'cs', 'fret', 'jc', 'noe', 'pre', 'rdc', 'rh']
        single_modes = ['jc']
        for mode in single_modes[::-1]:
            output_file = os.path.join(abs_output, "drksh3_unfolded_trades_%s.txt"%mode)
            main(exp_data, bc_data, epochs=1000, mode=mode, output_file=output_file, verbose=True)