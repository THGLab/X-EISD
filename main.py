import numpy as np
import pandas as pd
import os

abs_path = "../eisd-pkg"

# experimental data file names
relative_path = os.path.join(abs_path, "experimental_data")
EXP_DATA_FILENAME = {
    'rh'  : None,  #20.3
    'rdc' : os.path.join(relative_path, "drksh3_exp_rdcs.txt"),
    'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
    'noe' : os.path.join(relative_path, "8AAC_noes.txt"),
    'jc'  : os.path.join(relative_path, "drksh3_JC_exp_clean.txt"),
    'fret': None, # 0.55
    'cs'  : os.path.join(relative_path, "drksh3_CS_exp_mod.txt"),
    'saxs': os.path.join(relative_path, "unfolded_saxs_exp.txt")
}

# back calculated data file names
relative_path = os.path.join(abs_path, "back_calc_data/Trades")
TRADES_BC_DATA_FILENAME = {
    'rh'  : os.path.join(relative_path, "drksh3_RH_data.txt"),
    'rdc' : os.path.join(relative_path, "drksh3_RDCs_structRDCs.txt"),
    'pre' : os.path.join(relative_path, "drksh3_PREs_structdists.txt"),
    'noe' : os.path.join(relative_path, "drksh3_NOEs_structdists.txt"),
    'jc'  : os.path.join(relative_path, "drksh3_JC_alphas.txt"),
    'fret': os.path.join(relative_path, "drksh3_FRET_structEFF.txt"),
    'cs'  : os.path.join(relative_path, "drksh3_CS_structdata.txt"),
    'saxs': os.path.join(relative_path, "drksh3_SAXS_data.txt")
}

relative_path = os.path.join(abs_path, "back_calc_data/mixpool")
MIXED_BC_DATA_FILENAME = {
    'rh'  : os.path.join(relative_path, "drksh3_mixpool_RH.txt"),
    'rdc' : os.path.join(relative_path, "drksh3_mixpool_RDCs.txt"),
    'pre' : os.path.join(relative_path, "drksh3_mixpool_PREs.txt"),
    'noe' : os.path.join(relative_path, "drksh3_mixpool_NOEs.txt"),
    'jc'  : os.path.join(relative_path, "drksh3_mixpool_JC.txt"),
    'fret': os.path.join(relative_path, "drksh3_mixpool_FRET.txt"),
    'cs'  : os.path.join(relative_path, "drksh3_mixpool_CS.txt"),
    'saxs': os.path.join(relative_path, "drksh3_mixpool_SAXS.txt")
}

relative_path = os.path.join(abs_path, "back_calc_data/ENSEMBLE_PED")
ENSEMBLE_BC_DATA_FILENAME = {
    'rh'  : os.path.join(relative_path, "drksh3_ENSEMBLE_RH.txt"),
    'rdc' : os.path.join(relative_path, "drksh3_ENSEMBLE_RDCs.txt"),
    'pre' : os.path.join(relative_path, "drksh3_ENSEMBLE_PREs.txt"),
    'noe' : os.path.join(relative_path, "drksh3_ENSEMBLE_NOEs.txt"),
    'jc'  : os.path.join(relative_path, "drksh3_ENSEMBLE_JC.txt"),
    'fret': os.path.join(relative_path, "drksh3_ENSEMBLE_FRET.txt"),
    'cs'  : os.path.join(relative_path, "drksh3_ENSEMBLE_CS.txt"),
    'saxs': os.path.join(relative_path, "drksh3_ENSEMBLE_SAXS.txt")
}

class Stack():
    def __init__(self, name, data, sigma=None, mu=None):
        self.name = name
        self.data = data
        self.sigma = sigma
        self.mu = mu

    def get_idx(self):
        pass




def read_data(filenames, mode):
    """
    The main function to read all the back calculated files

    Parameters
    ----------
    filenames: dict
        This parameter is a dictionary of properties with their relative path to the data file.

    mode: str
        This parameter must be one of the following:
            - 'exp': experimental data
            - 'trades': back calculated data for TRADES pool
            - 'mixed': back calculated data for MIXED pool
            - 'ensemble': back calculated data for ENSEMBLE pool

    Returns
    -------
    dict: A dictionary of properties with their pandas data frame

    """
    if mode == 'exp':
        return {
            # property name: Stack(name, exp data, sigma, mu)
            'rh'  : Stack('rh', 20.3, 0.3, None),
            'rdc' : Stack('rdc', pd.read_csv(filenames['rdc']), None, None),
            'pre' : Stack('pre', pd.read_csv(filenames['pre']), None, None),
            'noe' : Stack('noe', pd.read_csv(filenames['noe']), None, None),
            'jc'  : Stack('jc', pd.read_csv(filenames['jc']), None, None),
            'fret': Stack('fret', 0.55, 0.02, None),
            'cs'  : Stack('cs', pd.read_csv(filenames['cs']), None, None),
            'saxs': Stack('saxs', pd.read_csv(filenames['saxs']), None, None)
        }

    elif mode == 'bc':
        pre_data = pd.read_csv(filenames['pre'], header=None, delim_whitespace=True, index_col=0)   # shape: None, 68
        pre_data.index = range(pre_data.shape[0]) # just wanted to keep the indices 0-indexed

        noe_data = pd.read_csv(filenames['noe'], header=None, delim_whitespace=True, index_col=0)   # shape: None, 93
        noe_data.index = range(noe_data.shape[0]) # just wanted to keep the indices 0-indexed

        jc_data = pd.read_csv(filenames['jc'], header=None, delim_whitespace=True, index_col=0)   # shape: None, 47
        jc_data.index = range(noe_data.shape[0]) # just wanted to keep the indices 0-indexed

        jc_data = pd.read_csv(filenames['jc'], header=None, delim_whitespace=True, index_col=0)   # shape: None, 47
        jc_data.index = range(noe_data.shape[0]) # just wanted to keep the indices 0-indexed

        return {
            # property name: Stack(name, back calc data, sigma, mu)
            'rh': Stack('rh', pd.read_csv(filenames['rh'], header=None),  # shape: None, 1
                        0.812, None),
            'rdc': Stack('rdc',
                         pd.read_csv(filenames['rdc'], header=None, delim_whitespace=True, index_col=0), # shape: None, 28
                         0.88, None),
            'pre': Stack('pre', pre_data, 0.0001, None),
            'noe': Stack('noe', noe_data, 0.0001, None),
            'jc': Stack('jc', jc_data,
                        (np.sqrt(0.14), np.sqrt(0.03), np.sqrt(0.08)),
                        (6.51, -1.76, 1.6)),
            'fret': Stack('fret', 0.55, 0.02, None),
            'cs': Stack('cs', pd.read_csv(filenames['cs']), None, None),
            'saxs': Stack('saxs', pd.read_csv(filenames['saxs']), None, None)

        }


