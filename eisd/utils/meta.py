import os

def meta_data(abs_path):
    """
    This function provides the path to the data that is used in our papers.
    Todo: replace with a download module after publication.

    Parameters
    ----------
    abs_path: str
        The path tho the directory with the main subdirectories.

    Returns
    -------
    dict: dictionary of dictionaries for file names to:
        - exp: eight experimental data
        - trades: eight back calc data for Trades pool
        - tades_uf: eight back calc data for Trades pool (same as trades just wanted to embarrassingly distinguish keys)
        - mixed: eight back calc data for Mixed pool
        - ensemble: eight back calc data for Ensemble pool

    """
    # experimental data file names
    relative_path = os.path.join(abs_path, "experimental_data")
    EXP_DATA_FILENAME = {
        'rh'  : None,
        'rdc' : os.path.join(relative_path, "drksh3_exp_rdcs.txt"),
        'pre' : os.path.join(relative_path, "drksh3_pres.txt"),
        'noe' : os.path.join(relative_path, "8AAC_noes.txt"),
        'jc'  : os.path.join(relative_path, "drksh3_JC_exp_clean.txt"),
        'fret': None,
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

    return {'exp':EXP_DATA_FILENAME, 'trades': TRADES_BC_DATA_FILENAME, 'trades_uf': TRADES_BC_DATA_FILENAME,
            'mixed': MIXED_BC_DATA_FILENAME, 'ensemble': ENSEMBLE_BC_DATA_FILENAME}
