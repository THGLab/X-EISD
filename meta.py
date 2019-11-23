import os


def meta_data():
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

    return EXP_DATA_FILENAME, TRADES_BC_DATA_FILENAME, MIXED_BC_DATA_FILENAME, ENSEMBLE_BC_DATA_FILENAME