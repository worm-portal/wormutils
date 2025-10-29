import numpy as np

def assign_worm_db_col_dtypes(
        df,
        float_cols=["G", "H", "S", "Cp",
                    "V", "a1.a", "a2.b",
                    "a3.c", "a4.d", "c1.e",
                    "c2.f", "omega.lambda", "z.T",
                    "azero", "neutral_ion_type", "regenerate_dissrxn",
                    "logK1", "logK2", "logK3", "logK4",
                    "logK5", "logK6", "logK7", "logK8",
                    "T1", "T2", "T3", "T4", "T5", "T6",
                    "T7", "T8"],
        str_cols=["name", "abbrv", "state", "formula",
                  "ref1", "ref2", "date",
                  "E_units", "tag", "dissrxn", "formula_ox",
                  "formula_modded", "formula_ox_modded", 
                  "P1", "P2", "P3", "P4", "P5", "P6",
                  "P7", "P8"],
        NA_string=""):

    df = df.map(lambda x: float("NaN") if x == NA_string else x)
    
    for col in float_cols:
        if col in df.columns:
            df[col] = df[col].astype(float)
    for col in str_cols:
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df