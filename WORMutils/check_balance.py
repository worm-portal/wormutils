from chemparse import parse_formula

def check_balance(formulas, stoich):
    """
    Check that a chemical reaction is balanced. If not, get missing composition.
    
    Parameters
    ----------
    formulas : list of str
        A list of species formulas that match the order of
        the stoichiometric reaction coefficients in the `stoich` parameter.
    
    stoich : list of numeric
        A list of stoichiometric reaction coefficients that match the order of
        the species formulas in the `formulas` parameter. Reactants are
        negative.
        
    Returns
    -------
    A printed warning and a dictionary of the missing composition if the
    reaction is unbalanced.
    """
    
    if len(formulas) != len(stoich):
        raise Exception("The number of species formulas does not match the "
              "number of stoichiometric coefficients in the reaction.")
    
    # sum all elements, +, and - by their reaction coefficient
    all_dict = {}
    for i,s in enumerate(formulas):
        s_dict = parse_formula(s)
        s_dict = {key: stoich[i]*s_dict[key] for key in s_dict.keys()}
        all_dict = {k: all_dict.get(k, 0) + s_dict.get(k, 0) for k in set(all_dict) | set(s_dict)}
    
    # sum + and - as Z (charge)
    if "+" not in list(all_dict.keys()):
        all_dict["+"] = 0
    if "-" not in list(all_dict.keys()):
        all_dict["-"] = 0
    all_dict["Z"] = all_dict["+"] - all_dict["-"]
    del all_dict["+"]
    del all_dict["-"]
    
    # delete all elements with a value of 0 (balanced)
    for key in list(all_dict.keys()):
        if all_dict[key] == 0:
            del all_dict[key]
    
    # print warnings, prepare missing composition dictionary
    if len(list(all_dict.keys())) > 0:
        missing_composition_dict = {k:[-all_dict[k]] for k in all_dict.keys()}
        print("Warning! The reaction is unbalanced. It is missing this composition:")
        print(pd.DataFrame(missing_composition_dict).to_string(index=False))
    else:
        missing_composition_dict = {}
    
    return missing_composition_dict