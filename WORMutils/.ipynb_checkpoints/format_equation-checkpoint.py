from .chemlabel import chemlabel

def format_equation(species, stoich, charge_sign_at_end=False):
    """
    Format a chemical equation to display in HTML
    (e.g., Plotly plots)
    
    Parameters
    ----------
    species : list of str
        List of species in the reaction
        
    stoich : list of numeric
        List of stoichiometric reaction coefficients (reactants are negative)
    
    charge_sign_at_end : bool, default False
        Display charge with sign after the number (e.g. SO4 2-)?
        
    
    Returns
    -------
    A formatted chemical formula string.
    """
    reactants_list = []
    products_list = []
    for i,s in enumerate(species):
        s_f = chemlabel(s, charge_sign_at_end=charge_sign_at_end)
        if stoich[i] < 0:
            if stoich[i] != -1:
                entry = str(abs(stoich[i])) + " " + s_f
            else:
                entry = s_f
            reactants_list.append(entry)
        elif stoich[i] > 0:
            if stoich[i] != 1:
                entry = str(stoich[i]) + " " + s_f
            else:
                entry = s_f
            products_list.append(entry)
    
    reactants_together = " + ".join(reactants_list)
    products_together = " + ".join(products_list)
    
    equation_str = " → ".join([reactants_together, products_together])
    
    return equation_str