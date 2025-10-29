import re

def parse_formula_ox(formula_ox_string):
    """
    Convert a WORM formula_ox string into a dictionary of element oxidation states and their quantities.
    For example, an input of "2Fe+3 Fe+2 4O-2" corresponding to magnetite would return the following:
    `{'Fe+3': 2.0, 'Fe+2': 1.0, 'O-2': 4.0}`.

    Parameters
    ----------
    formula_ox_string : str
        The formula_ox string, e.g, from the WORM thermodynamic database column called formula_ox

    Returns
    -------
    out : dict
        A dictionary where each key represents an element in a specific
        oxidation state, and its value is the number of that element in the
        chemical species' formula.
    """
    
    split_list = formula_ox_string.split()
    split_list_clean = [s.replace(" ", "") for s in split_list]
    
    try:
        elem_ox_names = [re.findall(r"^(?:\d+|)([A-Z].*$)", s)[0] for s in split_list_clean]
    except:
        elem_ox_names = []
    
    elem_ox_list = []
    for s in split_list:
        coeff = re.findall(r"(\d+)[A-Z]", s)
        if len(coeff) == 0:
            coeff = 1
        else:
            coeff = float(coeff[0])
        elem_ox_list.append(coeff)
    
    return {key:val for key,val in zip(elem_ox_names, elem_ox_list)}