import re

def chemlabel(name, charge_sign_at_end=False):
    
    """
    Format a chemical formula to display subscripts and superscripts in HTML
    (e.g., Plotly plots)
    Example, "CH3COO-" becomes "CH<sub>3</sub>COO<sup>-</sup>"
    
    Parameters
    ----------
    name : str
        A chemical formula.
    
    charge_sign_at_end : bool, default False
        Display charge with sign after the number (e.g. SO4 2-)?
        
    
    Returns
    -------
    A formatted chemical formula string.
    """
    
    # format only the first part of the name if it has "_(input)"
    if len(name.split("_(input)"))==2:
        if name.split("_(input)")[1] == '':
            name = name.split("_(input)")[0]
            input_flag=True
    else:
        input_flag = False
    
    name = _html_chemname_format(name, charge_sign_at_end=charge_sign_at_end)
    
    # add " (input)" to the end of the name
    if input_flag:
        name = name+" (input)"
    
    return(name)


def _html_chemname_format(name, charge_sign_at_end=False):
    
    """
    Function duplicated from pyCHNOSZ
    """
    
    p = re.compile(r'(?P<sp>[-+]\d*?$)')
    name = p.sub(r'<sup>\g<sp></sup>', name)
    charge = re.search(r'<.*$', name)

    name_no_charge = re.match(r'(?:(?!<|$).)*', name).group(0)
    mapping = {"0": "<sub>0</sub>", "1": "<sub>1</sub>", "2": "<sub>2</sub>",
               "3": "<sub>3</sub>", "4": "<sub>4</sub>", "5": "<sub>5</sub>",
               "6": "<sub>6</sub>", "7": "<sub>7</sub>", "8": "<sub>8</sub>",
               "9": "<sub>9</sub>", ".":"<sub>.</sub>"}
    name_no_charge_formatted = "".join([mapping.get(x) or x
                                        for x in list(name_no_charge)])

    if charge != None:
        name = name_no_charge_formatted + charge.group(0)
    else:
        name = name_no_charge_formatted

    if charge_sign_at_end:
        if "<sup>-" in name:
            name = name.replace("<sup>-", "<sup>")
            name = name.replace("</sup>", "-</sup>")
        if "<sup>+" in name:
            name = name.replace("<sup>+", "<sup>")
            name = name.replace("</sup>", "+</sup>")

    return(name)