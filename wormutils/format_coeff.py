import math

def format_coeff(coeff):
    
    """
    Format a reaction coefficient for Plotly/html display.
    """
    if coeff == 1 or coeff == -1:
        coeff = ""
    elif coeff.is_integer() and coeff < 0:
        coeff = str(-int(coeff))
    elif coeff.is_integer() and coeff > 0:
        coeff = str(int(coeff))
    else:
        if coeff < 0:
            coeff = _float_to_formatted_fraction(-coeff)
        else:
            coeff = _float_to_formatted_fraction(coeff)

    if coeff != "":
        coeff = coeff + " "

    return coeff

    
def _float_to_formatted_fraction(x, error=0.000001):
    
    """
    Format a fraction for html.
    """
    f = _float_to_fraction(x, error=error)
    
    whole_number_float = int((f[0]-(f[0]%f[1]))/f[1])
    remainder_tuple = (f[0]%f[1], f[1])
    
    if remainder_tuple[0] == 0:
        return str(whole_number_float)
    else:
        if whole_number_float == 0:
            whole_number_float = ""
        return "{0}<sup>{1}</sup>&frasl;<sub>{2}</sub>".format(
                whole_number_float, remainder_tuple[0], remainder_tuple[1])


def _float_to_fraction (x, error=0.000001):
    
    """
    Convert a float into a fraction. Works with floats like 2.66666666.
    Solution from https://stackoverflow.com/a/5128558/8406195
    """
    n = int(math.floor(x))
    x -= n
    if x < error:
        return (n, 1)
    elif 1 - error < x:
        return (n+1, 1)

    # The lower fraction is 0/1
    lower_n = 0
    lower_d = 1
    # The upper fraction is 1/1
    upper_n = 1
    upper_d = 1
    while True:
        # The middle fraction is (lower_n + upper_n) / (lower_d + upper_d)
        middle_n = lower_n + upper_n
        middle_d = lower_d + upper_d
        # If x + error < middle
        if middle_d * (x + error) < middle_n:
            # middle is our new upper
            upper_n = middle_n
            upper_d = middle_d
        # Else If middle < x - error
        elif middle_n < (x - error) * middle_d:
            # middle is our new lower
            lower_n = middle_n
            lower_d = middle_d
        # Else middle is our best fraction
        else:
            return (n * middle_d + middle_n, middle_d)