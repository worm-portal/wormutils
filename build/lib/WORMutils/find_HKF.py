import pandas as pd

def find_HKF(Gh=float('NaN'), V=float('NaN'), Cp=float('NaN'),
             Gf=float('NaN'), Hf=float('NaN'), Saq=float('NaN'),
             Z=float('NaN'), a1=float('NaN'), a2=float('NaN'), a3=float('NaN'),
             a4=float('NaN'), c1=float('NaN'), c2=float('NaN'),
             omega=float('NaN'), organic=False, volatile=False,
             HKF_scale=True, DEW=False, phase_TrPr=None, aq_complex=False,
             print_eq=False):
    
    """
    Estimate HKF parameters from standard state thermodynamic properties of an
    aqueous organic molecule.
    
    Parameters
    ----------
    Gh : numeric
        Standard state partial molal Gibbs free energy of hydration in cal/mol.
    
    V : numeric
        Standard state partial molal volume in cm3/mol.
    
    Cp : numeric
        Standard state partial molal heat capacity in cal/mol/K.
    
    Gf : numeric
        Standard state partial molal Gibbs free energy of formation in cal/mol.
    
    Hf : numeric
        Standard state partial molal enthalpy of formation in cal/mol.
    
    Saq : numeric
        Standard state partial molal third law entropy in cal/mol/K.
    
    Z : numeric
        The net charge of the molecule.

    a1, a2, a3, a4, c1, c2, omega : numeric, optional
        Parameters for the revised Helgeson Kirkham Flowers (HKF) equation of
        state. If these are not provided, they will be estimated using published
        correlation methods.

    organic : bool, default True
        Is this molecule organic? If so, correlations from Shock and Helgeson
        1990 will be used to estimate the value of the HKF parameter omega. If
        the molecule is not organic, correlations from Shock and Helgeson 1988
        will be used to obtain omega instead.

    volatile : bool, default False
        Is this molecule volatile? If volatile=True and organic=True,
        equation 60 from Shock and Helgeson 1990 will be used to estimate the
        HKF parameter omega. If volatile=False and organic=True, equation 61
        will be used instead.

    HKF_scale : bool, default True
        Should the output contain scaled HKF parameters according to the common
        convention?

    DEW : bool, default False
        Estimate HKF parameters according to Sverjensky et al. 2014? If
        so, provides compatibility with the Deep Earth Water (DEW) model.

    phase_TrPr : str, optional
        Required for estimating HKF equation of state parameters for neutral
        species using the DEW model. What is the phase of the species at 25 °C
        and 1 bar when not dissolved in water? Can be "cr", "gas", or "liq".

    aq_complex : bool, default False
        Determines whether the estimated a1 parameter will be representative of
        an aqueous complex for the sake of the Deep Earth Water (DEW) model. If
        True, equation 129 from Sverjensky 2019 will be used. If False, equation
        8 in Appendix 1 of Sverjensky et al. 2014 will be used.
        
    Returns
    ----------
    out_dict : dict
        A dictonary of properties and parameters.
    """

    eqs_used = [] # stores a list of strings that describes which equations were used in what order

    Tr = 298.15 # K
    theta = 228 # K
    
    # Born constants
    Y = -5.802E-05 # 1/K
    Q = 5.903E-07 # (1/bar)
    X = -3.09*10**-7 # 1/K**2
    
    pfunk = 2601
    conv = 41.8393

    eta = 1.66027*10**5 # angstroms*cal/mol
    eqs_used.append("eta = {} angstroms*cal/mol, Y = {} 1/K, Q = {} (1/bar), X = {} 1/K**2".format(eta, Y, Q, X))

    # define abs_protonBorn, mentioned in text after Eq 47 in Shock and Helgeson 1988
    abs_protonBorn = 0.5387 * 10**5
    eqs_used.append("abs_protonBorn = 0.5387 * 10**5 cal/mol, mentioned in text after Eq 47 in Shock and Helgeson 1988")

    # estimate omega if it isn't already supplied by the user
    if pd.isnull(omega):

        if not pd.isnull(Gh) and Z == 0:
            eqs_used.append("Gh is provided and charge equals zero so estimate omega from Plyasunov and Shock 2001...")

            Gh = Gh*4.184/1000 # convert Gh to kJ/mol temporarily due to the convention of using Joules in Plyasunov and Shock 2001
            omega = (2.61+(324.1/(Gh-90.6)))/10**-5
            eqs_used.append("omega = {} J/mol = (2.61+(324.1/(Gh-90.6)))/10**-5, Eq 8 in Plyasunov and Shock 2001".format("{0:.5g}".format(omega)))
            omega = omega/4.184 # convert to calorie-based units
            Gh = Gh*1000/4.184 # convert to calorie-based units
    
        elif Z == 0 and not DEW:
            eqs_used.append("Gh is not provided and charge equals zero so estimate omega for neutral solutes from Shock and Helgeson 1990...")

            if organic:
                if volatile:
                    omega = -1514.4*Saq
                    eqs_used.append("omega = {} cal/mol = -1514.4*Saq, Eq 60 in Shock and Helgeson 1990".format("{0:.5g}".format(omega)))
                else:
                    omega = -1514.4*Saq + 0.34*10**5
                    eqs_used.append("omega = {} cal/mol = -1514.4*Saq + 0.34*10**5, Eq 61 in Shock and Helgeson 1990".format("{0:.5g}".format(omega)))
            else:
                # TODO: this assumption is for metal complexes. What about inorganic neutral species that are not complexes?
                omega = -0.038E5
                eqs_used.append("omega = -0.038E5 cal/mol, Eq 42 in Sverjensky et al 1997".format("{0:.5g}".format(omega)))
        
        elif DEW:
            if Z == 0 and not isinstance(phase_TrPr, str):
                msg = ("Warning: please specify the phase of this species at 25 "
                       "degrees C and 1 bar when not dissolved in water by "
                       "setting the phase_TrPr parameter to either 'cr',"
                       "'gas', or 'liq'. Assuming the phase is 'cr' for now...")
                print(msg)
                phase_TrPr = "cr"
            
            if Z == 0 and phase_TrPr in ["c", "cr"]:
                omega = 0.3E5 # eq 125 in Sverjensky 2019
            elif Z == 0 and phase_TrPr in ["g", "gas", "l", "liq"]:
                omega = -0.3E5 # eq 126 in Sverjensky 2019
            else:
                if Z > 0:
                    Bz = 0.544E5*Z # Eq 123 in Sverjensky 2019
                else:
                    Bz = 1.62E5*Z # Eq 124 in Sverjensky 2019
                
                omega = -1514.4*Saq + Bz # Eq 122 in Sverjensky 2019, Eq 58 in Shock and Helgeson 1990
        
        elif Z != 0 and not DEW:
            eqs_used.append("Gh is not provided and charge does not equal zero so estimate omega for ionic species from Shock and Helgeson 1990...")
                
            # define alphaZ (described in text after Eq 59 in Shock and Helgeson 1990)
            if (abs(Z) == 1):
                alphaZ = 72
            elif (abs(Z) == 2):
                alphaZ = 141
            elif (abs(Z) == 3):
                alphaZ = 211
            elif (abs(Z) == 4):
                alphaZ = 286
            else:
                alphaZ = float('NaN')
            if print_eq and alphaZ != float('NaN'):
                eqs_used.append("alphaZ = {} because charge = {}, described in text after Eq 59 in Shock and Helgeson 1990".format(alphaZ, Z))
                
            # define BZ
            BZ = ((-alphaZ*eta)/(Y*eta - 100)) - Z * abs_protonBorn  # Eq 55 in Shock and Helgeson 1990
            eqs_used.append("BZ = {} = ((-alphaZ*eta)/(Y*eta - 100)) - Z * abs_protonBorn, Eq 55 in Shock and Helgeson 1990".format("{0:.5g}".format(BZ)))
    
            if organic:
                omega = -1514.4*Saq + BZ  # Eq 58 in Shock and Helgeson 1990
                eqs_used.append("omega = {} cal/mol = -1514.4*Saq + BZ, Eq 58 in Shock and Helgeson 1990".format("{0:.5g}".format(omega)))
    
            else:
                ### METHOD FOR INORGANIC AQUEOUS ELECTROLYTES USING SHOCK AND HELGESON 1988:
                rej = (Z**2 *(eta * Y - 100)/(Saq - alphaZ))  # Eqs 46+56+57 in Shock and Helgeson 1988
                eqs_used.append("rej = {} = (Z**2 *(eta * Y - 100)/(Saq - alphaZ)), Eqs 46+56+57 in Shock and Helgeson 1988".format("{0:.5g}".format(rej)))
        
                #find ion absolute omega*10**-5
                omega_abs_ion = (eta*(Z**2))/rej # Eq 45 in Shock and Helgeson 1988
                eqs_used.append("omega_abs_ion = {} cal/mol = (eta*(charge**2))/rej, Eq 45 in Shock and Helgeson 1988".format("{0:.5g}".format(omega_abs_ion)))
        
                #find ion omega
                omega = omega_abs_ion-(Z*abs_protonBorn) # Eq 47 in Shock and Helgeson 1988
                eqs_used.append("omega = {} cal/mol = omega_abs_ion-(Z*abs_protonBorn), Eq 47 in Shock and Helgeson 1988".format("{0:.5g}".format(omega)))
    
        else:
            omega = float('NaN')

    # find delta V solvation (cm3/mol)
    Vs = -omega*Q*conv
    eqs_used.append("Vs = {} cm3/mol = -omega*Q*conv, Eq 5 in Shock and Helgeson 1988, delta V solvation".format("{0:.5g}".format(Vs)))

    Vn = V - Vs
    eqs_used.append("Vn cm3/mol = {} cm3/mol = V - Vs, Eq 4 in Shock and Helgeson 1988, delta V nonsolvation".format("{0:.5g}".format(Vn)))

    sigma = 1.11*Vn + 1.8
    eqs_used.append("sigma cm3/mol = {} cm3/mol = 1.11*Vn + 1.8, Eq 87 in Shock and Helgeson".format("{0:.5g}".format(sigma)))
        
    Cps = omega*Tr*X
    eqs_used.append("Cps = {} cal/mol/K= omega*Tr*X, Eq 35 in Shock and Helgeson 1988, delta Cp solvation".format("{0:.5g}".format(Cps)))
        
    # find delta Cp nonsolvation (cal/mol*K)
    Cpn = Cp - Cps  # Eq 29 in Shock and Helgeson 1988
    eqs_used.append("Cpn = {} cal/mol/K = Cp - Cps, Eq 29 in Shock and Helgeson 1988, delta Cp nonsolvation".format("{0:.5g}".format(Cpn)))

    # calculate a1-a4
    if not pd.isnull(Gh) and Z == 0 and not DEW:
        eqs_used.append("Gh is provided and charge is neutral, so estimate a1, a2, and a4 from Plysunov and Shock 2001")

        # temporarily convert Gh into units of kJ/mol as per the convention of Plyasunov and Shock 2001
        Gh = Gh*4.184/1000

        if pd.isnull(a1):
            a1 = ((0.820-(1.858*10**-3*Gh))*V)/10
            eqs_used.append("a1 = {} J/mol/bar = ((0.820-(1.858*10**-3*Gh))*V)/10, Eq 10 in Plyasunov and Shock 2001".format("{0:.5g}".format(a1)))
        else:
            a1 = a1*4.184
            eqs_used.append("a1 = {} J/mol/bar, supplied by user".format("{0:.5g}".format(a1)))

        if pd.isnull(a2):
            a2 = ((0.648+((0.00481)*(Gh)))*V)/1E-2
            eqs_used.append("a2 = {} J/mol = ((0.648+((0.00481)*(Gh)))*V)/1E-2, Eq 11 in Plyasunov and Shock 2001".format("{0:.5g}".format(a2)))
        else:
            a2 = a2*4.184
            eqs_used.append("a2 = {} J/mol, supplied by user".format("{0:.5g}".format(a2)))

        if pd.isnull(a4):
            a4 = (8.10-(0.746*a2*1E-2)+(0.219*Gh))/1E-4
            eqs_used.append("a4 = {} (J*K)/mol = (8.10-(0.746*a2*1E-2)+(0.219*Gh))/1E-4, Eq 12 in Plyasunov and Shock 2001".format("{0:.5g}".format(a4)))
        else:
            a4 = a4*4.184
            eqs_used.append("a4 = {} (J*K)/mol, supplied by user".format("{0:.5g}".format(a4)))
        
        # convert Gh, a1, a2, and a4 into calorie-based units
        Gh = Gh*1000/4.184
        a1 = a1/4.184
        a2 = a2/4.184
        a4 = a4/4.184
    
    else:
        eqs_used.append("Gh is unavailable and/or charge is not 0")

        if DEW:
            Vs = -omega*(0.05903*10**-5*conv)

            eqs_used.append("Vs = {} cm3/mol = -omega*(0.05903*10**-5*conv), Eq 120 in Sverjensky 2019".format("{0:.5g}".format(Vs)))
            
            Vn = V - Vs
            
            eqs_used.append("Vn = {} cm3/mol = V - Vs, nonsolvation contribution to volume".format("{0:.5g}".format(Vn)))

            sigma = 1.11*Vn + 1.8
            eqs_used.append("sigma = {} cm3/mol = 1.11*Vn + 1.8, Eq 130 in Sverjensky 2019".format("{0:.5g}".format(sigma)))

            if pd.isnull(a1):
                if Z == 0 or aq_complex:
                    a1 = (0.1942*Vn + 1.5204)/10
                    eqs_used.append("a1 = {} cal/mol/bar = 0.1942*Vn + 1.5204, sign-corrected version of Eq 12 in Sverjensky 2019 (sign of y-intercept was flipped)".format("{0:.5g}".format(a1)))
                
                else:
                    a1 = ((0.1304*abs(Z) - 0.0217)*Vn  + (1.4567*abs(Z)+0.6187))/10
                    eqs_used.append("a1 = {} cal/mol/bar = ((0.1304*abs(Z) - 0.0217)*Vn  + (1.4567*abs(Z)+0.6187))/10, Eq 8 in Appendix 1 of Sverjensky et al 2014".format("{0:.5g}".format(a1)))
            else:
                eqs_used.append("a1 = {} cal/mol/bar, supplied by user".format("{0:.5g}".format(a1)))
        
        else:
            if pd.isnull(a1):
                a1 = 1.3684E-2*Vn+0.1765
                eqs_used.append("a1 = {} cal/mol/bar = 1.3684E-2*Vn+0.1765, Eq 85 in Shock and Helgeson 1988".format("{0:.5g}".format(a1)))
            else:
                eqs_used.append("a1 = {} cal/mol/bar, supplied by user".format("{0:.5g}".format(a1)))

        if pd.isnull(a2):
            a2 = (sigma/conv-a1)*pfunk
            eqs_used.append("a2 = {} cal/mol = (sigma/conv-a1)*pfunk, Eq 8 in Shock and Helgeson 1988, rearranged to solve for a2".format("{0:.5g}".format(a2)))
        else:
            eqs_used.append("a2 = {} cal/mol, supplied by user".format("{0:.5g}".format(a2)))

        if pd.isnull(a4):
            a4 = -4.134*a2-27790
            eqs_used.append("a4 = {} (cal*K)/mol = -4.134*a2-27790, Eq 88 in Shock and Helgeson 1988".format("{0:.5g}".format(a4)))
        else:
            eqs_used.append("a4 = {} (cal*K)/mol, supplied by user".format("{0:.5g}".format(a4)))

    # calculate c2
    if not pd.isnull(Gh) and Z == 0:

        # temporarily convert Gh into kJ/mol
        Gh = Gh*4.184/1000

        if pd.isnull(c2):
            c2 = (21.4 + 0.849*Gh)/1E-4
            eqs_used.append("c2 = {} (J*K)/mol = 21.4+(0.849*Gh), Eq 14 in Plyasunov and Shock 2001".format("{0:.5g}".format(c2)))
        else:
            c2 = c2*4.184
            eqs_used.append("c2 = {} (J*K)/mol, supplied by user".format("{0:.5g}".format(c2)))
        
        # convert Gh and c2 into calorie-based units
        Gh = Gh*1000/4.184
        c2 = c2/4.184
        
    else:
        if pd.isnull(c2):
            c2 = (0.2037*Cp - 3.0346)/10**-4
            eqs_used.append("c2 = {} (cal*K)/mol = (0.2037*Cp - 3.0346)*10**-4, Eq 89 in Shock and Helgeson 1988".format("{0:.5g}".format(c2)))
        else:
            eqs_used.append("c2 = {} (cal*K)/mol, supplied by user".format("{0:.5g}".format(c2)))

    if pd.isnull(c1):
        c1 = Cpn-(c2*(1/(Tr-theta))**2)
        eqs_used.append("c1 = {} cal/mol/K = Cpn-((c2*(1/(Tr-theta))**2), Eq 31 in Shock and Helgeson 1988, rearranged to solve for c1".format("{0:.5g}".format(c1)))
    else:
        eqs_used.append("c1 = {} cal/mol/K, supplied by user".format("{0:.5g}".format(c1)))

    if pd.isnull(a3):
        a3 = ((Vn/conv)-a1-a2/pfunk)*(Tr-theta)-(a4/pfunk)
        eqs_used.append("a3 = {} (cal*K)/mol/bar = ((Vn/conv)-a1-a2/pfunk)*(Tr-theta)-(a4/pfunk), after Eq 11 in Shock and Helgeson 1988, rearranged to solve for a3.".format("{0:.5g}".format(a3)))
    else:
        eqs_used.append("a3 = {} (cal*K)/mol/bar, supplied by user".format("{0:.5g}".format(a3)))

    if HKF_scale:
        a1 = a1*10
        a2 = a2*10**-2
        a3 = a3
        a4 = a4*10**-4
        c1 = c1
        c2 = c2*10**-4
        omega = omega*10**-5

    out_dict = {
        "G": Gf,
        "H": Hf,
        "S": Saq,
        "Cp": Cp,
        "V": V,
        "a1": a1,
        "a2": a2,
        "a3": a3,
        "a4": a4,
        "c1": c1,
        "c2": c2,
        "omega": omega,
        "Z": Z,
        "Vs": Vs,
        "Vn": Vn,
        "sigma": sigma}

    return out_dict, eqs_used