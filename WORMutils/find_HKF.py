import pandas as pd
from datetime import datetime

def find_HKF(Gh=float('NaN'), V=float('NaN'), Cp=float('NaN'),
             Gf=float('NaN'), Hf=float('NaN'), Saq=float('NaN'),
             Z=float('NaN'), a1=float('NaN'), a2=float('NaN'), a3=float('NaN'),
             a4=float('NaN'), c1=float('NaN'), c2=float('NaN'),
             omega=float('NaN'), organic=False, organic_acid=False, volatile=False,
             HKF_scale=True, DEW=False, phase_TrPr=None, aq_complex=False,
             name=None, abbrv=None, formula=None, azero=None, formula_ox=None,
             dissrxn=None, tag=None, neutral_ion_type=0,
             wrm_data_output=False, print_eq=False):
    
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

    organic : bool, default False
        Is this molecule organic? If so, correlations from Shock and Helgeson
        1990 or Plyasunov and Shock 2001 will be used to estimate certain HKF
        parameters. If the molecule is not organic, correlations from Shock and
        Helgeson 1988 will be used to obtain parameters instead.

    organic_acid : bool, default False
        Is this molecule an organic acid or acid anion? If so, correlations
        from Shock 1995 will be used to estimate certain HKF parameters, unless
        a Gibbs free energy of hydration is provided, in which Plyasunov and
        Shock 2001 will be used.

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

    name : str, optional
        Name of the compound. Used when `wrm_data_output`is True. If this
        compound is going to be used in conjunction with the WORM database, or
        be used with the AqEquil package, then ensure that the name you choose
        does not contain spaces (like "acetic-acid" or "my-custom-compound"). 
        
    abbrv : str, optional
        Abbreviation of the compound (e.g., "HexOOH" for hexanoic acid). Used
        when `wrm_data_output`is True.
    
    formula : str, optional
        Chemical formula for the compound. Used when `wrm_data_output`is True.
    
    azero : str, optional
         The azero parameter of the aqueous compound. Used when
         `wrm_data_output` is True. If no azero parameter is defined, then one
         will be estimated based on charge.
    
    formula_ox : str, optional
        Quantities of elements and their oxidation states in the compound. For
        example, methane's formula_ox would be 'C-4 4H+' and hexanoic acid's
        would be 'C-3 4C-2 O-2 C- H+'. Used when `wrm_data_output`is True.

    dissrxn : str, optional
        A dissociation reaction compatible with the 'dissrxn' column in the WORM
        database. This can be blank if this compound is meant to be a basis
        species. Used when `wrm_data_output`is True.
    
    tag : str, optional
        A tag compatible with the 'tag' column of the WORM database. Will be
        blank by default, representing a nonbasis species. Used when
        `wrm_data_output`is True.
    
    neutral_ion_type : int, default 0
        A neutral ion type compatible with the 'neutral_ion_type' column in the
        WORM database. The default of 0 means that it is not treated specially
        if the molecule has a charge of 0. Used when `wrm_data_output`is True.
    
    wrm_data_output : bool, default False
        Output the results as a dataframe that is compatible with the WORM
        database and is therefore importable into AqEquil and pyCHNOSZ packages?
        
    print_eq : bool, default False
        Print equations used in estimation? Equations are printed in the order
        they are calculated.
        
    Returns
    ----------
    hkf : dict or pandas.DataFrame
        If `wrm_data_output` is False, returns a dictonary of properties and
        parameters. Otherwise, returns a Pandas dataframe with a format that
        can be imported as thermodynamic data into AqEquil or pyCHNOSZ packages.

    eq : list of strings
        A list of steps used to perform the estimation. 
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

            if organic_acid and Z==0:
                omega = 661.98*Saq - 58740
                eqs_used.append("omega = {} cal/mol = 661.98*Saq - 58740, Eq 22 in Shock 1995".format("{0:.5g}".format(omega)))
            elif organic:
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
                
            if organic:
                eqs_used.append("Gh is not provided, charge does not equal zero, and species is organic so estimate omega for ionic species from Shock and Helgeson 1990...")
                
                # define BZ
                BZ = ((-alphaZ*eta)/(Y*eta - 100)) - Z * abs_protonBorn  # Eq 55 in Shock and Helgeson 1990
                eqs_used.append("BZ = {} = ((-alphaZ*eta)/(Y*eta - 100)) - Z * abs_protonBorn, Eq 55 in Shock and Helgeson 1990".format("{0:.5g}".format(BZ)))

                omega = -1514.4*Saq + BZ  # Eq 58 in Shock and Helgeson 1990
                eqs_used.append("omega = {} cal/mol = -1514.4*Saq + BZ, Eq 58 in Shock and Helgeson 1990".format("{0:.5g}".format(omega)))
    
            else:
                eqs_used.append("Gh is not provided, charge does not equal zero, and species is inorganic so estimate omega for ionic species from Shock and Helgeson 1988...")

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

    Cps = omega*Tr*X
    eqs_used.append("Cps = {} cal/mol/K= omega*Tr*X, Eq 35 in Shock and Helgeson 1988, delta Cp solvation".format("{0:.5g}".format(Cps)))
        
    # find delta Cp nonsolvation (cal/mol*K)
    Cpn = Cp - Cps  # Eq 29 in Shock and Helgeson 1988
    eqs_used.append("Cpn = {} cal/mol/K = Cp - Cps, Eq 29 in Shock and Helgeson 1988, delta Cp nonsolvation".format("{0:.5g}".format(Cpn)))

    # calculate a1-a4
    if not pd.isnull(Gh) and Z == 0 and not DEW:
        eqs_used.append("Gh is provided and charge is neutral, so estimate a1, a2, and a4 from Plysunov and Shock 2001")

        sigma = float("NaN")
        eqs_used.append("estimation of sigma is not required using this method from Plyasunov and Shock 2001...")
        
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
            eqs_used.append("sigma = {} cm3/mol = 1.11*Vn + 1.8, Eq 130 in Sverjensky 2019 (Eq 87 in Shock and Helgeson 1988)".format("{0:.5g}".format(sigma)))

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
            if organic_acid:
                sigma = 1.07143*Vn + 3.0
                eqs_used.append("sigma cm3/mol = {} cm3/mol = 1.07143*Vn + 3.0, Eq 23 in Shock 1995".format("{0:.5g}".format(sigma)))
                
            elif organic and Z==0:
                sigma = 1.0125*Vn
                eqs_used.append("sigma cm3/mol = {} cm3/mol = 1.0125*Vn, Eq 63 in Shock and Helgeson 1990".format("{0:.5g}".format(sigma)))
    
            else:
                sigma = 1.11*Vn + 1.8
                eqs_used.append("sigma cm3/mol = {} cm3/mol = 1.11*Vn + 1.8, Eq 87 in Shock and Helgeson 1988".format("{0:.5g}".format(sigma)))
        
            
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
            if organic_acid and Z==0:
                c2 = (0.0988*Cp - 4.961)/10**-4
                eqs_used.append("c2 = {} (cal*K)/mol = (0.0988*Cp - 4.961)/10**-4, Eq 28 in Shock 1995".format("{0:.5g}".format(c2)))
            elif organic_acid and Z != 0:
                c2 = (0.01212*Cp - 4.106)/10**-4
                eqs_used.append("c2 = {} (cal*K)/mol = (0.01212*Cp - 4.106)/10**-4, Eq 29 in Shock 1995".format("{0:.5g}".format(c2)))
            elif organic and Z==0:
                c2 = (0.0676*Cp - 4.054)/10**-4
                eqs_used.append("c2 = {} (cal*K)/mol = (0.0676*Cp - 4.054)/10**-4, Eq 64 in Shock and Helgeson 1990".format("{0:.5g}".format(c2)))
            else:
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

    hkf = {
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
        "sigma": sigma,
        "organic": organic,
        "organic_acid": organic_acid}



    if not wrm_data_output:
        return hkf, eqs_used
        
    else:
        if name != None:
            hkf["name"] = name
        if abbrv != None:
            hkf["abbrv"] = abbrv
        if formula != None:
            hkf["formula"] = formula
        if formula_ox != None:
            hkf["formula_ox"] = formula_ox
        if azero != None:
            hkf["azero"] = azero
        if dissrxn != None:
            hkf["dissrxn"] = dissrxn
        if tag != None:
            hkf["tag"] = tag
        if neutral_ion_type != None:
            hkf["neutral_ion_type"] = neutral_ion_type
        
        this_date = datetime.today().strftime('%Y%m%d') 
        
        def _est_azero(Z):
            azero = 4
            if Z == 2:
                azero = 6
            elif Z == 3:
                azero = 9
            elif Z == 4:
                azero = 11

            eqs_used.append("azero estimated as {}".format(azero)) # todo: track down where this comes from
            return azero
        
        if 'azero' in list(hkf.keys()):
            if hkf['azero'] == None:
                azero = _est_azero(hkf['Z'])
            else:
                azero = hkf['azero']
        else:
            azero = _est_azero(hkf['Z'])
        
        if hkf['organic'] or hkf['organic_acid']:
            cat_1 = "organic_aq"
        else:
            cat_1 = "inorganic_aq"
        
        if 'dissrxn' in list(hkf.keys()):
            if hkf['dissrxn'] != None:
                dissrxn = hkf['dissrxn']
            else:
                dissrxn = ""
        else:
            dissrxn = ""
        
        if 'tag' in list(hkf.keys()):
            if hkf['tag'] != None:
                tag = hkf['tag']
            else:
                tag = ""
        else:
            tag = ""
        
        if 'formula_ox' in list(hkf.keys()):
            if hkf['formula_ox'] != None:
                formula_ox = hkf['formula_ox']
            else:
                formula_ox = ""
        else:
            formula_ox = ""
        
        data={
            "name": [hkf['name']],
            "abbrv": [hkf['abbrv']],
            "formula": [hkf['formula']],
            "state": ["aq"],
            "ref1": ["findHKF"],
            "ref2": [""],
            "date": [this_date],
            "E_units": ["cal"],
            "G": [hkf['G']],
            "H": [hkf['H']],
            "S": [hkf['S']],
            "Cp": [hkf['Cp']],
            "V": [hkf['V']],
            "a1.a": [hkf['a1']],
            "a2.b": [hkf['a2']],
            "a3.c": [hkf['a3']],
            "a4.d": [hkf['a4']],
            "c1.e": [hkf['c1']],
            "c2.f": [hkf['c2']],
            "omega.lambda": [hkf['omega']],
            "z.T":[hkf['Z']],
            "azero":[azero],
            "neutral_ion_type":[hkf['neutral_ion_type']],
            "dissrxn":[dissrxn],
            "tag":[tag],
            "formula_ox":[formula_ox],
            "category_1":[cat_1],
            "category_2":[""],
        }
        df = pd.DataFrame(data)
        return df, eqs_used

def find_HKF_test(print_eq=False):
    
    """
    Test the HKF estimation function by regenerating published values.
    
    Parameters
    ----------
    print_eq : bool, default False
        Print equations used in estimation?
    """

    print("SELECT ITEMS FROM SHOCK AND HELGESON 1988, TABLE 12\n---------------------------------------------\n")

    print("Be+2\n---------")
    print("Input parameters:")
    print("find_HKF(Gf=-83500, Hf=-91500, Saq=-55.7, Cp=-1.3, V=-25.4, Z=2, organic=False)\n")
    out, eq = find_HKF(Gf=-83500, Hf=-91500, Saq=-55.7, Cp=-1.3, V=-25.4, Z=2, organic=False)
    pub = {"omega":"1.9007", "a1":"-1.0684", "a2":"-10.3901",
           "a3":"9.8338", "a4":"-2.3495", "c1":"22.9152", "c2":"-3.2994"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"], 4)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"], 4)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"], 4)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 4)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"], 4)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"], 4)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"], 4)))
    print("")
    for e in eq:
        print(e)
    print("")

    print("S2O6-2\n---------")
    print("Input parameters:")
    print("find_HKF(Gf=--231000, Hf=-280400, Saq=30, Cp=-46.5, V=43.3, Z=-2, organic=False)\n")
    out, eq = find_HKF(Gf=--231000, Hf=-280400, Saq=30, Cp=-46.5, V=43.3, Z=-2, organic=False)
    pub = {"omega":"2.7587", "a1":"8.6225", "a2":"13.2724",
           "a3":"0.5334", "a4":"-3.3277", "c1":"4.3301", "c2":"-12.5066"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"], 4)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"], 4)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"], 4)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 4)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"], 4)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"], 4)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"], 4)))
    print("")
    for e in eq:
        print(e)
    print("")

    print("SELECT ITEMS FROM SHOCK AND HELGESON 1990, TABLE 6\n---------------------------------------------\n")
    
    print("1-hexanamine\n---------")
    print("Input parameters:")
    print("find_HKF(Gf=14860, Hf=-46320, Saq=60.2, Cp=144, V=121.6, Z=0, organic=True)\n")
    out, eq = find_HKF(Gf=14860, Hf=-46320, Saq=60.2, Cp=144, V=121.6, Z=0, organic=True)
    pub = {"omega":"-0.5717", "a1":"18.2115", "a2":"28.2822",
           "a3":"12.6611", "a4":"-3.9481", "c1":"127.1903", "c2":"5.6804"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"], 4)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"], 4)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"], 4)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 4)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"], 4)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"], 4)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"], 4)))
    print("")
    for e in eq:
        print(e)
    print("")

    print("n-hexylbenzene\n---------")
    print("Input parameters:")
    print("find_HKF(Gf=40390, Hf=-25590, Saq=76.2, Cp=208.4, V=177, Z=0, organic=True)\n")
    out, eq = find_HKF(Gf=40390, Hf=-25590, Saq=76.2, Cp=208.4, V=177, Z=0, organic=True)
    pub = {"omega":"-0.8140", "a1":"25.7106", "a2":"43.2732",
           "a3":"13.8894", "a4":"-4.5678", "c1":"180.5115", "c2":"10.0338"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"], 4)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"], 4)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"], 4)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 4)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"], 4)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"], 4)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"], 4)))
    print("")
    for e in eq:
        print(e)
    print("")
    
    print("SELECT ITEMS FROM SHOCK 1995, TABLE 4\n---------------------------------------------\n")

    print("hexanoic acid\n---------")
    print("Input parameters:")
    print("find_HKF(Gf=-87120, Hf=-139290, Saq=69.5, Cp=125.0, V=116.55, Z=0, organic=True, organic_acid=True)\n")
    out, eq = find_HKF(Gf=-87120, Hf=-139290, Saq=69.5, Cp=125.0, V=116.55, Z=0, organic=True, organic_acid=True)
    pub = {"omega":"-0.1266", "a1":"17.6709", "a2":"33.3251",
           "a3":"-2.9700", "a4":"-4.1566", "c1":"108.8183", "c2":"7.3890"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"], 4)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"], 4)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"], 4)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 4)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"], 4)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"], 4)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"], 4)))
    print("")
    for e in eq:
        print(e)
    print("")
    
    print("hexanoate\n---------")
    print("Input parameters:")
    print("ind_HKF(Gf=-80490, Hf=-139870, Saq=45.3, Cp=90, V=102.21, Z=-1, organic=True, organic_acid=True)\n")
    out, eq = find_HKF(Gf=-80490, Hf=-139870, Saq=45.3, Cp=90, V=102.21, Z=-1, organic=True, organic_acid=True)
    pub = {"omega":"0.9427", "a1":"16.0700", "a2":"29.6995",
           "a3":"-2.1530", "a4":"-4.0067", "c1":"104.8115", "c2":"-3.0151"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"], 2)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"], 4)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"], 4)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 4)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"], 4)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"], 4)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"], 4)))
    print("")
    for e in eq:
        print(e)
    print("")

    
    print("PLYASUNOV AND SHOCK 2001, TABLE 4\n---------------------------------------------")
    print("Input parameters published in the table are converted from kJ or J into calorie-based units for find_HKF()\n")
    
    print("SO2\n---------")
    print("Input parameters:")
    print("Gh=-0.51/4.184*1000, V=39.0, Cp=146/4.184, Z=0, organic=True\n")
    out, eq = find_HKF(Gh=-0.51/4.184*1000, V=39.0, Cp=146/4.184, Z=0, organic=True)
    pub = {"omega":"-0.95", "a1":"32.02", "a2":"25.17",
           "a3":"18.71", "a4":"-10.79", "c1":"93.2", "c2":"20.97"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"]*4.184, 1)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"]*4.184, 2)))
    print("")
    for e in eq:
        print(e)
    print("")
    
    print("Pyridine\n---------")
    print("Input parameters:")
    print("Gh=-11.7/4.184*1000, V=77.1, Cp=306/4.184, Z=0, organic=True\n")
    out, eq = find_HKF(Gh=-11.7/4.184*1000, V=77.1, Cp=306/4.184, Z=0)
    pub = {"omega":"-0.56", "a1":"64.89", "a2":"45.62",
           "a3":"69.94", "a4":"-28.50", "c1":"278.1", "c2":"11.47"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"], 2)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"]*4.184, 1)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"]*4.184, 2)))
    print("")
    for e in eq:
        print(e)
    print("")
    
    print("1,4-Butanediol\n---------")
    print("Input parameters:")
    print("Gh=-37.7/4.184*1000, V=88.23, Cp=347/4.184, Z=0, organic=True\n")
    out, eq = find_HKF(Gh=-37.7/4.184*1000, V=88.23, Cp=347/4.184, Z=0, organic=True)
    pub = {"omega":"0.08", "a1":"78.50", "a2":"41.17",
           "a3":"76.32", "a4":"-30.87", "c1":"369.2", "c2":"-10.61"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"]*4.184, 1)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"]*4.184, 2)))
    print("")
    for e in eq:
        print(e)
    print("")
    
    print("beta-alanine\n---------")
    print("Input parameters:")
    print("Gh=-74/4.184*1000, V=58.7, Cp=76/4.184, Z=0, organic=True\n")
    out, eq = find_HKF(Gh=-74/4.184*1000, V=58.7, Cp=76/4.184, Z=0, organic=True)
    pub = {"omega":"0.64", "a1":"56.17", "a2":"17.14",
           "a3":"54.55", "a4":"-20.90", "c1":"165.5", "c2":"-41.43"}
    print("Published: {}, \tCalculated: {}, \tomega*10**-5".format(pub["omega"], round(out["omega"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta1*10".format(pub["a1"], round(out["a1"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta2*10**-2".format(pub["a2"], round(out["a2"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta3".format(pub["a3"], round(out["a3"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \ta4*10**-4".format(pub["a4"], round(out["a4"]*4.184, 2)))
    print("Published: {}, \tCalculated: {}, \tc1".format(pub["c1"], round(out["c1"]*4.184, 1)))
    print("Published: {}, \tCalculated: {}, \tc2*10**-4".format(pub["c2"], round(out["c2"]*4.184, 2)))
    print("")
    for e in eq:
        print(e)
    print("")