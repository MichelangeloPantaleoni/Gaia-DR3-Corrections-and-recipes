import numpy as np
import pandas as pd
import astropy.units as u
from scipy import interpolate
from astropy.coordinates import SkyCoord

def correct_parallax(df, mode = 'Apellaniz'):
    '''
    Computes de parallax zero point bias in the Gaia DR3 and uses it to correct
    the parallax, using Maíz Apellaniz 2022 [2022A%26A...657A.130M] recipe or the
    one in Lindegren et al. 2021 [2021A%26A...649A...4L].

    The "mode" parameter allows to select the recipe for the correction;
    'Apellaniz' (default) or 'Lindegren'.

    Inputs a pandas dataframe (df) with these Gaia DR3 parameters:
        phot_g_mean_mag
        parallax
        nu_eff_used_in_astrometry
        pseudocolour
        ecl_lat
        astrometric_params_solved
    Outputs the same dataframe but including these columns:
        'Plx_zero'     => The zero point parallax bias in milliarcseconds
        'Plx_corr'     => The corrected parallax in milliarcseconds

    If phot_g_mean_mag is not present then it assumes the star is a faint one (G = 21.0 mag)
    and computes the parallax zero point with that information.

    Values for phot_g_mean_mag outside the range 6.0 - 21.0 mag are computed as the closest
    boundary of the range, and values inside the range are linearly interpolated.

    Objects with 2-parameter solutions (astrometric_params_solved = 3) have no parallax
    measurements and thus the function retrieves nan for the zero point and the corrected
    parallax.

    Objects with 5-parameter solutions (astrometric_params_solved = 31) are the ones that
    show differences in regards to the choosen recipe.

    Objects with 6-parameter solutions (astrometric_params_solved = 95) follow the same
    recipe in Lindegren and in Apellaniz.
    '''
    gmag, lat, npar = df[['phot_g_mean_mag', 'ecl_lat', 'astrometric_params_solved']].values.T
    gmag, lat = gmag.astype(float), lat.astype(float)
    gmag[np.where(np.isnan(gmag))[0]] = 21.0 # Stars with non-valid G are assumed to be faint
    nueff = (df['nu_eff_used_in_astrometry'].fillna(0)+df['pseudocolour'].fillna(0)).values

    # For objects with 6-parameter solutions (Lindegren and Apellaniz do the same)
    gcut6 = [  6.0 ,  10.8 ,  11.2 ,  11.8 ,  12.2 ,  12.9 ,  13.1 ,  15.9 ,  16.1 ,  17.5 ,  19.0 ,  20.0 ,  21.0 ]
    q600  = [-27.85, -28.91, -26.72, -29.04, -12.39, -18.99, -38.29, -36.83, -28.37, -24.68, -15.32, -13.73, -29.53]
    q601  = [ -7.78,  -3.57,  -8.74,  -9.69,  -2.16,  -1.93,  +2.59,  +4.20,  +1.99,  -1.37,  +4.01, -10.92, -20.34]
    q602  = [+27.47, +22.92,  +9.36, +13.63, +10.23, +15.90, +16.20, +15.76,  +9.28,  +3.52,  -6.03,  -8.30, -18.74]
    q610  = [-32.1 ,  +7.7 , -30.3 , -49.4 , -92.6 , -57.2 , -10.5 , +22.3 , +50.4 , +86.8 , +29.2 , -74.4 , -39.5 ]
    q611  = [+14.4 , +12.6 ,  +5.6 , +36.3 , +19.8 ,  -8.0 ,  +1.4 , +11.1 , +17.2 , +19.8 , +14.1 ,+196.4 ,+326.8 ]
    q612  = [ +9.5 ,  +1.6 , +17.2 , +17.7 , +27.6 , +19.9 ,  +0.4 , +10.0 , +13.7 , +21.3 ,  +0.4 , -42.0 ,-262.3 ]
    q620  = [  -67.,  -572., -1104., -1129.,  -365.,  -554.,  -960., -1367., -1351., -1380.,  -563.,  +536., +1598.]

    # For objects with 5-parameter solutions:
    if mode == 'Apellaniz':
        # Maíz Apellaniz 2022 [2022A%26A...657A.130M]
        gcut5 = [  6.0 ,   7.4 ,   9.2 ,  10.8 ,  11.2 ,  11.8 ,  12.2 ,  12.9 ,  13.1 ,  15.9 ,  16.1 ,  17.5 ,  18.0 ,  19.0 ,  20.0 ,  21.0 ]
        q500  = [-54.33,  -8.17, -27.11, -20.03, -33.38, -36.95, -16.99, -26.49, -37.51, -32.82, -33.18, -25.02, -22.88, -18.40, -12.65, -18.22]
        q501  = [-11.97, -10.06,  -7.60,  -2.78, -12.21, -11.55,  -3.67, -10.63,  +3.33,  +5.15,  -1.31,  +5.83,   0.40,  +5.98,  -4.57, -15.24]
        q502  = [+25.39, +24.12, +22.48, +10.48,  +5.51,  -1.65, +15.81, +21.77, +20.41,  +6.50,  +5.48,  +6.57,  -5.10,  -6.46,  -7.46, -18.54]
        q510  = [-31.1 , -13.4 ,  +9.3 , +11.6 ,-132.3 ,-158.7 ,-109.9 , -76.0 ,  -2.9 ,  -9.10, -56.8 , -39.2 , -46.5 ,   0   ,   0   ,   0   ]
        q511  = [+19.1 , +23.7 , +29.6 , +34.8 , -10.2 , +13.2 , +63.2 , +43.2 , +29.6 , +12.2 , -38.1 , -29.1 , -35.4 ,  +5.5 , +97.9 ,+128.2 ]
        q520  = [-2529., -2529., -2529., -2529., -2529., -2529., -3625., -4353., -1675., -1341., -1705., -1284.,  -896.,   0   ,   0   ,   0   ]
        q530  = [  0   ,   0   ,   0   ,   0   ,   0   ,   0   ,   0   ,   0   , +32.1 ,+168.0 ,+112.1 ,+196.3 ,+126.5 ,   0   ,   0   ,   0   ]
        q540  = [-358.1,-358.1 ,-358.1 , -78.6 ,+203.1 ,-155.3 ,-144.2 , +23.6 , +99.5 ,+129.3 ,+153.1 ,+218.0 ,+190.2 ,+276.6 ,   0   ,   0   ]
    elif mode == 'Lindegren':
        # Lindegren et al. 2021 [2021A%26A...649A...4L]
        gcut5 = [  6.0 ,  10.8 ,  11.2 ,   11.8 ,   12.2 ,  12.9 ,   13.1 ,   15.9 ,   16.1 ,   17.5 ,   19.0 ,  20.0 ,   21.0 ]
        q500  = [-26.98, -27.23, -30.33,  -33.54,  -13.65, -19.53,  -37.99,  -38.33,  -31.05,  -29.18,  -18.40, -12.65,  -18.22]
        q501  = [ -9.62,  -3.07,  -9.23,  -10.08,   -0.07,  -1.64,   +2.63,   +5.61,   +2.83,   -0.09,   +5.98,  -4.57,  -15.24]
        q502  = [+27.40, +23.04,  +9.08,  +13.28,   +9.35, +15.86,  +16.14,  +15.42,   +8.59,   +2.41,   -6.46,  -7.46,  -18.54]
        q510  = [-25.1 , +35.3 , -88.4 , -126.7 , -111.4 , -66.8 ,   -5.7 ,    0   ,    0   ,    0   ,    0   ,   0   ,    0   ]
        q511  = [ -0.0 , +15.7 , -11.8 ,  +11.6 ,  +40.6 , +20.6 ,  +14.0 ,  +18.7 ,  +15.5 ,  +24.5 ,   +5.5 , +97.9 , +128.2 ]
        q520  = [-1257., -1257., -1257.,  -1257.,  -1257., -1257.,  -1257.,  -1189.,  -1404.,  -1165.,    0   ,   0   ,    0   ]
        q530  = [  0   ,   0   ,   0   ,    0   ,    0   ,   0   , +107.9 , +243.8 , +105.5 , +189.7 ,    0   ,   0   ,    0   ]
        q540  = [  0   ,   0   ,   0   ,    0   ,    0   ,   0   , +104.3 , +155.2 , +170.7 , +325.0 , +276.6 ,   0   ,    0   ]
    else:
        raise ValueError("  Wrong 'mode': You have to choose either 'Apellaniz' or 'Lindegren'.  ")

    c1 = -0.24*(nueff <= 1.24).astype(int)+(nueff-1.48)*((nueff > 1.24) & (nueff <= 1.72)).astype(int)+0.24*(nueff > 1.72).astype(int)
    c2 = +0.24**3*(nueff <= 1.24).astype(int)+(1.48-nueff)**3*((nueff > 1.24) & (nueff <= 1.48)).astype(int)
    c3 = (nueff-1.24)*(nueff <= 1.24).astype(int)
    c4 = (nueff-1.72)*(nueff > 1.72).astype(int)
    b1 = np.sin(np.deg2rad(lat))
    b2 = np.sin(np.deg2rad(lat))*np.sin(np.deg2rad(lat))-1/3
    qq500 = np.interp(gmag, gcut5, q500)
    qq501 = np.interp(gmag, gcut5, q501)
    qq502 = np.interp(gmag, gcut5, q502)
    qq510 = np.interp(gmag, gcut5, q510)
    qq511 = np.interp(gmag, gcut5, q511)
    qq520 = np.interp(gmag, gcut5, q520)
    qq530 = np.interp(gmag, gcut5, q530)
    qq540 = np.interp(gmag, gcut5, q540)
    qq600 = np.interp(gmag, gcut6, q600)
    qq601 = np.interp(gmag, gcut6, q601)
    qq602 = np.interp(gmag, gcut6, q602)
    qq610 = np.interp(gmag, gcut6, q610)
    qq611 = np.interp(gmag, gcut6, q611)
    qq612 = np.interp(gmag, gcut6, q612)
    qq620 = np.interp(gmag, gcut6, q620)
    z5 = qq500 + qq501*b1 + qq502*b2 + qq510*c1 + qq511*c1*b1 + qq520*c2 + qq530*c3 + qq540*c4
    z6 = qq600 + qq601*b1 + qq602*b2 + qq610*c1 + qq611*c1*b1 + qq612*c1*b2 + qq620*c2
    zpvals = z5*(npar == 31).astype(int)+z6*(npar == 95).astype(int)+999999*(npar == 3).astype(int) # Returns 999999 for npar = 3
    df['Plx_zero'] = zpvals/1000 # Change from microarcsec to milliarcsec
    df['Plx_zero'] = df['Plx_zero'].replace(999999/1000, np.nan) # No parallax <=> npar = 3 => NaNs
    df['Plx_corr'] = df['parallax']-df['Plx_zero'] # Correct parallax with the zero point bias
    return df


def correct_parallax_error(df):
    '''
    Corrects parallax uncertainties in Gaia DR3.
    Adapted from Maíz Apellaniz 2022 [2022A%26A...657A.130M], Appendix A.

    Inputs a pandas dataframe (df) with these Gaia DR3 parameters:
        phot_g_mean_mag
        parallax_error
        ruwe
        astrometric_params_solved
    Outputs the same dataframe but including the column:
        'Plx_err_corr' => The corrected parallax uncertainty

    If phot_g_mean_mag is not present then it assumes the star is a faint one (G = 21.0 mag)
    and computes the parallax zero point with that information.
    '''
    plx_err, gmag, ruwe, npar = df[['parallax_error', 'phot_g_mean_mag', 'ruwe', 'astrometric_params_solved']].values.T
    gmag, plx_err = gmag.astype(float), plx_err.astype(float)
    gmag[np.where(np.isnan(gmag))[0]] = 21.0 # Stars with non-valid G are assumed to be faint
    gref = np.array([ 6.50,  7.50,  8.50,  9.50, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75,
                     13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25, 16.75, 17.25, 17.75])
    kref = np.array([2.62, 2.38, 2.06, 1.66, 1.22, 1.41, 1.76, 1.76, 1.90, 1.92,
                     1.61, 1.50, 1.39, 1.35, 1.24, 1.20, 1.19, 1.18, 1.18, 1.14])
    reduced_gmag = (gmag < 6)*6+(gmag > 18)*18+((gmag <= 18) & (gmag >= 6))*gmag
    tck = interpolate.splrep(gref, kref, k = 1)
    k = interpolate.splev(reduced_gmag, tck, der = 0)
    geref = [6.00, 12.50, 13.50, 14.50, 15.50, 16.50, 17.50]
    keref = [0.50,  0.50,  1.01,  1.28,  1.38,  1.44,  1.32]
    new_tck = interpolate.splrep(geref, keref, k = 1)
    new_k = interpolate.splev(reduced_gmag, new_tck, der = 0)
    k = k*(1+((new_k <= 0.5)*0.5+(new_k > 0.5)*new_k)*(ruwe > 1.4))
    k = k*(1+(npar == 95*np.ones(len(npar)))*0.25)
    df['Plx_err_corr'] = np.sqrt((k*plx_err)**2+0.0103**2)
    return df


def calculate_angular_covariance(theta):
    '''
    Calculates the angular covariance of the parallax in Gaia DR3, between two
    sources separated by an angle theta in the sky.
    Adapted from Maíz, Pantaleoni & Barbá 2021 [2021A%26A...649A..13M].
    '''
    if isinstance(theta, float) or isinstance(theta, int):
        theta = np.array([theta])
        not_a_list = True
    else:
        not_a_list = False
    a, b, lmda, phi = 0.6, 0.94, 1.05, -(5/18)*np.pi
    V_LMC_0 = 46.2 # μas^2
    V_LMC = V_LMC_0*(a*np.exp(-theta/1.4)+(1-a)*(np.cos((2*np.pi*theta/lmda)+phi)/np.cos(phi))*(b*np.exp(-(theta/0.35)**0.8)+(1-b)))
    V_QSO = np.zeros(len(theta))
    for i in range(len(theta)):
        if theta[i] < 3.0:
            V_QSO[i] = 60.0
        elif theta[i] >= 3.0 and theta[i] <= 40.0:
            V_QSO[i] = 1.62*(40-theta[i])
        elif theta[i] > 40.0 and theta[i] <= 140.0:
            V_QSO[i] = 0
        else:
            V_QSO[i] = 1.62*(140-theta[i])
    V = V_LMC+V_QSO
    if not_a_list:
        V = V[0]
    return V


def calculate_group_parallax(ra, dec, Plx_corr, Plx_err_corr):
    '''
    Calculates the weighted average parallax of a group of stars.
    Adapted from Maíz, Pantaleoni & Barbá 2021 [2021A%26A...649A..13M].

    Inputs:
        ra, dec; Coordinates in degrees.
        Plx_corr, Plx_err_corr; Corrected parallax and its uncertainty in mas.
    Outputs:
        Plx_group, Plx_err_group; Corrected group parallax and its uncertainty in mas.

    TO-DO list (work-in-progress):
    - Make it work for a Pandas DataFrame structure instead of the single values.
    - Make sure it does the same as the IDL implementation by Jesús (in particular
      in compact groups and mainly for the the group parallax uncertainties).
    '''
    # Calculate Group Parallax
    weights = (1/Plx_err_corr**2)/np.sum(1/Plx_err_corr**2)
    Plx_group = np.sum(weights*Plx_corr)
    # Calculate Group Parallax Uncertainty
    pairwise_term = 0
    coords = SkyCoord(ra = ra*u.deg, dec = dec*u.deg, frame = 'icrs')
    for i in range(len(coords)-1):
        for j in range(i+1, len(coords)):
            sep = coords[i].separation(coords[j]).deg
            pairwise_term += weights[i]*weights[j]*calculate_angular_covariance(sep)
    Plx_err_group = np.sqrt(np.sum((weights**2)*(Plx_err_corr**2))+2*pairwise_term)
    return [Plx_group, Plx_err_group]


def correct_proper_motion(pmra, pmdec, ra, dec, phot_g_mean_mag):
    '''
    Corrects the proper motions in Gaia DR3, for bright (G <= 13) stars.
    Adapted from the code in Cantat-Gaudin & Brandt [2021A%26A...649A.124C] (Appendix A).

    Input : source position , coordinates and G magnitude from Gaia EDR3.
    Output: corrected proper motions.
    '''
    if np.isnan(phot_g_mean_mag) or (phot_g_mean_mag >= 13.0):
        return pmra, pmdec
    else:
        table = [[  0.00,   9.00,   9.50,  10.00, 10.50, 11.00, 11.50, 11.75, 12.00, 12.25, 12.50, 12.75],
                 [  9.00,   9.50,  10.00,  10.50, 11.00, 11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00],
                 [ 18.40,  14.00,  12.80,  13.60, 16.20, 19.40, 21.80, 17.70, 21.30, 25.70, 27.30, 34.90],
                 [ 33.80,  30.70,  31.40,  35.70, 50.00, 59.90, 64.20, 65.60, 74.80, 73.60, 76.60, 68.90],
                 [-11.30, -19.40, -11.80, -10.50,  2.10,  0.20,  1.00, -1.90,  2.10,  1.00,  0.50, -2.90]]
        table = np.array(table)
        Gmin = table[0]
        Gmax = table[1]
        omegaX = table[2][(Gmin <= phot_g_mean_mag) & (Gmax > phot_g_mean_mag)][0]
        omegaY = table[3][(Gmin <= phot_g_mean_mag) & (Gmax > phot_g_mean_mag)][0]
        omegaZ = table[4][(Gmin <= phot_g_mean_mag) & (Gmax > phot_g_mean_mag)][0]
        sin_ra,  cos_ra  = np.sin(np.deg2rad(ra)),  np.cos(np.deg2rad(ra))
        sin_dec, cos_dec = np.sin(np.deg2rad(dec)), np.cos(np.deg2rad(dec))
        pmra_corr = -sin_dec*cos_ra*omegaX-sin_dec*sin_ra*omegaY+cos_dec*omegaZ
        pmdec_corr = sin_ra*omegaX-cos_ra*omegaY
        return pmra-pmra_corr/1000., pmdec-pmdec_corr/1000.


def correct_proper_motion_error(pmra_error, pmdec_error, phot_g_mean_mag, ruwe,
                                astrometric_params_solved, only_internal = False,
                                only_external = False):
    '''
    Corrects the uncertainties in the proper motions from Gaia DR3 (random + systematic).
    It uses the same recipe as for the parallax uncertainty corrections in Maíz Apellaniz 2022
    [2022A%26A...657A.130M], Appendix A, but using a value of 23 µas/a (from Lindegren
    et al. 2021 [2021A%26A...649A...2L]) for the systematic uncertainty.

    If only_internal = True, then it only calculates the correction with the internal (random) uncertainties.
    If only_external = True, then it only calculates the correction with the external (systematic) uncertainties.

    If a transformation from equatorial to Galactic coordinates is expected, the procedure
    should start by using this function with only_internal = True, then make the coordinate
    transformation, and finally use this function again but with only_external = True (which
    adds the systematic uncertainty in quadrature).
    '''
    if only_internal and only_external:
        print("ERROR: only_internal and only_external can't both be True at the same time")
        return None
    if only_external:
        pmra_error_corrected = pmra_error**2+0.023**2
        pmdec_error_corrected = pmdec_error**2+0.023**2
    phot_g_mean_mag = np.where(np.isnan(phot_g_mean_mag), 21.0, phot_g_mean_mag) # Stars with non-valid G are assumed to be faint
    gref = np.array([ 6.50,  7.50,  8.50,  9.50, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75,
                     13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25, 16.75, 17.25, 17.75])
    kref = np.array([2.62, 2.38, 2.06, 1.66, 1.22, 1.41, 1.76, 1.76, 1.90, 1.92,
                     1.61, 1.50, 1.39, 1.35, 1.24, 1.20, 1.19, 1.18, 1.18, 1.14])
    reduced_gmag = (phot_g_mean_mag < 6)*6+(phot_g_mean_mag > 18)*18+((phot_g_mean_mag <= 18) & (phot_g_mean_mag >= 6))*phot_g_mean_mag
    tck = interpolate.splrep(gref, kref, k = 1)
    k = interpolate.splev(reduced_gmag, tck, der = 0)
    geref = [6.00, 12.50, 13.50, 14.50, 15.50, 16.50, 17.50]
    keref = [0.50,  0.50,  1.01,  1.28,  1.38,  1.44,  1.32]
    new_tck = interpolate.splrep(geref, keref, k = 1)
    new_k = interpolate.splev(reduced_gmag, new_tck, der = 0)
    k = k*(1+((new_k <= 0.5)*0.5+(new_k > 0.5)*new_k)*(ruwe > 1.4))
    k = k*(1+(astrometric_params_solved == 95*np.ones(len(astrometric_params_solved)))*0.25)
    if only_internal:
        pmra_error_corrected = k*pmra_error
        pmdec_error_corrected = k*pmdec_error
    else:
        pmra_error_corrected = np.sqrt((k*pmra_error)**2+0.0230**2)
        pmdec_error_corrected = np.sqrt((k*pmdec_error)**2+0.0230**2)
    return pmra_error_corrected, pmdec_error_corrected


def calculate_group_proper_motion(pmra, pmdec, pmra_error, pmdec_error, pmra_pmdec_corr = 0):
    '''
    Calculates the weighted average proper motion of a group of stars.

    TO-DO list (work-in-progress):
    - Implement a calculation for the uncertainty in the group's proper motion,
      possibly taking into account some angular covariance that might be present
      in Gaia DR3.
    - Implement the calculation to obtain the covariance matrix between the group's
      proper motion components.
    - Make it work for a Pandas DataFrame structure instead of the single values.
    '''
    # Calculate Group Proper Motion
    pmra_weights = (1/pmra_error**2)/np.sum(1/pmra_error**2)
    pmdec_weights = (1/pmdec_error**2)/np.sum(1/pmdec_error**2)
    pmra_group = np.sum(pmra_weights*pmra)
    pmdec_group = np.sum(pmdec_weights*pmdec)
    # Calculate Group Proper Motion Uncertainties
    pmra_error_group, pmdec_error_group = [], []
    # Calculate Group Proper Motion Correlations
    pmra_pmdec_corr_group = []
    return [pmra_group, pmdec_group, pmra_error_group, pmdec_error_group, pmra_pmdec_corr_group]


def full_galactic_transform(ra, dec, ra_error, dec_error, parallax, parallax_error,
                            pmra, pmdec, pmra_error, pmdec_error, ra_dec_corr = 0,
                            ra_parallax_corr = 0, ra_pmra_corr = 0, ra_pmdec_corr = 0,
                            dec_parallax_corr = 0, dec_pmra_corr = 0, dec_pmdec_corr = 0,
                            parallax_pmra_corr = 0, parallax_pmdec_corr = 0, pmra_pmdec_corr = 0):
    '''
    Transforms equatorial coordinates and proper motions to Galactic coordinates
    and proper motions by considering the uncertainties and the full covariance
    matrix (with the 10 correlations given by Gaia) between all the astrometric
    parameters. If the correlations are not given it assumes they are zero.

    Created by Michelangelo Pantaleoni González in 2022 (Python) and 2025 (IDL)
    Based on Butkevich & Lindegren 2022; Gaia DR3 Documentation, Section 4.1.7.
    https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html
    '''
    # ICRS coordinates of the north Galactic pole and the longitude of intersecting planes
    ra_pole, dec_pole, l_omega = np.deg2rad(192.85948), np.deg2rad(27.12825), np.deg2rad(32.93192)

    # Rotation matrices
    Rz_lomg = np.array([[ np.cos(-l_omega),  np.sin(-l_omega), 0],
                        [-np.sin(-l_omega),  np.cos(-l_omega), 0],
                        [                0,                 0, 1]])
    Rx_90dec = np.array([[1,                         0,                         0],
                         [0,  np.cos(np.pi/2-dec_pole),  np.sin(np.pi/2-dec_pole)],
                         [0, -np.sin(np.pi/2-dec_pole),  np.cos(np.pi/2-dec_pole)]])
    Rz_a90 = np.array([[ np.cos(ra_pole+np.pi/2),  np.sin(ra_pole+np.pi/2), 0],
                       [-np.sin(ra_pole+np.pi/2),  np.cos(ra_pole+np.pi/2), 0],
                       [                       0,                        0, 1]])
    AgT = Rz_lomg @ Rx_90dec @ Rz_a90 # Eq. 4.61

    # Galactic coordinates
    ra_rad, dec_rad = np.deg2rad(ra), np.deg2rad(dec)
    r_icrs = np.array([np.cos(ra_rad)*np.cos(dec_rad), np.sin(ra_rad)*np.cos(dec_rad), np.sin(dec_rad)]) # Eq. 4.58
    r_gal = AgT.dot(r_icrs) # Eq. 4.60
    l_rad = np.arctan2(r_gal[1], r_gal[0]) # Eq. 4.64
    b_rad = np.arctan2(r_gal[2], np.hypot(r_gal[1], r_gal[0])) # Eq. 4.64
    l = (np.rad2deg(l_rad)+360)%360 # Between 0º and 360º
    b = np.rad2deg(b_rad) # Between -90º and 90º

    # Galactic proper motion
    p_icrs = np.array([-np.sin(ra_rad), np.cos(ra_rad), 0]) # Eq. 4.65
    q_icrs = np.array([-np.cos(ra_rad)*np.sin(dec_rad), -np.sin(ra_rad)*np.sin(dec_rad), np.cos(dec_rad)]) # Eq. 4.65
    pm_icrs = p_icrs*pmra+q_icrs*pmdec # Eq. 4.67
    p_gal = np.array([-np.sin(l_rad), np.cos(l_rad), 0]) # Eq. 4.66
    q_gal = np.array([-np.cos(l_rad)*np.sin(b_rad), -np.sin(l_rad)*np.sin(b_rad), np.cos(b_rad)]) # Eq. 4.66
    pm_gal = AgT.dot(pm_icrs) # Eq. 4.69
    pml = p_gal.dot(pm_gal) # Eq. 4.71
    pmb = q_gal.dot(pm_gal) # Eq. 4.71

    # Build the covariance matrix
    ra_dec_cov = ra_error*dec_error*ra_dec_corr
    ra_parallax_cov = ra_error*parallax_error*ra_parallax_corr
    ra_pmra_cov = ra_error*pmra_error*ra_pmra_corr
    ra_pmdec_cov = ra_error*pmdec_error*ra_pmdec_corr
    dec_parallax_cov = dec_error*parallax_error*dec_parallax_corr
    dec_pmra_cov = dec_error*pmra_error*dec_pmra_corr
    dec_pmdec_cov = dec_error*pmdec_error*dec_pmdec_corr
    parallax_pmra_cov = parallax_error*pmra_error*parallax_pmra_corr
    parallax_pmdec_cov = parallax_error*pmdec_error*parallax_pmdec_corr
    pmra_pmdec_cov = pmra_error*pmdec_error*pmra_pmdec_corr
    C = np.array([[ra_error**2    , ra_dec_cov      , ra_parallax_cov   , ra_pmra_cov      , ra_pmdec_cov      ],
                  [ra_dec_cov     , dec_error**2    , dec_parallax_cov  , dec_pmra_cov     , dec_pmdec_cov     ],
                  [ra_parallax_cov, dec_parallax_cov, parallax_error**2 , parallax_pmra_cov, parallax_pmdec_cov],
                  [ra_pmra_cov    , dec_pmra_cov    , parallax_pmra_cov , pmra_error**2    , pmra_pmdec_cov    ],
                  [ra_pmdec_cov   , dec_pmdec_cov   , parallax_pmdec_cov, pmra_pmdec_cov   , pmdec_error**2    ]]) # Eq. 4.75

    # Transform the covariance matrix
    G = np.array([p_gal, q_gal]).dot(AgT.dot(np.array([p_icrs, q_icrs]).T)) # Eq. 4.81
    J = np.array([[G[0][0], G[0][1], 0,       0,       0],
                  [G[1][0], G[1][1], 0,       0,       0],
                  [      0,       0, 1,       0,       0],
                  [      0,       0, 0, G[0][0], G[0][1]],
                  [      0,       0, 0, G[1][0], G[1][1]]]) # Eq. 4.80
    C_gal = J @ C @ J.T # Eq. 4.79

    # Transformed uncertainties and correlations
    l_error, b_error, parallax_error, pml_error, pmb_error = np.sqrt(np.diag(C_gal))
    l_b_corr = C_gal[0][1]/(l_error*b_error)
    l_parallax_corr = C_gal[0][2]/(l_error*parallax_error)
    l_pml_corr = C_gal[0][3]/(l_error*pml_error)
    l_pmb_corr = C_gal[0][4]/(l_error*pmb_error)
    b_parallax_corr = C_gal[1][2]/(b_error*parallax_error)
    b_pml_corr = C_gal[1][3]/(b_error*pml_error)
    b_pmb_corr = C_gal[1][4]/(b_error*pmb_error)
    parallax_pml_corr = C_gal[2][3]/(parallax_error*pml_error)
    parallax_pmb_corr = C_gal[2][4]/(parallax_error*pmb_error)
    pml_pmb_corr = C_gal[3][4]/(pml_error*pmb_error)

    # Pack results together
    vals = [l, b, parallax, pml, pmb]
    errors = [l_error, b_error, parallax_error, pml_error, pmb_error]
    corr = [l_b_corr, l_parallax_corr, l_pml_corr, l_pmb_corr, b_parallax_corr,
            b_pml_corr, b_pmb_corr, parallax_pml_corr, parallax_pmb_corr, pml_pmb_corr]
    return [vals, errors, corr]
