import numpy as np
import pandas as pd

def correct_edr3_gband(df, columns_in = ['bp_rp', 'astrometric_params_solved', 'phot_g_mean_mag', 'phot_g_mean_flux'], columns_out = ['Gmag_corrected', 'Gflux_corrected']):
    '''
    Corrects the Gaia EDR3 (but not Gaia DR3; see the note) G-band photometry (both
    magnitude and flux).
    Adapted from Brown et al. 2021 [2021A%26A...650C...3G], which is the corrigendum
    of Riello et al. 2021 [2021A%26A...649A...3R].

    Input columns:
        bp_rp
        astrometric_params_solved
        phot_g_mean_mag
        phot_g_mean_flux
    Output columns:
        Gmag_corrected
        Gflux_corrected

    Note:
    This is valid only for Gaia EDR3 and is no longer needed for Gaia DR3, since
    these corrections where applied for DR3 by the DPAC, as stated in the source
    for the code: https://github.com/agabrown/gaiaedr3-6p-gband-correction.
    To give an example where DR3 changed the EDR3 photometry check Gaia (E)DR3 5333451634704851328.
    '''
    gmag_corrected, gflux_corrected = [], []
    for i in range(len(df)):
        br, ap, gmag, gflux = df.loc[i, columns_in].values
        if np.isnan(br) or np.isnan(ap) or np.isnan(gmag) or np.isnan(gflux):
            return np.nan, np.nan
        if np.isscalar(br) or np.isscalar(ap) or np.isscalar(gmag) or np.isscalar(gflux):
            br = np.float64(br)
            ap = np.int64(ap)
            gmag = np.float64(gmag)
            gflux = np.float64(gflux)
        if not (br.shape == ap.shape == gmag.shape == gflux.shape):
            raise ValueError('Function parameters must be of the same shape!')
        do_not_correct = np.isnan(br) | (gmag < 13) | (ap == 31)
        bright_correct = np.logical_not(do_not_correct) & (gmag >= 13) & (gmag <= 16)
        faint_correct = np.logical_not(do_not_correct) & (gmag > 16)
        br_c = np.clip(br, 0.25, 3.0)
        correction_factor = np.ones_like(gmag)
        correction_factor[faint_correct] = 1.00525-0.02323*br_c[faint_correct]+0.01740*np.power(br_c[faint_correct],2)-0.00253*np.power(br_c[faint_correct], 3)
        correction_factor[bright_correct] = 1.00876-0.02540*br_c[bright_correct]+0.01747*np.power(br_c[bright_correct],2)-0.00277*np.power(br_c[bright_correct], 3)
        gmag_corrected.append(gmag-2.5*np.log10(correction_factor))
        gflux_corrected.append(gflux*correction_factor)
    df[columns_out[0]] = np.array(gmag_corrected)
    df[columns_out[1]] = np.array(gflux_corrected)
    return df


def correct_gband(df, column_in = 'phot_g_mean_mag', column_out = 'G_corr'):
    '''
    Corrects the Gaia DR3 G-band photometry (magnitude).
    Adapted from Weiler et al. (in preparation), wich will make this recipe public.

    Note:
    If working with Gaia EDR3 you must first implement the correction of the
    function correct_edr3_gband, but this is unnecesary when working with
    Gaia DR3, since Gaia DR3 already shows the corrected Gaia EDR3 photometry.
    '''
    # Load Weiler & Maíz table
    df_corr = pd.read_csv('Weiler_Corrections.csv')
    df_corr = df_corr.rename(columns = {df_corr.columns[0]:'G_mag', df_corr.columns[2]:'G_bias'})[['G_mag', 'G_bias']]
    # Apply Weiler & Maíz correction
    df[column_out] = df[column_in] + np.interp(df[column_in], df_corr['G_mag'], df_corr['G_bias'])
    return df


def photometry_error(df, mode = 'Apellaniz', columns_in = ['phot_g_mean_flux', 'phot_bp_mean_flux', 'phot_rp_mean_flux', 'phot_g_mean_flux_error', 'phot_bp_mean_flux_error', 'phot_rp_mean_flux_error'], columns_out = ['G_err', 'BP_err', 'RP_err']):
    '''
    Calculates the photometric uncertainty in each Gaia band using the photometric
    flux and its error.

    If mode == 'CDS' it follows the equations in Note (G1) from
    https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/I/350?format=html&tex=true

    If mode == 'Apellaniz', then it uses the systematic uncertainty values
    from Maíz Apellániz & Weiler (2024) [2024arXiv240721388M].

    Inputs a pandas dataframe (df) with these Gaia DR3 parameters:
        phot_g_mean_flux
        phot_bp_mean_flux
        phot_rp_mean_flux
        phot_g_mean_flux_error
        phot_bp_mean_flux_error
        phot_rp_mean_flux_error
    Outputs the same dataframe but including these columns:
        G_err  ==> Uncertainty in the G band
        BP_err ==> Uncertainty in the BP band
        RP_err ==> Uncertainty in the RP band
    '''
    # Extract photometry from dataframe
    gband, gflux, gfluxerr, bpband, bpflux, bpfluxerr, rpband, rpflux, rpfluxerr = df[columns_in].values

    # Calculate sistematic uncertainties
    if mode == 'CDS':
        sigG0  = np.ones(len(df))*2.7553202
        sigBP0 = np.ones(len(df))*2.7901700
        sigRP0 = np.ones(len(df))*3.7793818
    elif mode == 'Apellaniz':
        sigG0  = np.empty(len(df))
        sigBP0 = np.empty(len(df))
        sigRP0 = np.empty(len(df))
        for i in range(len(df)):
            if not (np.isnan(gband[i]) or gband[i] == ''):
                if gband[i] <= 13.0:
                    sigG0[i] = 4.8
                else:
                    sigG0[i] = 3.9
                if gband[i] <= 10.87:
                    sigBP0[i] = 3.8
                    sigRP0[i] = 5.7
                else:
                    sigBP0[i] = 4.0
                    sigRP0[i] = 3.4
            else:
                sigG0[i] = np.nan
                if not np.isnan(bpband[i]):
                    # There are 5309 stars in total in Gaia DR3
                    # having BP photometry but no G photometry.
                    sigBP0[i] = 4.0
                else:
                    sigBP0[i] = np.nan
                if not np.isnan(rpband[i]):
                    # There are 3848128 stars in total in Gaia DR3
                    # having RP photometry but no G photometry.
                    sigRP0[i] = 5.7
                else:
                    sigRP0[i] = np.nan
    else:
        raise ValueError("  ERROR: Wrong mode in photometry_error()\n Choose between 'CDS' and 'Apellaniz'")
    cnst = 2.5/np.log(10)
    # Calculate G uncertainties:
    df[columns_out[0]] = np.hypot(cnst*(gfluxerr/gflux),    sigG0/1000)
    # Calculate BP uncertainties:
    df[columns_out[1]] = np.hypot(cnst*(bpfluxerr/bpflux), sigBP0/1000)
    # Calculate RP uncertainties:
    df[columns_out[2]] = np.hypot(cnst*(rpfluxerr/rpflux), sigRP0/1000)
    return df
