PRO FULL_GALACTIC_TRANSFORM, ra, dec, ra_error, dec_error, plx, plx_error, pmra, $
                             pmdec, pmra_error, pmdec_error, ra_dec_corr=ra_dec_corr, $
                             ra_plx_corr=ra_plx_corr, ra_pmra_corr=ra_pmra_corr, $
                             ra_pmdec_corr=ra_pmdec_corr, dec_plx_corr=dec_plx_corr, $
                             dec_pmra_corr=dec_pmra_corr, dec_pmdec_corr=dec_pmdec_corr, $
                             plx_pmra_corr=plx_pmra_corr, plx_pmdec_corr=plx_pmdec_corr, $
                             pmra_pmdec_corr=pmra_pmdec_corr, vals=vals, errors=errors, corr=corr
; =========================================================================================================
;+
; NAME:
;   FULL_GALACTIC_TRANSFORM
;
; PURPOSE:
;   Transforms coordinates and proper motions from equatorial to Galactic, by
;   taking into account the full covariance matrix.
;
; CREDITS:
;   Created by Michelangelo Pantaleoni González in 2022 (Python) and 2025 (IDL).
;   Based on Butkevich & Lindegren (2022), Gaia DR3 Documentation §4.1.7.
;   https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html
;
; INPUTS:
;   Astrometric parameters in equatorial coordinates:
;       ra    - Right Ascension in degrees
;       dec   - Declination in degrees
;       plx   - Parallax in milliarcseconds (mas) | Coordinate independent
;       pmra  - Proper motion component in right ascension in milliarcseconds per year (mas/a)
;       pmdec - Proper motion component in declination in mas/a
;   Uncertainties in the astrometric parameters:
;       ra_error    - Uncertainty in ra, in degrees
;       dec_error   - Uncertainty in dec, in degrees
;       plx_error   - Uncertainty in plx, in mas
;       pmra_error  - Uncertainty in pmra, in mas/a
;       pmdec_error - Uncertainty in pmdec, in mas/a
;   Correlation coefficients between the astrometric parameters (optional; the
;   values that are not parsed are assumed to be zero):
;       ra_dec_corr     - Correlation between ra and dec
;       ra_plx_corr     - Correlation between ra and plx
;       ra_pmra_corr    - Correlation between ra and pmra
;       ra_pmdec_corr   - Correlation between ra and pmdec
;       dec_plx_corr    - Correlation between dec and plx
;       dec_pmra_corr   - Correlation between dec and pmra
;       dec_pmdec_corr  - Correlation between dec and pmdec
;       plx_pmra_corr   - Correlation between plx and pmra
;       plx_pmdec_corr  - Correlation between plx and pmdec
;       pmra_pmdec_corr - Correlation between pmra and pmdec
;
; OUTPUTS:
;   vals; an array containing the 5 astrometric values in Galactic coordinates
;       l   - Galactic longitude in degrees
;       b   - Galactic latitude in degrees
;       plx - Parallax in mas | Unchanged since inputed
;       pml - Proper motion component in Galactic Longitude in mas/a
;       pmb - Proper motion component in Galactic Latitude in mas/a
;   errors; an array containing the uncertainties for the 5 astrometric values
;       l_error   - Uncertainty in l, in degrees
;       b_error   - Uncertainty in b, in degrees
;       par_error - Uncertainty in plx, in mas
;       pml_error - Uncertainty in pml, in mas/a
;       pmb_error - Uncertainty in pmb, in mas/a
;   corr; an array containing the 10 correlations between astrometric values
;       l_b_corr     - Correlation between l and b
;       l_plx_corr   - Correlation between l and plx
;       l_pml_corr   - Correlation between l and pml
;       l_pmb_corr   - Correlation between l and pmb
;       b_plx_corr   - Correlation between b and plx
;       b_pml_corr   - Correlation between b and pml
;       b_pmb_corr   - Correlation between b and pmb
;       plx_pml_corr - Correlation between plx and pml
;       plx_pmb_corr - Correlation between plx and pmb
;       pml_pmb_corr - Correlation between pml and pmb
;
; EXAMPLE:
; Example case using star Gaia DR3 4598003837468318464.
;   ra = 266.41766795249,
;   ra_error = 0.0195,
;   dec = 31.50467082215,
;   dec_error = 0.0205,
;   parallax = 0.5263,
;   parallax_error = 0.0255,
;   pmra = 4.572,
;   pmra_error = 0.026,
;   pmdec = 1.966,
;   pmdec_error = 0.031,
;   ra_dec_corr = 0.1059,
;   ra_parallax_corr = -0.1897,
;   ra_pmra_corr = 0.1667,
;   ra_pmdec_corr = 0.0971,
;   dec_parallax_corr = -0.0670,
;   dec_pmra_corr = 0.0486,
;   dec_pmdec_corr = 0.0624,
;   parallax_pmra_corr = -0.1347,
;   parallax_pmdec_corr = -0.2701,
;   pmra_pmdec_corr = 0.1921
;
;   FULL_GALACTIC_TRANSFORM, ra, dec, ra_error, dec_error, plx, plx_error, pmra, $
;                            pmdec, pmra_error, pmdec_error, ra_dec_corr=ra_dec_corr, $
;                            ra_plx_corr=ra_plx_corr, ra_pmra_corr=ra_pmra_corr, $
;                            ra_pmdec_corr=ra_pmdec_corr, dec_plx_corr=dec_plx_corr, $
;                            dec_pmra_corr=dec_pmra_corr, dec_pmdec_corr=dec_pmdec_corr, $
;                            plx_pmra_corr=plx_pmra_corr, plx_pmdec_corr=plx_pmdec_corr, $
;                            pmra_pmdec_corr=pmra_pmdec_corr, vals=vals, errors=errors, corr=corr
;
;   print, vals
;   print, errors
;   print, corr
;-
; ========================================================================================================
    ; Correlations must be zero if not passed
    if n_elements(ra_dec_corr) eq 0 then ra_dec_corr = 0.0
    if n_elements(ra_plx_corr) eq 0 then ra_plx_corr = 0.0
    if n_elements(ra_pmra_corr) eq 0 then ra_pmra_corr = 0.0
    if n_elements(ra_pmdec_corr) eq 0 then ra_pmdec_corr = 0.0
    if n_elements(dec_plx_corr) eq 0 then dec_plx_corr = 0.0
    if n_elements(dec_pmra_corr) eq 0 then dec_pmra_corr = 0.0
    if n_elements(dec_pmdec_corr) eq 0 then dec_pmdec_corr = 0.0
    if n_elements(plx_pmra_corr) eq 0 then plx_pmra_corr = 0.0
    if n_elements(plx_pmdec_corr) eq 0 then plx_pmdec_corr = 0.0
    if n_elements(pmra_pmdec_corr) eq 0 then pmra_pmdec_corr = 0.0

    ; Convert inputs to double precission
    ra              = DOUBLE(ra)
    dec             = DOUBLE(dec)
    ra_error        = DOUBLE(ra_error)
    dec_error       = DOUBLE(dec_error)
    plx             = DOUBLE(plx)
    plx_error       = DOUBLE(plx_error)
    pmra            = DOUBLE(pmra)
    pmdec           = DOUBLE(pmdec)
    pmra_error      = DOUBLE(pmra_error)
    pmdec_error     = DOUBLE(pmdec_error)
    ra_dec_corr     = DOUBLE(ra_dec_corr)
    ra_plx_corr     = DOUBLE(ra_plx_corr)
    ra_pmra_corr    = DOUBLE(ra_pmra_corr)
    ra_pmdec_corr   = DOUBLE(ra_pmdec_corr)
    dec_plx_corr    = DOUBLE(dec_plx_corr)
    dec_pmra_corr   = DOUBLE(dec_pmra_corr)
    dec_pmdec_corr  = DOUBLE(dec_pmdec_corr)
    plx_pmra_corr   = DOUBLE(plx_pmra_corr)
    plx_pmdec_corr  = DOUBLE(plx_pmdec_corr)
    pmra_pmdec_corr = DOUBLE(pmra_pmdec_corr)

    ; Define conversion constants
    DPI = DOUBLE(!PI)
    D2R = DPI/180.0D
    R2D = 180.0D/DPI

    ; ICRS coordinates of the north Galactic pole and the longitude of intersecting planes
    ra_pole =  192.85948D * D2R
    dec_pole =  27.12825D * D2R
    l_omega  =  32.93192D * D2R

    ; Rotation matrices
    Rz_lomg = [[ COS(-l_omega),  SIN(-l_omega), 0.0], $
               [-SIN(-l_omega),  COS(-l_omega), 0.0], $
               [      0.0     ,       0.0     , 1.0]]
    Rx_90dec = [[1.0,                    0.0,                   0.0], $
                [0.0,  COS(DPI/2 - dec_pole), SIN(DPI/2 - dec_pole)], $
                [0.0, -SIN(DPI/2 - dec_pole), COS(DPI/2 - dec_pole)]]
    Rz_a90 = [[ COS(ra_pole + DPI/2),  SIN(ra_pole + DPI/2), 0.0], $
              [-SIN(ra_pole + DPI/2),  COS(ra_pole + DPI/2), 0.0], $
              [          0.0        ,           0.0        , 1.0]]
    AgT = Rz_lomg ## (Rx_90dec ## Rz_a90)                                         ; Eq. 4.61

    ; Convert input coordinates to radians
    ra_rad  = ra  * D2R
    dec_rad = dec * D2R

    ; Compute Galactic coordinates
    r_icrs = [COS(ra_rad)*COS(dec_rad), SIN(ra_rad)*COS(dec_rad), SIN(dec_rad)]   ; Eq. 4.58
    r_gal  = AgT ## r_icrs                                                        ; Eq. 4.60
    l_rad  = ATAN(r_gal[1], r_gal[0])                                             ; Eq. 4.64
    b_rad  = ATAN(r_gal[2], SQRT(r_gal[0]^2 + r_gal[1]^2))                        ; Eq. 4.64
    l      = l_rad * R2D
    IF l LT 0 THEN l = l + 360.0
    IF l GE 360.0 THEN l = l - 360.0
    b      = b_rad * R2D

    ; Compute Galactic proper motion
    p_icrs = [-SIN(ra_rad), COS(ra_rad), 0.0]                                     ; Eq. 4.65
    q_icrs = [-COS(ra_rad)*SIN(dec_rad), -SIN(ra_rad)*SIN(dec_rad), COS(dec_rad)] ; Eq. 4.65
    pm_icrs = p_icrs * pmra + q_icrs * pmdec                                      ; Eq. 4.67
    p_gal = [-SIN(l_rad), COS(l_rad), 0.0]                                        ; Eq. 4.66
    q_gal = [-COS(l_rad)*SIN(b_rad), -SIN(l_rad)*SIN(b_rad), COS(b_rad)]          ; Eq. 4.66
    pm_gal = AgT ## pm_icrs                                                       ; Eq. 4.69
    pml    = p_gal ## pm_gal                                                      ; Eq. 4.71
    pmb    = q_gal ## pm_gal                                                      ; Eq. 4.71

    ; Build the covariance matrix C
    ra_dec_cov     = ra_error   * dec_error   * ra_dec_corr
    ra_plx_cov     = ra_error   * plx_error   * ra_plx_corr
    ra_pmra_cov    = ra_error   * pmra_error  * ra_pmra_corr
    ra_pmdec_cov   = ra_error   * pmdec_error * ra_pmdec_corr
    dec_plx_cov    = dec_error  * plx_error   * dec_plx_corr
    dec_pmra_cov   = dec_error  * pmra_error  * dec_pmra_corr
    dec_pmdec_cov  = dec_error  * pmdec_error * dec_pmdec_corr
    plx_pmra_cov   = plx_error  * pmra_error  * plx_pmra_corr
    plx_pmdec_cov  = plx_error  * pmdec_error * plx_pmdec_corr
    pmra_pmdec_cov = pmra_error * pmdec_error * pmra_pmdec_corr
                                                                                   ; Eq. 4.75
    C = [[ra_error^2   , ra_dec_cov    , ra_plx_cov    , ra_pmra_cov    , ra_pmdec_cov  ], $
         [ra_dec_cov   , dec_error^2   , dec_plx_cov   , dec_pmra_cov   , dec_pmdec_cov ], $
         [ra_plx_cov   , dec_plx_cov   , plx_error^2   , plx_pmra_cov   , plx_pmdec_cov ], $
         [ra_pmra_cov  , dec_pmra_cov  , plx_pmra_cov  , pmra_error^2   , pmra_pmdec_cov], $
         [ra_pmdec_cov , dec_pmdec_cov , plx_pmdec_cov , pmra_pmdec_cov , pmdec_error^2 ]]

    ; Transform the covariance matrix by applying the rotation
    icrsmat = [[p_icrs[0], p_icrs[1], p_icrs[2]], [q_icrs[0], q_icrs[1], q_icrs[2]]]
    galmat  = [[p_gal[0] ,  p_gal[1],  p_gal[2]], [q_gal[0] ,  q_gal[1],  q_gal[2]]]
    G = galmat ## (AgT ## Transpose(icrsmat))                                     ; Eq. 4.81
    J = [[G[0,0] , G[1,0] , 0.0 , 0.0    ,    0.0], $
         [G[0,1] , G[1,1] , 0.0 , 0.0    ,    0.0], $
         [0.0    , 0.0    , 1.0 , 0.0    ,    0.0], $
         [0.0    , 0.0    , 0.0 , G[0,0] , G[1,0]], $
         [0.0    , 0.0    , 0.0 , G[0,1] , G[1,1]]]                               ; Eq. 4.80
    C_gal = J ## (C ## Transpose(J))                                              ; Eq. 4.79

    ; Compute transformed uncertainties
    l_error   = SQRT(C_gal[0, 0])
    b_error   = SQRT(C_gal[1, 1])
    par_error = SQRT(C_gal[2, 2])
    pml_error = SQRT(C_gal[3, 3])
    pmb_error = SQRT(C_gal[4, 4])

    ; Compute transformed correlations
    l_b_corr     = C_gal[0,1]/(l_error*b_error)
    l_plx_corr   = C_gal[0,2]/(l_error*par_error)
    l_pml_corr   = C_gal[0,3]/(l_error*pml_error)
    l_pmb_corr   = C_gal[0,4]/(l_error*pmb_error)
    b_plx_corr   = C_gal[1,2]/(b_error*par_error)
    b_pml_corr   = C_gal[1,3]/(b_error*pml_error)
    b_pmb_corr   = C_gal[1,4]/(b_error*pmb_error)
    plx_pml_corr = C_gal[2,3]/(par_error*pml_error)
    plx_pmb_corr = C_gal[2,4]/(par_error*pmb_error)
    pml_pmb_corr = C_gal[3,4]/(pml_error*pmb_error)

    ; Pack results together
    vals   = [l, b, plx, pml, pmb]
    errors = [l_error, b_error, par_error, pml_error, pmb_error]
    corr   = [l_b_corr, l_plx_corr, l_pml_corr, l_pmb_corr, b_plx_corr, $
              b_pml_corr, b_pmb_corr, plx_pml_corr, plx_pmb_corr, pml_pmb_corr]

end
