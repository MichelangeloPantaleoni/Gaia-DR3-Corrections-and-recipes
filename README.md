# Gaia DR3 Corrections and recipes
Gaia DR3 parallax, proper motion and photometric corrections aswell as distance estimation algorithms used in Pantaeloni et al. (in review).

## Parallax corrections
Astrometry in Gaia DR3 remained unchanged with respect to Gaia EDR3 ([Vallenari et al. 2023](https://ui.adsabs.harvard.edu/#abs/2023A%26A...674A...1G)), which means that the current corrections are valid for both data releases.

**Parallax zero-point:**

The parallax values, $\varpi$, in Gaia DR3 (as retrieved from [ESA’s Gaia Archive](https://gea.esac.esa.int/archive/) or [CDS Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3)) have associated to them a systematic bias, $\varpi_0$. Therefore the parallax must be corrected by simply removing this zero-point value:

$\varpi_c = \varpi-\varpi_0$

The algorithm that computes $\varpi_0$ is described in [Maíz Apellániz (2022)](https://ui.adsabs.harvard.edu/abs/2022A&A...657A.130M), which improved on the parallax zero-point recipe in [Lindegren et al. (2021b)](https://ui.adsabs.harvard.edu/#abs/2021A%26A...649A...4L) for the bright stars. The ```correct_parallax()``` function performs this correction following the recipe in Maíz Apellániz (2022) by default (the parameter ```mode``` allows the user to select also 'Lindegren' as a recipe).

The parallax bias is a function of the source brightness, ecliptic latitude, pseudocolor and the number of parameters for the astrometric solution, so in general, these values have to be provided (```phot_g_mean_mag```, ```ecl_lat```, ```pseudocolor``` and ```astrometric_params_solved```). If ```phot_g_mean_mag``` is not present then we assume the star to be faint (meaning $G = 21.0$ mag) and compute the parallax zero-point from there. Values for ```phot_g_mean_mag``` outside the $6.0 < G < 18.0$ range are computed as if they were at the closest boundary of the range, and values inside that range are linearly interpolated between the values described in table 4 in Maíz Apellániz (2022) and tables 9 and 10 in Lindegren et al. (2021b).

Objects with 2-parameter solutions (```astrometric_params_solved``` = 3) have no parallax measurements and thus the function retrieves ```nan``` for the zero point and the corrected parallax. Objects with 5-parameter solutions (```astrometric_params_solved``` = 31) are the ones that show differences in regards to the choosen recipe. Objects with 6-parameter solutions (```astrometric_params_solved``` = 95) give the same result using the recipe in Lindegren et al. (2021b) or that of Maíz Apellániz (2022).

**Parallax uncertainties:**

Gaia DR3 significantly underestimates the catalogue uncertainties for the parallax, $\sigma_{\text{int}}$ (uncertainty "internal" to the catalouge), which has prompted optimistically precise results in the recent literature. To correct these internal random uncertainties we scale them up by a constant factor $k$, whose value is calculated following the recipe in [Maíz Apellániz (2022)](https://ui.adsabs.harvard.edu/abs/2022A%2526A...657A.130M). Then, a previously unaccounted systematic uncertainty of $\sigma_{\text{sys}} = 10.3$ $\mu as$ ([Maíz Apellániz et al. 2021c](https://ui.adsabs.harvard.edu/#abs/2021A%26A...649A..13M)) is added in quadrature, following [Fabricius et al. (2021)](https://ui.adsabs.harvard.edu/#abs/2021A%26A...649A...5F), to yield an "external", total (random + systematic) uncertainty of

$\sigma_{\text{ext}} = \sqrt{(k\sigma_{\text{int}})^2+\sigma_{\text{sys}}^2}$

The calculations are handled by the ```correct_parallax_error()``` function, which computes the multiplicative factor $k$ as a function of the brightness of the source and the ```RUWE``` value. The dependence on the ```RUWE``` has been tested up to $3.0$, which is far beyond the $1.4$ threshold recommended for reliable single-source astrometric solutions ([Lindegren 2018](https://www.semanticscholar.org/paper/Re-normalising-the-astrometric-chi-square-in-Gaia/94f6f242b43ada2675fd46b811bc86584a906019)). ```RUWE``` values beyond 1.4 are thought to express a non-negligible probability of giving bad results for the parallax, as the astrometric solutions provided by the Gaia DPAC are not fitting correctly the single-source astrometry (probably due to unrecognized binarity). This prompts researchers to select sources with lower ```RUWE``` values, but in fact sources with larger ones can still be usefull to some extent. When above ```RUWE``` of 1.4 the $k$ factor simply increases rapidly (reproducing a reasonable increase in the parallax uncertainty).

The resulting $\varpi_c \pm \sigma_{\text{ext}}$ values improve the accuracy of Gaia DR3 $\varpi \pm \sigma_{\text{int}}$ values for single sources by a significant amount. Papers that directly make use of the Gaia DR3 values are prone to make overconfident assessments of parallax uncertainties and systematically overestimate distances for bright stars.
 

**Parallax for groups of stars**
When calculating distances

## Proper motion corrections

## Astrometric correlations when tranforming to Galactic coordinates

## Distances to OB stars

## Photometric corrections
