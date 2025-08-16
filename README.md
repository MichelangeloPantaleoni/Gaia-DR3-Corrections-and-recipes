# Gaia DR3 Corrections and recipes
This repository presents Python code for the Gaia DR3 parallax, proper motion and photometric corrections aswell as distance estimation algorithms used in Pantaeloni et al. (in review).

## Parallax corrections
Astrometry in Gaia DR3 remained unchanged with respect to Gaia EDR3 ([Vallenari et al. 2023](https://ui.adsabs.harvard.edu/#abs/2023A&A...674A...1G)), which means that the current corrections are valid for both data releases.

### Parallax zero-point correction:

The parallax values, $\varpi$, in Gaia DR3 (as retrieved from [ESA’s Gaia Archive](https://gea.esac.esa.int/archive/) or [CDS Vizier](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3)) have a systematic bias, $\varpi_0$, associated to them. The parallax bias was first uncovered in [Lindegren et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A&A...616A...2L) by taking Gaia DR2 parallaxes of quasars all over the sky (which are expected to have true parallaxes indistinguishable from zero) and showing that the distribution was indeed not centered at zero parallax. For Gaia DR3 the parallax bias was recalculated ([Lindegren et al. 2021a](https://ui.adsabs.harvard.edu/#abs/2021A&A...649A...4L)). Therefore the parallax for any star must be corrected by simply removing this zero-point value:

$$\varpi_c = \varpi-\varpi_0$$

The algorithm to compute $\varpi_0$ is described in [Maíz Apellániz (2022)](https://ui.adsabs.harvard.edu/abs/2022A&A...657A.130M), which improved on the parallax zero-point recipe in Lindegren et al. (2021a) for the bright stars. The ```correct_parallax()``` function performs this correction following the recipe in Maíz Apellániz (2022) by default (the parameter ```mode``` allows the user to select also 'Lindegren' as a recipe).

The parallax bias is a function of the source brightness, ecliptic latitude, pseudocolor and the number of parameters for the astrometric solution, so in general, these values have to be provided (```phot_g_mean_mag```, ```ecl_lat```, ```pseudocolor``` and ```astrometric_params_solved```). If ```phot_g_mean_mag``` is not present then we assume the star to be faint (meaning $G = 21.0$ mag) and compute the parallax zero-point from there. Values for ```phot_g_mean_mag``` outside the $6.0 < G < 18.0$ range are computed as if they were at the closest boundary of the range, and values inside that range are linearly interpolated between the values described in table 4 in Maíz Apellániz (2022) and tables 9 and 10 in Lindegren et al. (2021a).

Objects with 2-parameter solutions (```astrometric_params_solved``` = 3) have no parallax measurements and thus the function retrieves ```nan``` for the zero point and the corrected parallax. Objects with 5-parameter solutions (```astrometric_params_solved``` = 31) are the ones that show differences in regards to the choosen recipe. Objects with 6-parameter solutions (```astrometric_params_solved``` = 95) give the same result using the recipe in Lindegren et al. (2021a) or that of Maíz Apellániz (2022).

### Parallax uncertainties:

Gaia DR3 significantly underestimates the catalogue uncertainties for the parallax, $\sigma_{\text{int}}$ (uncertainty "internal" to the catalouge), which has prompted optimistically precise results in the recent literature. To correct these internal random uncertainties we scale them up by a constant factor $k$, whose value is calculated following the recipe in Maíz Apellániz (2022). Then, a previously unaccounted systematic uncertainty of $\sigma_{\text{sys}} = 10.3$ $\mu as$ ([Maíz Apellániz et al. 2021c](https://ui.adsabs.harvard.edu/#abs/2021A&A...649A..13M)) is added in quadrature, following [Fabricius et al. (2021)](https://ui.adsabs.harvard.edu/#abs/2021A&A...649A...5F), to yield an "external", total (random + systematic) uncertainty of

$$\sigma_{\text{ext}} = \sqrt{(k\sigma_{\text{int}})^2+\sigma_{\text{sys}}^2}$$

The calculations are handled by the ```correct_parallax_error()``` function, which computes the multiplicative factor $k$ as a function of the brightness of the source and the ```RUWE``` value. The dependence on the ```RUWE``` has been tested up to $3.0$, which is far beyond the $1.4$ threshold recommended for reliable single-source astrometric solutions ([Lindegren 2018](https://www.semanticscholar.org/paper/Re-normalising-the-astrometric-chi-square-in-Gaia/94f6f242b43ada2675fd46b811bc86584a906019)). ```RUWE``` values beyond $1.4$ are thought to express a non-negligible probability of giving bad results for the parallax, as the astrometric solutions provided by the Gaia DPAC are not fitting correctly the single-source astrometry (probably due to unrecognized binarity). This prompts researchers to select sources with lower ```RUWE``` values, but in fact sources with larger ones can still be usefull to some extent. When above ```RUWE``` of $1.4$ the $k$ factor simply increases rapidly (reproducing a reasonable increase in the parallax uncertainty).

The resulting $\varpi_c \pm \sigma_{\text{ext}}$ values improve the accuracy of Gaia DR3 $\varpi \pm \sigma_{\text{int}}$ values for single sources by a significant amount. Papers that directly make use of the Gaia DR3 values are prone to make overconfident assessments of parallax uncertainties and systematically overestimate distances for bright stars.
 

### Parallax for groups of stars:

When estimating the distance to a group of $n$ stars (all assumed to lie at the same distance; wich is something that can be done for cluster members, as the depth of the cluster is typically negligible compared to the distance to it), one might be tempted to take the average of all the individual distances, derived from the parallaxes for each individual star in the group. However, this approach is problematic as it amplifies biases in the distance estimation procedure. A more robust method is to first compute the average parallax and then convert that average into a unique distance. While both methods agree in the limit of negligible uncertainties, this is never the case for real data, and the second approach should therefore be preferred. An improvement on this method can be made by making sure the group's average parallax, $\varpi_g$, is a weighted mean, where each weight, $w_i$, is inversely related to the uncertainty of each individual previously-corrected parallax in the group, $\varpi_{c,i}$. This follows the procedure presented in [Campillay et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2137C):

$$\varpi_g = \sum_{i = 1}^{n} w_i \varpi_{c,i}$$
, where the weight for $i$-th star is
$$w_i = \frac{1/\sigma_{ext, i}^2}{\sum_{i = 1}^{n} 1/\sigma_{ext, i}^2}$$

Both the group parallax and the associated uncertainty of the group parallax (see the next subsection) are calculated by the ```calculate_group_parallax()``` function.

### Parallax uncertainties for groups of stars:

The aforementioned all-sky Gaia DR3 parallax zero-point baseline is not the full story, as parallax seems to vary also for samples of equally distant sources in narrow fields. For this, quasars don't suffice, as they are too separated for these bias short-scale variations to be noticed. To adress these systematic 2D oscillations in parallax, dense fields like the Magellanic Clouds are typically used ([Arenou et al. 2018](https://ui.adsabs.harvard.edu/abs/2018A&A...616A..17A) for Gaia DR2 and [Lindegren et al. 2021b](https://ui.adsabs.harvard.edu/abs/2021A&A...649A...2L/abstract) for Gaia DR3). This so-called "checkered pattern" is the result of the [mission's scanning law](https://www.youtube.com/watch?v=lntlPhXSRQY). Even if in Gaia DR3 the amplitude of this bias has diminished significantly with respect to Gaia DR2, it is still important to consider this, as for closely separated sources (like in the case of stellar clusters) the bias shared by the sources might counteract or contribute to the baseline systematic, $\varpi_0$.


## Proper motion corrections

## Astrometric correlations when tranforming to Galactic coordinates

## Distances to OB stars

## Photometric corrections
