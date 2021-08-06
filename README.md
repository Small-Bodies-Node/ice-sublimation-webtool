# ice-sublimation
Cometary ice sublimation model, based on Cowan & A'Hearn (1979, EM&P 21, 155). An online version of this tool is at https://pds-smallbodies.astro.umd.edu/tools/ma-evap/index.shtml.

## Overview
This work, originally developed by M. A'Hearn, provides a simple tool to calculate the sublimation of ices under various circumstances. The calculations are all based on the methods described by Cowan and A'Hearn (1979 Moon and Planets 21, 155-171), which in turn are based on earlier work by Delsemme and others as referenced in that paper. The calculations have not been substantially altered from the original publication and but updated vapor pressures and latent heats have been used. These changes to the input parameters lead to only small changes in the resultant sublimation rates. We note, as pointed out to us by W. Huebner and D. Boice, that the use of empirical results for heat of sublimation and vapor pressure lead the results to be inconsistent with the Clausius-Clapeyron equation, but this has negligible effect on the results for most (but not all) physical situations. Results are available for pure water, pure CO, pure CO2.

## fastrot
A rapidly rotating nucleus is one for which the thermal inertia is large enough that parallels of latitude become isotherms. The average sublimation over the nucleus is then a function of the obliquity, i.e., of the orientation of the rotation axis. The sublimation is relatively high if the nucleus is pole-on toward the sun (obliquity = 0°) and much smaller if the axis is perpendicular to the comet-sun line (obliquity = 90°). All calculations assume a spherical body. Note that a non-rotating comet is thus identical to a comet that is pole-on toward the sun. The approach given here calculates the sublimation separately at each latitude (and for a rapid rotator the sublimation is constant all the way around the parallel of latitude, even on the night side) and then calculates the appropriate average over the entire surface, i.e., the average over all 4*pi*R^2 of the surface, including areas where the actual sublimation is zero.

The rapidly rotating nucleus model is available as `cgifastrot.f`.  This code was used by the PDS Small Bodies Node for an online tool.  It requires `sublime.f`.  To compile:
```
gfortran cgifastrot.f sublime.f -o fastrot
```

Usage:
```
fastrot <species> <Av> <Air> <rh> <obliquity>
```
1. species - integer specifying the desired ice: (1) H2O, (2) H2O+CH4 clathrate, (3) CO2, (4) CO.
2. Av - visual albedo
3. Air - infrared albedo
4. rh - heliocentric distance in au
5. obliquity - 90 - angle between rotation axis and the solar direction

## Data
Tabulated results are available for the `fastrot` model and pole-on case, which is identical both to the non-rotating case and to the case of zero thermal inertia. The visual Bond albedo is 5% and the thermal emissivity is 100%. To calculate average sublimation for other obliquities (angles between equatorial plane and the comet-sun line), other values of the albedo and emissivity, or other heliocentric distances, this provided code may be used. The example output is given in the `data/` directory.
