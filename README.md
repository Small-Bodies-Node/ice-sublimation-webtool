# ice-sublimation-webtool
Cometary ice sublimation model, based on Cowan & A'Hearn (1979, EM&P 21, 155).
***
## Attibution

## Notice


***
## Overview
A rapidly rotating nucleus is one for which the thermal inertia is large enough that parallels of latitude become isotherms. The average sublimation over the nucleus is then a function of the obliquity, i.e., of the orientation of the rotation axis. The sublimation is relatively high if the nucleus is pole-on toward the sun (obliquity = 0°) and much smaller if the axis is perpendicular to the comet-sun line (obliquity = 90°). All calculations assume a spherical body. Note that a non-rotating comet is thus identical to a comet that is pole-on toward the sun. The approach given here calculates the sublimation separately at each latitude (and for a rapid rotator the sublimation is constant all the way around the parallel of latitude, even on the night side) and then calculates the appropriate average over the entire surface, i.e., the average over all 4*pi*R^2 of the surface, including areas where the actual sublimation is zero.

***

## Web Tool
### Setup
Flask is the only python dependency which is not a part of the standard library. It can be added from pip via `pip install Flask`.

### Usage
To launch the web tool, navigate to the project's top level directory, `ice-sublimation-webtool/`, and enter the command `python app.py`.

Your command line will then tell you which port the web tool is assessable from. By default, this will be http://127.0.0.1:5000/. Follow the link in your browser to access the project. 

***

## Python Scripts

## fastrot.py
The rapidly rotating nucleus model is available as `fastrot.py`. 

### Setup
This script does not have any dependencies outside the python3 standard library. 

### Usage
```
fastrot.py [-h] --Av visual albedo --Air infrared albedo --rh
	heliocentric_distance --obl obliquity [--temp temperature]
        [--verbosity verbosity]
        species

example:
 python fastrot.py 'H2O' --Av=0.00 --Air=0.05 --rh=3.98 --obl=90
```
1. species - string (or integer for backwards compatibility) specifying the desired ice: (1) H2O, (2) H2O+CH4 clathrate, (3) CO2, (4) CO.
2. Av - visual albedo
3. Air - infrared albedo
4. rh - heliocentric distance in au
5. obliquity - 90 - angle between rotation axis and the solar direction

## survey_fastrot.py
This script can be used to call `fastrot.py` over a parameter space with a single call.

### Set up
`survey_fastrot.py` does not introduce any additional packages outside the standard library. 
The only added requirement is that both `survey_fastrot.py` and `fastrot.py` are kept within the same directory.

### Usage
```
usage: survey_fastrot.py [-h] --species_list species [species ...] --Av_list visual albedo [visual albedo ...] --Air_list infrared albedo [infrared albedo ...]
                         --rh_list heliocentric_distance [heliocentric_distance ...] --obl_list obliquity [obliquity ...]
                        
example:
python survey_fastrot.py --species_list 'H2O' 'H2O-CH4' --Av_list 0.05 0.10 --Air_list 0 --rh_list 3.98 --obl_list 40 41 42
```
1. species_list species [species ...]
                        Desired ice species to be considered. 
                        The valid inputs are: H2O, H2O-CH4, CO2, CO
2. Av_list - visual albedo [visual albedo ...]
3. Air_list - infrared albedo [infrared albedo ...]
4. rh_list - heliocentric_distance [heliocentric_distance ...]
5. obl_list - obliquity [obliquity ...]

***

## Data
Tabulated results are also available for the `fastrot.py` model in the following [Google Drive folder](https://drive.google.com/drive/folders/1j3tAJtPspPRGPS2Jj_XYl1ZPzT5nkOnz?usp=sharing). This lookup table is too large to be stored directory on Github as the file for each species is roughly 80 MB.

The directory `example-data/` provides the pole-on case, which is identical both to the non-rotating case and to the case of zero thermal inertia. The visual Bond albedo is 5% and the thermal emissivity is 100%. This directory was included to provide a direct comparison between the updated model and the original FORTRAN codes. The entirety of the original version can be found within `oldVersion/` and the results to compare to `example-data/` are located within `oldVersion/data/`.

***

## Version History
The original iteration of this work was developed by M. A'Hearn. He created a simple tool to calculate the sublimation of ices under various circumstances. The calculations are all based on the methods described by Cowan and A'Hearn (1979 Moon and Planets 21, 155-171), which in turn are based on earlier work by Delsemme and others as referenced in that paper. The calculations have not been substantially altered from the original publication and but updated vapor pressures and latent heats have been used. These changes to the input parameters lead to only small changes in the resultant sublimation rates. We note, as pointed out to us by W. Huebner and D. Boice, that the use of empirical results for heat of sublimation and vapor pressure lead the results to be inconsistent with the Clausius-Clapeyron equation, but this has negligible effect on the results for most (but not all) physical situations. Results are available for pure water, pure CO, pure CO2.

In the summer of 2021, M. Van Selous rewrote the model using python and created the updated web tool and lookup table.