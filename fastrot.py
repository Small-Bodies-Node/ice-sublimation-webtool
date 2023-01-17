"""
Description
-----------
This program calculates the average sublimation per unit area for a rapidly
rotating cometary nucleus.  For a sufficiently rapid rotation, or equivalently
for sufficiently high thermal inertia, a parallel of latitude is an isotherm and
this is assumed by the program.

The method integrates an energy balance equation with the Newton-Raphson method
to derive an equilibrium temperature.  Simpson's rule is used to integrate over
latitude.

The properties of the ice are handled in the separate subroutine sublime which
provides the vapor pressure and latent heat, as well as their derivatives, for a
given temperature.


Modification History
--------------------
Michael S. P. Kelley January 2023
  - Make file output optional.
  - Improve error handing.
  - Clarify definition of obliquity.
  - Number of latitude steps added as a parameter.
  - Simplify the code a bit.
  - Force vapor pressure to 0 for very low CO2 and CO temperatures.

Mark Van Selous July 2021
  - Rewrote cgifastrot.f and sublime.f in python as fastrot.py
  - Increased the latitude step size from 41 to 181.
  - Thanks to this increased step size and having more accurate trig functions
    than were available in fortran77, there is very slight alteration in the
    final results.

B. Prager 06/24
  - Original fastrot.f renamed to cgifastrot.f to make a distinction.
  - Removed file i/o from cgifastrot.f to avoid permissions issues with cgi
    scripts.
  - Added command line input of parameters. Parameter 1: Species. Parameter 2:
    Visual Albedo. Parameter 3: Infared Albedo. Parameter 4: Heliocentric
    Distance. Parameter 5: inclination.
  - Initialized the parameters such that they can accept nine characters of
    data.
"""

import csv
import math
import sys
import json
import logging

speciesList = ["H2O", "H2O-CH4", "CO2", "CO"]

# Constants
sigma = 5.67e-5
f0 = 1.39e6
dyncm2 = 1.33322e3  # conversion from Torr to dyne/cm2
boltz = 1.38e-16
ergcal = 6.953e-17
proton = 1.67e-24

tstart = {
    "H2O": 190,
    "H2O-CH4": 190,
    "CO2": 100,
    "CO": 60,
}


def sublime(species, temperature):
    """
    Calculates the latent heat of sublimation and the vapor pressure of the
    solid for various ices and the derivatives thereof.


    Parameters
    ----------
    species : str
        Ice species to consider:
          - 'H2O'
          - 'H2O_CH4'
          - 'CO2'
          - 'CO'

    temperature: float
        Temperature (Kelvin).  If temperature <= 0 K, then the code defaults to
        a species dependent initial value:
          - H2O: 190 K
          - H2O-CH4: 190 K
          - CO2: 100 K
          - CO: 60 K


    Returns
    -------
    mass : float

    xlt : float

    xltprim : float

    press : float

    pprim : float

    temperature : float
        The output temperature is either the input temperature or the species
        dependent initial temperature if the input temperature is less than 0 K.

    """
    mass, xlt, xltprim, press, pprim = None, None, None, None, None

    # If temp <=0, then use the species dependent starting point.
    if temperature <= 0:
        temperature = tstart[species]

    t = temperature
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    t5 = t4 * t
    t6 = t3 * t3

    if species == "H2O":
        mass = 18.0
        xlt = 12420.0 - 4.8 * t
        xltprim = -4.8

        # from Marti & Mauersberger (1993 GRL 20, 363)
        press = -2663.5 / t + 12.537
        press = 10.0 * 10.0**press
        pprim = (2663.5 / t2) * press

    elif species == "H2O-CH4":
        mass = 18.0
        xlt = 12160.0 + 0.5 * t - 0.033 * t2
        xltprim = 0.5 - 0.066 * t

        # from Marti & Mauersberger (1993 GRL 20, 363)
        press = -2663.5 / t + 12.537
        press = 10.0 * 10.0**press
        pprim = (+2663.5 / t2) * press

    elif species == "CO2":
        mass = 44.0
        xlt = 6269.0 + 9.877 * t - 0.130997 * t2 + 6.2735e-4 * t3 - 1.2699e-6 * t4
        xltprim = 9.877 - 0.261994 * t + 1.88205e-3 * t2 - 5.0796e-6 * t3

        if t <= 20:
            logging.warn("CO2 temperature < 20 K")
            press = 0
            pprim = 0
        else:
            press = (
                21.3807649e0
                - 2570.647e0 / t
                - 7.78129489e4 / t2
                + 4.32506256e6 / t3
                - 1.20671368e8 / t4
                + 1.34966306e9 / t5
            )
            pprim = (
                2570.647e0 / t2
                + 1.556258978e5 / t3
                - 12.97518768e6 / t4
                + 4.82685472e8 / t5
                - 6.7483153e9 / t6
            )
            press = dyncm2 * 10.0**press
            pprim = pprim * press

    elif species == "CO":
        mass = 28
        if t > 68.127:
            sys.exit(f"error in CO temp, T = {t}")
        elif t > 61.544:
            xlt = 1855 + 3.253 * t - 0.06833 * t2
            xltprim = 3.253 - 0.13666 * t
            press = (
                16.8655152e0 - 748.151471e0 / t - 5.84330795e0 / t2 + 3.93853859e0 / t3
            )
            pprim = 748.15147e0 / t2 + 11.6866159e0 / t3 - 11.81561577e0 / t4
        elif t >= 14.0:
            xlt = (
                1893
                + 7.331 * t
                + 0.01096 * t2
                - 0.0060658 * t3
                + 1.166e-4 * t4
                - 7.8957e-7 * t5
            )
            xltprim = (
                7.331 + 0.02192 * t - 0.0181974 * t2 + 4.664e-4 * t3 - 3.94785e-6 * t4
            )

            press = (
                18.0741183e0
                - 769.842078e0 / t
                - 12148.7759e0 / t2
                + 2.7350095e5 / t3
                - 2.9087467e6 / t4
                + 1.20319418e7 / t5
            )
            pprim = (
                769.842078e0 / t2
                + 24297.5518 / t3
                - 820502.85e0 / t4
                + 11634986.8e0 / t5
                - 60159709.0e0 / t6
            )
            press = dyncm2 * 10.0**press
            pprim = pprim * press

        else:
            logging.warn("CO temperature < 14 K")
            xlt = (
                1893
                + 7.331 * t
                + 0.01096 * t2
                - 0.0060658 * t3
                + 1.166e-4 * t4
                - 7.8957e-7 * t5
            )
            xltprim = (
                7.331 + 0.02192 * t - 0.0181974 * t2 + 4.664e-4 * t3 - 3.94785e-6 * t4
            )

            press = 0
            pprim = 0

    mass = mass * proton
    xlt = xlt * ergcal
    xltprim = xltprim * ergcal

    return mass, xlt, xltprim, press, pprim, temperature


def run_model(species, Av, Air, rh, obliquity, nlat, temperature0=-1, verbosity=1):
    """
    A call of this function replicates the behavior of the original cgifastrot.f
    script. After reading validating the input parameters, run_model() will
    iterate through mainloop().


    Parameters
    ----------

    species : str
        Ice species to consider
          - H2O
          - H2O_CH4
          - CO2
          - CO

    Av : float (Av > 0)
        Visual albedo

    Air : float
        Infrared albedo

    rh : float
        Heliocentric distance (in au)

    obliquity : float
        Obliquity, angle between the object's rotational axis and its orbital
        axis.

    nlat : int
        Number of latitude bands to calculate.

    temperature0: float
        Initial temperature guess. If this parameter is not specified, a species
        dependent initial value will be used:
          - H2O: 190 K
          - H2O-CH4: 190 K
          - CO2: 100 K
          - CO: 60 K'

    verbosity: int
        Used to specify the logging level
          - verbosity = 0: Only the final results (or a fatal error) will be
            displayed).
          - verbosity = 1: The input parameters are also displayed.
          - Otherwise: Additional output will be displayed for debugging
            purposes.


    Returns
    -------
    species: str
        Inputted species

    obliquity: float
        Inputted obliquity

    rh: float
        Inputted heliocentric distance (au)

    rlog: float
        rlog = log10(rh)

    Av : float (Av > 0)
        Inputted visual albedo

    Air : float
        Inputted infrared albedo

    zbar: float
        Average sublimation (in molecules cm^-2 s^-1)

    zlog: float
        zlog = log10(zbar)

    """

    if species not in speciesList:
        logging.error(
            f'The inputted species of "{species}" is not one of {speciesList}'
        )
        raise ValueError("Invalid species.")

    if Av < 0:
        logging.error(
            f"A visual albedo of {Av} is not a valid input."
            " Please input a value greater than 0."
        )
        raise ValueError("Invalid visual albedo.")

    if verbosity == 0:
        logging.basicConfig(level="WARNING")
    elif verbosity == 1:
        logging.basicConfig(level="INFO")
    else:
        logging.basicConfig(level="DEBUG")

    logging.info("Input Parameters:")
    logging.info(
        f"Species = {species}, Avis = {Av}, Air = {Air}, r_H = {rh}, Obl = {obliquity}"
    )

    incl = (90 - obliquity) * math.pi / 180  # radians

    z = [0] * nlat  # sublimation rate as a function of latitude
    sin_latitude = [-1]
    delta_sin_latitude = 2.0 / (nlat - 1)  # sin(latitude) step size

    for idx in range(1, nlat):
        sin_latitude.append(sin_latitude[0] + idx * delta_sin_latitude)

    niter_total = 0  # total number of iterations for all latitudes
    for i in range(0, nlat):
        latitude = math.asin(sin_latitude[i])
        frac = 0  # insolation scale factor at latitude

        if latitude <= -incl:
            frac = 0
        elif latitude > incl:
            frac = sin_latitude[i] * math.cos(incl)
        else:
            x1 = (
                math.cos(incl)
                * sin_latitude[i]
                * (math.acos(-math.tan(latitude) * (1 / math.tan(incl))))
                / math.pi
            )
            x2 = (
                math.sin(incl)
                * math.cos(latitude)
                * math.sin(math.acos(-math.tan(latitude) / math.tan(incl)))
                / math.pi
            )
            frac = x1 + x2

        niter = 0  # number of iterations for this latitude
        temperature = temperature0
        if frac > 0:
            while niter < 100000:
                z[i], temperature, converged = main_loop(
                    species, Av, Air, rh, frac, temperature
                )
                niter += 1
                if converged:
                    break
            else:
                raise RuntimeError("Energy balance iteration did not converge.")

        logging.debug(
            "obliquity: %f, latitude: %f, z: %g, iterations: %d",
            obliquity,
            latitude * 180 / math.pi,
            z[i],
            niter,
        )
        niter_total += niter

    zbar = 0.0
    for i in range(0, nlat - 1):
        zbar = zbar + 0.5 * (z[i] + z[i + 1]) * delta_sin_latitude

    zbar = zbar / 2
    zlog = math.log10(zbar)
    rlog = math.log10(rh)

    output = {
        "species": species,
        "obliquity": obliquity,
        "r_H": rh,
        "rlog": rlog,
        "Av": Av,
        "Air": Air,
        "Zbar": zbar,
        "Zlog": zlog,
    }

    logging.info("Final Results:")
    logging.info(output)
    return output


def main_loop(species, Av, Air, rh, frac, temperature):
    """Calculate temperature and sublimation rate.


    Parameters
    ----------
    species: str
        Inputted species

    Av : float
        Visual albedo (Av > 0)

    Air : float
        Infrared albedo

    rh: float
        Heliocentric distance (au)

    frac : float
        Insolation scale factor of this obliquity and latitude

    temperature : float
        Estimated temperature at latitude (Kelvins).


    Returns
    -------
    z : float
        Sublimation rate.

    temperature : float
        Updated temperature estimate (Kelvins).

    converged : bool
        `True` if the energy equation is balanced.

    """

    mass, xlt, xltprim, press, pprim, temperature = sublime(species, temperature)
    root = 1 / math.sqrt(mass * 2 * math.pi * boltz)
    root_t = math.sqrt(temperature)
    sun = f0 * frac * (1.0 - Av) / rh**2
    radiat = (1 - Air) * sigma * temperature**4
    evap = root / root_t * press * xlt
    phi = radiat + evap - sun
    z = max(evap / xlt, 1e-30)

    drad = 4 * radiat / temperature
    x1 = pprim * xlt
    x2 = press * xltprim

    devap = root / root_t * (x1 + x2)
    phipri = drad + devap

    dt = math.copysign(min(10, abs(phi / phipri / 2)), phi / phipri)
    temperature -= dt

    converged = abs(phi / sun) < 1e-4 or abs(phi) < 1e-4

    return z, temperature, converged


############
description1 = (
    "This program calculates the average sublimation per unit area for a rapidly rotating cometary"
    " nucleus. For a sufficiently rapid rotation, or equivalently for sufficiently high thermal inertia,"
    " a parallel of latitude is an isotherm and this is assumed by the program."
)

description2 = (
    "The method integrates an energy balance equation with the Newton-Raphson method"
    " to derive an equilibrium temperature.  Simpson's rule is used to integrate over"
    " latitude."
)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="\n\n".join([description1, description2]),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("species", choices=speciesList, help="Ice species to consider.")
    parser.add_argument(
        "--Av",
        metavar="visual_albedo",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--Air",
        metavar="infrared_albedo",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--rh",
        metavar="heliocentric_distance",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--obl",
        metavar="obliquity",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--nlat", metavar="n", type=int, default=181, help="Number of latitude steps"
    )
    parser.add_argument(
        "--temp",
        metavar="temperature",
        type=float,
        default=-1,
        help="Not passing a starting temperature will default to a species dependent starting value:\n"
        "H2O: 190 K\n"
        "H2O-CH4: 190 K\n"
        "CO2: 100 K\n"
        "CO: 60 K\n",
    )
    parser.add_argument(
        "-o", metavar="filename", dest="filename", help="Save results to this file name"
    )
    parser.add_argument(
        "--format", choices=["json", "csv"], default="json", help="output file format"
    )
    parser.add_argument(
        "--verbosity",
        "-v",
        metavar="verbosity",
        type=int,
        default=0,
        help="By default (verbosity = 0), only the final result will be displayed in stdout."
        " A verbosity of 1 will output the logger messages as well.",
    )

    try:
        args = parser.parse_args()
        results = {
            "status": "success",
            "results": run_model(
                args.species,
                args.Av,
                args.Air,
                args.rh,
                args.obl,
                args.nlat,
                args.temp,
                args.verbosity,
            ),
        }
    except Exception as e:
        results = {"status": "failure", "message": str(e)}

    print(json.dumps(results))

    if results["status"] != "failure" and args.filename is not None:
        if args.format == "json":
            with open(args.filename, "w") as json_file:
                json.dump(results, json_file)
        else:
            fieldnames = list(results["results"].keys())
            with open(args.filename, "w") as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(results["results"])
