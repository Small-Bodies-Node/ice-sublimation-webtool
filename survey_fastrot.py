"""
    Description
    -----------
    survey_fastrot.py runs the main ice sublimation code (fastrot.py) over a desired parameter space.

    The code accepts both list or standalone values for each of the five input parameters. It then takes the product
    of these sets and computes the rate of sublimation for each combination of input parameters.


    The parameters listed below correspond are describing the elements of the five parameter sets
    (species_set, Av_set, Air_set, rh_set, obl_set)

    Parameters
    ----------
    species : str or int
        Desired ice species to be considered
        - 1: 'H2O'
        - 2: 'H20_CH4'
        - 3: 'CO2'
        - 4: 'CO'
    Av : float (Av > 0)
        Visual albedo
    Air : float
        Infrared albedo
    rh : float
        Heliocentric distance (in au)
    obliquity : float
        Obliquity - 90 - angle between rotation axis and the solar direction


    Returns
    -------
    output_json : dict
        Python dictionary containing the exact output which is also stored within `results/output.json`
        Its only key-value pair, output_json["results"], is a list of the output dictionaries from each call of fastrot.py.
        The keys for each of these inner dictionaries are:
            - "species" : str
            - "Av" : float
            - "Air" : float
            - "r_H" : float
            - "rlog" : float
            - "obliquity" : float
            - "Zbar" : float
            - "Zlog" : float
    'results/output.json` : .json file
    `results/output.csv` : .csv file
"""
import os
import csv
import fastrot
from itertools import product
from json import dump


def survey_fastrot(species_set, Av_set, Air_set, rh_set, obl_set, nlat):
    search_space = []
    arguments = locals()
    for param in [species_set, Av_set, Air_set, rh_set, obl_set, nlat]:
        search_space.append(param if isinstance(param, list) else [param])
    results = []
    for inputs in product(*search_space):
        results.append(fastrot.run_model(*inputs))

    output_json = {"results": results}
    json_path = os.path.join(os.getcwd(), "results", "output.json")
    with open(json_path, "w") as json_file:
        dump(output_json, json_file)

    fieldnames = results[0].keys()
    csv_path = os.path.join(os.getcwd(), "results", "output.csv")
    with open(csv_path, "w") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    return output_json


description = "Use this program to iterate `fastrot.py` over a desired parameter space."

if __name__ == "__main__":
    import argparse

    speciesList = ["H2O", "H2O-CH4", "CO2", "CO"]

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--species_set", 
        metavar="species", choices=speciesList, 
        nargs="+", 
        required=True, help="Ice species to consider.")
    parser.add_argument(
        "--Av_set",
        nargs="+",
        metavar="visual_albedo",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--Air_set",
        nargs="+",
        metavar="infrared_albedo",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--rh_set",
        nargs="+",
        metavar="heliocentric_distance",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--obl_set",
        nargs="+",
        metavar="obliquity",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--nlat", metavar="n", type=int, default=181, help="Number of latitude steps"
    )

    try:
        args = parser.parse_args()
        survey_fastrot(
            args.species_set, args.Av_set, args.Air_set, args.rh_set, args.obl_set, args.nlat
        )
    except Exception as e:
        print(e)
