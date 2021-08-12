"""
    Description
    -----------
    Runs the main script (fastrot.py) and saves the results in JSON.

    The advantage of calling survey_fastrot.py is that it
    iterates over a user defined search space in a single call.
"""
import os
import csv
import time
import math

def survey_fastrot(species_list, Av_list, Air_list, rh_list, obl_list):
    import fastrot
    from itertools import product
    from json import dump

    search_space = []
    arguments = locals()
    len_search_space = 1
    for param in [Av_list, Air_list, rh_list, obl_list]:
        search_space.append(param if isinstance(param, list) else [param])
        len_search_space = len_search_space*len(param)

    start_time = time.time()
    for species in species_list:
        file_name = '_'.join([species, 'output.csv'])
        csv_path = os.path.join(os.getcwd(), 'results', file_name)
        temp_results = []
        for idx, inputs in enumerate(product(*search_space)):
            try:
                next_result = fastrot.run_model(species, *inputs, verbosity=0)
            except Exception as e:
                next_result = {
			"species": species,
    			"obliquity": inputs[3],
        		"r_H": inputs[2],
        		"rlog": math.log10(inputs[2]),
        		"Av": inputs[0],
        		"Air": inputs[1],
        		"Zbar": e,
        		"Zlog": e,
    		}
            temp_results.append(next_result)

            if idx == 0:
                fieldnames = temp_results[0].keys()
                with open(csv_path, 'w') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                    writer.writeheader()
                    for row in temp_results:
                        writer.writerow(row)
                temp_results = []

            elif idx % 5000 == 1:
                percent_complete = idx/(len_search_space*4)
                print('Percent Complete:', round(percent_complete*100, 4))
                duration = time.time() - start_time
                print('Duration:', duration)
                remaining_estimate = duration/percent_complete - duration
                print('Will be finished in ', int(remaining_estimate), ' seconds')

                with open(csv_path, 'a') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                    for row in temp_results:
                        writer.writerow(row)
                temp_results = []

        percent_complete = idx / (len_search_space * 4)
        print('Percent Complete:', round(percent_complete * 100, 4))
        duration = time.time() - start_time
        print('Duration:', duration)
        remaining_estimate = duration / percent_complete - duration
        print('Will be finished in ', int(remaining_estimate), ' seconds')

        with open(csv_path, 'a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            for row in temp_results:
                writer.writerow(row)
        temp_results = []

    print('all done')

if __name__ == '__main__':
    speciesList = ['H2O', 'H2O-CH4', 'CO2', 'CO']

    from helperfunctions import linspace

    Av_list = linspace(0.5, .95, 19)
    Air_list = linspace(0.5, .95, 19)
    rh_list = linspace(1, 30, 30)
    obl_list = linspace(0, 90, 91)

    survey_fastrot(speciesList, Av_list, Air_list, rh_list, obl_list)

