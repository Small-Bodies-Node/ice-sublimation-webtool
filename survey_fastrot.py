"""
    Description
    -----------
    Runs the main script (fastrot.py) and saves the results in JSON.

    The advantage of calling survey_fastrot.py is that it
    iterates over a user defined search space in a single call.
"""
import os
import csv


def survey_fastrot(species_list, Av_list, Air_list, rh_list, obl_list):
    import fastrot
    from itertools import product
    from json import dump

    search_space = []
    arguments = locals()
    for param in [species_list, Av_list, Air_list, rh_list, obl_list]:
        search_space.append(param if isinstance(param, list) else [param])
    results = []
    for inputs in product(*search_space):
        results.append(fastrot.run_model(*inputs))

    output_json = {'results': results}
    json_path = os.path.join(os.getcwd(), 'results', 'output.json')
    with open(json_path, 'w') as json_file:
        dump(output_json, json_file)

    print('result[0]')
    print(results[0])
    fieldnames = results[0].keys()
    csv_path = os.path.join(os.getcwd(), 'results', 'output.csv')
    with open(csv_path, 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    return output_json


if __name__ == '__main__':
    survey_fastrot('H2O', [0.2, 0.3, 0.4], [0.05], 3.98, [45, 46, 47])
