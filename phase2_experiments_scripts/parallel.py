'''
 # @ Author: Ning An
 # @ Create Time: 2024-11-05 10:45:51
 # @ Modified by: Ning An
 # @ Modified time: 2024-11-06 10:33:51
 '''

from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import os
import sys
sys.path.append(os.path.abspath(".."))
from constants import *
from monte_carlo_utility import *
from monte_carlo_main import *

def process_case(parallel_param):
    """
    parallel_param type: tuple
    parallel_param format: (t, myact, index, k, output_dir, datapackage)
    """
    (t, act, index, k, output_dir, datapackage) = parallel_param
    print(f"Processing {t} in process: {os.getpid()}")

    simu = SimulationScript()

    if t == "static":
        perform_static(index, datapackage, output_dir, k, act, t)
    elif t == "pedigree":
        perform_stochastic(index, datapackage, output_dir, k, act, t)
    
    print(f"{t} simulation is done.")


if __name__ == "__main__":
    start_time = time.time()

    simu = SimulationScript()

    for param in COMBINED_PARAMETERS:
        print(f"{param[4]} simulation is running...")
        tech_matrix, bio_matrix, cf_matrix, activities = prepare_datapackage_matrices(param[0], param[1], param[5])
        parallel_params = []
        k = 0
        for t in DIST_TYPE:
            if t == "static":
                datapackage = create_static_datapackage(tech_matrix, bio_matrix, cf_matrix, activities)
                print("Datapackage are formatted.")
                for myact, index in param[2]:
                    k += 1
                    parallel_params.append((t, myact, index, k, param[3], datapackage))
            elif t == "pedigree":
                datapackage = create_stochastic_datapackage(tech_matrix, bio_matrix, cf_matrix, activities, param[4], param[5])
                print("Datapackage are formatted.")
                for myact, index in param[2]:
                    k += 1
                    parallel_params.append((t, myact, index + 1, k, param[3], datapackage))  # because add another column
            
        max_workers = os.cpu_count() if os.cpu_count() else 4
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_case, parallel_param): parallel_param for parallel_param in parallel_params}
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    print(f'Generated an exception: {exc}')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total execution time(s): {elapsed_time}")
