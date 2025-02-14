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
        static_lca(index, datapackage, output_dir, k, t, act)
    elif t == "pedigree":
        stochastic_lca(index, datapackage, output_dir, k, t, act)
    
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
                datapackage = create_static_datapackage(tech_matrix, bio_matrix, cf_matrix)
                print("Datapackage are formatted.")
                for myact, index in param[2]:
                    k += 1

                    # Pedigree: each case as different functional unit, becasue 1 column for foreground system, so the index has to add 1.
                    # parallel_params.append((t, myact, (index + 1), k, param[3], datapackage))
                    
                    # Pedigree: each case has the same functinoal unit, the index in the constants.py always set to 0, so don't have to add 1 in this case.
                    parallel_params.append((t, myact, index, k, param[3], datapackage))

                    manual_tech_Data = -tech_matrix
                    np.fill_diagonal(manual_tech_Data, -manual_tech_Data.diagonal())

                    # Pedigree: each case as different functional unit, becasue 1 column for foreground system, so the index has to add 1.
                    # lca_score_manual = simu.manual_lca(manual_tech_Data, bio_matrix, cf_matrix, index + 1)

                    # Pedigree: each case has the same functinoal unit, the index in the constants.py always set to 0, so don't have to add 1 in this case.
                    lca_score_manual = simu.manual_lca(manual_tech_Data, bio_matrix, cf_matrix, index)
                    print(f"Manually calculated lca score: {lca_score_manual, myact}")
            elif t == "pedigree":
                datapackage = create_stochastic_datapackage(tech_matrix, bio_matrix, cf_matrix, activities, param[4], param[5])
                print("Datapackage are formatted.")
                for myact, index in param[2]:
                    k += 1
                    # parallel_params.append((t, myact, (index + 1), k, param[3], datapackage))  # pedigree: because add another column
                    parallel_params.append((t, myact, index, k, param[3], datapackage))  # pedigree: same functional nuit, dno't have to add 1.

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
