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

def process_case(parallel_param):
    """
    parallel_param type: tuple
    parallel_param format: ((t, "0", myact, index, k, output_dir), matrices)
    """
    (t, u, myact, index, k, output_dir) = parallel_param[0]
    print(f"Processing {t}_{u} in process: {os.getpid()}")

    matrices = parallel_param[1]
    A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = matrices

    simu = SimulationScript()

    if t == "static":
        simu.perform_static(index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C, output_dir, k, t, myact)
        print(f"{t}_{u} simulation is done.")
    elif t == "uniform":
        dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
        simu.perform_simu(index, dp_stochastic, output_dir, k, myact, t, u)
        print(f"{t}_{u} simulation is done.")
    elif t == "log-normal":
        dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
        simu.perform_simu(index, dp_stochastic, output_dir, k, myact, t, u)
        print(f"{t}_{u} simulation is done.")


if __name__ == "__main__":
    start_time = time.time()

    simu = SimulationScript()

    for param in COMBINED_PARAMETERS:
        print(f"{param[4]} simulation is running...")
        A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = simu.build_bw_matrix(param[0], param[1])
        print("Matrices are formatted.")

        matrices = (A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)

        parallel_params = []
        k = 0
        for t in DIST_TYPE:
            if t == "static":
                for myact, index in param[2]:
                    k += 1
                    parallel_params.append((t, 0, myact, index, k, param[3]))
            elif t == "uniform":
                for u in U_UNIFORM:
                    for myact, index in param[2]:
                        k += 1
                        parallel_params.append((t, u, myact, index, k, param[3]))
            elif t == "log-normal":
                for u in U_LOG:
                    for myact, index in param[2]:
                        k += 1
                        parallel_params.append((t, u, myact, index, k, param[3]))
        
        max_workers = os.cpu_count() if os.cpu_count() else 4
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_case, (parallel_param, matrices)): (parallel_param, matrices) for parallel_param in parallel_params}
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    print(f'Generated an exception: {exc}')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total execution time(s): {elapsed_time}")
