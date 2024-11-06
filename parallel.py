'''
 # @ Author: Ning An
 # @ Create Time: 2024-11-05 10:45:51
 # @ Modified by: Ning An
 # @ Modified time: 2024-11-06 10:33:51
 '''

from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import os
from constants import *
from monte_carlo_utility import *


def process_case(parallel_param):
    """
    parallel_param type: tuple
    parallel_param format: ((t, "0", myact, index, k), matrices)
    """
    (t, u, myact, index, k) = parallel_param[0]
    print(f"Processing {t}_{u} in process: {os.getpid()}")

    matrices = parallel_param[1]
    A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = matrices

    simu = SimulationScript()

    if t == "baseline":
        simu.perform_baseline(index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C, SMALL_DIR_OUTPUT, k, t, myact)
        print(f"{t}_{u} simulation is done.")
    elif t == "uniform":
        dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
        simu.perform_simu(index, dp_stochastic, SMALL_DIR_OUTPUT, k, myact, t, u)
        print(f"{t}_{u} simulation is done.")
    elif t == "log-normal":
        dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
        simu.perform_simu(index, dp_stochastic, SMALL_DIR_OUTPUT, k, myact, t, u)
        print(f"{t}_{u} simulation is done.")


if __name__ == "__main__":
    start_time = time.time()

    simu = SimulationScript()

    A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = simu.build_bw_matrix(SMALL_A_FILE, SMALL_S_FILE)
    print("Matrices are formatted.")

    matrices = (A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)

    parallel_params = []
    k = 0
    for t in DIST_TYPE:
        if t == "baseline":
            for myact, index in SMALL_CHOSEN_ACT:
                k += 1
                parallel_params.append((t, 0, myact, index, k))
        elif t == "uniform":
            for u in U_UNIFORM:
                for myact, index in SMALL_CHOSEN_ACT:
                    k += 1
                    parallel_params.append((t, u, myact, index, k))
        elif t == "log-normal":
            for u in U_LOG:
                for myact, index in SMALL_CHOSEN_ACT:
                    k += 1
                    parallel_params.append((t, u, myact, index, k))
    
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
