from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
import time
import os
from constants import *
from monte_carlo_utility import *


def process_case(parallel_param):
    """
    parallel_param type: tuple
    parallel_param format: ((t, "0", myact, index, k), matrices)
    """
    print(f"Processing {t}_{u} in process: {os.getpid()}")

    (t, u, myact, index, k) = parallel_param[0]
    matrices = parallel_param[1]
    A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = matrices

    simu = SimulationScript()

    if t == "baseline":
        simu.perform_baseline(index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C, BIG_DIR_OUTPUT, k, t, myact)
        print(f"{t}_{u} simulation is done.")
    elif t == "uniform":
        dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
        simu.perform_simu(index, dp_stochastic, BIG_DIR_OUTPUT, k, myact, t, u)
        print(f"{t}_{u} simulation is done.")
    elif t == "log-normal":
        dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
        simu.perform_simu(index, dp_stochastic, BIG_DIR_OUTPUT, k, myact, t, u)
        print(f"{t}_{u} simulation is done.")


if __name__ == "__main__":
    start_time = time.time()

    simu = SimulationScript()

    A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = simu.build_bw_matrix(BIG_A_FILE, BIG_S_FILE)
    print("Matrices are formatted.")

    matrices = (A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)

    parallel_params = []
    k = 0
    for t in DIST_TYPE:
        if t == "baseline":
            for myact, index in BIG_CHOSEN_ACT:
                k += 1
                parallel_params.append((t, 0, myact, index, k))
        elif t == "uniform":
            for u in U_UNIFORM:
                for myact, index in BIG_CHOSEN_ACT:
                    k += 1
                    parallel_params.append((t, u, myact, index, k))
        elif t == "log-normal":
            for u in U_LOG:
                for myact, index in BIG_CHOSEN_ACT:
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
