from constants import *
import time
from monte_carlo_utility import *


if __name__ == "__main__":

    simu = SimulationScript()

    start_time = time.time()

    # ---------- RUN SIMULATION ---------- 
    
    # # Form matrices for bw
    # A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = simu.build_bw_matrix(BIG_A_FILE, BIG_S_FILE)
    # print("Matrices are formatted.")

    # # Run the simulation
    # k = 0
    # for t in DIST_TYPE:
    #     # This is the baseline case
    #     if t == "baseline":
    #         for myact, index in BIG_CHOSEN_ACT:
    #             k += 1
    #             print(f"----------- Starting CASE {k}, activity: {myact} -----------")
    #             lca = simu.perform_baseline(index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C, BIG_DIR_OUTPUT, t)

    #             print(f"CASE {k} simulation is done.")

    #     # This is the uniform case
    #     elif t == "uniform":
    #         for u in U_UNIFORM:
    #             for myact, index in BIG_CHOSEN_ACT:
    #                 k += 1
    #                 print(f"----------- Starting CASE {k}, activity: {myact} -----------")

    #                 # Add uncertainty
    #                 dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
    #                 print(f"Uncertainty added. ({t} distribution with uncertainty {u})")

    #                 # Perform lca
    #                 lca = simu.perform_simu(index, dp_stochastic, BIG_DIR_OUTPUT, k, myact, t, u)

    #                 print(f"CASE {k} simulation is done.")
    #     # This is the log-normal case
    #     elif t == "log-normal":
    #         for u in U_LOG:
    #             for myact, index in BIG_CHOSEN_ACT:
    #                 k += 1
    #                 print(f"----------- Starting CASE {k}, activity: {myact} -----------")

    #                 # Add uncertainty
    #                 dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
    #                 print(f"Uncertainty added. ({t} distribution with uncertainty {u})")

    #                 lca = simu.perform_simu(index, dp_stochastic, BIG_DIR_OUTPUT, k, myact, t, u)

    #                 print(f"CASE {k} simulation is done.")

    # print("All simulations completed.")

    # ---------- DRAW PLOTS ---------- 
    # generate plots for small dataset
    data = simu.get_plot(FOLDER_PATH, DB_SIZE[0])
    simu.draw_plot(data, PLT_COMP_TYPE[0], DB_SIZE[0], PLT_SAVE_PATH)
    simu.draw_plot(data, PLT_COMP_TYPE[1], DB_SIZE[0], PLT_SAVE_PATH)

    # generate plots for big dataset
    data = simu.get_plot(FOLDER_PATH, DB_SIZE[1])
    simu.draw_plot(data, PLT_COMP_TYPE[0], DB_SIZE[1], PLT_SAVE_PATH)
    simu.draw_plot(data, PLT_COMP_TYPE[1], DB_SIZE[1], PLT_SAVE_PATH)

    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_hour = round(elapsed_time / 60 / 60, 2)
    print(f"Total execution time(s): {elapsed_time}")
    # print(f"Total execution time(h): {elapsed_hour}")