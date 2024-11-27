from constants import *
import time
import pandas as pd
# from monte_carlo_utility import *
from new_script import *

simu = SimulationScript_test()

def data_prepare():
    # technosphere matrix
    tech_df = pd.read_table("exiobase_2022_small/A.txt", sep='\t', low_memory=False, header=None)
    raw_tech = tech_df.iloc[3:, 2:].astype('float').to_numpy()
    activities = simu.get_activities(tech_df.iloc[3:, 0], tech_df.iloc[3:, 1])
    tech_matrix = simu.form_tech_matrix(raw_tech)

    # biosphere matrix
    bio_df = pd.read_csv("exiobase_2022_small/satellite/S.txt", header=[0,1], index_col=[0], sep='\t', low_memory=False)
    bio_matrix = simu.form_bio_matrix(bio_df, GHG)
    if tech_matrix.shape[1] != bio_matrix.shape[1] or bio_matrix.shape[0] != len(GHG):
        print("The shape of matrices don't match!")

    # ATTENTION: have to make sure the file has the same order as biosphere emissions.
    # cf_matrix = simu.form_cf_matrix("EXIOBASE-ecoinvent-bio-bw-GHG.csv", METHOD)
    cf_matrix = cf_matrix = np.diagflat(OLD_CFS)

    simu.save_metadata(["extra"] + activities, "technosphere")
    simu.save_metadata(GHG, "biosphere")

    return tech_matrix, bio_matrix, cf_matrix

def data_prepare_addional_data():
    # technosphere matrix
    tech_df = pd.read_table("exiobase_2022_small/A.txt", sep='\t', low_memory=False, header=None)
    raw_tech = tech_df.iloc[3:, 2:].astype('float').to_numpy()
    activities = simu.get_activities(tech_df.iloc[3:, 0], tech_df.iloc[3:, 1])
    tech_matrix, a = simu.form_tech_matrix(raw_tech)

    extend_data = pd.read_csv("foreground_system_small_changed.csv").iloc[:, :2]
    tech_matrix_new = simu.extend_matrix(tech_matrix, extend_data, activities, is_technosphere=True)
    if not (tech_matrix_new.shape[0] == tech_matrix.shape[0]+1 and tech_matrix_new.shape[1] == tech_matrix.shape[1]+1):
        print("Add column and row to technosphere failed!")

    # biosphere matrix
    bio_df = pd.read_csv("exiobase_2022_small/satellite/S.txt", header=[0,1], index_col=[0], sep='\t', low_memory=False)
    bio_matrix = simu.form_bio_matrix(bio_df, GHG)

    extend_data = pd.DataFrame([{"Exiobase_big_col (matrix B)": "N2O - combustion - air",
                                 "Amount": 5.62,
                                 "Exchange uncertainty type": 1,
                                 "Exchange loc": 1.7263316639056,
                                 "Exchange scale": 0,
                                 "Exchange negative": False}]).iloc[:, :2]
    bio_matrix_new = simu.extend_matrix(bio_matrix, extend_data, GHG, is_technosphere=False)
    if not (bio_matrix_new.shape[0] == bio_matrix.shape[0] and bio_matrix_new.shape[1] == bio_matrix.shape[1]+1):
        print("Add column and row to biosphere failed!")

    # ATTENTION: have to make sure the file has the same order as biosphere emissions.
    # cf_matrix = simu.form_cf_matrix("EXIOBASE-ecoinvent-bio-bw-GHG.csv", METHOD)
    cf_matrix = CFS

    simu.save_metadata(activities, "technosphere")
    simu.save_metadata(GHG, "biosphere")

    return tech_matrix_new, bio_matrix_new, cf_matrix

def create_datapackages(tech_matrix, bio_matrix, cf_matrix):
    datapackage_data = simu.prepare_bw_matrix(tech_matrix, bio_matrix, cf_matrix)
    datapackage = simu.prepare_datapackage(datapackage_data)
                       
    return datapackage

def run_experiments(datapackage):
    for act in SMALL_CHOSEN_ACT:
        print(act[1])
        simu.perform_simulation(act[1], datapackage)


if __name__ == "__main__":
    # original baseline
    tech_matrix, bio_matrix, cf_matrix = data_prepare()
    datapackage = create_datapackages(tech_matrix, bio_matrix, cf_matrix)
    run_experiments(datapackage)

    # # additional data baseline
    # tech_matrix, bio_matrix, cf_matrix = data_prepare_addional_data()
    # datapackage = create_datapackages(tech_matrix, bio_matrix, cf_matrix)
    # run_experiments(datapackage)

    # # additional data MC simulation
    # tech_matrix, bio_matrix, cf_matrix = data_prepare_addional_data()
    # datapackage = create_datapackages(tech_matrix, bio_matrix, cf_matrix)
    # run_experiments(datapackage)


    # data_prepare_addional_data()
    # simu = SimulationScript()

    # start_time = time.time()

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

    # # ---------- DRAW PLOTS ---------- 
    # # generate plots for small dataset
    # data = simu.get_plot(FOLDER_PATH, DB_SIZE[0])
    # simu.draw_plot(data, PLT_COMP_TYPE[0], DB_SIZE[0], PLT_SAVE_PATH)
    # simu.draw_plot(data, PLT_COMP_TYPE[1], DB_SIZE[0], PLT_SAVE_PATH)

    # # generate plots for big dataset
    # data = simu.get_plot(FOLDER_PATH, DB_SIZE[1])
    # simu.draw_plot(data, PLT_COMP_TYPE[0], DB_SIZE[1], PLT_SAVE_PATH)
    # simu.draw_plot(data, PLT_COMP_TYPE[1], DB_SIZE[1], PLT_SAVE_PATH)

    # end_time = time.time()
    # elapsed_time = end_time - start_time
    # elapsed_hour = round(elapsed_time / 60 / 60, 2)
    # print(f"Total execution time(s): {elapsed_time}")
    # # print(f"Total execution time(h): {elapsed_hour}")