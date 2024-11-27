from constants import *
import pandas as pd
from monte_carlo_utility_new import *


simu = SimulationScript()

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

    # characterization factor matrix
    # ATTENTION: have to make sure the file has the same order as biosphere emissions.
    # cf_matrix = simu.form_cf_matrix("EXIOBASE-ecoinvent-bio-bw-GHG.csv", METHOD)
    cf_matrix = np.diagflat(OLD_CFS)

    simu.save_metadata(activities, "technosphere")
    simu.save_metadata(activities, "biosphere")

    return tech_matrix, bio_matrix, cf_matrix

def data_prepare_addional_data():
    # technosphere matrix
    tech_df = pd.read_table("exiobase_2022_small/A.txt", sep='\t', low_memory=False, header=None)
    raw_tech = tech_df.iloc[3:, 2:].astype('float').to_numpy()
    activities = simu.get_activities(tech_df.iloc[3:, 0], tech_df.iloc[3:, 1])
    tech_matrix = simu.form_tech_matrix(raw_tech)

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

    # characterization factor matrix
    # ATTENTION: have to make sure the file has the same order as biosphere emissions.
    # cf_matrix = simu.form_cf_matrix("EXIOBASE-ecoinvent-bio-bw-GHG.csv", METHOD)
    cf_matrix = np.diagflat(CFS)

    # add the additional column name to record the index
    simu.save_metadata(["extra"] + activities, "technosphere")
    simu.save_metadata(["extra"] + activities, "biosphere")

    return tech_matrix_new, bio_matrix_new, cf_matrix

def create_datapackages(tech_matrix, bio_matrix, cf_matrix):
    datapackage_data = simu.prepare_bw_matrix(tech_matrix, bio_matrix, cf_matrix)
    datapackage = simu.prepare_datapackage(datapackage_data)
                       
    return datapackage

def run_experiments(datapackage):
    for act in SMALL_CHOSEN_ACT:
        lca_score = simu.perform_simulation(act[1], datapackage)
        print(lca_score)


if __name__ == "__main__":
    # baseline, exiobase small
    tech_matrix, bio_matrix, cf_matrix = data_prepare()
    datapackage = create_datapackages(tech_matrix, bio_matrix, cf_matrix)
    run_experiments(datapackage)