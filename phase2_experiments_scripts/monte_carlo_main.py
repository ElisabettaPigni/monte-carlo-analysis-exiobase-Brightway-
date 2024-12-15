import pandas as pd
from monte_carlo_utility import *
import os
import sys
sys.path.append(os.path.abspath(".."))
from constants import *

simu = SimulationScript()

def prepare_datapackage_matrices(a_file, s_file, extend_file):
    # get all activities
    activities = simu.get_activities(a_file, delimiter='\t')
    activities = ["extra_column"] + activities  # activities with extra column.

    # background database technosphere matrix
    tech_df = pd.read_table(a_file, sep='\t', header=None, low_memory=False)
    raw_tech = tech_df.iloc[3:, 2:].astype('float').to_numpy()

    # add extra data to technosphere 
    extend_data_tech = pd.read_csv(extend_file, delimiter=";")
    extend_data_amount = extend_data_tech.iloc[:, :2]
    tech_matrix_extended = simu.extend_matrix(raw_tech, extend_data_amount, activities, is_technosphere=True)
    if not (raw_tech.shape[0]+1 == tech_matrix_extended.shape[0] and raw_tech.shape[1]+1 == tech_matrix_extended.shape[1]):
        print("Add column and row to technosphere failed!")

    tech_matrix = simu.form_tech_matrix(tech_matrix_extended)  # tech_data without extra column

    # biosphere matrix
    bio_df = pd.read_csv(s_file, header=[0,1], index_col=[0], sep='\t', low_memory=False)
    raw_bio = simu.form_bio_matrix(bio_df, GHG)

    # add extra data to biosphere
    extend_data_bio = pd.DataFrame([{"Exiobase_big_col (matrix B)": "N2O - combustion - air",
                                 "Amount": 5.62,
                                 "Exchange uncertainty type": 1,
                                 "Exchange loc": 1.7263316639056,
                                 "Exchange scale": 0,
                                 "Exchange negative": False}])
    extend_data_bio_amount = extend_data_bio.iloc[:, :2]
    bio_matrix = simu.extend_matrix(raw_bio, extend_data_bio_amount, GHG, is_technosphere=False)
    if not (raw_bio.shape[0] == bio_matrix.shape[0] and raw_bio.shape[1]+1 == bio_matrix.shape[1]):
        print("Add column and row to biosphere failed!")
    
    # characterization factor matrix
    cf_matrix = np.diagflat(CFS)

    return [tech_matrix, bio_matrix, cf_matrix, activities]

def create_static_datapackage(tech_matrix, bio_matrix, cf_matrix, activities):
    datapackage_data = simu.prepare_bw_matrix(tech_matrix, bio_matrix, cf_matrix, activities)
    datapackage = simu.prepare_datapackage(datapackage_data)
    
    return datapackage
    
def create_stochastic_datapackage(tech_matrix, bio_matrix, cf_matrix, activities, exiobase_type, extend_file):
    simu.metadata = [{activities.index(act): (act)} for act in activities]

    # prepare datapackage data
    datapackage_data = simu.prepare_bw_matrix(tech_matrix, bio_matrix, cf_matrix, ['extra_column'] + activities)
    tech_data, tech_indices, tech_flip = datapackage_data[0]
    bio_data, bio_indices = datapackage_data[1]

    # find and save pedigree undertainty
    if exiobase_type == "EXIOBASE":
        country_region, sector_seccat, region_sector_dfs = simu.map_pedigree_uncertainty(COUNTRY_FILE, SECTOR_FILE, GSD_FILE)
        for i in range(len(simu.metadata)):
            if 0 in simu.metadata[i]: # skip extra column
                continue
            
            for index, act in simu.metadata[i].items():
                gsd = simu.find_pedigree_uncertainty(act, country_region, sector_seccat, region_sector_dfs)
                simu.metadata[i][index] = (act, gsd)
    else:
        df = pd.read_csv(GSD_SMALL_FILE, delimiter=";")
        for i in range(len(simu.metadata)):
            if 0 in simu.metadata[i]: # skip extra column
                continue
            
            for index, act in simu.metadata[i].items():
                gsd = df[df[df.columns[0]] == act][df.columns[1]].iloc[0]
                simu.metadata[i][index] = (act, gsd)

    # find and save specific uncertainty
    uncertainty_spec_df = pd.read_csv(extend_file, delimiter=";")
    uncertainty_df = uncertainty_spec_df.loc[:, [uncertainty_spec_df.columns[0], uncertainty_spec_df.columns[3]]].copy()
    uncertainty_df = uncertainty_df.where(pd.notnull(uncertainty_df), 0)
    uncertainty_spec = np.zeros([len(activities)])
    for act in activities:
        if uncertainty_df[uncertainty_df.iloc[:, 0] == act][uncertainty_spec_df.columns[3]].empty:
            uncertainty_spec[activities.index(act)] = 0
        else:
            uncertainty_spec[activities.index(act)] = uncertainty_df[uncertainty_df.iloc[:, 0] == act][uncertainty_spec_df.columns[3]].iloc[0]
    for i in range(len(simu.metadata)):
        if 0 in simu.metadata[i]:
            for key, value in simu.metadata[i].items():
                simu.metadata[i][key] = (value, uncertainty_spec.tolist())

    # add pedigree and specific uncertainty
    tech_uncertainty = simu.add_uncertainty(tech_data, tech_indices, tech_flip)
    bio_uncertainty = simu.add_uncertainty(bio_data, bio_indices, None)
    
    # add multifunctionality negative
    extend_data_tech = pd.read_csv(extend_file, delimiter=";")
    extend_data_bio = pd.DataFrame([{"Exiobase_big_col (matrix B)": "N2O - combustion - air",
                                 "Amount": 5.62,
                                 "Exchange uncertainty type": 1,
                                 "Exchange loc": 1.7263316639056,
                                 "Exchange scale": 0,
                                 "Exchange negative": False}])
    tech_uncertainty = simu.add_multifunctionality_negative(extend_data_tech, extend_data_tech.columns[0], "Exchange negative", tech_uncertainty, tech_indices, activities)
    bio_uncertainty = simu.add_multifunctionality_negative(extend_data_bio, extend_data_bio.columns[0], "Exchange negative", bio_uncertainty, bio_indices, activities)

    datapackage = simu.prepare_datapackage(datapackage_data, uncertainty=[tech_uncertainty, bio_uncertainty, None])

    return datapackage

def perform_static(index, datapackage, directory, k, act, t):
    """
    Perform static simulation.
    """
    lca = bc.LCA(
        demand={index: 1},
        data_objs=[datapackage],
    )
    lca.lci()
    lca.lcia()

    os.makedirs(directory, exist_ok=True)
    filename = os.path.join(directory, f"CASE_{k}_{t}_MC_simulations_{act}.csv")

    with open(filename, "w") as file:
        file.write("kg CO2eq\n") # Write the header
        file.write(f"{lca.score}")
        print(f"Static LCA result saved to {filename}.")

def perform_stochastic(index, datapackage, directory, k, act, t):
    """
    Perform Monte Carlo simulation and save the lca score.
    """
    lca = bc.LCA(
        demand={index: 1},
        data_objs=[datapackage],
        use_distributions=True,
    )
    lca.lci()
    lca.lcia()

    print('LCA score (with uncertainty): ', lca.score)

    os.makedirs(directory, exist_ok=True)
    filename = os.path.join(directory, f"CASE_{k}_{t}_MC_simulations_{act}.csv")
    
    # TODO: these can be passed as parameters.
    batch_size = 50
    num_batches = 10

    with open(filename, "w") as file:
        file.write("kg CO2eq\n")
        for p in range(num_batches):
            batch_results = [lca.score for _ in zip(range(batch_size), lca)]
            df_batch = pd.DataFrame(batch_results, columns=["kg CO2eq"])
            df_batch.to_csv(file, header=False, index=False)
            print(f"Batch {p} saved to {filename}.")

    print(f"Results saved to {filename}.")

# TODO
def manual_lca(tech_data, bio_data, cf_matrix, act):
    manual_tech_Data = -tech_data
    np.fill_diagonal(manual_tech_Data, -manual_tech_Data.diagonal())
    # TODO: because all data is positive in dp, but is it the same when using matrix calculation.
    lca_score_manual = simu.manual_lca(manual_tech_Data, bio_data, cf_matrix, act[1])
    print(f"Manually calculated lca score: {lca_score_manual}")

# This is the static simulation case with only exiobase background system.
def lca_background_system(a_file, s_file, chosen_act):
    # technosphere matrix
    tech_df = pd.read_table(a_file, sep='\t', header=None, low_memory=False)
    raw_tech = tech_df.iloc[3:, 2:].astype('float').to_numpy()
    tech_matrix = simu.form_tech_matrix(raw_tech)

    # biosphere matrix
    bio_df = pd.read_csv(s_file, header=[0,1], index_col=[0], sep='\t', low_memory=False)
    bio_matrix = simu.form_bio_matrix(bio_df, GHG)
    if tech_matrix.shape[1] != bio_matrix.shape[1] or bio_matrix.shape[0] != len(GHG):
        print("The shape of matrices don't match!")

    # characterization factor matrix
    cf_matrix = np.diagflat(CFS)

    # create datapackage
    datapackage_data = simu.prepare_bw_matrix(tech_matrix, bio_matrix, cf_matrix)
    datapackage = simu.prepare_datapackage(datapackage_data)
    
    # LCA simulation
    for act in chosen_act:
        lca_score = simu.perform_simulation(act[1], datapackage)
        print(f"Brightway calculated lca score: {lca_score}")

    return lca_score
