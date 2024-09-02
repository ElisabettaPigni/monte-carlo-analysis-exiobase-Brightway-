import bw_processing as bwp #this is needed for brightway 2.5 package
import bw2calc as bc #this is needed for brightway 2.5 package
import bw2data as bd #this is needed for brightway 2.5 package
from scipy import sparse #This is necessary to create the sparse matrix, which is a lighter matrix in which zero values have been removed.
import pandas as pd
import numpy as np
from random import sample
import random
import os
from statistic_analysis import StatisticAnalysis


class SimulationScript:
    # Get all activities
    def get_activities(self, A_file_path):
        A_raw = pd.read_table(A_file_path)
        countries = A_raw['region'].drop_duplicates().iloc[2:].tolist()
        sectors = list(A_raw.iloc[2:,1].drop_duplicates())
        activities = [ x + '-' + y for x in countries for y in sectors]

        return activities


    # Choose activities for experiments
    def choose_activities(self, activities, amount):
        chosen_activities = []
        indices = sample(range(len(activities) - 1), amount)

        chosen_activities = [(activities[i], i) for i in indices]

        return chosen_activities
    

    # Build matrix adapt to bw
    def build_bw_matrix(self, A_file_path, S_file_path):
        activities = self.get_activities(A_file_path) # Get all activities

        A_raw = pd.read_table(A_file_path)
        A_IO = A_raw.iloc[2:,2:].astype('float').values
        I = np.identity(len(A_IO))
        A = A_IO - I
        np.fill_diagonal(A, -A.diagonal()) # then change back again, but only the diagonal

        Asparse = sparse.coo_array(A) # technology matrix A as sparse object and then coordinates
        a_data = Asparse.data # amounts or values
        a_indices = np.array([tuple(coord) for coord in np.transpose(Asparse.nonzero())], dtype=bwp.INDICES_DTYPE) # indices of each exchange
        a_flip = np.array([False if i[0] == i[1] else True for i in a_indices ]) # Numerical sign of the inputs needs to be flipped negative

        # import environemntal extensions
        S_raw = pd.read_table(S_file_path, header=[0,1], index_col=[0])

        GHG_rows = ["CO2 - combustion - air",
                    "CO2 - non combustion - Cement production - air",
                    "CO2 - non combustion - Lime production - air",
                    "CO2 - waste - biogenic - air", 
                    "CO2 - waste - fossil - air",
                    "CO2 - agriculture - peat decay - air", 
                    "CH4 - agriculture - air",
                    "CH4 - waste - air",
                    "CH4 - combustion - air",
                    "CH4 - non combustion - Extraction/production of (natural) gas - air",
                    "CH4 - non combustion - Extraction/production of crude oil - air",
                    "CH4 - non combustion - Mining of antracite - air",
                    "CH4 - non combustion - Mining of bituminous coal - air",
                    "CH4 - non combustion - Mining of coking coal - air",
                    "CH4 - non combustion - Mining of lignite (brown coal) - air",
                    "CH4 - non combustion - Mining of sub-bituminous coal - air",
                    "CH4 - non combustion - Oil refinery - air",
                    "N2O - combustion - air",
                    "N2O - agriculture - air",
                    "HFC - air",
                    "SF6 - air"]

        S = S_raw.loc[GHG_rows]

        B = S.to_numpy()

        Bsparse = sparse.coo_array(B) # Intervention matrix B as sparse object and then coordinates
        b_data = Bsparse.data # amounts or values
        b_indices_remap = [[i[0] + len(activities),i[1]] for i in np.transpose(Bsparse.nonzero())] # need to make sure biosphere indices are different from technosphere
        b_indices = np.array([tuple(coord) for coord in b_indices_remap], dtype=bwp.INDICES_DTYPE)

        # matrix of charachterisation factors (method E.F. 3.1 - Climate change): CFs = [bw.ef.co2, bw.sox, ..., ]
        CFs = [1., 29.8, 273., 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 1., 1., 25200., 14600., 29.8, 1., 273., 29.8, 1., 1.]

        C = np.matrix(np.zeros((len(CFs), len(CFs))))
        C_diag = np.matrix(CFs)
        np.fill_diagonal(C, C_diag)

        Csparse = sparse.coo_array(C) # Sparse C  matrix of characterisation factors
        c_data =  Csparse.data 
        c_indices_remap = [[i[1] + len(activities),i[1]+ len(activities)] for i in np.transpose(Csparse.nonzero())] # same indices as in B
        c_indices = np.array([tuple(coord) for coord in c_indices_remap], dtype=bwp.INDICES_DTYPE) 

        return A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip
    

    # Perform baseline simulation
    def perform_baseline(self, myact, index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C):
        # Creating a datapackage
        dp_static = bwp.create_datapackage()
        dp_static.add_persistent_vector(
            matrix='technosphere_matrix',
            indices_array=a_indices,
            data_array=a_data,
            flip_array=a_flip,
        )
        dp_static.add_persistent_vector(
            matrix='biosphere_matrix',
            indices_array=b_indices,
            data_array=b_data,
        )
        dp_static.add_persistent_vector(
            matrix='characterization_matrix',
            indices_array=c_indices,
            data_array=c_data,
        )

        lca = bc.LCA(
            demand={activities.index(myact): 1}, # using index because the argument of FU is an INTEGER
            data_objs=[dp_static],
        )
        lca.lci()
        lca.lcia()

        # check with normal matrix operation
        f = np.zeros(len(A))
        f[index] = 1 # functional unit

        print('ioscore',np.sum(C.dot(B.dot((np.linalg.inv(A_)).dot(f))))) # matrix operation
        print('LCA score: ',lca.score)  # similar values, when rounded. All is good

        return lca


    # Ready the matries for simulation
    def add_uncertainty(self, t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip):
        if t == "uniform":
            results_a = []
            results_b = []

            #set the uncertainty for matrix A
            min_val_a = a_data - (a_data * u)
            max_val_a = a_data + (a_data * u)

            # Iterate over indices and check corresponding values in a_flip
            for i in range(len(a_data)):
                if not a_flip[i]:  # If a_flip at index i is False
                    parameters_a = (0, a_data[i], np.NaN, np.NaN, np.NaN, np.NaN, False)
                else:  # If a_flip at index i is True
                    parameters_a = (4, np.NaN, np.NaN, np.NaN, min_val_a[i], max_val_a[i], False)
                results_a.append(parameters_a)
            
            a_uncertainty = np.array(results_a, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )

            min_val_b = b_data - (b_data * u)
            max_val_b = b_data + (b_data * u)

            for i in range(len(b_data)):
                parameters_b = (4, np.NaN, np.NaN, np.NaN, min_val_b[i], max_val_b[i], False)
                results_b.append(parameters_b)
            
            b_uncertainty = np.array(results_b, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )
        elif t == "log-normal":
            results_a = []
            results_b = []

            #add uncertainty for A
            mu_a = np.log(a_data)
            sigma = np.log(u)

            # Iterate over indices and check corresponding values in a_flip
            for i in range(len(a_data)):
                if not a_flip[i]:  # If a_flip at index i is False
                    parameters_a = (0, a_data[i], np.NaN, np.NaN, np.NaN, np.NaN, False)
                else:  # If a_flip at index i is True
                    parameters_a = (2, mu_a[i], sigma, np.NaN, np.NaN, np.NaN, False)
                results_a.append(parameters_a)

            a_uncertainty = np.array(results_a, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )

            #add uncertainty for B
            mu_b = np.log(b_data)

            # Iterate over indices and check corresponding values in a_flip
            for i in range(len(b_data)):
                parameters_b = (2, mu_b[i], sigma, np.NaN, np.NaN, np.NaN, False)
                results_b.append(parameters_b)

            b_uncertainty = np.array(results_b, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )

        dp_stochastic = bwp.create_datapackage()

        dp_stochastic.add_persistent_vector(
            matrix='technosphere_matrix',
            indices_array=a_indices,
            data_array=a_data,
            flip_array=a_flip,
            distributions_array=a_uncertainty,
        )
        dp_stochastic.add_persistent_vector(
            matrix='biosphere_matrix',
            indices_array=b_indices,
            data_array=b_data,
            distributions_array=b_uncertainty,
        )
        dp_stochastic.add_persistent_vector(
            matrix='characterization_matrix',
            indices_array=c_indices,
            data_array=c_data,
        )

        return dp_stochastic
    

    # Perform Monte Carlo simulation
    def perform_simu(self, myact, dp_stochastic, activities):
        lca = bc.LCA(
            demand={activities.index(myact): 1},
            data_objs=[dp_stochastic],
            use_distributions=True,
        )
        lca.lci()
        lca.lcia()

        print('LCA score: ',lca.score)
        
        return lca
    

    # Save results to csv
    def save_result(self, directory, lca, k, myact):
        os.makedirs(directory, exist_ok=True)  # Create the directory if it does not exist
        filename = os.path.join(directory, f"CASE_{k}_MC_simulations_{myact}.csv") # Define filename for saving results
        
        # Define simulation parameters
        batch_size = 50
        num_batches = 10

        with open(filename, "w") as file:
            file.write("kg CO2eq\n") # Write the header
            for p in range(num_batches):
                # Run simulations for the current batch
                batch_results = [lca.score for _ in zip(range(batch_size), lca)]
                
                # Convert to DataFrame
                df_batch = pd.DataFrame(batch_results, columns=["kg CO2eq"])
                
                # Save batch results to CSV, appending to the file
                df_batch.to_csv(file, header=False, index=False)
            
                # Print progress
                print(f"Batch {p} saved to {filename}.")

        print(f"Results saved to {filename}.")


if __name__ == "__main__":
    # --------------------- Configuration --------------------- 
    dist_type = ["baseline", "uniform", "log-normal"] # Define the types of distribution
    u_uniform = [0.1, 0.2, 0.3, 0.5] # Define the uncertainty for uniform distribution
    u_log = [1.01, 1.1, 2] # Define the uncertainty for log distribution
    amount = 4 # This is the amount of activities for 1 CASE

    # File paths
    dir_input = f"{os.getcwd()}/exiobase_2022_small"
    dir_output = f"{os.getcwd()}/output"
    A_file_path = f"{os.getcwd()}/exiobase_2022_small/A.txt"
    S_file_path = f"{os.getcwd()}/exiobase_2022_small/satellite/S.txt"
    
    # Make sure execute in the right directory
    if os.getcwd() != dir_input:
        os.chdir(dir_input)
    print("Directed to the right directory.")


    # --------------------- Simulation --------------------- 
    simu = SimulationScript()

    # Choose activities (If you want to change the chosen activities, please uncomment the following code.)
    # activities = simu.get_activities(A_file_path)
    # chosen_activities = simu.choose_activities(activities, amount)
    # print("The following activities are chosen:")
    # for myact, _ in chosen_activities:
    #     print(myact, end=",")
    chosen_activities = ["RoW-Services", "EU28-Biodiesels", "EU28-Agriculture-Forestry-Fishing", "EU28-Basic iron and steel and of ferro-alloys and first products thereof"]
    
    # Adapt matrices for bw
    A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = simu.build_bw_matrix(A_file_path, S_file_path)
    print("Matrices are formatted.")

    # Run the simulation
    k = 0
    for t in dist_type:
        # This is the baseline case
        if t == "baseline":
            for myact, index in chosen_activities:
                k += 1
                print(f"----------- Starting CASE {k}, activity: {myact} -----------")
                lca = simu.perform_baseline(myact, index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C)
                simu.save_result(dir_output, lca, k, myact)

                print(f"CASE {k} simulation is done.")

        # This is the uniform case
        elif t == "uniform":
            for u in u_uniform:
                for myact, index in chosen_activities:
                    k += 1
                    print(f"----------- Starting CASE {k}, activity: {myact} -----------")

                    # Add uncertainty
                    dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
                    print(f"Uncertainty added to CASE {k}")

                    # Perform lca
                    lca = simu.perform_simu(myact, dp_stochastic, activities)
                    simu.save_result(dir_output, lca, k, myact)

                    print(f"CASE {k} simulation is done.")
        # This is the log-normal case
        elif t == "log-normal":
            for u in u_log:
                for myact, index in chosen_activities:
                    k += 1
                    print(f"----------- Starting CASE {k}, activity: {myact} -----------")
                    # Add uncertainty
                    dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
                    print(f"Uncertainty added. ({t} distribution with uncertainty {u})")

                    lca = simu.perform_simu(myact, dp_stochastic, activities)
                    simu.save_result(dir_output, lca, k, myact)

                    print(f"CASE {k} simulation is done.")


    print("All simulations completed.")


    # --------------------- Statistic Analysis --------------------- 
    # stat = StatisticAnalysis()
    # stat.simu_plot()
    # stat.matrix_plot(A, A_IO, B)