import bw_processing as bwp #this is needed for brightway 2.5 package
import bw2calc as bc #this is needed for brightway 2.5 package
import bw2data as bd #this is needed for brightway 2.5 package
from scipy import sparse #This is necessary to create the sparse matrix, which is a lighter matrix in which zero values have been removed.
import pandas as pd
import numpy as np
import os
from statistic_analysis import StatisticAnalysis


class SimulationScript:
    def __init__(self):
        pass


    # Build matrix adapt to bw
    def build_bw_matrix(self, A_file_path, S_file_path):
        A_raw = pd.read_table(A_file_path)

        countries = A_raw['region'].drop_duplicates().iloc[2:].tolist()
        sectors = list(A_raw.iloc[2:,1].drop_duplicates())
        activities = [ x + '-' + y for x in countries for y in sectors]
        units = (['MEUR' for i in range(0,len(activities))])
        extensions = ["349b29d1-3e58-4c66-98b9-9d1a076efd2e", "0795345f-c7ae-410c-ad25-1845784c75f5", "20185046-64bb-4c09-a8e7-e8a9e144ca98", "0795345f-c7ae-410c-ad25-1845784c75f5",
        "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5",
        "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "349b29d1-3e58-4c66-98b9-9d1a076efd2e",
        "349b29d1-3e58-4c66-98b9-9d1a076efd2e", "35d1dff5-b535-4628-9826-4a8fce08a1f2", "e5ba3517-a93f-422e-9ce2-f1007fcf6c06", "0795345f-c7ae-410c-ad25-1845784c75f5",
        "349b29d1-3e58-4c66-98b9-9d1a076efd2e", "20185046-64bb-4c09-a8e7-e8a9e144ca98", "0795345f-c7ae-410c-ad25-1845784c75f5", "349b29d1-3e58-4c66-98b9-9d1a076efd2e", "349b29d1-3e58-4c66-98b9-9d1a076efd2e"
        ]

        A_IO = A_raw.iloc[2:,2:].astype('float').values
        I = np.identity(len(A_IO))
        A_ = I - A_IO
        A = -A_
        np.fill_diagonal(A, -A.diagonal()) # then change back again, but only the diagonal

        Asparse = sparse.coo_array(A) # technology matrix A as sparse object and then coordinates
        a_data = Asparse.data # amounts or values
        a_indices = np.array([tuple(coord) for coord in np.transpose(Asparse.nonzero())], dtype=bwp.INDICES_DTYPE) # indices of each exchange
        a_flip = np.array([False if i[0] == i[1] else True for i in a_indices ]) # Numerical sign of the inputs needs to be flipped negative

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

        # matrix of charachterisation factors (method E.F. 3.1 - Climate change)
        CFs = [1., 29.8, 273., 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 1., 1., 25200., 14600., 29.8, 1., 273., 29.8, 1., 1.]

        C = np.matrix(np.zeros((len(CFs), len(CFs))))
        C_diag = np.matrix(CFs)
        np.fill_diagonal(C, C_diag)

        Csparse = sparse.coo_array(C) # Sparse C  matrix of characterisation factors
        c_data =  Csparse.data 
        c_indices_remap = [[i[1] + len(activities),i[1]+ len(activities)] for i in np.transpose(Csparse.nonzero())] # same indices as in B
        c_indices = np.array([tuple(coord) for coord in c_indices_remap], dtype=bwp.INDICES_DTYPE) 

        return activities, A, A_IO, B, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip
    
    
    # Perform baseline simulation
    def normal_lca(self, activities, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C):
        # Creating the datapackage
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

        myact = activities[1] # to find the desired code for the functional unit
        print(myact)

        lca = bc.LCA(
            demand={activities.index(myact): 1}, # using index because the argument of FU is an INTEGER
            data_objs=[dp_static],
        )
        lca.lci()
        lca.lcia()
        lca.score

        f = np.zeros(len(A)) # check with normal matrix operation
        f[1] = 1 # functional unit

        print('ioscore',np.sum(C.dot(B.dot((np.linalg.inv(A_)).dot(f))))) # matrix operation
        print('lca score',lca.score)  # similar values, when rounded. All is good


    # -------------------------------------- STOCHASTIC LCA ----------------------------------------------
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
    def perform_simu(self, dp_stochastic, activities):
        myact = activities[1] # to find the desired code for the functional unit
        print(f"myact: {myact}")

        lca = bc.LCA(
            demand={activities.index(myact): 1},
            data_objs=[dp_stochastic],
            use_distributions=True,
        )
        lca.lci()
        lca.lcia()
        
        return lca


    # Save results to csv
    def save_result(self, directory, lca, k):
        os.makedirs(directory, exist_ok=True)  # Create the directory if it does not exist
        filename = os.path.join(directory, f"CASE_{k}_MC_simulations.csv") # Define filename for saving results
        
        # Define simulation parameters
        batch_size = 50
        num_batches = 10

        file.write("kg CO2eq\n") # Write the header

        with open(filename, "w") as file:
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
    # Monte Carlo simulation Description:
        # CASE 1: uniform distribution with 10% uncertainty
        # CASE 2: uniform distribution with 20% uncertainty
        # CASE 3: uniform distribution with 30% uncertainty
        # CASE 4: log-normal distribution with 1.01 uncertainty
        # CASE 5: log-normal distribution with 1.1 uncertainty 
        # CASE 6: log-normal distribution with 2 uncertainty 

    u_uniform = [0.1, 0.2, 0.3, 0.5]
    u_log = [1.01, 1.1, 2]
    dist_type = ["uniform", "log-normal", "baseline"]

    # File paths
    dir_input = "/Users/bp45th/Documents/github/monte-carlo-analysis-exiobase-Brightway-/exiobase_2022_small"
    dir_output = "/Users/bp45th/Documents/github/monte-carlo-analysis-exiobase-Brightway-/output"
    A_file_path = "/Users/bp45th/Documents/github/monte-carlo-analysis-exiobase-Brightway-/exiobase_2022_small/A.txt"
    S_file_path = "/Users/bp45th/Documents/github/monte-carlo-analysis-exiobase-Brightway-/exiobase_2022_small/satellite/S.txt"
    
    # Make sure execute in the right directory
    if os.getcwd() != dir_input:
        os.chdir(dir_input)
    
    print("Directed to the right directory.")

    simu = SimulationScript()

    # Adapt matrices for bw
    activities, A, A_IO, B, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip = simu.build_bw_matrix(A_file_path, S_file_path)
    print("Matrices is formatted.")


    # --------------------- Simulation --------------------- 
    k = 0
    for t in dist_type:
        # This is the uniform case
        if t == "uniform":
            for u in u_uniform:
                k += 1
                print(f"Starting CASE {k}.")

                # Add uncertainty
                dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
                print(f"Uncertainty added to CASE {k}")

                # Perform lca
                lca = simu.perform_simu(dp_stochastic, activities)
                print(f"CASE {k} simulation is done.")

                simu.save_result(dir_output, lca, k)

        # This is the log-normal case
        elif t == "log-normal":
            for u in u_log:
                k += 1
                print(f"Starting CASE {k}.")

                # Add uncertainty
                dp_stochastic = simu.add_uncertainty(t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip)
                print(f"Uncertainty added. ({t} distribution with uncertainty {u})")

                lca = simu.perform_simu(dp_stochastic, activities)
                print(f"CASE {k} simulation is done.")

                simu.save_result(dir_output, lca, k)

        elif t == "normal_lca":
            pass

    print("All simulations completed.")


    # --------------------- Statistic Analysis --------------------- 
    stat = StatisticAnalysis()
    stat.simu_plot()
    stat.matrix_plot(A, A_IO, B)